import numpy as np
import pandas as pd
import emcee
import corner
import ev
from ev.run import run_montecarlo
import radvel
import scipy.stats as ss
import modules_KOBEsim.SimulTime_KOBEsim as SimulTime
import os
import math
from astropy.time import Time
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch


#=====================================
# Model and functions to run the MCMC
#=====================================


"""
Parameters of the RV model:
---------------------------
- Vsys: systemic velocity (m/s)
- P: orbital period (d)
- K: RV semi-amplitude (m/s)
- t0: JD of conjuction (d)
- e: eccentricity
- w: argument of the periastron (rad)
- m: parameter of the linear trend
- q: parameter of the quadratic trend
"""


def model_RV(Vsys, P, K, t0, e, w, m, q, t, planet, wh):
    RV = Vsys
    if planet:
        tperi = radvel.orbit.timetrans_to_timeperi(t0, P, e, w)
        keplerian = radvel.kepler.rv_drive(np.array(t), [P, tperi, e, w, K])
        RV += keplerian
        if wh:
            RV += m*(t - t0) + q*(t - t0)**2   # whitening

    return RV


def log_likelihood(theta, t, rv, erv, planet, wh):
    shape_theta = theta.shape
    if len(shape_theta) == 1:  # this case is for emcee. Bayev works with multiple dims
        theta = theta.reshape(1,-1)
        shape_theta = theta.shape
    log_like = np.zeros(shape_theta[0])
    for s in range(shape_theta[0]):
        if planet:
            if wh:
                Vsys, P, K, t0, e, w, jitter, m, q = theta[s]
            else:
                Vsys, P, K, t0, e, w, jitter = theta[s]
                m, q = 0, 0
        else:
            Vsys, jitter = theta[s]
            P, K, t0, e , w, m, q = 0, 0, 0, 0, np.pi/2, 0, 0
        model = model_RV(Vsys, P, K, t0, e, w, m, q, t, planet, wh)
        sigma2 = erv ** 2 + jitter ** 2
        log_like[s] = -0.5 * (len(t)*np.log(2*np.pi) + np.sum((rv - model) ** 2 / sigma2 + np.log(sigma2)))

    return log_like


def log_prior(theta, Priors, prior_type, param_names):
    shape_theta = theta.shape
    if len(shape_theta) == 1:    # this case is for emcee
        theta = theta.reshape(1,-1)
        shape_theta = theta.shape
    log_pr = np.zeros(shape_theta[0])
    for s in range(shape_theta[0]):
        for ind_p,p_name in enumerate(param_names):
            param = theta[s][ind_p]
            if prior_type[p_name] == 'u' and Priors[p_name][0] < param < Priors[p_name][1]:
                log_pr[s] += np.log(1.0/(Priors[p_name][1] - Priors[p_name][0]))
            elif prior_type[p_name] == 'g' and param > 0:
                log_pr[s] += np.log(1.0/(np.sqrt(2.0*np.pi)*Priors[p_name][1])) - 0.5*(param-Priors[p_name][0])**2/Priors[p_name][1]**2
            else:
                log_pr[s] += -np.inf

    return log_pr


def log_probability(theta, t, rv, erv, Priors,  prior_type, param_names, planet, wh):
    lprior = log_prior(theta, Priors, prior_type, param_names)
    if not np.isfinite(lprior):
        return -np.inf, -np.inf

    return log_likelihood(theta, t, rv, erv, planet, wh) + lprior, lprior


#======================
#        MCMC
#======================


def fitMCMC(n_steps, t, rv, erv, Priors, prior_type, planet, wh = False):
    fraction_samples = 0.15
    if planet:
        param_names = ['Vsys', 'P', 'K', 't0', 'e', 'w', 'jitter']
        if wh:
            param_names.append('m')
            param_names.append('q')
    else:
        param_names = ['Vsys', 'jitter']
    ndim  = len(param_names)
    nwalkers = 4*ndim

    for p_name in param_names:
        if prior_type[p_name] == 'u':
            p0_n = [np.random.uniform(Priors[p_name][0], Priors[p_name][1], nwalkers)]
        elif prior_type[p_name] == 'g':
            p0_n = [np.random.normal(loc = Priors[p_name][0], scale = 0.2, size = nwalkers)]
        try:
            p0 = np.concatenate((p0, p0_n), axis = 0)
        except:
            p0 = p0_n
    p0 = p0.transpose()

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args = ([t, rv, erv, Priors,  prior_type, param_names, planet, wh]))
    state = sampler.run_mcmc(p0, nsteps = n_steps, progress = True)
    sampler.reset()
    sampler.run_mcmc(state.coords[np.argmax(state.log_prob)] + 1e-2 * np.random.randn(nwalkers, ndim), nsteps = int(n_steps/2), progress = True)

    flat_samples = sampler.get_chain(flat = True)

    method_args = {'nbins':100, 'nsamples': int(fraction_samples*n_steps/2), 'densityestimation': 'histogram'}
    ev = run_montecarlo(flat_samples, log_likelihood, log_prior, ([t, rv, erv, planet, wh]), ([Priors, prior_type, param_names]), method_args, estimator='perrakis', nmc=50)

    if wh:
        rv_wh = do_whitening(flat_samples, rv, t)
    else:
        rv_wh = False

    return flat_samples, ev, rv_wh, param_names


def do_whitening(flatsamples_wh, rv_prewh, t):
    Vsys, P, K, t0, e, w, jitter, m, q = np.median(flatsamples_wh, axis=0)
    rv_wh = rv_prewh - model_RV(0, P, 0, t0, e, w, m, q, t, planet = True, wh = True)

    return rv_wh


def plot_fitMCMC(t, rv, erv, star, flatsamples, param_names, wh):
    median_parameters_H1 = np.median(flatsamples, axis = 0)
    t_plot = np.linspace(t[0], t[-1], 2*int(t[-1]-t[0]))
    Vsys, P, K, t0, e, w, jitter = median_parameters_H1[:7]
    m, q = 0, 0

    phase = ((np.array(t) - t0) % P) / P
    RV_fit = model_RV(Vsys, P, K, t0, e, w, m, q, t_plot, True, False)

    plt.figure(figsize = (12, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios = [3, 1])
    gs.update(wspace = 0, hspace = 0)
    ax1 = plt.subplot(gs[0])
    ax1.scatter(t, rv, c = 'm', s = 60, label = 'RV data')
    ax1.errorbar(t, rv, yerr = erv, c = 'm', linestyle = "none")
    ax1.plot(t_plot, RV_fit, linewidth = 2, c = 'black', label = '$H_1$')
    ax1.plot(t_plot, np.full(shape = len(t_plot), fill_value = Vsys), linewidth = 2, linestyle = 'dashdot', c = 'black', label = '$H_0$')
    ax1.set_ylabel('RV (m/s)', fontsize = 20)
    ax1.legend(loc = 4, fontsize = 17)
    ax1.set_xticklabels([])

    # plot confidence intervals
    y_models = []
    inds = np.random.randint(len(flatsamples), size = 1000)
    for ind2 in inds:
        Vsys, P, K, t0, e, w, jitter = median_parameters_H1[:7]
        y_model = model_RV(Vsys, P, K, t0, e, w, m, q, t_plot, True, False)
        y_models.append(y_model)
    y_models = np.array(y_models)
    mu,sigma1 = y_models.mean(0), y_models.std(0)
    ax1.fill_between(t_plot, mu - 2 * sigma1, mu + 2 * sigma1, alpha = 0.1, color = 'grey')
    ax1.fill_between(t_plot, mu - sigma1, mu + sigma1, alpha = 0.2, color = 'grey')

    # residuals
    ind = [np.where(np.round(t_plot) == np.round(t)[i])[0][0] for i in range(len(t))]
    residuals = rv - RV_fit[ind]
    ax2 = plt.subplot(gs[1])
    ax2.scatter(t, residuals, s = 15, color = 'black')
    ax2.hlines(0, t[0], t[-1], linestyle = 'dotted', color = 'black')
    ax2.errorbar(t, residuals, yerr = erv, c = 'black', linestyle = "none")
    ax2.set_xlabel(r'Time (JD)', fontsize = 20)
    ax2.set_ylabel('Residuals', fontsize = 20)

    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 15)
    ax1.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 15)
    ax2.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
    plt.savefig(f'Output/Figure/fit_{star}_wh{wh}_Obs#{len(rv)}.pdf', dpi = 300, bbox_inches = 'tight', pad_inches = 0.2)

    fig = corner.corner(flatsamples, labels = param_names)
    plt.savefig(f'Output/Figure/Corner_posteriors_{star}_wh{wh}_Obs#{len(rv)}.pdf', dpi = 300, bbox_inches = 'tight', pad_inches = 0.2)


#======================
#     Find best lBF
#======================


# funtion to reduce the time between observations (i.e. trade-off between number of observations and timespan)
def beta_difdays(lBF_init, lBF, dif_days, beta_param):
    dif_lBF_weig = lBF - lBF_init
    for i in range(len(lBF)):
        dif_lBF_weig[i] = dif_lBF_weig[i] * ss.beta.pdf(x = dif_days[i]/max(dif_days), a = beta_param[0], b = beta_param[1])

    return dif_lBF_weig


def best_lBF(n_steps, median_parameters_H1, lBF_init, sBF_init, schedule_JD, t, rv, erv, Priors, prior_type, min_alt, t_exp, obs, star, Nph, beta, beta_param, max_days_apart):
    Vsys, P, K, t0, e, w, jitter =  median_parameters_H1[:7]
    m, q  = 0, 0

    phase_array, t_cand_array = SimulTime.time_sim(min_alt, t_exp, obs, star, schedule_JD, t0, P, t, Nph, max_days_apart)
    t_before = t[-1]

    rv_cand_array = np.array([])
    erv_cand_array = np.array([])
    ev_H1_array = np.array([])
    sev_H1_array = np.array([])
    ev_H0_array = np.array([])
    sev_H0_array = np.array([])
    number_phase = 0

    for ind, ph in enumerate(phase_array):

        if ph == 0:  # not going too far condition in SimulTime_KOBEsim
            rv_cand_array = np.append(rv_cand_array, np.nan)
            erv_cand_array = np.append(erv_cand_array, np.nan)
            ev_H1_array = np.append(ev_H1_array, np.nan)
            sev_H1_array = np.append(sev_H1_array, np.nan)
            ev_H0_array= np.append(ev_H0_array, np.nan)
            sev_H0_array = np.append(sev_H0_array, np.nan)
            continue

        t_new = t_cand_array[ind]
        number_phase += 1
        rv_new = model_RV(Vsys, P, K, t0, e, w, m, q, t = [t_new], planet = True, wh = False)
        erv_new = np.median(erv)
        rv_cand_array = np.append(rv_cand_array, rv_new)
        erv_cand_array = np.append(erv_cand_array, erv_new)

        t = np.append(t, t_new)
        rv = np.append(rv, rv_new)
        erv =  np.append(erv, erv_new)

        print(f'{star}: KOBEsim testing orbital phase {number_phase}/{len(t_cand_array[t_cand_array!= 0])}')
        flatsamples_H1_new, ev_H1_new, _, _ = fitMCMC(n_steps, t, rv, erv, Priors, prior_type, planet = True, wh = False)
        ev_H1_array = np.append(ev_H1_array, np.median(ev_H1_new))
        sev_H1_array = np.append(sev_H1_array, np.std(ev_H1_new))

        flatsamples_H0_new, ev_H0_new, _, _ = fitMCMC(n_steps, t, rv, erv, Priors, prior_type, planet = False, wh = False)
        ev_H0_array = np.append(ev_H0_array, np.median(ev_H0_new))
        sev_H0_array = np.append(sev_H0_array, np.std(ev_H0_new))

        t = t = np.delete(t, -1)
        rv = np.delete(rv, -1)
        erv = np.delete(erv, -1)

    lBF = ev_H1_array - ev_H0_array
    slBF = np.sqrt(sev_H1_array**2 + sev_H0_array**2)
    dif_lBF_original = lBF - lBF_init
    s_dif_lBF_original = np.sqrt(sBF_init**2 + slBF**2)

    # Consequence of not going too far condition in SimulTime_KOBEsim (remove NaNs)
    ind_nan = np.argwhere(np.isnan(ev_H1_array))
    lBF = np.delete(lBF, ind_nan)
    slBF = np.delete(slBF, ind_nan)
    dif_lBF_original = np.delete(dif_lBF_original,ind_nan)
    s_dif_lBF_original = np.delete(s_dif_lBF_original,ind_nan)
    phase_array = np.delete(phase_array,ind_nan)
    t_cand_array = np.delete(t_cand_array,ind_nan)
    rv_cand_array = np.delete(rv_cand_array, ind_nan)
    erv_cand_array = np.delete(erv_cand_array, ind_nan)

    if beta:
        dif_lBF_weig = beta_difdays(lBF_init, lBF, t_cand_array - t_before, beta_param)
        ind_best = np.where(dif_lBF_weig == np.nanmax(dif_lBF_weig))[0][0]
    else:
        ind_best = np.where(dif_lBF_weig == np.nanmax(dif_lBF_weig))[0][0]

    best_phase = phase_array[ind_best]
    best_t = t_cand_array[ind_best]

    # Save file
    JD_obsnight = [math.floor(jd) for jd in t_cand_array]
    cday_long = Time(JD_obsnight, format='jd', scale='utc').isot
    cday = np.array([t[:10] for t in cday_long])
    ph_reshape =  np.array(phase_array).reshape(-1,phase_array.shape[0])[0]
    df = pd.DataFrame({'Calendar_day':cday[np.argsort(dif_lBF_weig)[::-1]], 'JD':t_cand_array[np.argsort(dif_lBF_weig)[::-1]], 'phase': ph_reshape[np.argsort(dif_lBF_weig)[::-1]],
    'lBF': np.round(lBF[np.argsort(dif_lBF_weig)[::-1]],3), 'sigma_lBF': np.round(slBF[np.argsort(dif_lBF_weig)[::-1]],3),
    'delta_lBF': np.round(dif_lBF_original[np.argsort(dif_lBF_weig)[::-1]],3), 'sigma_delta_lBF': np.round(s_dif_lBF_original[np.argsort(dif_lBF_weig)[::-1]],3)})
    df.to_csv(f'Output/File/{star}_KOBEsim_Obs#{len(rv)+1}.csv', index=False)
    priority = df.index.values

    return phase_array, lBF, slBF, best_phase, best_t, priority


def plot_bestlBF(phase_array, lBF, incert_lBF, best_lBF, dif_lBF, sigma_dif_lBF, calendar_day_best, priority, star, rv, wh):
    ind_best = np.where(lBF == best_lBF)[0][0]
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 6))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 15)
    ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
    plt.axhline(best_lBF, color = 'grey', linewidth = 2, linestyle = 'dashed')
    plt.axhline(best_lBF - dif_lBF[np.where(lBF == best_lBF)[0][0]], color = 'grey' , linewidth = 2, linestyle = 'dashed')
    plt.axvline(phase_array[np.where(lBF == best_lBF)[0][0]], color = 'grey', alpha = 0.6, linewidth = 2, linestyle = 'dashed')

    color_map = 'winter'

    sc = plt.scatter(phase_array, lBF, c = priority[::-1], cmap = color_map, s = 90, edgecolor = 'k')
    cbar = fig.colorbar(sc, ticks = [min(priority) + 0.5, max(priority) - 0.5])
    cbar.ax.set_yticklabels(['Low', 'High'], fontsize = 17)
    cbar.set_label('Priority', rotation = 270, fontsize = 17)
    cNorm = mpl.colors.Normalize(vmin = min(priority), vmax=max(priority))
    mapper = cm.ScalarMappable(norm = cNorm, cmap=color_map)
    colorerr = [(mapper.to_rgba(p)) for p in priority[::-1]]
    for i in range(len(lBF)):
        ax.errorbar(phase_array[i], lBF[i], yerr = incert_lBF[i], c = colorerr[i], linestyle = "none", lw = 2)
    minlBF = min([min(lBF), best_lBF-dif_lBF[np.where(lBF == best_lBF)[0][0]]])
    ymin = minlBF - (0.2*(abs(minlBF-max(lBF))))
    ymax = max(lBF) + (0.2*(abs(max(lBF)-min(lBF))))
    delta_ylim = abs(ymin-ymax)
    plt.ylim(ymin, ymax)

    if phase_array[ind_best] < 0.4:
        xtext = 0.6
        xtext2 = 0.91
    else:
        xtext = 0.042
        xtext2 = 0.01
    if dif_lBF[ind_best] < 0.4*delta_ylim:
        ytext = ymax - 0.26*delta_ylim
        ytext2 = ymax - 0.29*delta_ylim
    else:
        ytext = best_lBF - 0.36*delta_ylim
        ytext2 = best_lBF - 0.39*delta_ylim
    plt.text(phase_array[np.where(lBF == best_lBF)[0][0]] + 0.01, ytext, fr'$\phi$ = {np.round(phase_array[np.where(lBF == best_lBF)[0][0]], 2)}', color = 'gray',
    fontsize = 15, rotation = 90)
    plt.text(phase_array[np.where(lBF == best_lBF)[0][0]] - 0.03, ytext2, f'{calendar_day_best}', color = 'gray', fontsize = 15, rotation = 90)
    xyA = (xtext2 + 0.03, best_lBF)
    xyB = (xtext2 + 0.03, best_lBF - dif_lBF[np.where(lBF == best_lBF)[0][0]])
    coordsA = "data"
    coordsB = "data"
    con1 = ConnectionPatch(xyA, xyB, coordsA, coordsB, arrowstyle="<|-|>", shrinkA = 5, shrinkB = 5, mutation_scale = 20, fc = "gray", color = 'gray')
    ax.add_artist(con1)
    plt.text(xtext2 - 0.01, (2*best_lBF - dif_lBF[ind_best])/2 - 0.22*delta_ylim, f'$\Delta$ $\ln$(B$_{{10}}$) = {round(dif_lBF[ind_best],2)} $\pm$ {round(sigma_dif_lBF[ind_best],2)}', color = 'gray', fontsize = 15, rotation = 90)
    ax.set_xlabel(r'$\phi$', fontsize = 17)
    ax.set_ylabel(f'$\ln$(B$_{{10}}$)', fontsize = 17)
    plt.xlim(-0.05, 1.05)
    plt.savefig(f'Output/Figure/lBF_KOBEsim_{star}_wh{wh}_Obs#{len(rv) + 1}.pdf', dpi = 300, bbox_inches = 'tight', pad_inches = 0.2)
