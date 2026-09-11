import numpy as np
import pandas as pd
import emcee
import corner
import importlib
import bayev.lib as lib
import bayev.perrakis as perrakis
from bayev.run import run_montecarlo
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
from pathlib import Path
import re


#=====================================
# Model and functions to run the MCMC
#=====================================


#Parameters of the RV model:
#---------------------------
#- Vsys: systemic velocity (m/s)
#- P: orbital period (d)
#- K: RV semi-amplitude (m/s)
#- t0: JD of conjuction (d)
#- e: eccentricity
#- w: argument of the periastron (rad)
#- m: parameter of the linear trend
#- q: parameter of the quadratic trend
#-----------------------------
#We use a parametrization with 1) secosw = np.sqrt(e)*np.cos(w), and 2) sesinw = np.sqrt(e)*np.sin(w) to avoid problems in the sampling since the priors of these new parameters can be gaussian centered in 0.


def model_RV(Vsys, P, K, t0, secosw, sesinw, m, q, t, planet):
    Vsys = np.atleast_1d(Vsys)
    P = np.atleast_1d(P)
    K = np.atleast_1d(K)
    t0 = np.atleast_1d(t0)
    secosw = np.atleast_1d(secosw)
    sesinw = np.atleast_1d(sesinw)
    t = np.asarray(t)

    RV_array = np.zeros((len(Vsys), len(t)))
    t_reference = t0[0]

    for i in range(len(Vsys)):
        RV = Vsys[i] + m * (t - t_reference) + q * (t - t_reference)**2

        if planet:
            params = radvel.model.Parameters(
                num_planets=1,
                basis="per tc secosw sesinw k",
            )
            params["per1"].value = P[i]
            params["k1"].value = K[i]
            params["tc1"].value = t0[i]
            params["secosw1"].value = secosw[i]
            params["sesinw1"].value = sesinw[i]

            keplerian = radvel.model._standard_rv_calc(
                t,
                params,
                radvel.model.Vector(params),
            )
            RV += keplerian

        RV_array[i, :] = RV

    return RV_array


def log_likelihood(theta, t, rv, erv, planet, wh):
    shape_theta = theta.shape
    if len(shape_theta) == 1:  # this case is for emcee. Bayev works with multiple dims
        theta = theta.reshape(1,-1)
        shape_theta = theta.shape
    log_like = np.zeros(shape_theta[0])
    for s in range(shape_theta[0]):
        if planet:
            if wh:
                Vsys, P, K, t0, secosw, sesinw, jitter, m, q = theta[s]
            else:
                Vsys, P, K, t0, secosw, sesinw, jitter = theta[s]
                m, q = 0, 0
        else:
            Vsys, jitter = theta[s]
            P, K, t0, secosw, sesinw, m, q = 0, 0, 0, 0, 0, 0, 0
        model = model_RV([Vsys], [P], [K], [t0], [secosw], [sesinw], m, q, t, planet)
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
            elif prior_type[p_name] == 'gt' and Priors[p_name][2] < param < Priors[p_name][3]:
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


def fitMCMC(n_steps, mult_nw, t, rv, erv, Priors, prior_type, merit_f, planet, wh = False):
    fraction_samples = 0.15
    if planet:
        param_names = ['Vsys', 'P', 'K', 't0', 'secosw', 'sesinw', 'jitter']
        if wh:
            param_names.append('m')
            param_names.append('q')
    else:
        param_names = ['Vsys', 'jitter']
    ndim  = len(param_names)
    nwalkers =  mult_nw * ndim

    for p_name in param_names:
        if prior_type[p_name] == 'u':
            p0_n = [np.random.uniform(Priors[p_name][0], Priors[p_name][1], nwalkers)]
        elif prior_type[p_name] == 'g':
            p0_n = [np.random.normal(loc = Priors[p_name][0], scale = 0.2, size = nwalkers)]
        elif prior_type[p_name] == 'gt':
            param_val = [np.random.normal(loc = Priors[p_name][0], scale = 0.2, size = nwalkers)]
            while np.any(np.abs(param_val) > 1):
                param_val = [np.random.normal(loc = Priors[p_name][0], scale = 0.2, size = nwalkers)]
            p0_n = param_val
        try:
            p0 = np.concatenate((p0, p0_n), axis = 0)
        except:
            p0 = p0_n
    p0 = p0.transpose()

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args = ([t, rv, erv, Priors,  prior_type, param_names, planet, wh]), blobs_dtype = np.dtype('object'))
    state = sampler.run_mcmc(p0, nsteps = n_steps, progress = True)
    sampler.reset()
    sampler.run_mcmc(state.coords[np.argmax(state.log_prob)] + 1e-2 * np.random.randn(nwalkers, ndim), nsteps = int(n_steps/2), progress = True)

    flat_samples = sampler.get_chain(flat = True)

    if merit_f == 'BF':
        method_args = {'nbins':200, 'nsamples': int(fraction_samples * n_steps/2), 'densityestimation': 'histogram'}
        ev = run_montecarlo(flat_samples, log_likelihood, log_prior, ([t, rv, erv, planet, wh]), ([Priors, prior_type, param_names]), method_args, estimator='perrakis', nmc=1000)
        merit_f_val = ev

    elif merit_f == 'eK':
        k_index = param_names.index("K")
        k_samples = flat_samples[:, k_index]
        k_low, k_high = np.percentile(k_samples, [15.865, 84.135])
        sigma_k = 0.5 * (k_high - k_low)
        merit_f_val = sigma_k

    if wh:
        rv_wh = do_whitening(flat_samples, rv, t)
    else:
        rv_wh = False

    return flat_samples, merit_f_val, rv_wh, param_names


def do_whitening(flatsamples_wh, rv_prewh, t):
    _, _, _, _, _, _, _, m, q = np.median(flatsamples_wh, axis=0)
    rv_wh = rv_prewh - model_RV(0, [0], 0, [0], [0], [0], m, q, t, False)

    return rv_wh


def output_stem(star, n_observations, whitening):
    target = re.sub(r"[^a-z0-9]+", "-", star.lower()).strip("-")
    return f"{target}_nobs-{n_observations}_wh-{str(whitening).lower()}"


def plot_fitMCMC(t, rv, erv, star, flatsamples, param_names, wh):
    t = np.asarray(t).reshape(-1)
    rv = np.asarray(rv).reshape(-1)
    erv = np.asarray(erv).reshape(-1)

    if not (len(t) == len(rv) == len(erv)):
        raise ValueError(
            f"Inconsistent RV input sizes: "
            f"len(t)={len(t)}, len(rv)={len(rv)}, len(erv)={len(erv)}."
        )

    t_plot = np.linspace(t[0], t[-1], 2 * int(t[-1] - t[0]))
    Vsys, P, K, t0, secosw, sesinw, jitter = flatsamples[:, :7].T
    m, q = 0, 0

    phase = ((np.array(t) - np.median(t0)) % np.median(P)) / np.median(P)
    RV_fit = model_RV(Vsys, P, K, t0, secosw, sesinw, m, q, t_plot, True)

    plt.figure(figsize = (12, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios = [3, 1])
    gs.update(wspace = 0, hspace = 0)
    ax1 = plt.subplot(gs[0])
    ax1.scatter(t, rv, c = 'm', edgecolor = 'w', lw = 1, s = 60, label = 'RV data')
    ax1.errorbar(t, rv, yerr = erv, c = 'm', linestyle = "none")

    muH1, sigmaH1 = RV_fit.mean(0), RV_fit.std(0)
    muH0, sigmaH0 = Vsys.mean(0), Vsys.std(0)
    ax1.plot(t_plot, muH1, linewidth = 2, c = 'black', label = '$H_1$')
    ax1.plot(t_plot, np.full(shape = len(t_plot), fill_value = muH0), linewidth = 2, linestyle = 'dashdot', c = 'black', label = '$H_0$')
    ax1.set_ylabel('RV (m/s)', fontsize = 20)
    ax1.legend(loc = 4, fontsize = 17)
    ax1.set_xticklabels([])
    ax1.fill_between(t_plot, muH1 - 2 * sigmaH1, muH1 + 2 * sigmaH1, alpha = 0.1, color = 'grey')
    ax1.fill_between(t_plot, muH1 - sigmaH1, muH1 + sigmaH1, alpha = 0.2, color = 'grey')

    # residuals
    ind = [np.where(np.round(t_plot) == np.round(t)[i])[0][0] for i in range(len(t))]
    residuals = rv - muH1[ind]
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

    figure_dir = Path("outputs/figures")
    figure_dir.mkdir(parents = True, exist_ok = True)
    stem = output_stem(star, len(rv), wh)

    plt.savefig(
        figure_dir / f"{stem}_rv-fit.pdf",
        dpi = 300,
        bbox_inches = "tight",
        pad_inches = 0.2)

    fig = corner.corner(flatsamples, labels = param_names)
    plt.savefig(
    figure_dir / f"{stem}_posterior-corner.pdf",
    dpi = 300,
    bbox_inches = "tight",
    pad_inches = 0.2)


#======================
#     Find best lBF
#======================


# funtion to reduce the time between observations (i.e. trade-off between number of observations and timespan)
def beta_difdays(lBF_init, lBF, dif_days, beta_param):
    dif_lBF_weig = lBF - lBF_init
    for i in range(len(lBF)):
        dif_lBF_weig[i] = dif_lBF_weig[i] * ss.beta.pdf(x = dif_days[i]/max(dif_days), a = beta_param[0], b = beta_param[1])

    return dif_lBF_weig


def best_lBF(n_steps, mult_nw, flatsamples_H1, lBF_init, sBF_init, schedule_JD, t, rv, erv, Priors, prior_type, min_alt, t_exp, obs, star, Nph, beta, beta_param, max_days_apart, wh):
    stem = output_stem(star, len(rv), wh)
    Vsys, P, K, t0, secosw, sesinw, jitter = flatsamples_H1[:, :7].T
    m, q  = 0, 0

    phase_array, t_cand_array = SimulTime.time_sim(min_alt, t_exp, obs, star, schedule_JD, np.median(t0), np.median(P), t, Nph, max_days_apart)
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

        rv_new_array = model_RV(Vsys, P, K, t0, secosw, sesinw, m, q, np.array([t_new]), True)
        rv_new = rv_new_array.mean(0)

        erv_mcmc = rv_new_array.std(0)
        erv_gaussian = np.random.normal(loc = np.median(erv), scale = np.std(erv))
        erv_new = np.sqrt(erv_mcmc**2 + erv_gaussian**2)

        rv_cand_array = np.append(rv_cand_array, rv_new)
        erv_cand_array = np.append(erv_cand_array, erv_new)

        t = np.append(t, t_new)
        rv = np.append(rv, rv_new)
        erv =  np.append(erv, erv_new)

        print(f'{star}: KOBEsim testing orbital phase {number_phase}/{len(t_cand_array[t_cand_array!= 0])}')
        _, ev_H1_new, _, _ = fitMCMC(n_steps, mult_nw, t, rv, erv, Priors, prior_type, 'BF', planet = True, wh = False)
        ev_H1_array = np.append(ev_H1_array, np.median(ev_H1_new))
        sev_H1_array = np.append(sev_H1_array, np.std(ev_H1_new))

        _, ev_H0_new, _, _ = fitMCMC(n_steps, mult_nw, t, rv, erv, Priors, prior_type, 'BF', planet = False, wh = False)
        ev_H0_array = np.append(ev_H0_array, np.median(ev_H0_new))
        sev_H0_array = np.append(sev_H0_array, np.std(ev_H0_new))

        t = t = np.delete(t, -1)
        rv = np.delete(rv, -1)
        erv = np.delete(erv, -1)

    lBF = ev_H1_array - ev_H0_array
    slBF = np.sqrt(sev_H1_array ** 2 + sev_H0_array ** 2)
    dif_lBF_original = lBF - lBF_init
    s_dif_lBF_original = np.sqrt(sBF_init ** 2 + slBF ** 2)

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
    else:
        dif_lBF_weig = dif_lBF_original
    ind_best = np.where(dif_lBF_weig == np.nanmax(dif_lBF_weig))[0][0]


    best_phase = phase_array[ind_best]
    best_t = t_cand_array[ind_best]

    # Save file
    sort_idx = np.argsort(dif_lBF_weig)[::-1]
    JD_obsnight = [math.floor(jd) for jd in t_cand_array]
    cday_long = Time(JD_obsnight, format='jd', scale='utc').isot
    cday = np.array([t[:10] for t in cday_long])
    ph_reshape =  np.array(phase_array).reshape(-1,phase_array.shape[0])[0]
    df = pd.DataFrame({
        'Calendar_day': cday[sort_idx],
        'JD': t_cand_array[sort_idx],
        'phase': np.round(ph_reshape[sort_idx], 3),
        'lBF': np.round(lBF[sort_idx], 3),
        'sigma_lBF': np.round(slBF[sort_idx], 3),
        'delta_lBF': np.round(dif_lBF_original[sort_idx], 3),
        'sigma_delta_lBF': np.round(s_dif_lBF_original[sort_idx], 3),
        'dif_lBF_weighted': dif_lBF_weig[sort_idx]
    })

    table_dir = Path("outputs/tables")
    table_dir.mkdir(parents = True, exist_ok = True)
    df.to_csv(table_dir / f"{stem}_next-observation-candidates_lBF.csv", index = False)

    priority = np.empty(len(sort_idx), dtype = int)
    priority[sort_idx] = np.arange(1, len(sort_idx) + 1)[::-1]
    priority = np.asarray(priority)

    return phase_array, lBF, slBF, best_phase, best_t, priority


def plot_bestlBF(phase_array, lBF, incert_lBF, lBF_init, dif_lBF, sigma_dif_lBF, calendar_day_best, priority, star, rv, wh):
    
    # Convert inputs to numpy arrays
    phase_array = np.asarray(phase_array)
    lBF = np.asarray(lBF)
    incert_lBF = np.asarray(incert_lBF)
    dif_lBF = np.asarray(dif_lBF)
    sigma_dif_lBF = np.asarray(sigma_dif_lBF)
    priority = np.asarray(priority)

    # Identify the best observing date
    ind_best = np.where(priority == np.max(priority))[0][0]

    best_phase = phase_array[ind_best]
    best_lBF_value = lBF[ind_best]
    best_dif = dif_lBF[ind_best]
    best_sigma_dif = sigma_dif_lBF[ind_best]

    # Figure
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 6))

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 15)
    ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')

    # Expected log Bayes Factor after the observation
    ax.axhline(best_lBF_value, color = 'grey', linewidth = 2, linestyle = 'dashed')
    # Current best log Bayes Factor
    ax.axhline(lBF_init, color = 'grey', linewidth = 2, linestyle = 'dashed')

    # Best observing phase
    ax.axvline(best_phase, color = 'grey', alpha = 0.6, linewidth = 2, linestyle = 'dashed')

    # Scatter plot
    color_map = 'winter'
    sc = ax.scatter(phase_array, lBF, c = priority, cmap = color_map, s = 90, edgecolor = 'k')

    cbar = fig.colorbar(sc, ax = ax, ticks = [np.min(priority), np.max(priority)])
    cbar.ax.set_yticklabels(['Low', 'High'], fontsize = 17)
    cbar.set_label('Priority', rotation = 270, fontsize = 17)

    #cNorm = mpl.colors.Normalize(vmin = np.min(priority), vmax = np.max(priority))
    #mapper = cm.ScalarMappable(norm = cNorm, cmap = color_map)
    #colorerr = [mapper.to_rgba(p) for p in priority]
    # for i in range(len(lBF)):
    #     ax.errorbar(phase_array[i], lBF[i], yerr = incert_lBF[i], color = colorerr[i], linestyle='none', lw=2)

    lower = np.min(lBF) #- incert_lBF)
    upper = np.max(lBF) #+ incert_lBF)

    min_lBF = min(lower, best_lBF_value, best_lBF_value - best_dif)
    max_lBF = max(upper, best_lBF_value, best_lBF_value - best_dif)
    y_range = max_lBF - min_lBF

    if y_range == 0:
        y_range = 1.0

    ymin = min_lBF - 0.1 * y_range
    ymax = max_lBF + 0.1 * y_range

    delta_ylim = ymax - ymin

    ax.set_ylim(ymin, ymax)

    # annotations
    if best_phase < 0.4:
        xtext2 = 0.91
    else:
        xtext2 = 0.03

    if best_dif < 0.4 * delta_ylim:
        ytext = ymax - 0.26 * delta_ylim
        ytext2 = ymax - 0.29 * delta_ylim
    else:
        ytext = best_lBF_value - 0.36 * delta_ylim
        ytext2 = best_lBF_value - 0.39 * delta_ylim

    # Best phase and date annotations
    ax.text(best_phase + 0.01, ytext, rf'$\phi$ = {np.round(best_phase, 2)}', color = 'gray', fontsize = 15, rotation = 90)
    ax.text(best_phase - 0.03, ytext2, f'{calendar_day_best}', color = 'gray', fontsize = 15, rotation = 90)

    xyA = (xtext2 + 0.03, best_lBF_value)
    xyB = (xtext2 + 0.03, best_lBF_value - best_dif)

    con1 = ConnectionPatch(xyA, xyB, "data", "data", arrowstyle = "<|-|>", shrinkA = 5, shrinkB = 5, mutation_scale = 20, fc = "gray", color = "gray")
    ax.add_artist(con1)

    ax.text(xtext2, 0.5,
            (rf'$\Delta$ $\ln$(B$_{{10}}$) = 'f'{round(best_dif, 2)} 'rf'$\pm$ {round(best_sigma_dif, 2)}'),
            transform = ax.transAxes, va = 'center', rotation = 90,
            color = 'gray', fontsize = 15)

    ax.set_xlabel(r'$\phi$', fontsize = 17)
    ax.set_ylabel(r'$\ln$(B$_{10}$)', fontsize = 17)
    ax.set_xlim(-0.05, 1.05)

    figure_dir = Path("outputs/figures")
    figure_dir.mkdir(parents = True, exist_ok = True)
    stem = output_stem(star, len(rv), wh)

    fig.savefig(figure_dir / f"{stem}_bayes-factor-by-phase.pdf", dpi = 300, bbox_inches = "tight", pad_inches = 0.2)

    plt.close(fig)


#======================
#     Find best eK
#======================

def best_eK(n_steps, mult_nw, flatsamples_H1, eK_init, schedule_JD, t, rv, erv, Priors, prior_type, min_alt, t_exp, obs, star, Nph, beta, beta_param, max_days_apart, wh):
    stem = output_stem(star, len(rv), wh)
    Vsys, P, K, t0, secosw, sesinw, jitter = flatsamples_H1[:, :7].T
    m, q  = 0, 0

    phase_array, t_cand_array = SimulTime.time_sim(min_alt, t_exp, obs, star, schedule_JD, np.median(t0), np.median(P), t, Nph, max_days_apart)
    t_before = t[-1]

    rv_cand_array = np.array([])
    erv_cand_array = np.array([])
    eK_array = np.array([])
    number_phase = 0

    for ind, ph in enumerate(phase_array):

        if ph == 0:  # not going too far condition in SimulTime_KOBEsim
            rv_cand_array = np.append(rv_cand_array, np.nan)
            erv_cand_array = np.append(erv_cand_array, np.nan)
            eK_array = np.append(eK_array, np.nan)
            continue

        t_new = t_cand_array[ind]
        number_phase += 1

        rv_new_array = model_RV(Vsys, P, K, t0, secosw, sesinw, m, q, np.array([t_new]), True)
        rv_new = rv_new_array.mean(0)

        erv_mcmc = rv_new_array.std(0)
        erv_gaussian = np.random.normal(loc = np.median(erv), scale = np.std(erv))
        erv_new = np.sqrt(erv_mcmc**2 + erv_gaussian**2)

        rv_cand_array = np.append(rv_cand_array, rv_new)
        erv_cand_array = np.append(erv_cand_array, erv_new)

        t = np.append(t, t_new)
        rv = np.append(rv, rv_new)
        erv =  np.append(erv, erv_new)

        print(f'{star}: KOBEsim testing orbital phase {number_phase}/{len(t_cand_array[t_cand_array!= 0])}')
        _, eK_new, _, _ = fitMCMC(n_steps, mult_nw, t, rv, erv, Priors, prior_type, 'eK', planet = True, wh = False)
        eK_array = np.append(eK_array, eK_new)
       
        t = t = np.delete(t, -1)
        rv = np.delete(rv, -1)
        erv = np.delete(erv, -1)

    # Consequence of not going too far condition in SimulTime_KOBEsim (remove NaNs)
    ind_nan = np.argwhere(np.isnan(eK_array))
    eK_array = np.delete(eK_array, ind_nan)
    phase_array = np.delete(phase_array,ind_nan)
    t_cand_array = np.delete(t_cand_array,ind_nan)
    rv_cand_array = np.delete(rv_cand_array, ind_nan)
    erv_cand_array = np.delete(erv_cand_array, ind_nan)

    dif_eK_original = eK_init - eK_array
    if beta:
        dif_days = t_cand_array - t_before
        dif_eK = dif_eK_original * ss.beta.pdf(x = dif_days/max(dif_days), a = beta_param[0], b = beta_param[1])
    else:
        dif_eK = dif_eK_original.copy()

    ind_best = np.where(dif_eK == np.nanmax(dif_eK))[0][0]

    best_phase = phase_array[ind_best]
    best_t = t_cand_array[ind_best]

    # Save file
    JD_obsnight = [math.floor(jd) for jd in t_cand_array]
    cday_long = Time(JD_obsnight, format='jd', scale='utc').isot
    cday = np.array([t[:10] for t in cday_long])
    ph_reshape =  np.array(phase_array).reshape(-1,phase_array.shape[0])[0]
    df = pd.DataFrame({'Calendar_day':cday[np.argsort(dif_eK)[::-1]], 
                     'JD':t_cand_array[np.argsort(dif_eK)[::-1]], 
                     'phase': ph_reshape[np.argsort(dif_eK)[::-1]],
                     'eK': np.round(eK_array[np.argsort(dif_eK)[::-1]],3), 
                     'delta_eK': np.round(dif_eK_original[np.argsort(dif_eK)[::-1]],3)})

    table_dir = Path("outputs/tables")
    table_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(
        table_dir / f"{stem}_next-observation-candidates_eK.csv",
        index = False)

    priority = df.index.values

    return phase_array, eK_array, best_phase, best_t, priority