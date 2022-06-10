import modules_KOBEsim.inputs as inputs
import modules_KOBEsim.extract_data as ext_data
import modules_KOBEsim.run_MCMC as run_MCMC
import numpy as np
from astropy.time import Time
import math


"""
KOBEsim finds the next optimum date to observe a given target to maximize the efficiency of the RV observations based on Bayesian evidence (O. Balsalobre-Ruza et al., 2022).

INPUT Parameters
================
- obs*: observatory coordinates and altitude -> latitude(deg) longitude(deg) height(m)
- star*: name of the star to search for coordinates with simbad
- file*: if ASCII -> columns: jd, rv, erv; If FITS -> columns: data['OBJ_DATE_BJD'], data['SPECTRO_CCF_RV'], data['SPECTRO_CCF_RV_ERR']. Units RV and ERV in m/s.
- sch: ASCII file with the schedule of the program (unique columns with JD). If not given, every date can be a candidate.
- P: period of the hypothetical planet in days. If not given, maximum power from the periodogram.
- t0: JD of the conjunction.
- minalt: minimum altitude to observe the star. Default: 20 deg.
- texp: exposure time. Default: 700 s.
- Nph: number of orbital sub-phases to compare. Default: 20.
- beta: prioritize observing in closest dates by means of weighting the increase in the Bayes factor with a beta function. True as default.
- ab: Parameters of the beta function. Default: Beta(a = 1, b = 5).
- wh: whitening. True to substract a linear and a quadratic trend on the data. Default: False.
- max_da: maximum days apart searching for the next optimum observing date. Default: 90 d.
- n: Number of steps per walker for the emcee warm-up phase. Default: 20 000

*: Mandatory inputs.

Returns
=======
- CSV file with the candidate dates to be the next observation sorted by priority.
- Plot ln(lB_10) vs. orbital phase.

Example running the script from terminal
========================================
python run_KOBEsim.py -obs 37.22 -2.55 2168 -star "HD 14737" -file "./RV_data_HD14737.csv" -sch "./schedule_granted_dates.csv" -P 86.5

NOTE: hereafter lBF is ln(B_10) = ln(evidence_planet_hypothesys) - ln(evidence_null_hypothesis)
"""


#            Inputs
#---------------------------------
obs, star, path_rv, path_sch, P_peak, t0_input, min_alt, t_exp, Nph, beta, beta_param, wh, max_days_apart, n_steps = inputs.get()
schedule_JD = ext_data.extract_schedule(path_sch)
jd, rv, erv = ext_data.extract_rv(path_rv)


#            Priors
#---------------------------------
# u -> uniform prior -> [from, to]
# g -> gaussian prior -> [mean, sd]
def Prior_def(t, t0_input, P_peak):
    if t0_input != None:
        t0_prior_type = 'g'
        t0_prior_value = [t0_input, 2]
    else:
        t0_prior_type = 'u'
        t0_prior_value = [t[0], t[0] + P_peak]

    prior_type = {'Vsys': 'u', 'P': 'g', 'K': 'u', 't0': t0_prior_type, 'e': 'u', 'w': 'u', 'jitter': 'u', 'm': 'u', 'q': 'u'}
    Priors = {'Vsys':[-1e9, 1e9], 'P':[P_peak, 2], 'K':[0, 1e4], 't0': t0_prior_value, 'e':[0, 1], 'w':[0, np.pi], 'jitter':[0, 100], 'm':[-100, 100], 'q':[-100, 100]}
    return Priors, prior_type


#==================================
#             MAIN
#==================================


if P_peak == None:
    Priors, prior_type = Prior_def(jd, t0_input, P_peak = 250) # temporary priors to run the prewhitening
    P_peak, rv = ext_data.periodogram_Ppeak(n_steps, jd, rv, erv, Priors, prior_type, wh)
    P_peak_from_periodogram = True
else:
    P_peak_from_periodogram = False

Priors, prior_type = Prior_def(jd, t0_input, P_peak)



# Initial lBF
print(f'Fit H0 and H1 models for {star} data')
flatsamples_H1, ev_H1_init, rv_wh, paramnamesH1 = run_MCMC.fitMCMC(n_steps, jd, rv, erv, Priors, prior_type, planet = True, wh = wh)
if wh and not P_peak_from_periodogram:   # If P_peak_from_periodogram, prewhitening already done
    rv = rv_wh
    flatsamples_H1, ev_H1_init, _, paramnamesH1 = run_MCMC.fitMCMC(n_steps, jd, rv, erv, Priors, prior_type, planet = True)

flatsamples_H0, ev_H0_init, _, _ = run_MCMC.fitMCMC(n_steps, jd, rv, erv, Priors, prior_type, planet = False)
median_parameters_H1 = np.median(flatsamples_H1, axis = 0)

# Fit and Corner plots
run_MCMC.plot_fitMCMC(jd, rv, erv, star, flatsamples_H1, paramnamesH1, wh)

lBF_init = np.median(ev_H1_init) - np.median(ev_H0_init)
slBF_init = np.sqrt(np.std(ev_H1_init)**2 + np.std(ev_H0_init)**2)
print(f'Initial ln(B01) = {round(lBF_init,3)} +- {round(slBF_init,3)}')


# Select best next phase to observe
n_phase_cand, lBF_cand, slBF_cand, best_phase, best_t, priority =  run_MCMC.best_lBF(n_steps, median_parameters_H1, lBF_init, slBF_init, schedule_JD, jd, rv, erv, Priors, prior_type, min_alt, t_exp, obs, star, Nph, beta, beta_param, max_days_apart)
ind_best = np.where(n_phase_cand == best_phase)[0][0]
cday = Time(math.floor(best_t), format = 'jd', scale = 'utc').isot
run_MCMC.plot_bestlBF(n_phase_cand, lBF_cand, slBF_cand, lBF_cand[ind_best], lBF_cand - lBF_init, np.sqrt(slBF_cand**2 + slBF_init**2), cday, priority, star, rv, wh)

print(f'TARGET {star}')
print('-------------------------')
print(f'Optimal phase = {round(best_phase,3)}')
print(f'Optimal next observing date around {round(best_t,3)} JD -> {cday[:10]}')
print(f'Predicted Deltaln(B01) = {round(lBF_cand[ind_best] - lBF_init,3)} +- {round(np.sqrt(slBF_cand[ind_best]**2 + slBF_init**2),3)}')
print(f'Predicted ln(B01) = {round(lBF_cand[ind_best],3)} +- {round(slBF_cand[ind_best],3)}')
