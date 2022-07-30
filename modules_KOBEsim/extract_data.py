import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import radvel
from astropy.timeseries import LombScargle
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
from astroplan import Observer
import modules_KOBEsim.run_MCMC as run_MCMC
import modules_KOBEsim.SimulTime_KOBEsim as SimulTime
from matplotlib.ticker import AutoMinorLocator


def extract_rv(path_rv):
    try: # FITS
        file = fits.open(path_file)
        jd = file[1].data['OBJ_DATE_BJD'] + 2400000
        rv = file[1].data['SPECTRO_CCF_RV']
        erv =  file[1].data['SPECTRO_CCF_ERV']
    except: # ASCII
        df_data = pd.read_csv(path_rv, names = ['jd','rv','erv'], header = 0, delimiter = " ", index_col = False, comment = '#')
        jd = df_data.jd.values
        rv = df_data.rv.values
        erv = df_data.erv.values
    return jd, rv, erv


def extract_schedule(path_sch):
    if path_sch != None:
        df_sch = pd.read_csv(path_sch, names = ['day_JD'], delimiter = "\t", index_col = False, comment = '#')
        schedule_JD = [Time(d, format = 'jd') for d in df_sch.day_JD.values]
    else:
        schedule_JD = None
    return schedule_JD


def periodogram_Ppeak(n_steps, mult_nw, t, rv, erv, Priors, prior_type, whitening):
    if whitening:
        rv_prewh = rv
        flatsamples_wh, ev, rv = run_MCMC.fitMCMC(n_steps, mult_nw, t, rv, erv, Priors, prior_type, planet = True, wh = whitening)

        # Plot before and after the prewhitening
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.93, 5.5))
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 15)
        ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
        plt.scatter(t, rv_prewh, color = 'crimson', s = 150, edgecolor = 'k', alpha = 0.8, label = 'Before prewhitening')
        plt.errorbar(t, rv_prewh, yerr = erv, color = 'crimson', lw = 2, linestyle = "none", alpha = 0.8)
        plt.scatter(t, rv, color = 'cornflowerblue', s = 150, edgecolor = 'k', alpha = 0.8, label = 'After prewhitening')
        plt.errorbar(t, rv, yerr = erv, color = 'cornflowerblue', lw = 2, linestyle = "none", alpha = 0.8)
        plt.legend(loc = 'upper left', fontsize = 25)
        plt.xlabel('Date (JD)', fontsize = 20)
        plt.ylabel('RV (m/s)', fontsize = 20)
        plt.tick_params(axis='both', labelsize = 16)
        plt.xscale('log')
        plt.show()

    # Periodogram
    frequency, power = LombScargle(t, rv, erv).autopower()
    periodicity = 1/frequency
    P_peak = periodicity[np.argmax(power)]
    FAP01, FAP005, FAP001 = LombScargle(t, rv, erv).false_alarm_level([0.1, 0.05, 0.01])
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 7), sharex = True, sharey = True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major', labelsize = 15)
    ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
    plt.plot(1/frequency, power/FAP01, c = 'k', lw = 2)

    plt.scatter(periodicity[np.argmax(power)], max(power)/FAP01, color = 'crimson', s = 50, label = f'$P_{{peak}}$ = {round(periodicity[np.argmax(power)],1)} d')
    plt.xscale('log')
    try:
        plt.ylim(0, np.max(power)/FAP01 + 0.1)
    except:
        pass
    plt.xlim(min(periodicity), max(periodicity))
    plt.legend(loc = 'upper right', fontsize = 17)
    plt.ylabel('GLS / FAP 0.1', fontsize = 18)
    plt.xlabel(r'P (d)', fontsize = 18)
    plt.tick_params(axis = 'both', labelsize = 16)
    plt.show()

    return P_peak, rv
