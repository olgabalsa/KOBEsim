import os
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
from astroplan import Observer


# Function to simulate the candidate dates to observe at each phase of the RV curve
# Considerations: granted nights if a schedule is given, target's visibility during exposure time, and twilights
def time_sim(min_alt, t_exp, obs, star, schedule_JD, t0, P, t, N_phases, max_days_apart):

    Obs_loc = EarthLocation(lat = obs[0]*u.deg, lon = obs[1]*u.deg, height = obs[2]*u.m)
    Observatory_class = Observer(Obs_loc)
    t_exp = t_exp/(24*3600)       # from s to jd
    tar = Simbad.query_object(star)
    RA, DEC = tar['RA'][0], tar['DEC'][0]
    coord = SkyCoord(f'{RA} {DEC}', unit = (u.hourangle, u.deg))

    # Array of phases to cover all posibilities
    boundaries_phase = np.linspace(0, 1, N_phases + 1)
    phase_array = np.zeros(N_phases)
    t_cand = np.zeros([N_phases])

    t_now = Time(Time.now(), scale = 'utc').jd
    t_start = max(t_now, t[-1]) + 1
    day = int(t_start)

    while ((not np.all(phase_array)) & (day < t_start + max_days_apart)) | (sum(phase_array != 0) < 1): # Not going too far condition: continue searching for candidate dates if 1) we are not max_days_apart days further away and still orbital phases posibilities empty, or 2) just 1 date as possible even if we are more than max_days_apart days away
        if schedule_JD != None:
            if (day < schedule_JD[-1]) and (day not in schedule_JD):  # If the candidate day is before the ending of the schedule, check if it is assigned
                day += 1
                continue

        else:
            tw_ev =  Observatory_class.twilight_evening_astronomical(Time(day, format = 'jd', scale = 'utc')).value  # Beginning of the night in JD
            tw_mo = Observatory_class.twilight_morning_astronomical(Time(day + 1, format = 'jd', scale = 'utc')).value # Night is over

            phase_day = ((tw_ev - t0) % P) / P
            index_phase = np.where(boundaries_phase < phase_day)[0][-1]
            if phase_array[index_phase] == 0:
                t_night = np.arange(tw_ev, tw_mo, step = t_exp)
                alt_night = coord.transform_to(AltAz(location = Obs_loc, obstime = Time(t_night, format = 'jd'))).alt.value
                index_mask_minalt = np.where(alt_night > min_alt)[0]
                if np.any([index_mask_minalt[i+1] - index_mask_minalt[i] for i in range(len(index_mask_minalt)-1)]) == 1: # During night the target has an altitude > min_alt during a time longer than t_exp
                    phase_array[index_phase] = phase_day
                    t_cand[index_phase] = tw_ev
            day += 1

    return phase_array, t_cand
