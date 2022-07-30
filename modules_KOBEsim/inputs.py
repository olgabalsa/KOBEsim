import argparse
import sys

def get():

    parser = argparse.ArgumentParser()
    parser.add_argument("-obs", "--observatory_coordinates", nargs=3, type = float, default = None, metavar = ('lat','long','alt'), help = "List: [Latitude (deg), Longitude (deg), Height (m)]")
    parser.add_argument("-obs_n", "--observatory_name", default = None, help = "Name of the observatory to get the coordinates from EarthLocation (astropy).")
    parser.add_argument("-star", "--star_name", default = None, help = "Name of the star to search for coordinates with simbad")
    parser.add_argument("-file", "--file_path", default = None, help = "Path of your File with the RV Data gathered so far")
    parser.add_argument("-sch", "--schedule_path", default = None, help = "Path of your File with the Schedule")
    parser.add_argument("-P", "--period", type = float, default = None, help = "Period to be targeted with KOBEsim")
    parser.add_argument("-t0", "--time_transit", default = None, type = float, help = "Inferior conjunction JD")
    parser.add_argument("-minalt", "--minimum_altitude", default = 20, help = "Minimum altitude to observe the star (deg)")
    parser.add_argument("-texp", "--t_exp", default = 700, help = "Exposure time (s)")
    parser.add_argument("-Nph", "--N_phases", default = 20, help = "Number of orbital sub-phases")
    parser.add_argument("-beta", "--beta_function_bool", default = True, help = "True to use a beta function to reduce the time gap between observations")
    parser.add_argument("-ab", "--parameters_beta_function", default = [1, 5], help = "Parameters of the Beta(a, b). List: [a, b]")
    parser.add_argument("-wh", "--whitening", default = False, help = "True to substract a linear or quadratic trend on the data before extract T with the periodogram")
    parser.add_argument("-max_da", "--max_days_apart", default = 90, type = float, help = "Maximum days apart searching for the next optimum observing date")
    parser.add_argument("-n", "--n_steps", default = 20000, type = float, help = "Number of steps for the emcee")
    parser.add_argument("-nw", "--mult_n_walkers", default = 4, type = float, help = "Multiple of the number of parameters for the number of walkers for the emcee (Number of walkers = nw * number of parameters)")

    args = parser.parse_args()


    obs_coords = args.observatory_coordinates
    obs_name = args.observatory_name
    if obs_coords != None:
        lat, long, height = args.observatory_coordinates
    elif (obs_coords == None) and (obs_name != None):
        try:
            if (obs_name == 'CAHA') or (obs_name == 'caha'):
                lat, long, height = 37.22, -2.55, 2168
            else:
                coords = EarthLocation.of_site(obs_name)
                lat, long, height = coords.lat.value, coords.lon.value, coords.height.value
        except:
            print('Observatoy name not resolved. Please, check the available observatory names at astropy.coordinates.EarthLocation.get_site_names or provide the observatory coordinates as "-obs lat long height" in units deg, deg, m, respectively.')
            sys.exit()
    else:
        print('Observatoy coordinates or name are mandatory. Provide them when calling the script as "-obs LAT LON HEIGHT" in units deg, deg, m, respectively, or "-obs_n NAME".')
        sys.exit()
    star = args.star_name
    if star == None:
        print(r'Star name is mandatory. Provide it when calling the script as "-star "NAME"".')
        sys.exit()
    path_rv = args.file_path
    if path_rv == None:
        print(r'File with RV timeseries is mandatory. Provide it when calling the script as "-file "path"".')
        sys.exit()

    path_sch = args.schedule_path
    P_peak = args.period
    t0_input = args.time_transit
    min_alt = args.minimum_altitude
    t_exp = args.t_exp
    Nph = int(args.N_phases)
    beta = args.beta_function_bool
    a, b = args.parameters_beta_function
    wh = args.whitening
    max_days_apart = args.max_days_apart
    n_steps = int(args.n_steps)
    mult_nw = int(args.mult_n_walkers)

    def bool_str(val):
        if val == 'False' or val == False:
            return False
        else:
            return True

    beta = bool_str(beta)
    wh = bool_str(wh)

    return [lat, long, height], star, path_rv, path_sch, P_peak, t0_input, min_alt, t_exp, Nph, beta, [a, b], wh, max_days_apart, n_steps, mult_nw
