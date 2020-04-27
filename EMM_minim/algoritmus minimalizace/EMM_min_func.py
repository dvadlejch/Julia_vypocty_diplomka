# 27.4.2020
# file contains functions used in EMM minim process

## packages
import numpy as np
from scipy.optimize import fsolve
from scipy.special import j0, j1


## defining a function calculating the fluorescence signal for given parameters

def get_hist_fit(fotkor, voltages, t_res, t_measurement, background_photocounts):
    # function returns DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma

    # pomocne promene
    fotkor_shape = np.shape(fotkor)
    # casova skala foton-kor. dat
    t_scale = np.array(range(0, fotkor_shape[0])) * t_res

    return 0

