# 27.4.2020
# file contains functions used in EMM minim process

## packages
import numpy as np
from scipy.optimize import fsolve
from scipy.special import j0, j1
from scipy.optimize import minimize


## defining a function calculating the fluorescence signal for given parameters

def get_hist_fit(fotkor, voltages, t_res, t_measure, background_photocounts):
    # function returns DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma

    # pomocne promene
    fotkor_shape = np.shape(fotkor)
    # casova skala foton-kor. dat
    t_scale = np.array(range(0, fotkor_shape[0])) * t_res

    # ----- odecet pozadi
    bg_ph_sum = background_photocounts * t_measure  # celkovy pocet fotonu pozadi za cas mereni
    last_bin_ratio = fotkor[fotkor_shape[0] - 2, :] / fotkor[fotkor_shape[0] - 3, :]  # pomer mezi county v poslednim/predposlednim binu
    bg_ph_per_bin = bg_ph_sum / (fotkor_shape[0] - 2 + last_bin_ratio)

    fotkor[:fotkor_shape[0] - 2, :] = fotkor[:fotkor_shape[0] - 2, :] - bg_ph_per_bin
    fotkor[fotkor_shape[0] - 2, :] = fotkor[fotkor_shape[0] - 2, :] - bg_ph_per_bin * last_bin_ratio
    #-----------------------------

    #------- odhad RF frekvence i s nejistotou
    # odhad periody triggeru
    T_trig = (fotkor[fotkor_shape[0] - 2, :] / fotkor[fotkor_shape[0] - 3, :]) * t_res + t_scale[fotkor_shape[0] - 2]

    T_trig_sigma = t_res
    # frekvence buzeni pasti

    drive_freq = 1 / T_trig
    drive_freq_sigma = 1 / T_trig ** 2 * T_trig_sigma
    Omega = 2 * np.pi * np.mean( drive_freq )
    Omega_sigma = 2 * np.pi * np.sqrt( sum( drive_freq_sigma**2) / fotkor_shape[1] )

    ####### definice fce, ktera vraci likehood, pomoci ktereho budu fitovat

    def likehood_transform(x, Omega, S, time_step, sigma):
        # definuju funkci vracejici logaritmus pravdepodobnosti, ze z distrubuce dane sinusovkou, co fituji vyberu pozorovane body
        # predpokladam, ze kazdy bod je normalne rozdelen kolem sinusovky

        len_S = len(S)
        # print(len_S)
        S_fit = x[0] * (1 + x[1] * np.cos(Omega * time_step * np.arange(0, len_S) + x[2]))

        sum_term = ((S - S_fit) / sigma) ** 2
        log_term = np.log(np.ones(len_S) * sigma * np.sqrt(2 * np.pi))

        return (0.5 * np.sum(sum_term) + np.sum(log_term))  # vraci -log( likehood)

    def likehood_transform_jac(x, Omega, S, time_step, sigma):
        # vektor jacob. likehood fce
        len_S = len(S)
        sum_term0 = 2 / sigma ** 2 * (1 + x[1] * np.cos(Omega * time_step * np.arange(0, len_S) + x[2])) * (
                    x[0] * x[1] * np.cos(Omega * time_step * np.arange(0, len_S) + x[2]) + x[0] - S)

        sum_term1 = 2 / sigma ** 2 * x[0] * np.cos(Omega * time_step * np.arange(0, len_S) + x[2]) * (
                    x[0] * x[1] * np.cos(Omega * time_step * np.arange(0, len_S) + x[2]) + x[0] - S)

        sum_term2 = (-2 / sigma ** 2) * x[0] * x[1] * (
                    x[0] * x[1] * np.cos(Omega * time_step * np.arange(0, len_S) + x[2]) + x[0] - S) * np.sin(
            Omega * time_step * np.arange(0, len_S) + x[2])

        return (0.5 * np.array([np.sum(sum_term0), np.sum(sum_term1), np.sum(sum_term2)]))
    ##############


    return 0

