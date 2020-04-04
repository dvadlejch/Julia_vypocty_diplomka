# 19.9.2019
# file contains useful functions related to photon-correlation method

## packages
import numpy as np
from scipy.optimize import fsolve
from scipy.special import j0, j1


## defining a function calculating the fluorescence signal for given parameters

def get_A(decay_rate, detun):
    # returns driven oscillator amplitude in point given by detun param.
    return 1 / 2 * (decay_rate / 2 - detun * 1j) / (detun ** 2 + (decay_rate / 2) ** 2)


def fl_signal(beta, laser_detun, Omega, decay_rate):
    # function calculates deltaS/S0 photon-correlation signal
    # input: beta, laser detuning, RF drive freq, decay rate
    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    numer = 2 * j0(beta) * j1(beta) * np.absolute(np.conj(A) * A_plus - A * np.conj(A_minus))
    denom = j0(beta) ** 2 * np.absolute(A) ** 2 + j1(beta) ** 2 * (np.absolute(A_plus) ** 2 + np.absolute(A_minus) ** 2)

    return numer / denom


def signal_phase(laser_detun, Omega, decay_rate):
    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    return np.angle(np.conj(A) * A_plus - A * np.conj(A_minus))


def get_beta(Omega, decay_rate, laser_detun, norm_mod_amp):
    # function calculates the corresponding beta for given known parameters
    # input: drive frequency, decay_rate, laser detuning, normalized modulation amplitude given by ph_corr_signal
    # output: float value of beta
    def root_func(beta, laser_detun, Omega, decay_rate, norm_mod_amp):
        return fl_signal(beta, laser_detun, Omega, decay_rate) - norm_mod_amp

    sol = fsolve(root_func, np.array([0]), args=(laser_detun, Omega, decay_rate, norm_mod_amp))
    return float(sol)

# function fl_signal contains approximation of the second order term in the fluorescence
# I'm gonna define similar function without the approximation
def fl_signal_second_order(beta, laser_detun, Omega, decay_rate, time, phi, phi_prime, scale):
    # function calculates photon-correlation signal
    # input: beta, laser detuning, RF drive freq, decay rate, time vector, phase of omega term, phase of 2*omega term
    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    Delta_S = 2 * j0(beta) * j1(beta) * np.absolute(np.conj(A) * A_plus - A * np.conj(A_minus))
    S_0 = j0(beta) ** 2 * np.absolute(A) ** 2 + j1(beta) ** 2 * (np.absolute(A_plus) ** 2 + np.absolute(A_minus) ** 2)
    Delta_S_2 = 2*j1(beta)**2 * np.absolute(A_plus*np.conj(A_minus))

    return scale*(S_0 + Delta_S*np.cos(Omega*time + phi) + Delta_S_2*np.cos(2*Omega*time + phi_prime) )

## defining fluorescence func
def fl_signal_second_order_mean(beta, laser_detun, Omega, decay_rate, time, phi, phi_prime, scale):
    # function calculates photon-correlation signal
    # input: beta, laser detuning, RF drive freq, decay rate, time vector, phase of omega term, phase of 2*omega term
    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    Delta_S = 2 * j0(beta) * j1(beta) * np.absolute(np.conj(A) * A_plus - A * np.conj(A_minus))
    #S_0 = j0(beta) ** 2 * np.absolute(A) ** 2 + j1(beta) ** 2 * (np.absolute(A_plus) ** 2 + np.absolute(A_minus) ** 2)
    Delta_S_2 = 2*j1(beta)**2 * np.absolute(A_plus*np.conj(A_minus))

    return scale*( Delta_S*np.cos(Omega*time + phi) + Delta_S_2*np.cos(2*Omega*time + phi_prime) )

## in order to fit second order term I must tweak get_beta function

def get_beta_second(Omega, decay_rate, laser_detun, norm_mod_amp_second):
    # function calculates the corresponding beta for given known parameters
    # input: drive frequency, decay_rate, laser detuning, normalized modulation amplitude given by ph_corr_signal
    # output: float value of beta


    def root_func(beta, laser_detun, Omega, decay_rate, norm_mod_amp_second):
        A_minus = get_A(decay_rate, laser_detun - Omega)
        A_plus = get_A(decay_rate, laser_detun + Omega)
        A = get_A(decay_rate, laser_detun)

        Delta_S = 2 * j0(beta) * j1(beta) * np.absolute(np.conj(A) * A_plus - A * np.conj(A_minus))
        # S_0 = j0(beta) ** 2 * np.absolute(A) ** 2 + j1(beta) ** 2 * (np.absolute(A_plus) ** 2 + np.absolute(A_minus) ** 2)
        Delta_S_2 = 2 * j1(beta) ** 2 * np.absolute(A_plus * np.conj(A_minus))
        return Delta_S_2/Delta_S - norm_mod_amp_second

    sol = fsolve(root_func, np.array([0]), args=(laser_detun, Omega, decay_rate, norm_mod_amp_second))
    return float(sol)