# 27.4.2020
# file contains functions used in EMM minim process

## packages
import numpy as np
# from scipy.optimize import fsolve
# from scipy.special import j0, j1
from scipy.optimize import minimize
from scipy.optimize import least_squares


## defining a function calculating the fluorescence signal for given parameters

def get_hist_fit(fotkor, voltages, t_res, t_measure, background_photocounts, hist_sigma, max_phi_unc = 0.3, phi0 = 1.1,
                 sign_DeltaS = False):
    # function returns DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma, hist_sigma

    nu = (voltages[:,0] - voltages[:,1]) / (voltages[:,0] + voltages[:,1])
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
    ############## fitovani

    # cyklus fitujici vsechny foton-korelacni data

    x = np.zeros((3, fotkor_shape[1]))
    DeltaS_S_ratio = []
    sigmas = np.zeros((3, fotkor_shape[1]))
    Delta_S_S_ratio_sigma = []
    fot_phi = []
    fot_phi_sigma = []



    for i in range(fotkor_shape[1]):

        ##############------------- cast kodu maximalizujici likehood
        # -----
        # podminky urcujici prijimuti reseni
        # max_phi_unc = 0.3
        # phi0 = 1.1

        # --- zde budu zkouset postupne ruzne pocatecni body tak, aby minimalizace vybrala globalni minimum
        x0 = [fotkor[:fotkor_shape[0] - 2, i].mean(), 0.5 * (fotkor[:fotkor_shape[0] - 2, i].max()
                                                                     - fotkor[:fotkor_shape[0] - 2,
                                                                       i].min()) / fotkor[:fotkor_shape[0] - 2,
                                                                                   i].mean(), phi0]
        #     x0 = [fotkor[:fotkor_shape[0]-2,i].mean(), 0.5* ( fotkor[:fotkor_shape[0]-2,i].max()
        #         -fotkor[:fotkor_shape[0]-2,i].min() )/fotkor[:fotkor_shape[0]-2,i].mean(), phi0, 100]
        #     fit = minimize(likehood_transform, x0, args=(Omega, fotkor[:fotkor_shape[0]-2,i], t_res, hist_sigma[i] ), tol=1e-10 )
        fit = minimize(likehood_transform, x0, args=(Omega, fotkor[:fotkor_shape[0] - 2, i], t_res, hist_sigma[i]),
                       tol=1e-10,
                       jac=likehood_transform_jac)
        #     fit = minimize(likehood_transform_sigma, x0, args=(Omega, fotkor[:fotkor_shape[0]-2,i], t_res), tol=1e-10)
        #     print(fit)
        #     print('\n')
        if (np.sqrt(fit.hess_inv[2, 2]) > max_phi_unc) or (fit.x[1] < 0):
            x0 = [fotkor[:fotkor_shape[0] - 2, i].mean(), 0.5 * (fotkor[:fotkor_shape[0] - 2, i].max()
                                                                         - fotkor[:fotkor_shape[0] - 2,
                                                                           i].min()) / fotkor[
                                                                                       :fotkor_shape[0] - 2, i].mean(),
                  -phi0]
            #         x0 = [fotkor[:fotkor_shape[0]-2,i].mean(), 0.5* ( fotkor[:fotkor_shape[0]-2,i].max()
            #             -fotkor[:fotkor_shape[0]-2,i].min() )/fotkor[:fotkor_shape[0]-2,i].mean(), -phi0, 100]
            fit = minimize(likehood_transform, x0,
                           args=(Omega, fotkor[:fotkor_shape[0] - 2, i], t_res, hist_sigma[i]), tol=1e-10,
                           jac=likehood_transform_jac)
        #         fit = minimize(likehood_transform_sigma, x0, args=(Omega, fotkor[:fotkor_shape[0]-2,i], t_res), tol=1e-10)
        #         print(fit)
        #         print('\n')

        # ---- odhad nejistot parametru----
        #     C = fit.hess_inv  # variancni-kovariancni matice
        C = fit.hess_inv
        # -------------------------------
        x[:, i] = fit.x
        sigmas[:, i] = np.sqrt(np.diagonal(C))
        DeltaS_S_ratio.append(x[1, i])

        # ---- faze fot-kor signalu
        fot_phi.append(np.angle(DeltaS_S_ratio[i] * np.exp(1j * x[2, i])))

        # ---- sigma delta s ku s
        Delta_S_S_ratio_sigma.append(np.sqrt(C[1, 1]))

        # ---- sigma fot_phi
        fot_phi_sigma.append(np.sqrt(C[2, 2]))

    # S_0 = x[0, :]  # parametry S_0
    # DeltaS = x[1, :] * S_0  # delta S
    DeltaS_S_ratio = np.array(DeltaS_S_ratio)
    # doplneni zamenek, pokud skoci faze fotkor signalu o vice nez pi/2
    if sign_DeltaS:
        signchange_ind = np.argwhere(np.abs( np.angle( np.exp(1j*np.array(fot_phi) ) * np.exp(-1j * fot_phi[0]) ) )> np.pi/2)
        try:
            DeltaS_S_ratio[signchange_ind] = - DeltaS_S_ratio[signchange_ind]
        except:
            pass

    return( DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma, x, Omega, Omega_sigma, nu )

## fce vracejici koeficienty fitu zavislosti deltaS_S_ratio na \nu, prusecik s nulou, interval napeti pro dalsi iteraci
def get_DeltaS_S_nu_fit(DeltaS_S_ratio, nu, U_avg=500, iter_coef=0.25):
    # input: amplitudy modulace, prislusejici nu, hodnota napeti na axialnich el. kolem ktere hledam interval,
#           koef. pro dalsi iteraci napeti

    def MM_resid(x, deltaS_S, nu):
        return (deltaS_S - x[0] - x[1] * nu)

    def MM_line(x, nu):
        return (x[0] + x[1] * nu)

    x0 = [0.06, +0.1]
    fit = least_squares(MM_resid, x0, args=(DeltaS_S_ratio, nu),
                        ftol=1e-10, xtol=1e-10)
    linfit = fit.x

    nu_MM_zero = - linfit[0] / linfit[1]  # expected MM zero

    #--- vypocet napeti, kde budu hledat v pristi iteraci
    # naleznu interval ve kterem lezi minimum a urcim napeti, pro ktere bych mel merit v dalsi iteraci
    min_inverv_leng = np.abs(min(nu) - nu_MM_zero) * iter_coef
    min_interv = np.array([nu_MM_zero - min_inverv_leng, nu_MM_zero + min_inverv_leng])

    # prepocet na napeti
    U_5 = U_avg * (1 + min_interv)
    U_6 = U_avg * (1 - min_interv)

    # plot fitu
    nu_fit = np.linspace( nu_MM_zero - min_inverv_leng/iter_coef, nu_MM_zero + min_inverv_leng/iter_coef, 200)
    DeltaS_S_fit_nu = MM_line(linfit, nu_fit)

    return( U_5, U_6, min_interv, linfit, nu_fit, DeltaS_S_fit_nu )
##
def get_DeltaS_S_xz_fit(DeltaS_S_ratio, fot_phi, U_komp_x, DeltaS_S_min_z, fot_phi_min_z, gamma, iter_coef=0.25):
    # input: modulace, napeti na kompenzacni el., iteracni koef.
    # pouze pro dva body
    # funkce vraci napeti pro dalsi iteraci, fit

    def MM_resid(x, deltaS_S, nu):
        return (deltaS_S - x[0] - x[1] * nu)

    def MM_line(x, nu):
        return (x[0] + x[1] * nu)

    x0 = [0.06, +0.1]


    fit = least_squares(MM_resid, x0, args=(DeltaS_S_ratio, U_komp_x),
                        ftol=1e-10, xtol=1e-10)
    linfit = fit.x

    # ----- hledani bodu s odpovidajici pozadovanou hodnotou modulace  a faze
    # gamma = 45 / 180 * np.pi  # uhel mezi smerem z a svazkem Sxz

    DeltaS_S_ratio_xz_teor = DeltaS_S_min_z * np.cos(gamma)
    fot_phi_xz_teor = fot_phi_min_z

    U_komp_x_mozne_res = np.array(
        [(DeltaS_S_ratio_xz_teor - linfit[0]) / linfit[1], (-DeltaS_S_ratio_xz_teor - linfit[0]) / linfit[1]])

    sign_of_points_res = np.sign(MM_line(linfit, U_komp_x_mozne_res))
    sign_of_data_points = np.sign(MM_line(linfit, U_komp_x))  # zde jsem zjistil, na ktere strane se nachazi hledany bod

    # ---- obema moznym resenim priradim komplexni cisla podle toho, na ktere strane od nuly jsou
    # pak spocitam rozdil mezi timto prirazenym uhlem a pozadovanym uhlem a vyberu z moznych reseni nejlepsi schodu

    phase_dif = np.abs(
        np.angle(sign_of_data_points * sign_of_points_res * np.exp(1j *( np.array(fot_phi) - fot_phi_xz_teor) )) )

    # rozdil fazi mezi moznym resenim a pozadovanou fazi
    U_komp_x_res = U_komp_x_mozne_res[np.argmin(phase_dif)]

    # vyberu interval, pro dalsi iteraci
    U_komp_x_interval = [U_komp_x_res - np.abs(U_komp_x_res - U_komp_x).min() * iter_coef,
                         U_komp_x_res + np.abs(U_komp_x_res - U_komp_x).min() * iter_coef]

    # fit ke vraceni
    U_komp_x_fit = np.linspace(min(U_komp_x_mozne_res) - np.abs(U_komp_x - U_komp_x_res).max(),
                               max(U_komp_x_mozne_res) + np.abs(U_komp_x - U_komp_x_res).max(), 200)
    DeltaS_S_fit = MM_line(linfit, U_komp_x_fit)

    return( U_komp_x_interval, linfit, U_komp_x_fit, DeltaS_S_fit, DeltaS_S_ratio_xz_teor, fot_phi_xz_teor )
##
def get_DeltaS_S_xy_fit(DeltaS_S_ratio, fot_phi, U_komp_y, DeltaS_S_min_z,
                        DeltaS_S_min_xz, fot_phi_min_z, fot_phi_min_xz, gamma, epsilon, iter_coef=0.25):
    # input: modulace, napeti na kompenzacni el., iteracni koef.
    # pouze pro dva body
    # funkce vraci napeti pro dalsi iteraci, fit

    def MM_resid(x, deltaS_S, nu):
        return (deltaS_S - x[0] - x[1] * nu)

    def MM_line(x, nu):
        return (x[0] + x[1] * nu)

    x0 = [0.06, +0.1]


    fit = least_squares(MM_resid, x0, args=(DeltaS_S_ratio, U_komp_y),
                        ftol=1e-10, xtol=1e-10)
    linfit = fit.x

    # ----- hledani bodu s odpovidajici pozadovanou hodnotou modulace  a faze
    # gamma = 45 / 180 * np.pi  # uhel mezi smerem z a svazkem Sxz
    DeltaS_S_ratio_xy_teor_complex = np.sin(epsilon)/np.sin(gamma) * ( DeltaS_S_min_xz*np.exp(1j*fot_phi_min_xz) -
                                                                       DeltaS_S_min_z*np.cos(gamma)*np.exp(1j*fot_phi_min_z))
    DeltaS_S_ratio_xy_teor = np.abs(DeltaS_S_ratio_xy_teor_complex)
    fot_phi_xy_teor = np.angle(DeltaS_S_ratio_xy_teor_complex)

    U_komp_y_mozne_res = np.array(
        [(DeltaS_S_ratio_xy_teor - linfit[0]) / linfit[1], (-DeltaS_S_ratio_xy_teor - linfit[0]) / linfit[1]])

    sign_of_points_res = np.sign(MM_line(linfit, U_komp_y_mozne_res))
    sign_of_data_points = np.sign(MM_line(linfit, U_komp_y))  # zde jsem zjistil, na ktere strane se nachazi hledany bod

    # ---- obema moznym resenim priradim komplexni cisla podle toho, na ktere strane od nuly jsou
    # pak spocitam rozdil mezi timto prirazenym uhlem a pozadovanym uhlem a vyberu z moznych reseni nejlepsi schodu

    phase_dif = np.abs(
        np.angle(sign_of_data_points * sign_of_points_res * np.exp(1j *( np.array(fot_phi) - fot_phi_xy_teor) )) )

    # rozdil fazi mezi moznym resenim a pozadovanou fazi
    U_komp_y_res = U_komp_y_mozne_res[np.argmin(phase_dif)]

    # vyberu interval, pro dalsi iteraci
    U_komp_y_interval = [U_komp_y_res - np.abs(U_komp_y_res - U_komp_y).min() * iter_coef,
                         U_komp_y_res + np.abs(U_komp_y_res - U_komp_y).min() * iter_coef]

    # fit ke vraceni
    U_komp_y_fit = np.linspace(min(U_komp_y_mozne_res) - np.abs(U_komp_y - U_komp_y_res).max(),
                               max(U_komp_y_mozne_res) + np.abs(U_komp_y - U_komp_y_res).max(), 200)
    DeltaS_S_fit = MM_line(linfit, U_komp_y_fit)

    return( U_komp_y_interval, linfit, U_komp_y_fit, DeltaS_S_fit, DeltaS_S_ratio_xy_teor, fot_phi_xy_teor )

## funkce vracejici bety v ortogonalni bazi xyz
# muze byt dosazeno i DeltaS_S
def get_beta_xyz_phi_xyz(variables):
    # input: beta = [beta_z, beta_xz, beta_xy]  variables = [beta, fot_phi, gamma, epsilon]
    #       fot_phi = [phi_z, phi_xz, phi_xy]

    # output: beta_xyz = [beta_x, beta_y, beta_z]
    beta = variables[:3]
    fot_phi = variables[3:6]
    gamma = variables[6]
    epsilon = variables[7]

    beta_x_comp = 1 / np.sin(gamma) * (
                beta[1] * np.exp(1j * fot_phi[1]) - beta[0] * np.cos(gamma) * np.exp(1j * fot_phi[0]))
    beta_y_comp = 1 / np.cos(epsilon) * (np.sin(epsilon) * beta_x_comp - beta[2] * np.exp(1j * fot_phi[2]))

    return (np.array([np.abs(beta_x_comp), np.abs(beta_y_comp)]),
            np.array([np.angle(beta_x_comp), np.angle(beta_y_comp)])
            )