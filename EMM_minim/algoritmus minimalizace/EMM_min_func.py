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

#     # ----- odecet pozadi
#     bg_ph_sum = background_photocounts * t_measure  # celkovy pocet fotonu pozadi za cas mereni
#     last_bin_ratio = fotkor[fotkor_shape[0] - 2, :] / fotkor[fotkor_shape[0] - 3, :]  # pomer mezi county v poslednim/predposlednim binu
#     bg_ph_per_bin = bg_ph_sum / (fotkor_shape[0] - 2 + last_bin_ratio)

#     fotkor[:fotkor_shape[0] - 2, :] = fotkor[:fotkor_shape[0] - 2, :] - bg_ph_per_bin
#     fotkor[fotkor_shape[0] - 2, :] = fotkor[fotkor_shape[0] - 2, :] - bg_ph_per_bin * last_bin_ratio
#     #-----------------------------

    #----- pro vice ruznych casu mereni a pozadi
    # ----- odecet pozadi
    for i in range(fotkor_shape[1]):
        bg_ph_sum = background_photocounts[i] * t_measure[i]  # celkovy pocet fotonu pozadi za cas mereni
        last_bin_ratio = fotkor[fotkor_shape[0] - 2, i] / fotkor[fotkor_shape[0] - 3, i]  # pomer mezi county v poslednim/predposlednim binu
        bg_ph_per_bin = bg_ph_sum / (fotkor_shape[0] - 2 + last_bin_ratio)

        fotkor[:fotkor_shape[0] - 2, i] = fotkor[:fotkor_shape[0] - 2, i] - bg_ph_per_bin
        fotkor[fotkor_shape[0] - 2, i] = fotkor[fotkor_shape[0] - 2, i] - bg_ph_per_bin * last_bin_ratio
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

    return( DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma, x, Omega, Omega_sigma, nu, t_scale )


###### fitovaci funkce uzpusobena mereni v ose xy
def get_hist_Sxy_fit(fotkor, voltages, t_res, t_measure, iont_photocounts, max_phi_unc = 0.3, phi0 = 1.1,
                 sign_DeltaS = False):
    # function returns DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma, hist_sigma

    nu = (voltages[:,0] - voltages[:,1]) / (voltages[:,0] + voltages[:,1])
    # pomocne promene
    fotkor_shape = np.shape(fotkor)
    # casova skala foton-kor. dat
    t_scale = np.array(range(0, fotkor_shape[0])) * t_res

#     # ----- odecet pozadi
#     bg_ph_sum = background_photocounts * t_measure  # celkovy pocet fotonu pozadi za cas mereni
#     last_bin_ratio = fotkor[fotkor_shape[0] - 2, :] / fotkor[fotkor_shape[0] - 3, :]  # pomer mezi county v poslednim/predposlednim binu
#     bg_ph_per_bin = bg_ph_sum / (fotkor_shape[0] - 2 + last_bin_ratio)

#     fotkor[:fotkor_shape[0] - 2, :] = fotkor[:fotkor_shape[0] - 2, :] - bg_ph_per_bin
#     fotkor[fotkor_shape[0] - 2, :] = fotkor[fotkor_shape[0] - 2, :] - bg_ph_per_bin * last_bin_ratio
#     #-----------------------------

    #----- pro vice ruznych casu mereni a pozadi
    # ----- odecet pozadi
    all_photons_sum = np.sum(fotkor, axis=0) # celkovy pocet fotonu = background + odraz Sxy + iont
    # chci odecist background + odraz = celkovy - iont
    for i in range(fotkor_shape[1]):
        bg_ph_sum = all_photons_sum[i] - iont_photocounts[i] * t_measure[i]  # celkovy pocet fotonu pozadi za cas mereni
        last_bin_ratio = fotkor[fotkor_shape[0] - 2, i] / fotkor[fotkor_shape[0] - 3, i]  # pomer mezi county v poslednim/predposlednim binu
        bg_ph_per_bin = bg_ph_sum / (fotkor_shape[0] - 2 + last_bin_ratio)

        fotkor[:fotkor_shape[0] - 2, i] = fotkor[:fotkor_shape[0] - 2, i] - bg_ph_per_bin
        fotkor[fotkor_shape[0] - 2, i] = fotkor[fotkor_shape[0] - 2, i] - bg_ph_per_bin * last_bin_ratio
    #-----------------------------
    
    #------ hist sigma podle poctu fotonu od iontu:
    photon_sum = np.sum(fotkor, axis=0)
    hist_sigma = 0.09088658 * np.sqrt( photon_sum )
    print("hist sigma =",hist_sigma)
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

    return( DeltaS_S_ratio, Delta_S_S_ratio_sigma, fot_phi, fot_phi_sigma, x, Omega, Omega_sigma, nu, t_scale )
#################



## fce vracejici hodnoty pro vykresleni fitu histogramu
def get_hist_fit_values(t_scale, x, Omega):
    # input: casova osa histogramu, x=[S_0, DeltaS_S, fot_phi], Omega
    def fit_func(x, Omega, time_points):
        return x[0]*( 1 + x[1] * np.cos(Omega * time_points + x[2]) )

    time_fit = np.linspace(0, t_scale.max(), 200)
    fotkor_fit = fit_func(x, Omega, time_fit)
    return(time_fit, fotkor_fit)

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

# upravena verze hledajici minimum v ose xy -- je zde nenulovy prumet do osy z, ktery neni zanedbatelny
def get_DeltaS_S_xy_z_proj_fit(DeltaS_S_ratio, fot_phi, U_komp_y, DeltaS_S_min_z,
                        DeltaS_S_min_xz, fot_phi_min_z, fot_phi_min_xz, gamma, Sxy_vec, iter_coef=0.25):
    # input: modulace, napeti na kompenzacni el., iteracni koef.
    # Sxy_vec odpovida jednotkovemu vektoru ve smeru laseru Sxy
    # pouze pro dva body
    # funkce vraci napeti pro dalsi iteraci, fit
    
    
    #### overeno
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
    DeltaS_S_ratio_xy_teor_complex = Sxy_vec[0]/np.sin(gamma) * DeltaS_S_min_xz * np.exp(1j * fot_phi_min_xz) - (Sxy_vec[0]/np.tan(gamma) - Sxy_vec[2]) * DeltaS_S_min_z * np.exp(1j*fot_phi_min_z)
    
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

# def get_beta_xyz_phi_xyz(variables):
#     # input: beta = [beta_z, beta_xz, beta_xy]  variables = [beta, fot_phi, gamma, epsilon]
#     #       fot_phi = [phi_z, phi_xz, phi_xy]

#     # output: beta_xyz = [beta_x, beta_y, beta_z]
    
#     beta = variables[:3]
#     fot_phi = variables[3:6]
#     gamma = variables[6]
#     epsilon = variables[7]

#     beta_x_comp = 1 / np.sin(gamma) * (
#                 beta[1] * np.exp(1j * fot_phi[1]) - beta[0] * np.cos(gamma) * np.exp(1j * fot_phi[0]))
#     beta_y_comp = 1 / np.cos(epsilon) * (np.sin(epsilon) * beta_x_comp - beta[2] * np.exp(1j * fot_phi[2]))

#     return (np.array([np.abs(beta_x_comp), np.abs(beta_y_comp)]),
#             np.array([np.angle(beta_x_comp), np.angle(beta_y_comp)])
#             )

# korekce funkce vracejici prumety beta pro svazek xy mirici lib. smerem
def get_beta_xyz_phi_xyz_Sxy_z_proj(variables):
    # input: beta = [beta_z, beta_xz, beta_xy]  variables = [beta, fot_phi, gamma, a, c]
    #       fot_phi = [phi_z, phi_xz, phi_xy]
    # jednotkovy vektor svazku Sxy = [a,b,c]

    # output: beta_xyz = [beta_x, beta_y, beta_z]
    
    #### overeno
    beta = variables[:3]
    fot_phi = variables[3:6]
    gamma = variables[6]
    a = variables[7]
    c = variables[8]
    b = -np.sqrt(1 - c**2 - a**2)  # zde pozor na znamenko, je treba ho zadat rucne
    
    beta_x_comp = 1 / np.sin(gamma) * (
                beta[1] * np.exp(1j * fot_phi[1]) - beta[0] * np.cos(gamma) * np.exp(1j * fot_phi[0]))
    beta_y_comp = 1 / b * ( beta[2] * np.exp(1j * fot_phi[2]) - a*beta_x_comp - c*beta[0]*np.exp(1j * fot_phi[0]) )

    return (np.array([np.abs(beta_x_comp), np.abs(beta_y_comp)]),
            np.array([np.angle(beta_x_comp), np.angle(beta_y_comp)])
            )

### prepocet nu na z funkce:
def get_z_given_nu(nu, delta_z_ax):
    # funkce vraci axialni polohu iontu 
    # input: nu, [delta_z_ax_5, delta_z_ax_6]
    a = 0.000357087248516796 *1e6
    b = 0.0000614272209845667 *1e6
    c = 0.000214572720416111 *1e6
    
    return( a* nu + b* nu**3 + c *nu**5 + sum(delta_z_ax) * 0.5 )

# funkce vracejici amplitudu rf pole v zavislosti na nu a na ostatnich parametrech
def E_rf_asym_amp_nu(Vrf, phi, delta_z_ax, nu, f_interp):
    # Vrf = [Vrf_1, Vrf_3, Vrf_5, Vrf_6] # defaultne mam Vrf24 = 0 => asym drive
    # phi = [phi_1 = 0, phi_2, phi_56]  # phi_1 = 0 -- volba, dale pak phi_5=phi_6
    # delta_z_ax - [vychyleni ax_5, vychyleni ax_6]
    # nu = [] body, ve kterych chci fci vyhodnotit
    # f_interp = (E_field_rad, E_field_ax5, E_field_ax6) 
    
    # out: amplituda E_pole
    
    rad_amp = sum( Vrf[0:2] * np.exp(phi[0:2] * 1j ) )
    ax_5_amp = Vrf[2] * np.exp(phi[2] * 1j )
    ax_6_amp = Vrf[3] * np.exp(phi[2] * 1j )
    
    z = get_z_given_nu(nu, delta_z_ax)
    E_rf_complex = rad_amp * f_interp[0](z) + ax_5_amp * f_interp[1](z - delta_z_ax[0]) + ax_6_amp * f_interp[2](z - delta_z_ax[1])

    # amplituda a faze
    return np.abs(E_rf_complex)

# funkce vracejici fazi rf pole v zavislosti na nu a na ostatnich parametrech
def E_rf_asym_phase_nu(Vrf, phi, delta_z_ax, nu, f_interp):
    # Vrf = [Vrf_1, Vrf_3, Vrf_5, Vrf_6] # defaultne mam Vrf24 = 0 => asym drive
    # phi = [phi_1 = 0, phi_2, phi_56]  # phi_1 = 0 -- volba, dale pak phi_5=phi_6
    # delta_z_ax - [vychyleni ax_5, vychyleni ax_6]
    # nu = [] body, ve kterych chci fci vyhodnotit
    # f_interp = (E_field_rad, E_field_ax5, E_field_ax6) 
    
    # out: amplituda E_pole
    
    rad_amp = sum( Vrf[0:2] * np.exp(phi[0:2] * 1j ) )
    ax_5_amp = Vrf[2] * np.exp(phi[2] * 1j )
    ax_6_amp = Vrf[3] * np.exp(phi[2] * 1j )
    
    z = get_z_given_nu(nu, delta_z_ax)
    E_rf_complex = rad_amp * f_interp[0](z) + ax_5_amp * f_interp[1](z - delta_z_ax[0]) + ax_6_amp * f_interp[2](z - delta_z_ax[1])

    # amplituda a faze
    return np.angle( E_rf_complex )

def get_axial_EMM_fit_minim_alg(nu, E_rf, E_rf_sigma, E_field_rad_jedna, E_field_ax_5, E_field_ax_6):
    from scipy.optimize import least_squares

    #--------------
    def fit_resid_E_amp_weight(x, nu_data, E_amp_data, weight):
        Vrf_5 = 0.84*x[2]
        Vrf_6 = x[2]

        Vrf = [x[0], x[0], Vrf_5, Vrf_6]

        # fixni parametry
        phi_2 = 0
        phi_56 = x[1]

        phi = np.array( [0, phi_2, phi_56] )
        
        delta_z_ax = [-4.34665041e+01, -2.28337142e+01 ]
    
    
        return( np.sqrt(weight)* (E_rf_asym_amp_nu(Vrf, phi, delta_z_ax,nu_data,(E_field_rad_jedna , E_field_ax_5, E_field_ax_6) ) - E_amp_data) )
    #------------------
    x0 = np.array( [200, 0.1, 1] )

    bbounds = ([0, 0, 0],[900, 2*np.pi, 100] )

    fit = least_squares(fit_resid_E_amp_weight, x0, args=(nu, E_rf, 1/E_rf_sigma**2) , ftol=1e-10, xtol=1e-10,
                       bounds=bbounds)
    # nejistoty
    jac = fit.jac
    C = np.linalg.inv( np.transpose(jac) @ np.diag(1/E_rf_sigma**2) @ jac )
    sigmas_params = np.sqrt( np.diagonal(C) )
    
    return(fit.x, sigmas_params)