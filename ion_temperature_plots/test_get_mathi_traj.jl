# test vykresleni trajektorie

using Plots
using PhysicalConstants.CODATA2018
using Unitful
using Statistics
using LaTeXStrings

include("get_mathi_traj.jl")

# parametry pasti
Vrf = 500  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

T = 0.5e-3 # teplota iontu
E_ext = [1,0,0]
delta_phi = [0,0,0] # fazovy rozdil protejsich radialnich elektrod [x, y, 0]
phi = [pi/2,0,0]

# casovy rozsah
#tspan = range(0, 0.995e-5, length=401)
tspan = range(0, 6.145235e-7, length=601)
# analyticke reseni
(u_sec, u_IMM, u_EMM,u_EMM_phase, Per_sec) = get_mathi_traj(Vrf, Udc, Ω, T, E_ext, delta_phi, phi, tspan, div=true, sym_type=false)

u = u_sec + u_IMM + u_EMM # celkovy pohyb iontu

(E_kin_sec, T_kin_sec) = get_E_kin_1D(u_sec[:,1])
(E_kin_IMM, T_kin_IMM) = get_E_kin_1D(u_IMM[:,1])
(E_kin_EMM, T_kin_EMM) = get_E_kin_1D(u_EMM[:,1])



#pyplot()
gr()
plot(tspan/Per_sec[1], T_kin_EMM * 1e3, label="EMM", dpi=200)
plot!(tspan/Per_sec[1], T_kin_IMM * 1e3, label="IMM")
xlabel!(L"$ t/T_{\rm{sec}} $")
ylabel!(L" \rm{Teplota \,\,[mK]}")
#plot!(tspan, T_kin_EMM)
hline!( [mean(T_kin_EMM)*1e3] , linestyle=:dash, label=L"$ \left< T_{\rm{EMM}} \right> $")
hline!([ mean(T_kin_IMM)*1e3], linestyle=:dash, label=L"$ \left< T_{\rm{IMM}} \right> $")

#savefig("D:\\Onedrive_vut\\ÚPT\\diplomka\\Julia_vypocty_diplomka\\ion_temperature_plots\\temp_vs_t_sec_IMM.pdf")
#savefig("/home/dan/diplomka_winfiles/Julia_vypocty_diplomka/ion_temperature_plots/temp_vs_t_EMM_IMM.pdf")
