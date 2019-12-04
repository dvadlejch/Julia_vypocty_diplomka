# energie celkoveho pohybu

using Plots
using PhysicalConstants.CODATA2018
using Unitful
using Statistics
#using LaTeXStrings

include("get_mathi_traj.jl")

# parametry pasti
Vrf = 500  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

T = 0.5e-3 # teplota iontu
E_ext = [1,1,0]
delta_phi = [2e-4,2e-4,0] # fazovy rozdil protejsich radialnich elektrod [x, y, 0]
phi = [0,0,0]

# casovy rozsah
#tspan = range(0, 2*pi/Ω * 3, length=601)  # drive freq.
tspan = range(0, 3*6.145235e-7, length=601)   # sekularni freq.
# analyticke reseni
(u, Per_sec, E_kin_avg_an) = get_mathi_traj(Vrf, Udc, Ω, T, E_ext,
    delta_phi, phi, tspan, div=false, sym_type=true)

(E_kin_x, T_kin_x) = get_E_kin_1D(u[:,1])
(E_kin_y, T_kin_y) = get_E_kin_1D(u[:,2])
(E_kin_z, T_kin_z) = get_E_kin_1D(u[:,3])

println("-------------------------------")
println("E_kin_x = ", mean(E_kin_x))
println("E_kin_y = ", mean(E_kin_y))
println("E_kin_z = ", mean(E_kin_z))
println("-------------------------------")



#plots

gr()
plot(tspan/Per_sec[1], T_kin_z * 1e3)

plot(tspan, u[:,6])
