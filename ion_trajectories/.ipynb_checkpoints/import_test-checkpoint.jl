# vykresleni trajektorie iontu
# 27.11.2019

using Plots
include("ion_traj.jl")

# parametry pasti
Vrf = 400  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

# pocatecni podminky
u0 = [0, 0, 0, 1e-6,1e-7,1e-6] # v metrech

# externi DC pole
E_ext = [3000, 0, 0]
tspan = (0.0, 1.0e-6)  # casovy rozsah reseni

traj = get_ion_traj(Vrf, Udc, Ω, E_ext, u0, tspan) # trajektorie iontu

gr()
plot(traj, vars=(4,6))
