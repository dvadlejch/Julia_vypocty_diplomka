# test vykresleni trajektorie

using Plots
using PhysicalConstants.CODATA2018
using Unitful
using Statistics
include("get_mathi_traj.jl")

# parametry pasti
Vrf = 350  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

T = 5e-3 # teplota iontu
E_ext = [0,0,100]
phi = [0,0,0]

tspan = range(0, 1e-5, length=401)

u = get_mathi_traj(Vrf, Udc, Ω, T, E_ext, phi, tspan)
vz = u[:,3]
E_kin_z = get_E_kin_1D(vz)
T_kin_z = convert(Float64, ElementaryCharge/(1u"C")) * E_kin_z / convert(Float64, BoltzmannConstant/(1u"J*K^-1") )

gr()
plot(tspan, T_kin_z)
hline!([mean(T_kin_z)])
