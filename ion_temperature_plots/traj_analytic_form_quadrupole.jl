# teplota ionty grafy
# 27.11.2019

using Plots
using PhysicalConstants.CODATA2018
using Unitful

# parametry pasti
Vrf = 400  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

# konstanty
m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
e = convert(Float64, ElementaryCharge / 1u"C") # naboj iontu
z0 = 2.25e-3  # vzdalenost axialnich elektrod od stredu pasti [m]
r0 = 0.6167e-3 # vzdalenost radialnich elektrod od stredu pasti [m]
κ = 0.0597

# parametry mathieovy rovnice
a_x = -4* e*κ*Udc / (m*z0^2 * Ω^2)
a_y = a_x
a_z = - 2 * a_x
a = [a_x, a_y, a_z]

q_x = 2* e*Vrf / (m*r0^2 * Ω^2)
q_y = -q_x
q_z = 0
q = [q_x, q_y, q_z]


# amplituda sekul. pohybu
ω = 1/2 * Ω * sqrt.(a + 1/2 * q.^2)
T = 5e-4 # teplota iontu
u0 = sqrt.(2*BoltzmannConstant/(u"J*K^-1")*T./(m*ω.^2)) # amplitudy sekularniho pohybu

# externi elektricke pole
E_ext = [0, 0, 0]
us = e./(m*ω.^2) .* E_ext

# pocatecni faze v jednotlivych osach
phi = [0, 0, 0]

# funkce vracejici polohu a rychlost v case t
u(t) = (us + u0.*cos.(ω*t + phi)) .* (1 .+ q/2 .* cos(Ω*t))
vu(t) = -u0.*ω.*sin.(ω*t + phi) - 1/2*cos(Ω*t)*q.*ω.*u0.*sin.(ω*t+phi) - 1/2*Ω*
    sin(Ω*t)*u0.*q.*cos.(ω*t + phi) - 1/2*Ω*sin(Ω*t)*q.*us

### srovnani numer. reseni a anal. reseni
tspan = range(0,1e-6, length=201)
tspan_ode = (0, 1e-6)
u0_ode = [0,0,0, u0[1], u0[2], u0[3]]
# analiticke reseni
u_an = u.(tspan)
vu_an = vu.(tspan)

u_an_res = zeros(length(u_an),3)
# reshape
for i in 1:length(u_an)
    for j in 1:3
        u_an_res[i,j] = u_an[i][j]
    end
end

# ode solver
include("D:\\Onedrive_vut\\ÚPT\\diplomka\\Julia_vypocty_diplomka\\ion_trajectories\\ion_traj.jl")
u_ode = get_ion_traj(Vrf, Udc, Ω, E_ext, u0_ode, tspan_ode)

# plot srovnani
gr()
plot(tspan,u_an_res[:,1])
plot!(u_ode, vars=(0,4))
