# reseni trajektorie iontu ve kvadrupolovem poli - mathieova rovnice
# 27.11.2019

## importy
using DifferentialEquations
using PhysicalConstants.CODATA2018
using Unitful
using Plots


# konstanty
m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
e = convert(Float64, ElementaryCharge / 1u"C") # naboj iontu
z0 = 2.25e-3  # vzdalenost axialnich elektrod od stredu pasti [m]
r0 = 0.6167e-3 # vzdalenost radialnich elektrod od stredu pasti [m]
κ = 0.0597

# parametry pasti
Vrf = 400  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

# parametry mathieovy rovnice
a_x = -4* e*κ*Udc / (m*z0^2 * Ω^2)
a_y = a_x
a_z = - 2 * a_x

q_x = 2* e*Vrf / (m*r0^2 * Ω^2)
q_y = -q_x
q_z = 0

## reseni mathieovy rovnice

# funkce pro ODE solver
function mathi!(du, u, p, ξ)
    # p = [ax, ay, az, qx, qy, qz]
    # u = [vx, vy, vz, x, y, z]
    du[1] = - (p[1] - 2*p[4]*cos(2*ξ)) * u[4] # Dvx
    du[2] = - (p[2] - 2*p[5]*cos(2*ξ)) * u[5]
    du[3] = - (p[3] - 2*p[6]*cos(2*ξ)) * u[6]

    du[4] = u[1] # Dx
    du[5] = u[2]
    du[6] = u[3]
end

# poc podminky
u0 = [0, 0, 0, 1e-6, 1e-6, 1e-6]
p = [a_x, a_y, a_z, q_x, q_y, q_z]

# casovy rozsah
ξspan = (0, 35*pi)

# reseni ODEProblem
prob = ODEProblem(mathi!, u0, ξspan, p) # definice problemu pro ODE solver

sol = solve(prob, abstol=1e-8,reltol=1e-8) # reseni ODEProblem
#println(sol)

# plot
gr(); gui()
plot(sol, vars=(4,6))
