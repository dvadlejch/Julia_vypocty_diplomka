# reseni trajektorie iontu ve kvadrupolovem poli
# 26.11.2019

## importy
using DifferentialEquations
using PhysicalConstants.CODATA2018
using ForwardDiff
using Plots


# konstanty
m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
e = convert(Float64, ElementaryCharge / 1u"C") # naboj iontu
z0 = 2.25e-3  # vzdalenost axialnich elektrod od stredu pasti [m]
r0 = 0.6167e-3 # vzdalenost radialnich elektrod od stredu pasti [m]
κ = 0.0597
charge_mass_ratio = e/m

# parametry pasti
Vrf = 400  # napeti radialnich elektrod [V]
Udc = 1300  # napeti axialnich elektrod [V]
Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

# vypocet potencialu
#phi(X) = Vrf/(2*r0^2) * ( X[2]^2 - X[1]^2 ) * cos(Ω*X[4]) + κ*Udc/(2*z0^2) * ( 2*X[3]^2 - X[1]^2 - X[2]^2) # X = [x,y,z,t]
    ## zatim se mi jevi nejjednodussi analyticky vyraz symbolicky derivovat, abych dostal el. pole
#E = -1*[ (-κ*Udc/z0^2 - Vrf/r0^2 * Cos(Ω*t) )*x, (-κ*Udc/z0^2 + Vrf/r0^2 * Cos(Ω*t) )*y, 2*κ*Udc/z0^2 * z]

# fukce pro ODE solver
function ion_motion!(du, u, p, t)
    # funkce pocita zmeny souradnic u podle pohybove rovnice, je treba 6 rovnic prvniho radu
    # u = [vx, vy, vz, x, y, z]
    Vrf, Udc, Ω, z0, r0, κ, m, e, charge_mass_ratio, E_ext = p # parametry vypoctu

    c_x = (κ*Udc/z0^2 + Vrf/r0^2 * cos(Ω*t) ) # konstanta pro vypocet E_x
    c_y = (κ*Udc/z0^2 - Vrf/r0^2 * cos(Ω*t) ) # konstanta pro vypocet E_y
    c_z = -2*κ*Udc/z0^2 # konstanta pro vypocet E_z

    du[1] = charge_mass_ratio * (c_x*u[4] + E_ext[1]) # dvx
    du[2] = charge_mass_ratio * (c_y*u[5] + E_ext[2]) # dvy
    du[3] = charge_mass_ratio * (c_z*u[6] + E_ext[3]) # dvz

    du[4] = u[1] # dx
    du[5] = u[2] # dy
    du[6] = u[3] # dz
end

# pocatecni podminky
u0 = [0, 0, 0, 1e-6,1e-7,0] # v metrech

# externi DC pole
E_ext = [3000, 0, 0]
p = [Vrf, Udc, Ω, z0, r0, κ, m, e, charge_mass_ratio, E_ext] # parametry ODE

tspan = (0.0, 1.0e-6)  # casovy rozsah reseni
prob = ODEProblem(ion_motion!, u0, tspan, p) # definice problemu pro ODE solver

sol = solve(prob, reltol=1e-8, abstol=1e-8) # reseni ODEProblem

# plot
gr()
plot(sol, vars=(4,5))
