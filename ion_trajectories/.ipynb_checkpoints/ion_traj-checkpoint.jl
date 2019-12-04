# funkce vracejici trajektorii iontu v kvadrupolovem priblizeni
using DifferentialEquations
using PhysicalConstants.CODATA2018

function get_ion_traj(Vrf, Udc, Ω, E_ext, u0, tspan)
    # input: param. pasti, externi DC pole, poc. podm. - [vx, vy, vz, x, y, z], casovy rozsah
    # output: solve output : solveoutput.t - cas, solveoutput.u - [vx, vy, vz, x, y, z]

    # konstanty
    m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
    e = convert(Float64, ElementaryCharge / 1u"C") # naboj iontu
    z0 = 2.25e-3  # vzdalenost axialnich elektrod od stredu pasti [m]
    r0 = 0.6167e-3 # vzdalenost radialnich elektrod od stredu pasti [m]
    κ = 0.0597
    charge_mass_ratio = e/m

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

    # ode solver
    p = [Vrf, Udc, Ω, z0, r0, κ, m, e, charge_mass_ratio, E_ext] # parametry ODE

    prob = ODEProblem(ion_motion!, u0, tspan, p) # definice problemu pro ODE solver
    return solve(prob, reltol=1e-9, abstol=1e-9) # reseni ODEProblem
end
