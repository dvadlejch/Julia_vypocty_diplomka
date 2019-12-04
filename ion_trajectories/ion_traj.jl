# funkce vracejici trajektorii iontu v kvadrupolovem priblizeni
using DifferentialEquations
using PhysicalConstants.CODATA2018
using Unitful

function get_ion_traj(Vrf, Udc, Ω, E_ext, delta_phi, u0, tspan; sym_type=false)
    # input: param. pasti, externi DC pole, poc. podm. - [vx, vy, vz, x, y, z], casovy rozsah
    # output: solve output : solveoutput.t - cas, solveoutput.u - [vx, vy, vz, x, y, z]

    # konstanty
    m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
    e = convert(Float64, ElementaryCharge / 1u"C") # naboj iontu
    z0 = 2.25e-3  # vzdalenost axialnich elektrod od stredu pasti [m]
    r0 = 0.6167e-3 # vzdalenost radialnich elektrod od stredu pasti [m]
    κ = 0.0597
    α = 0.7627
    charge_mass_ratio = e/m

    # konstaty v ODE funkci
    c1 = κ*Udc/z0^2
    c2 = Vrf/r0^2
    c3_x = α * Vrf / (2*r0) * delta_phi[1] # faktor urcujici pole zpusobene fazovym rozdilem prot. el.
    c3_y = α * Vrf / (2*r0) * delta_phi[2] # faktor urcujici pole zpusobene fazovym rozdilem prot. el.

    # fukce pro ODE solver asym verze
    function ion_motion_asym!(du, u, p, t)
        # funkce pocita zmeny souradnic u podle pohybove rovnice, je treba 6 rovnic prvniho radu
        # u = [vx, vy, vz, x, y, z]
        c1, c2, c3_y, charge_mass_ratio, E_ext = p # parametry vypoctu

        c_x = (c1 + c2 * cos(Ω*t) ) # konstanta pro vypocet E_x
        c_y = (c1 - c2 * cos(Ω*t) ) # konstanta pro vypocet E_y
        c_z = -2*c1 # konstanta pro vypocet E_z



        du[1] = charge_mass_ratio * (c_x*u[4] + E_ext[1]) # dvx
        du[2] = charge_mass_ratio * (c_y*u[5] + E_ext[2] + c3_y*sin(Ω*t)) # dvy
        du[3] = charge_mass_ratio * (c_z*u[6] + E_ext[3]) # dvz

        du[4] = u[1] # dx
        du[5] = u[2] # dy
        du[6] = u[3] # dz
    end

    # fukce pro ODE solver sym verze
    function ion_motion_sym!(du, u, p, t)
        # funkce pocita zmeny souradnic u podle pohybove rovnice, je treba 6 rovnic prvniho radu
        # u = [vx, vy, vz, x, y, z]
        c1, c2, c3_x, c3_y, charge_mass_ratio, E_ext = p # parametry vypoctu

        c_x = (c1 + c2 * cos(Ω*t) ) # konstanta pro vypocet E_x
        c_y = (c1 - c2 * cos(Ω*t) ) # konstanta pro vypocet E_y
        c_z = -2*c1 # konstanta pro vypocet E_z

        # phase diff

        du[1] = charge_mass_ratio * (c_x*u[4] + E_ext[1] - c3_x/2*sin(Ω*t)) # dvx
        du[2] = charge_mass_ratio * (c_y*u[5] + E_ext[2] + c3_y/2*sin(Ω*t)) # dvy
        du[3] = charge_mass_ratio * (c_z*u[6] + E_ext[3]) # dvz

        du[4] = u[1] # dx
        du[5] = u[2] # dy
        du[6] = u[3] # dz
    end

    if sym_type
        # ode solver
        p = [c1, c2, c3_x, c3_y, charge_mass_ratio, E_ext] # parametry ODE

        prob = ODEProblem(ion_motion_sym!, u0, tspan, p) # definice problemu pro ODE solver
        return solve(prob, reltol=1e-9, abstol=1e-9) # reseni ODEProblem
    else
        # ode solver
        p = [c1, c2, c3_y, charge_mass_ratio, E_ext] # parametry ODE

        prob = ODEProblem(ion_motion_asym!, u0, tspan, p) # definice problemu pro ODE solver
        return solve(prob, reltol=1e-9, abstol=1e-9) # reseni ODEProblem
    end
end
