# funkce pro vypocet analitickeho vyjadreni trajektorii
# ve kvadrupolovem potencialu
# 28.11.2019

using PhysicalConstants.CODATA2018
using Unitful

function get_mathi_traj(Vrf, Udc, Ω, T, E_ext, phi, tspan, div=false)
    # funkce vraci matici rychlosti a poloh iontu [vx,vy,vz,x,y,z] pro dane casy
    # input: E_ext - externi DC pole, phi - poc. faze trajektorii, tspan - casovy vektor

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

    q_x = -2* e*Vrf / (m*r0^2 * Ω^2)
    q_y = -q_x
    q_z = 0
    q = [q_x, q_y, q_z]

    # amplituda sekul. pohybu
    ω = 1/2 * Ω * sqrt.(a + 1/2 * q.^2)
    u0 = sqrt.(2*BoltzmannConstant/(u"J*K^-1")*T./(m*ω.^2)) # amplitudy sekularniho pohybu
    us = e./(m*ω.^2) .* E_ext

    if div
        # funkce vracejici polohu a rychlost sec. poh, IMM, EMM
        u_sec(t) = u0.*cos.(ω*t + phi) + us # sekularni pohyb
        u_IMM(t) = u0.*cos.(ω*t + phi) .* q/2 .* cos(Ω*t) # IMM - pozor na rozdeleni IMM a EMM pro velka us
        u_EMM(t) = us .* q/2 .* cos(Ω*t) # EMM
        vu_sec(t) = -u0.*ω.*sin.(ω*t + phi)
        vu_IMM(t) = - 1/2*cos(Ω*t)*q.*ω.*u0.*sin.(ω*t+phi) - 1/2*Ω*
            sin(Ω*t)*u0.*q.*cos.(ω*t + phi)
        vu_EMM(t) = - 1/2*Ω*sin(Ω*t)*q.*us

        # reseni
        

    else
        # funkce vracejici polohu a rychlost v case t
        u(t) = (us + u0.*cos.(ω*t + phi)) .* (1 .+ q/2 .* cos(Ω*t))
        vu(t) = -u0.*ω.*sin.(ω*t + phi) - 1/2*cos(Ω*t)*q.*ω.*u0.*sin.(ω*t+phi) - 1/2*Ω*
            sin(Ω*t)*u0.*q.*cos.(ω*t + phi) - 1/2*Ω*sin(Ω*t)*q.*us

        # analiticke reseni
        u_an = u.(tspan)
        vu_an = vu.(tspan)

        # finalni usporadani matice u
        u_return = zeros(length(tspan), 6)

        for i in 1:length(tspan)
            for j in 1:3
                u_return[i, j] = vu_an[i][j]
                u_return[i, j+3] = u_an[i][j]
            end
        end

        return u_return
    end

end

function get_E_kin_1D(u::Array{Float64,1})

    # konstanty
    m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
    e = convert(Float64, ElementaryCharge / 1u"C")
    return 1/2 * m * (u.^2) / e
end
