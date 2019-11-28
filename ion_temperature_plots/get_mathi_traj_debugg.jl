using PhysicalConstants.CODATA2018
using Unitful

function get_mathi_traj_debug()
    m = 40 * convert(Float64,AtomicMassConstant / (1u"kg")) # hmotnost iontu
    e = convert(Float64, ElementaryCharge / 1u"C") # naboj iontu
    z0 = 2.25e-3  # vzdalenost axialnich elektrod od stredu pasti [m]
    r0 = 0.6167e-3 # vzdalenost radialnich elektrod od stredu pasti [m]
    κ = 0.0597
    α = 0.75 # dipol. koef. protejsich elektrod

    # parametry pasti
    Vrf = 500  # napeti radialnich elektrod [V]
    Udc = 1300  # napeti axialnich elektrod [V]
    Ω = 2*pi * 30e6 # budici frekvence pasti [Hz]

    T = 0.5e-3 # teplota iontu
    E_ext = [1,0,0]
    delta_phi_au = [1,1,0] # fazovy rozdil protejsich radialnich elektrod [x, y, 0]
    phi = [pi/2,0,0]

    # casovy rozsah
    #tspan = range(0, 0.995e-5, length=401)
    tspan = range(0, 6.145235e-7, length=601)
    # analyticke reseni

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
    # fazovy rozdil protejsich elektrod
    #delta_phi_au = copy(delta_phi)
    #push!(delta_phi_au, 0)

    # debugg
    div=true
    sym_type=false
    neg_sym_type = !sym_type

    if div
        print("div = true")
        if sym_type
            print("div = true & sym_type = true")
        else
            print("div = true & sym_type = false")
        end

    else
        print("div = false")
        if sym_type
            print("div = false & sym_type = true")
        else
            print("div = false & sym_type = false")
        end
    end


    if div
        # funkce vracejici polohu a rychlost sec. poh, IMM, EMM
        u_sec(t) = u0.*cos.(ω*t + phi) + us # sekularni pohyb
        u_IMM(t) = u0.*cos.(ω*t + phi) .* q/2 .* cos(Ω*t) # IMM - pozor na rozdeleni IMM a EMM pro velka us
        u_EMM(t) = us .* q/2 .* cos(Ω*t) # EMM
        vu_sec(t) = -u0.*ω.*sin.(ω*t + phi)
        vu_IMM(t) = - 1/2*cos(Ω*t)*q.*ω.*u0.*sin.(ω*t+phi) - 1/2*Ω*
            sin(Ω*t)*u0.*q.*cos.(ω*t + phi)
        vu_EMM(t) = - 1/2*Ω*sin(Ω*t)*q.*us

        if sym_type
            u_phase(t) =  1/8*r0*α*sin(Ω*t)*q.*delta_phi_au .* [-1,0,0] +
                1/8*r0*α*sin(Ω*t)*q.*delta_phi_au .* [0,1,0]

            vu_phase(t) = 1/8*r0*α*Ω*cos(Ω*t)*q.*delta_phi_au .* [-1,0,0] +
                1/8*r0*α*Ω*cos(Ω*t)*q.*delta_phi_au .* [0,1,0]

            u_EMM_phase_an = u_phase.(tspan)
            vu_EMM_phase_an = vu_phase.(tspan)
        end

        if neg_sym_type
            a = 1
            print("vlezlo to semka")

            u_phase1(t) = 1/4 * r0 * α * sin(Ω*t)*q.*delta_phi_au .* [0,1,0]

            fun2(t) = 1/4*r0*α*Ω*cos(Ω*t)*q.*delta_phi_au .* [0,1,0]
            u_EMM_phase_an = u_phase1.(tspan)
            vu_EMM_phase_an = fun2.(tspan)


        end

        # reseni
        u_sec_an = u_sec.(tspan)
        u_IMM_an = u_IMM.(tspan)
        u_EMM_an = u_EMM.(tspan)
        vu_sec_an = vu_sec.(tspan)
        vu_IMM_an = vu_IMM.(tspan)
        vu_EMM_an = vu_EMM.(tspan)


        # finalni usporadani matic u_sec, u_IMM, u_EMM
        u_sec_return = zeros(length(tspan), 6)
        u_IMM_return = zeros(length(tspan), 6)
        u_EMM_return = zeros(length(tspan), 6)
        u_EMM_phase_return = zeros(length(tspan), 6)

        #println(1/4 *m* u0[1]^2 * 1/8 * q[1]^2*Ω^2 / ElementaryCharge)
        #println(4/m *( e * q[1]*E_ext[1]/( (2*a[1] + q[1]^2)*Ω) )^2 / ElementaryCharge )
        #println("T = ", 2*pi/ω[1])

        for i in 1:length(tspan)
            for j in 1:3
                u_sec_return[i, j] = vu_sec_an[i][j]
                u_sec_return[i, j+3] = u_sec_an[i][j]
                u_IMM_return[i, j] = vu_IMM_an[i][j]
                u_IMM_return[i, j+3] = u_IMM_an[i][j]
                u_EMM_return[i, j] = vu_EMM_an[i][j]
                u_EMM_return[i, j+3] = u_EMM_an[i][j]
                u_EMM_phase_return[i, j] = vu_EMM_phase_an[i][j]
                u_EMM_phase_return[i, j+3] = u_EMM_phase_an[i][j]
            end
        end


        return (u_sec_return, u_IMM_return, u_EMM_return, u_EMM_phase_return, 2*pi./ω) # plus vraci vektor sekularnich period

    else

        # funkce vracejici polohu a rychlost v case t
        if sym_type

            u(t) = (us + u0.*cos.(ω*t + phi)) .* (1 .+ q/2 .* cos(Ω*t)) +
                1/8*r0*α*sin(Ω*t)*q.*delta_phi_au .* [-1,0,0] + 1/8*r0*α*sin(Ω*t)*
                q.*delta_phi_au .* [0,1,0]

            vu(t) = -u0.*ω.*sin.(ω*t + phi) - 1/2*cos(Ω*t)*q.*ω.*u0.*sin.(ω*t+phi) -
                1/2*Ω*sin(Ω*t)*u0.*q.*cos.(ω*t + phi) - 1/2*Ω*sin(Ω*t)*q.*us + 1/8*
                r0*α*Ω*cos(Ω*t)*q.*delta_phi_au .* [-1,0,0] + 1/8*r0*α*Ω*cos(Ω*t)*
                q.*delta_phi_au .* [0,1,0]

        else
            u(t) = (us + u0.*cos.(ω*t + phi)) .* (1 .+ q/2 .* cos(Ω*t)) +
                1/4 * r0*α*sin(Ω*t)*q.*delta_phi_au .* [0,1,0]

            vu(t) = -u0.*ω.*sin.(ω*t + phi) - 1/2*cos(Ω*t)*q.*ω.*u0.*sin.(ω*t+phi) -
                1/2*Ω*sin(Ω*t)*u0.*q.*cos.(ω*t + phi) - 1/2*Ω*sin(Ω*t)*q.*us +
                1/4*r0*α*Ω*cos(Ω*t)*q.*delta_phi_au .* [0,1,0]
        end

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

        return (u_return, 2*pi./ω) # plus vraci vektor sekularnich period
    end
end
