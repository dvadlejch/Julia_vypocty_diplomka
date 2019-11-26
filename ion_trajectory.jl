# reseni trajektorie iontu ve kvadrupolovem poli
# 26.11.2019

## importy
using DifferentialEquations
using PhysicalConstants.CODATA2018
using ForwardDiff


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

# vypocet potencialu
#phi(X) = Vrf/(2*r0^2) * ( X[2]^2 - X[1]^2 ) * cos(Ω*X[4]) + κ*Udc/(2*z0^2) * ( 2*X[3]^2 - X[1]^2 - X[2]^2) # X = [x,y,z,t]
