module comparison

include("propwing.jl")
import CCBlade
using PyPlot
close("all")



Vinf = 50.0  # freestream speed
J = 0.85  # prop advance ratio
etaprop_center = [0.3, 0.6]  # spanwise location of prop (center)
cw = 1.0  # 1 if clockwise, -1 if counter-clockwise (when viewed from upstream)
xtowing = 0.2  # distance from prop to wing (for wake contraction) - assuming it is the same for all props
gamma = 2*pi/180  # flight path angle
alpha_wing_relative_to_fuselage = 2*pi/180 # angle of attack of wing relative to fuselage reference line
alpha = gamma + alpha_wing_relative_to_fuselage  # angle of attack
rho = 1.225  # air density  (irrelevant in normalization)

# --------- propeller --------
# geometry
Rhub = 0.035/2
Rtip = 0.236/2

# chord_prop = Rtip*[0.158788, 0.15083, 0.146158, 0.145793, 0.148978, 0.153175, 0.154064, 0.145548, 0.119747]
# twist_prop = pi/180*[54.9242, 49.6058, 44.6923, 40.1836, 36.0798, 32.3809, 29.0869, 26.1977, 23.7134]
# r = linspace(Rhub+0.05*Rtip, 0.95*Rtip, length(chord_prop))
# r = Rtip*[2.31e-01, 4.00e-01, 5.11e-01, 6.15e-01, 7.20e-01, 7.89e-01, 8.60e-01, 9.30e-01, 9.65e-01, 9.99e-01]
# chord_prop = Rtip*[1.56e-01, 1.46e-01, 1.47e-01, 1.52e-01, 1.54e-01, 1.52e-01, 1.38e-01, 1.10e-01, 9.13e-02, 6.87e-02]
# twist_prop = pi/180*[5.27e+01, 4.39e+01, 3.85e+01, 3.35e+01, 2.91e+01, 2.69e+01, 2.47e+01, 2.35e+01, 2.25e+01, 2.19e+01]
r = 1e-3*[30, 45, 60, 75, 90, 105]
chord_prop = 1e-3*[11.88, 15.59, 18.81, 19.55, 18.32, 3.96]
twist_prop = pi/180*[32.5, 26.5, 23.5, 19, 16.5, 14]

twist_75 = propwing.interp1(r, twist_prop, [0.75*Rtip])[1]
pitch = 25*pi/180 - twist_75
twist_prop += pitch + 6.6*pi/180

B = 4  # number of blades

aftype = CCBlade.af_from_aerodynfile("NACA64_A17.dat")
n = length(r)
af = Array{CCBlade.AirfoilData}(n)
for i = 1:n
    af[i] = aftype
end

Dprop = Rtip*2  # propeller diameter
Omega = Vinf/(J*Dprop)*2*pi  # rotation rate

uwake, vwake, T, Q, reff = propwing.propanalysis(Rhub, Rtip, r, chord_prop, twist_prop, B, af, Vinf, Omega, rho, xtowing)

nprop = length(etaprop_center)
println("T = ", T*nprop)
n = Omega/(2*pi)
CT = T / (rho * n^2 * Dprop^4)



# --------- wing geometry -------
AR = 5.33
sectionspan = [0.64]
twist_wing = [0, 0]*pi/180
b = sectionspan[1]*2
chord_wing = [b/AR, b/AR]
tc = [0.12, 0.12]  # irrelevant for the aerodynamics
sweep = [0]*pi/180
dihedral = [0]
N = [100]


alpha = 13.5*pi/180  # near stall
L, Di, CP, cl, cllocal, Vinfeff, alphaeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, reff, uwake, vwake, etaprop_center, cw)

semispan = sum(sectionspan)
figure()
plot(CP.y/semispan, cllocal)
plot([0, 1], [1.2, 1.2], "--")

figure()
plot(CP.y/semispan, alphaeff*180/pi)

q = 0.5*rho*Vinf^2
S = 2*sum(CP.chord.*CP.ds)
CLmax = L/(q*S)
println("Vinf = ", Vinf)
println("L = ", L)
println("Di = ", Di)


# -----  comparison case with 4 props --------
etaprop_center = [0.2, 0.4, 0.6, 0.8]  # spanwise location of prop (center)

# iterate on Vinf until I match lift.
Vinf = 48.7

Dprop = sqrt(2*Dprop^2*50^2/4/Vinf^2)  # I know what the diameter needs to be for constant thrust because J and solidity are constant.
# resize prop for new diameter
r *= (Dprop/2)/Rtip
Rhub *= (Dprop/2)/Rtip
chord_prop *= (Dprop/2)/Rtip
Rtip = Dprop/2
# compute new rotation rate (for same J)
Omega = Vinf/(J*Dprop)*2*pi



uwake, vwake, T, Q, reff = propwing.propanalysis(Rhub, Rtip, r, chord_prop, twist_prop, B, af, Vinf, Omega, rho, xtowing)

nprop = length(etaprop_center)
println("T = ", T*nprop)

alpha = 13.5*pi/180
L, Di, CP, cl, cllocal, Vinfeff, alphaeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, reff, uwake, vwake, etaprop_center, cw)

semispan = sum(sectionspan)
figure()
plot(CP.y/semispan, cllocal)
plot([0, 1], [1.2, 1.2], "--")

q = 0.5*rho*Vinf^2
S = 2*sum(CP.chord.*CP.ds)
CLmax = L/(q*S)
println("Vinf = ", Vinf)
println("L = ", L)
println("Di = ", Di)

end
