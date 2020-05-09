module example

include("propwing.jl")
import CCBlade
using PyPlot
close("all")

# ------- common parameters (standard metric units) ------
rho = 1.225  # air density
Vinf = 10.0  # freestream speed
J = 0.5  # prop advance ratio
xtowing = 18.0  # distance from prop to wing (for wake contraction) - assuming it is the same for all props
etaprop_center = [40.0/94, 70.0/94]  # spanwise location of prop (center)
cw = 1.0  # 1 if clockwise, -1 if counter-clockwise (when viewed from upstream)
gamma = 3*pi/180  # flight path angle
alpha_wing_relative_to_fuselage = 2*pi/180 # angle of attack of wing relative to fuselage reference line

alpha = gamma + alpha_wing_relative_to_fuselage  # angle of attack


# --------- propeller --------
# geometry
Rhub = 3*.5
Rtip = 3*3.0

r = 3*[0.526, 0.628, 0.729, 0.831, 0.9132, 0.9586, 1.0332,
     1.1128, 1.1925, 1.2722, 1.3519, 1.4316, 1.5114, 1.5911,
     1.6708, 1.7505, 1.8302, 1.9099, 1.9896, 2.0693, 2.1490, 2.2287,
     2.3084, 2.3881, 2.4678, 2.5475, 2.6273, 2.7070, 2.7867, 2.8661, 2.9410]
chord_prop = 3*[0.6270, 0.6255, 0.6231, 0.6199, 0.6165, 0.6125, 0.6054, 0.5973, 0.5887,
          0.5794, 0.5695, 0.5590, 0.5479, 0.5362, 0.5240, 0.5111, 0.4977,
          0.4836, 0.4689, 0.4537, 0.4379, 0.4214, 0.4044, 0.3867, 0.3685,
          0.3497, 0.3303, 0.3103, 0.2897, 0.2618, 0.1920]

twist_prop = pi/180.0*[40.2273, 38.7657, 37.3913, 36.0981, 34.8803, 33.5899, 31.6400,
                   29.7730, 28.0952, 26.5833, 25.2155, 23.9736, 22.8421, 21.8075,
                   20.8586, 19.9855, 19.1800, 18.4347, 17.7434, 17.1005, 16.5013,
                   15.9417, 15.4179, 14.9266, 14.4650, 14.0306, 13.6210, 13.2343,
                   12.8685, 12.5233, 12.2138]
B = 3  # number of blades

aftype = CCBlade.af_from_aerodynfile("NACA64_A17.dat")
n = length(r)
af = Array{CCBlade.AirfoilData}(n)
for i = 1:n
    af[i] = aftype
end

Dprop = Rtip*2  # propeller diameter
Omega = Vinf/(J*Dprop)*2*pi  # rotation rate

uwake, vwake, T, Q, reff = propwing.propanalysis(Rhub, Rtip, r, chord_prop, twist_prop, B, af, Vinf, Omega, rho, xtowing)


# --------- wing geometry -------
sectionspan = [29, 65]
chord_wing = [43, 26, 11]
twist_wing = [0, 0, 0]*pi/180
tc = [0.13, 0.12, 0.11]
sweep = [25, 30]*pi/180
dihedral = [0, 0]
N = [29, 65]

# L, Di, CP, cl, cllocal, Vinfeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, reff, uwake, vwake, etaprop_center, cw)
# println(L)
# println(Di)
#
# figure()
# plot(CP.y, cl)
# plot(CP.y, cllocal)
#
# # ------ viscous drag -----
# # TODO
# Dp = 0.0
# D = Di + Dp
#
# # ------- add prop thrust to lift --------
# L += T*sin(gamma)  # add prop thrust to lift
# thrust_ratio = T*cos(gamma)/D  # must be larger than 1


# -------- diameter checking --------
nd = 100
dratio = linspace(0.02, 0.8, nd)
CLvec = zeros(nd)
CDivec = zeros(nd)

etaprop_center = [40.5828/94]  # spanwise location of prop (center)

for i = 1:nd
    rprop = reff/reff[end]*dratio[i]/2*sum(sectionspan)
    L, Di, CP, cl, cllocal, Vinfeff, alphaeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, rprop, uwake, vwake, etaprop_center, cw)
    q = 0.5*rho*Vinf^2
    S = 2*sum(CP.chord.*CP.ds)
    CLvec[i] = L/(q*S)
    CDivec[i] = Di/(q*S)
end

# println(CLvec)
# println(CDivec)

CLbad = [0.382905, 0.386424, 0.386774, 0.386081, 0.389872, 0.390167, 0.392373, 0.39475, 0.394911, 0.397192, 0.397608, 0.397071, 0.400073, 0.400052, 0.399881, 0.402262, 0.40314, 0.405263, 0.406519, 0.406577, 0.409258, 0.410152, 0.409989, 0.41318, 0.414047, 0.415491, 0.417441, 0.417875, 0.420247, 0.421134, 0.42098, 0.42363, 0.424618, 0.425658, 0.428337, 0.429073, 0.430862, 0.432741, 0.433258, 0.436002, 0.437109, 0.437732, 0.44078, 0.441917, 0.443226, 0.445345, 0.445991, 0.448085, 0.449637, 0.450339, 0.453482, 0.454775, 0.456139, 0.458766, 0.459777, 0.461669, 0.463841, 0.464765, 0.467362, 0.46903, 0.470428, 0.473236, 0.474219, 0.475866, 0.478265, 0.479343, 0.481874, 0.483983, 0.485173, 0.488129, 0.489759, 0.491443, 0.494142, 0.495331, 0.497667, 0.500125, 0.501265, 0.504214, 0.506247, 0.507659, 0.51031, 0.511974, 0.514132, 0.516882, 0.51815, 0.520937, 0.523536, 0.525125, 0.52798, 0.530019, 0.53199, 0.534845, 0.536503, 0.539055, 0.541797, 0.543271, 0.546195, 0.548757, 0.550623, 0.553543]
CDibad = [0.00588463, 0.00587591, 0.00587341, 0.00590396, 0.00588026, 0.00587782, 0.00588183, 0.00585698, 0.00585025, 0.0058566, 0.00585454, 0.005875, 0.00586265, 0.00586892, 0.00589654, 0.00588845, 0.00588124, 0.00588864, 0.00588322, 0.00589091, 0.00588995, 0.00588725, 0.00590671, 0.00589429, 0.00588971, 0.00589854, 0.00589212, 0.00589244, 0.00589664, 0.00589763, 0.00591354, 0.00591262, 0.00591137, 0.00592136, 0.00591176, 0.00591071, 0.00591893, 0.00591281, 0.00591794, 0.00591655, 0.00591542, 0.00592514, 0.00591515, 0.0059112, 0.00592028, 0.00591486, 0.00591883, 0.00592354, 0.0059209, 0.00592647, 0.00591972, 0.00591554, 0.00591953, 0.00591197, 0.00591032, 0.00591382, 0.0059059, 0.00590737, 0.00590384, 0.00589778, 0.00589761, 0.00588931, 0.00588846, 0.00589156, 0.00588417, 0.00588266, 0.00587735, 0.00586829, 0.00586781, 0.00585689, 0.00584926, 0.00584727, 0.00583654, 0.00583216, 0.00582569, 0.00581349, 0.00580967, 0.00579741, 0.00578652, 0.00578435, 0.00577362, 0.0057639, 0.00575583, 0.00574088, 0.00573432, 0.00572027, 0.00570297, 0.00569357, 0.00567768, 0.00566346, 0.00565299, 0.0056353, 0.00562337, 0.00560788, 0.00558971, 0.00557922, 0.00555981, 0.00554011, 0.00552628, 0.00550529]



rc("figure", figsize=(3.5, 2.6))
rc("font", size=10.0)  #, family="CMU Serif")
rc("lines", linewidth=1.5)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.16, top=0.97, right=0.95)
rc("axes", color_cycle=["348ABD", "A60628", "009E73"])

figure()
plot(dratio, CLbad)
plot(dratio, CLvec, linewidth=1)
xlabel("diameter/semispan")
ylabel(L"C_L")
xlim([0, 0.3])
ylim([0.38, 0.44])
legend(["original", "smooth"], loc="upper left")
# savefig("../MDAO_paper/figures/smoothCL.pdf", transparent=true)

figure()
plot(dratio, CDibad*10000)
plot(dratio, CDivec*10000, linewidth=1)
xlabel("diameter/semispan")
ylabel(L"$C_{Di}$ (counts)")
xlim([0, 0.8])
legend(["original", "smooth"])
# savefig("../MDAO_paper/figures/smoothCDi.pdf", transparent=true)





end
