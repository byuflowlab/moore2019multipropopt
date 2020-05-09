module example

include("propwing.jl")
import CCBlade
using PyPlot
close("all")

Vinf = 50.0  # freestream speed
J = 0.85  # prop advance ratio
etaprop_center = [0.3/0.64]  # spanwise location of prop (center)
cw = 1.0  # 1 if clockwise, -1 if counter-clockwise (when viewed from upstream) (symmetric wing is on the right viewed from upstream, so clockwise is inboard up)
alpha = 4*pi/180
xtowing = 0.2  # distance from prop to wing (for wake contraction) - assuming it is the same for all props

rho = 1.225  # air density  (irrelevant in normalization)



# Pprop = 5.5e3
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
println(pitch*180/pi)
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

#TODO: check thrust coefficient make sure it is 0.168
n = Omega/(2*pi)
CT = T / (rho * n^2 * Dprop^4)
# println(CT)
# println(Q*Omega)
# println(T*Vinf)


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

alpha = 4.0*pi/180
L, Di, CP, cl1, cllocal, Vinfeff, alphaeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, reff, uwake, vwake, etaprop_center, cw)

alpha = 0.0
L, Di, CP, cl2, cllocal, Vinfeff, alphaeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, reff, uwake, vwake, etaprop_center, cw)


expdata = [
9.772047626521212 0.36487438017230045;
17.456214395697152 0.3699581337937061;
25.141017429406453 0.3776505720018155;
28.36038871521223 0.3770728438055406;
31.51982388197718 0.4307570275417708;
34.66386144703521 0.42131104427976984;
37.734219579132215 0.4097793858774725;
40.93806601032445 0.3455497537656257;
52.987516489855814 0.2482967198442424;
56.206378764034916 0.24563204397860466;
59.276864149038595 0.23462212249364806;
62.42090171409665 0.22517613923164698;
65.57779182272824 0.2684255846210621;
68.72221114650627 0.2605448121110834;
74.86521796302029 0.24687275981862233;
81.22608175575076 0.2264143100135313;
87.36476197343806 0.19500320253148462;
92.03952475281122 0.16153059796140856
]


expdata0=[
9.619921363040625 0.029354838709677367;
17.404980340760133 0.032258064516128976;
25.111402359108773 0.04193548387096771;
28.256880733944946 0.04290322580645156;
31.402359108781134 0.09612903225806448;
34.469200524246375 0.08903225806451609;
37.53604193971165 0.08032258064516125;
40.76015727391874 0.04999999999999995;
52.87024901703797 -0.025161290322580687;
56.01572739187416 -0.05774193548387102;
59.161205766710346 -0.07806451612903226;
62.30668414154652 -0.08354838709677423;
65.37352555701176 -0.04548387096774195;
68.59764089121887 -0.021290322580645227;
74.73132372214937 -0.018709677419354864;
81.02228047182172 -0.013548387096774195;
87.07732634338137 -0.010967741935483832;
91.79554390563564 -0.008387096774193525;
]


rc("figure", figsize=(3.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.19, bottom=0.18, top=0.97, right=0.96)


semispan = sectionspan[1]
figure()
plot(CP.y/semispan, cl1)
plot(expdata[:, 1]/100, expdata[:, 2], "o")
# plot(CP.y/semispan, cllocal)

ax = gca()
ax[:set_color_cycle](nothing)

# figure()
plot(CP.y/semispan, cl2)
plot(expdata0[:, 1]/100, expdata0[:, 2], "o")

xlabel(L"2y/b")
ylabel(L"c_l")
legend(["VLM", "Exper."], loc="upper right")
# savefig("../MDAO_paper/figures/epemaval1.pdf", transparent=true)

xtowing = 0.2  # distance from prop to wing (for wake contraction) - assuming it is the same for all props

Vinf = 19.0  # freestream speed
J = 0.695  # prop advance ratio
etaprop_center = [0.35]  # spanwise location of prop (center)
cw = 1.0  # 1 if clockwise, -1 if counter-clockwise (when viewed from upstream)
alpha = 4*pi/180

rho = 1.225  # air density  (irrelevant in normalization)



# Pprop = 5.5e3
# --------- propeller --------
# geometry
Rhub = 0.084/2
Rtip = 0.41/2

r = Rtip*[2.31e-01, 4.00e-01, 5.11e-01, 6.15e-01, 7.20e-01, 7.89e-01, 8.60e-01, 9.30e-01, 9.65e-01, 9.99e-01]
chord_prop = Rtip*[1.56e-01, 1.46e-01, 1.47e-01, 1.52e-01, 1.54e-01, 1.52e-01, 1.38e-01, 1.10e-01, 9.13e-02, 6.87e-02]
twist_prop = pi/180*[5.27e+01, 4.39e+01, 3.85e+01, 3.35e+01, 2.91e+01, 2.69e+01, 2.47e+01, 2.35e+01, 2.25e+01, 2.19e+01]

twist_70 = propwing.interp1(r, twist_prop, [0.7*Rtip])[1]
pitch = 30*pi/180 - twist_70
twist_prop += pitch + 6.0*pi/180

B = 6  # number of blades

aftype = CCBlade.af_from_aerodynfile("naca642-015a-Re300-fixed.dat")
n = length(r)
af = Array{CCBlade.AirfoilData}(n)
for i = 1:n
    af[i] = aftype
end

Dprop = Rtip*2  # propeller diameter
Omega = Vinf/(J*Dprop)*2*pi  # rotation rate

uwake, vwake, T, Q, reff = propwing.propanalysis(Rhub, Rtip, r, chord_prop, twist_prop, B, af, Vinf, Omega, rho, xtowing)

#TODO: check thrust coefficient make sure it is 0.168
n = Omega/(2*pi)
CT = T / (rho * n^2 * Dprop^4)
println(CT)
# println(Q*Omega)
# println(T*Vinf)


# --------- wing geometry -------
# AR = 5.33
semispan = 2.58/2
sectionspan = [0.35*semispan, 0.65*semispan]
twist_wing = [0, 0, 0]*pi/180
b = semispan*2
chord_wing = [0.279, 0.279, 0.161]
tc = [0.12, 0.12, 0.12]  # irrelevant for the aerodynamics
sweep = [0, -5]*pi/180
dihedral = [0, 0]
N = [35, 65]
cmac = 0.24

L, Di, CP, cl, cllocal, Vinfeff, alphaeff = propwing.winganalysis(sectionspan, chord_wing, twist_wing, tc, sweep, dihedral, N, alpha, rho, Vinf, reff, uwake, vwake, etaprop_center, cw)

exp_x = [4.61E-02, 1.39E-01, 2.44E-01, 2.87E-01, 4.00E-01, 4.38E-01, 5.62E-01, 6.87E-01, 8.09E-01]
# exp_x2 = [4.57E-02, 1.39E-01, 2.44E-01, 2.87E-01, 3.99E-01, 4.38E-01, 5.62E-01, 6.86E-01, 8.09E-01]
# exp_x3 = [4.61E-02, 1.39E-01, 2.44E-01, 2.88E-01, 4.00E-01, 4.38E-01, 5.62E-01, 6.86E-01, 8.08E-01]
exp_y1 = [5.12E-01, 5.79E-01, 7.54E-01, 6.86E-01, 3.31E-01, -5.93E-02, 3.20E-01, 3.15E-01, 2.83E-01]
exp_y2 = [4.97E-01, 5.63E-01, 7.09E-01, 5.03E-01, 1.44E-01, -1.01E-01, 3.06E-01, 3.03E-01, 2.70E-01]
exp_y3 = [5.27E-01, 5.97E-01, 7.98E-01, 8.70E-01, 5.14E-01, -1.47E-02, 3.31E-01, 3.29E-01, 2.95E-01]

other_vlm_x = [9.66E-03, 2.95E-02, 4.82E-02, 6.69E-02, 8.47E-02, 1.04E-01, 1.23E-01, 1.42E-01, 1.61E-01, 1.81E-01, 2.00E-01, 2.19E-01, 2.38E-01, 2.56E-01, 2.75E-01, 2.95E-01, 3.14E-01, 3.33E-01, 3.61E-01, 3.97E-01, 4.34E-01, 4.70E-01, 5.08E-01, 5.43E-01, 5.80E-01, 6.16E-01, 6.54E-01, 6.89E-01, 7.26E-01, 7.62E-01, 7.99E-01, 8.35E-01, 8.72E-01, 9.08E-01, 9.45E-01, 9.81E-01]
other_vlm_y = [5.18E-01, 5.20E-01, 5.16E-01, 5.19E-01, 5.22E-01, 5.30E-01, 5.41E-01, 5.55E-01, 5.73E-01, 6.02E-01, 6.75E-01, 8.13E-01, 8.74E-01, 8.83E-01, 8.44E-01, 7.54E-01, 5.29E-01, 4.39E-01, 4.19E-01, 1.67E-01, 1.12E-01, 1.63E-01, 2.70E-01, 2.95E-01, 3.04E-01, 3.06E-01, 3.03E-01, 2.97E-01, 2.89E-01, 2.78E-01, 2.65E-01, 2.50E-01, 2.32E-01, 2.09E-01, 1.78E-01, 1.29E-01]


rc("figure", figsize=(3.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.16, bottom=0.18, top=0.97, right=0.96)

figure()
plot(CP.y/semispan, cl.*CP.chord/cmac)
errorbar(exp_x, exp_y1, yerr=[exp_y1-exp_y2, exp_y3-exp_y1], fmt="o")
plot(other_vlm_x, other_vlm_y, "k--")
xlim([0 , 1])
xlabel(L"2y/b")
ylabel(L"c_l\ c / c_{mac}")
legend(["Our VLM", "Epema VLM", "Experimental"])
# savefig("../MDAO_paper/figures/epemaval2.pdf", transparent=true)



end
