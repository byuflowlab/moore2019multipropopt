# module comparison

path,_ = splitdir(@__FILE__)
# using PyPlot
#Packages I haven't developed
import Interpolations
using PyPlot
close("all")
rc("figure", figsize=(3.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.19, bottom=0.18, top=0.97, right=0.95)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
import JLD
import CSV
import Dierckx
import StatsBase

# Packages I played a very minor role
import PreComp
import Composites
import VTKtools

# Packages I played a minor role
import OptParams
import Atmosphere
import CCBlade #new CCBlade package w/o ND splined af data
import BeamFEA
import Snopt
import Gradients

# Packages I played a major role in
import BrentMin
import MotorPower
import Akima
# import VLM
include("oldvlm.jl")
import BPM
import AirfoilPrep

printiter = 1

include("$path/objective.jl")
include("$path/propstructmisc.jl")
include("$path/propwingmisc.jl")
include("$path/setup.jl")
include("$path/save_vtk.jl")

plyprops = plyproperties()

ctlparams = controlparams(
objective=:none, #:eta_total_cruise #:mass,
printfreq = 50,
savevtk = 500,
printobjective=true,
printdetailedoutput=true,
printconstraintviolations=true,
)

# Define general constraints
constraints = [
:p_materialfailure,
:p_localbuckling,
:w_materialfailure,
:w_localbuckling,

:totalmass,
:machtip,
:Lift,
:Drag,
:noise,
# :p_chordlamthick,
# :range,
:takeoff_dist,
:stall,
:capacity,
]

# Prop inputs
p_n_akima = 7
Dprop = 0.236
# p_radii =  1e-3*[17.5, 30, 45, 60, 75, 90, 105]/(Dprop/2)
p_radii = [0.148, 0.254237, 0.381356, 0.508475, 0.635593, 0.762712, 0.889831, 1.0]
p_chord_pts = 1e-3*[9.88, 11.88, 15.59, 18.81, 19.55, 18.32, 13.96, 0.01]
p_twist_pts = pi/180*[35.0, 32.5, 26.5, 23.5, 19, 16.5, 14.0, 10.0] #+ p_pitch
J = 0.85  # prop advance ratio
velocity = 50.0
p_rpm = velocity/(J*Dprop)*60.0  # rotation rate
etaprop_center = [-0.3/0.64,0.3/0.64]

#check thrust coefficient make sure it is 0.168
twist_75 = interp1(p_radii, p_twist_pts, [0.75])[1]
p_pitch = 25*pi/180 - twist_75
# println(p_pitch*180/pi)
p_twist_pts += p_pitch + 4.9*pi/180

#Wing Inputs
w_n_linear = 5
w_nondim_halfspany = [0.0,1.0]
AR = 5.33
w_halfspan = 0.64
b = w_halfspan*2
 w_chord =[b/AR, b/AR]
w_aoa = 4.0*pi/180





table_af_wing_blown = JLD.load("$path/airfoils/af_prop_E212.jld")
table_3D_wing_blown = table_af_wing_blown["NDtable"]
spline_3D_blown = AirfoilPrep.SplineND_from_tableND(table_3D_wing_blown)


table_af_wing_free = JLD.load("$path/airfoils/af_prop_E212.jld")
table_3D_wing_free = table_af_wing_free["NDtable"]
spline_3D_freestream = AirfoilPrep.SplineND_from_tableND(table_3D_wing_free)

# # Load the prop airfoil ND table
# ClarkY_non_extrap = JLD.load("$(path)/airfoils/af_prop_E212.jld")
# ClarkY_non_extrap = ClarkY_non_extrap["NDtable"]
# omega = p_rpm*pi/30
# TSR = omega * Dprop/2 / velocity
#
# grid_alphas=[i for i in -180:1.0:180]
# Rtip = Dprop/2
# chord_prop_pt = p_chord_pts
# rpt = p_radii * Rtip
# r = linspace(p_radii[2],p_radii[end-1],p_n_akima) * Rtip
# chord_prop = Akima.interp_noderiv(rpt, chord_prop_pt, r)
# af = Array{Any}(p_n_akima)
# for i = 1:p_n_akima
#
#    r_over_R = r[i]
#    c_over_r = chord_prop[i]/r[i]
#    NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(ClarkY_non_extrap,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)
#    afspl = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)
#    af[i] = afspl
# end
# JLD.save("$(path)/airfoils/af_veldhius.jld", "af", af)

af = JLD.load("$(path)/airfoils/af_veldhius.jld")
af = JLD.load("$(path)/airfoils/af.jld")
af = af["af"]
af = fill(af[1],p_n_akima)

af_xy = CSV.read("$(path)/airfoils/e212-il.csv",header = true, delim = ',',datarow = 10,allowmissing=:none)
xaf1 = reverse(af_xy[:,1])/100.0
yaf1 = reverse(af_xy[:,2])/100.0

# nameout = "$(path)/AirfoilData/e212.dat"
# write(nameout,"x\ty\n")
# open(nameout,"a") do x
#     for i = 1:length(xaf1)
#         write(x,"$(xaf1[i])\t$(yaf1[i])\n")
#     end
# end

#Filter airfoil to get points where strain is evaluated
n_strainpts = 10
xafstrain = zeros(n_strainpts)
yafstrain = zeros(n_strainpts)
step_simplify = round(Int,length(xaf1)/n_strainpts)
idx = 1
for i = 1:n_strainpts
    xafstrain[i] = xaf1[idx[1]]
    yafstrain[i] = yaf1[idx[1]]
    idx += step_simplify
end

# Fill airfoil data for each section
p_xaf = zeros(p_n_akima,length(xaf1))
p_yaf = zeros(p_n_akima,length(yaf1))
p_xafstrain = zeros(p_n_akima,length(xafstrain))
p_yafstrain = zeros(p_n_akima,length(yafstrain))
for i = 1:p_n_akima
    p_xaf[i,:] = copy(xaf1) #break links since using in multiple areas
    p_yaf[i,:] = copy(yaf1)
    p_xafstrain[i,:] = copy(xafstrain)
    p_yafstrain[i,:] = copy(yafstrain)
end

#TODO: different airfoil for wing
w_xaf = zeros(w_n_linear,length(xaf1))
w_yaf = zeros(w_n_linear,length(yaf1))
w_xafstrain = zeros(w_n_linear,length(xafstrain))
w_yafstrain = zeros(w_n_linear,length(yafstrain))
for i = 1:w_n_linear
    w_xaf[i,:] = copy(xaf1)
    w_yaf[i,:] = copy(yaf1)
    w_xafstrain[i,:] = copy(xafstrain)
    w_yafstrain[i,:] = copy(yafstrain)
end


altitude_cruise = 3000.0
design4 = propwingdesign(;
    #Design Variables at this scope
    batt_mass = 0.0,
    w_aoa = w_aoa,
    p_rpm =p_rpm,
    p_pitch = 1.0,
    velocity = velocity,

    #Design Parameters at this scope
    constraints = constraints,
    af = af,
    p_xaf = p_xaf,
    p_yaf = p_yaf,
    p_xafstrain = p_xafstrain,
    p_yafstrain = p_yafstrain,
    w_xaf = w_xaf,
    w_yaf = w_yaf,
    w_xafstrain = w_xafstrain,
    w_yafstrain = w_yafstrain,
    spline_3D_freestream = spline_3D_freestream,
    spline_3D_blown = spline_3D_blown,

    # Design variables unique to the system function
    p_lam_tin = [ones(p_radii)*0.001*10;ones(p_radii)*0.001*85],
    p_orientation = [0.0,45.0],
    p_chord_pts = p_chord_pts,
    p_twist_pts = p_twist_pts,
    p_Rtip = Dprop/2,
    kv = 100.0,
    i0 = 2.0,
    w_chord_pts = w_chord,
    w_twist_pts = [0.0,0.0]+w_aoa,
    w_halfspan = w_halfspan,
    w_sweep_pts = zeros(w_chord),
    w_dihedral_pts = zeros(w_chord),
    w_lam_tin = [
        ones(w_nondim_halfspany)*0.001*20;#t1 highmodulus_uni
        ones(w_nondim_halfspany)*0.001*20;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.01*2;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.001*20;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.01*2;#t2 highmodulus_weave
        ],
    w_orientation = [0.0,45.0,45.0],
    w_airfoilthickness = [0.0,0.0], #Does nothing

    #Design variables unique to the system function
    etaprop_center = etaprop_center,
    p_n_akima = p_n_akima,
    p_radii =  p_radii,
    altitude = 0.0,
    numprops = 2,
    p_x_offset = ones(p_radii)*0.25,
    p_aerocenter = ones(p_radii)*0.275,
    p_usedmaterials = ["highmodulus_uni","highmodulus_weave"],
    p_webloc =[0.999],
    plyprops = plyprops,
    blades = 4,
    w_n_linear = w_n_linear,
    w_x_offset = ones(w_nondim_halfspany)*0.25,
    w_aerocenter = ones(w_nondim_halfspany)*0.27,
    w_nondim_halfspany = w_nondim_halfspany,
    w_usedmaterials = ["highmodulus_uni","highmodulus_weave","taylor_foam","highmodulus_weave","taylor_foam"],
    w_webloc = [0.5],
    prop_tilt = 0.0*pi/180,

    #Control parameters
    detailed = true,
    VTKfilename = "Veldhius.vtk",
    savefinalvtk = true,
    verification = false,
    plots = true,
    n_CP = 200,

)

#Call the PropWing System
detailed_output4 = system(ctlparams,design4)

# println(detailed_output4.CT)

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


# 0 degree aoa
w_aoa = 0.0*pi/180
design0 = propwingdesign(;
    #Design Variables at this scope
    batt_mass = 0.0,
    w_aoa = w_aoa,
    p_rpm =p_rpm,
    p_pitch = 1.0,
    velocity = velocity,

    #Design Parameters at this scope
    constraints = constraints,
    af = af,
    p_xaf = p_xaf,
    p_yaf = p_yaf,
    p_xafstrain = p_xafstrain,
    p_yafstrain = p_yafstrain,
    w_xaf = w_xaf,
    w_yaf = w_yaf,
    w_xafstrain = w_xafstrain,
    w_yafstrain = w_yafstrain,
    spline_3D_freestream = spline_3D_freestream,
    spline_3D_blown = spline_3D_blown,

    # Design variables unique to the system function
    p_lam_tin = [ones(p_radii)*0.001*10;ones(p_radii)*0.001*85],
    p_orientation = [0.0,45.0],
    p_chord_pts = p_chord_pts,
    p_twist_pts = p_twist_pts,
    p_Rtip = Dprop/2,
    kv = 100.0,
    i0 = 2.0,
    w_chord_pts = w_chord,
    w_twist_pts = [0.0,0.0]+w_aoa,
    w_halfspan = w_halfspan,
    w_sweep_pts = zeros(w_chord),
    w_dihedral_pts = zeros(w_chord),
    w_lam_tin = [
        ones(w_nondim_halfspany)*0.001*20;#t1 highmodulus_uni
        ones(w_nondim_halfspany)*0.001*20;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.01*2;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.001*20;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.01*2;#t2 highmodulus_weave
        ],
    w_orientation = [0.0,45.0,45.0],
    w_airfoilthickness = [0.0,0.0], #Does nothing

    #Design variables unique to the system function
    etaprop_center = etaprop_center,
    p_n_akima = p_n_akima,
    p_radii =  p_radii,
    altitude = 0.0,
    numprops = 2,
    p_x_offset = ones(p_radii)*0.25,
    p_aerocenter = ones(p_radii)*0.275,
    p_usedmaterials = ["highmodulus_uni","highmodulus_weave"],
    p_webloc =[0.999],
    plyprops = plyprops,
    blades = 4,
    w_n_linear = w_n_linear,
    w_x_offset = ones(w_nondim_halfspany)*0.25,
    w_aerocenter = ones(w_nondim_halfspany)*0.27,
    w_nondim_halfspany = w_nondim_halfspany,
    w_usedmaterials = ["highmodulus_uni","highmodulus_weave","taylor_foam","highmodulus_weave","taylor_foam"],
    w_webloc = [0.5],
    prop_tilt = 0.0*pi/180,

    #Control parameters
    detailed = true,
    VTKfilename = "Veldhius.vtk",
    savefinalvtk = true,
    verification = false,
    plots = true,
    n_CP = 200,

)

#Call the PropWing System
detailed_output0 = system(ctlparams,design0)


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


# chord = [p.chord for p in detailed_output0.panels]
# N = length(detailed_output0.panels)
# CPy = zeros(N)
# for i = 1:N  # CP
#     rcp = VLM.control_point(detailed_output0.panels[i])
#     CPy[i] = rcp[2]
# end

chord = detailed_output4.panels.chord
CPy = detailed_output4.panels.y

figname = "Veldhius"
figure(figname)
plot(CPy/design4.w_halfspan, detailed_output4.cl)
plot(expdata[:, 1]/100, expdata[:, 2], "o")
# plot(CP.y/semispan, cllocal)

ax = gca()
ax[:set_color_cycle](nothing)

# figure()
plot(CPy/design0.w_halfspan, detailed_output0.cl)
plot(expdata0[:, 1]/100, expdata0[:, 2], "o")

xlabel("Halfspan Location")
ylabel("Local Lift Coefficient")
legend(["VLM", "Exper."], loc="upper right")
savefig("../MDAO_paper/figures/$figname.pdf", transparent=true)


# Error Calculation

spl_VLM0deg = Dierckx.Spline1D(CPy/design0.w_halfspan, detailed_output0.cl)
spl_VLM4deg = Dierckx.Spline1D(CPy/design4.w_halfspan, detailed_output4.cl)


error0 = 0.0
for i = 1:length(expdata0[:,1])
    error0 = error0+(spl_VLM0deg(expdata0[i, 1]/100)-expdata0[i, 2]).^2
end
sumofsquareserror0 = sqrt(error0/length(expdata[:, 1]))
println("sumofsquareserror0 $sumofsquareserror0")

error4 = 0.0
for i = 1:length(expdata[:,1])
    error4 = error4+(spl_VLM4deg(expdata[i, 1]/100)-expdata[i, 2]).^2
end
sumofsquareserror4 = sqrt(error4/length(expdata[:, 1]))
println("sumofsquareserror4 $sumofsquareserror4")
