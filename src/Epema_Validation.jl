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
p_n_akima = 6
Dprop = 0.406
# p_radii =  1e-3*[17.5, 30, 45, 60, 75, 90, 105]/(Dprop/2)
p_radii = [0.205, .231, .400, .511, .615, .720, .789, .860, .930, .965, .999, 1.0]
p_chord_pts = Dprop/2*[1.56e-01, 1.56e-01, 1.46e-01, 1.47e-01, 1.52e-01, 1.54e-01, 1.52e-01, 1.38e-01, 1.10e-01, 9.13e-02, 6.87e-02, 6.0e-02]
p_twist_pts = pi/180*[5.27e+01, 5.27e+01, 4.39e+01, 3.85e+01, 3.35e+01, 2.91e+01, 2.69e+01, 2.47e+01, 2.35e+01, 2.25e+01, 2.19e+01, 2.09e+01] #+ p_pitch
J = 0.695  # prop advance ratio
velocity = 19.0
p_rpm = velocity/(J*Dprop)*60.0  # rotation rate
etaprop_center = [-0.35,0.35]

#check thrust coefficient make sure it is 0.168
twist_75 = interp1(p_radii, p_twist_pts, [0.75])[1]
println(twist_75*180/pi)
p_pitch = 30*pi/180 - twist_75
# println(p_pitch*180/pi)
p_twist_pts += p_pitch - 2.0*pi/180 +2.0*pi/180

#Wing Inputs
w_n_linear = 30
w_nondim_halfspany = [0.0,0.35,1.0]
w_halfspan = 2.58/2
b = w_halfspan*2
w_chord =[0.279, 0.279, 0.161]
w_sweep = [0.0,-5,0.0]*pi/180
w_aoa = 4.0*pi/180
cmac = 0.24




table_af_wing_blown = JLD.load("$path/airfoils/af_prop_E212.jld")
table_3D_wing_blown = table_af_wing_blown["NDtable"]
spline_3D_blown = AirfoilPrep.SplineND_from_tableND(table_3D_wing_blown)


table_af_wing_free = JLD.load("$path/airfoils/af_prop_E212.jld")
table_3D_wing_free = table_af_wing_free["NDtable"]
spline_3D_freestream = AirfoilPrep.SplineND_from_tableND(table_3D_wing_free)

# # Load the prop airfoil ND table
# E212_non_extrap = JLD.load("$(path)/airfoils/af_prop_E212.jld")
# E212_non_extrap = E212_non_extrap["NDtable"]
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
#    NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(E212_non_extrap,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)
#    afspl = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)
#    af[i] = afspl
# end
# JLD.save("$(path)/airfoils/af_epema.jld", "af", af)

af = JLD.load("$(path)/airfoils/af_epema.jld")
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
designE = propwingdesign(;
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
    w_twist_pts = zeros(w_chord)+w_aoa,
    w_halfspan = w_halfspan,
    w_sweep_pts = w_sweep,
    w_dihedral_pts = zeros(w_chord),
    w_lam_tin = [
        ones(w_nondim_halfspany)*0.001*20;#t1 highmodulus_uni
        ones(w_nondim_halfspany)*0.001*20;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.01*2;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.001*20;#t2 highmodulus_weave
        ones(w_nondim_halfspany)*0.01*2;#t2 highmodulus_weave
        ],
    w_orientation = [0.0,45.0,45.0],
    w_airfoilthickness = zeros(w_chord), #Does nothing

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
    blades = 6,
    w_n_linear = w_n_linear,
    w_x_offset = ones(w_nondim_halfspany)*0.25,
    w_aerocenter = ones(w_nondim_halfspany)*0.27,
    w_nondim_halfspany = w_nondim_halfspany,
    w_usedmaterials = ["highmodulus_uni","highmodulus_weave","taylor_foam","highmodulus_weave","taylor_foam"],
    w_webloc = [0.5],
    prop_tilt = 0.0*pi/180,

    #Control parameters
    detailed = true,
    VTKfilename = "Epema.vtk",
    savefinalvtk = true,
    verification = false,
    plots = false,
    n_CP = 200,

)

#Call the PropWing System
detailed_outputE = system(ctlparams,designE)

println(detailed_outputE.CT)

exp_x = [4.61E-02, 1.39E-01, 2.44E-01, 2.87E-01, 4.00E-01, 4.38E-01, 5.62E-01, 6.87E-01, 8.09E-01]
# exp_x2 = [4.57E-02, 1.39E-01, 2.44E-01, 2.87E-01, 3.99E-01, 4.38E-01, 5.62E-01, 6.86E-01, 8.09E-01]
# exp_x3 = [4.61E-02, 1.39E-01, 2.44E-01, 2.88E-01, 4.00E-01, 4.38E-01, 5.62E-01, 6.86E-01, 8.08E-01]
exp_y1 = [5.12E-01, 5.79E-01, 7.54E-01, 6.86E-01, 3.31E-01, -5.93E-02, 3.20E-01, 3.15E-01, 2.83E-01]
exp_y2 = [4.97E-01, 5.63E-01, 7.09E-01, 5.03E-01, 1.44E-01, -1.01E-01, 3.06E-01, 3.03E-01, 2.70E-01]
exp_y3 = [5.27E-01, 5.97E-01, 7.98E-01, 8.70E-01, 5.14E-01, -1.47E-02, 3.31E-01, 3.29E-01, 2.95E-01]

other_vlm_x = [9.66E-03, 2.95E-02, 4.82E-02, 6.69E-02, 8.47E-02, 1.04E-01, 1.23E-01, 1.42E-01, 1.61E-01, 1.81E-01, 2.00E-01, 2.19E-01, 2.38E-01, 2.56E-01, 2.75E-01, 2.95E-01, 3.14E-01, 3.33E-01, 3.61E-01, 3.97E-01, 4.34E-01, 4.70E-01, 5.08E-01, 5.43E-01, 5.80E-01, 6.16E-01, 6.54E-01, 6.89E-01, 7.26E-01, 7.62E-01, 7.99E-01, 8.35E-01, 8.72E-01, 9.08E-01, 9.45E-01, 9.81E-01]
other_vlm_y = [5.18E-01, 5.20E-01, 5.16E-01, 5.19E-01, 5.22E-01, 5.30E-01, 5.41E-01, 5.55E-01, 5.73E-01, 6.02E-01, 6.75E-01, 8.13E-01, 8.74E-01, 8.83E-01, 8.44E-01, 7.54E-01, 5.29E-01, 4.39E-01, 4.19E-01, 1.67E-01, 1.12E-01, 1.63E-01, 2.70E-01, 2.95E-01, 3.04E-01, 3.06E-01, 3.03E-01, 2.97E-01, 2.89E-01, 2.78E-01, 2.65E-01, 2.50E-01, 2.32E-01, 2.09E-01, 1.78E-01, 1.29E-01]

# chord = [p.chord for p in detailed_outputE.panels]
# N = length(detailed_outputE.panels)
# CPy = zeros(N)
# for i = 1:N  # CP
#     rcp = VLM.control_point(detailed_outputE.panels[i])
#     CPy[i] = rcp[2]
# end

chord = detailed_outputE.panels.chord
CPy = detailed_outputE.panels.y

figname = "Epema"
figure(figname)
plot(CPy/designE.w_halfspan, detailed_outputE.cl.*chord/cmac)
errorbar(exp_x, exp_y1, yerr=[exp_y1-exp_y2, exp_y3-exp_y1], fmt="o")
plot(other_vlm_x, other_vlm_y, "k--")
xlim([0 , 1])
xlabel("Halfspan Location")
ylabel("Local Lift Coefficient")#L"c_l\ c / c_{mac}")
legend(["Our VLM", "Epema VLM", "Epema Exper."])
savefig("../MDAO_paper/figures/$figname.pdf", transparent=true)


# Error Calculation

spl_VLM0deg = Dierckx.Spline1D(CPy/designE.w_halfspan, detailed_outputE.cl.*chord/cmac)



error0 = 0.0
for i = 1:length(exp_x)
    error0 = error0+(spl_VLM0deg(exp_x[i])-exp_y1[i]).^2
end
sumofsquareserror0 = sqrt(error0/length(exp_x))
println("sumofsquareserror0 $sumofsquareserror0")
