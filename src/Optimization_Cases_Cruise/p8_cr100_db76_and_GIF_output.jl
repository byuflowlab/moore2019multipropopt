# module moore2018multipropopt

path,_ = splitdir(@__FILE__); path = "$path/../"
# using PyPlot
#Packages I haven't developed
import Interpolations
# using PyPlot
import JLD
import CSV
import Dierckx

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

# Packages I played a major role
import BrentMin
import MotorPower
import Akima
# import VLM
include("$path/oldvlm.jl")
import BPM
import AirfoilPrep

printiter = 37171
#TODO verify this works with parallel


include("$path/objective.jl")
include("$path/propstructmisc.jl")
include("$path/propwingmisc.jl")
include("$path/setup.jl")
include("$path/save_vtk.jl")

# # include("$path/../../BeamFEA.jl/src/BeamFEA.jl")
# rc("figure", figsize=(5.8, 2.6))
# rc("font", size=10.0)
# rc("lines", linewidth=1.5)
# rc("lines", markersize=3.0)
# rc("legend", frameon=false)
# rc("axes.spines", right=false, top=false)
# rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.69)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
# close("all")

# function optimize(numprops,to_factor,db_factor)

    # define run controls
    ctlparams = controlparams(
    objective= :takeoff, #:range, #:mass,#:range #:ETA_total_cruise,
    printfreq = 70,
    savevtk = 70,
    printobjective=true,
    printdetailedoutput=true,
    printconstraintviolations=true,
    )

    designvars =  [
    # :p_orientation,
    :p_lam_t, #1-12
    :p_chord, #13-18
    :p_twist_d, #19-24
    :p_pitch_d_accel, #25
    :p_pitch_d_climb, #25
    :p_pitch_d_cruise, #25

    :velocity_accel,
    :velocity_climb,
    # :velocity_cruise,
    :p_Rtip,
    :p_rpm_accel,
    :p_rpm_climb,
    :p_rpm_cruise,
    # :number_blades,
    :kv, #28
    :i0, #29


    # :w_orientation,
    # :w_chord, #30-32
    :w_aoa_d_accel,   #multi
    :w_aoa_d_climb,   #multi
    :w_aoa_d_cruise,   #multi
    # :w_twist_d, #33-35
    # :w_halfspan,
    # :w_sweep_d,
    # :w_dihedral_d,
    # :w_lam_t, #36-end
    :batt_mass_accel,
    :batt_mass_climb,
    :batt_mass_cruise,
    ]

    # Define general constraints
    constraints = [
    :p_materialfailure,
    :p_localbuckling,
    :p_chordlamthick,
    # :w_materialfailure,
    # :w_localbuckling,
    # :w_chordlamthick,

    :totalmass,
    :machtip,
    :Lift,
    :Drag,
    :noise,

    :range,
    # :takeoff_dist,
    :cruisespeed,
    :radius,
    :flightpathangle,
    :stall,
    :capacity,
    :p_alphas,
    :p_separation,
    # :span,
    ]



    n_CP = 120.0
    payload = 684.0 + 245.0 #airframe + payload
    max_mass = 1230.0
    gravity = 9.81
    h_set = 15.24  #meters or 50 feet
    specif_energy = 300.0  #wh/kg
    required_range = 50.0*1000  #km
    altitude_accel = 0.0
    altitude_climb = 0.0
    altitude_cruise = 155.0
    cruise_factor = 100.0
    required_cruise = cruise_factor * 0.447 #MPH TO m/s


    # designvars, constraints, ctlparams = getoptparams()
    numprops = 8
    to_factor = 1
    db_factor = 1

    noise_factor = 76.0
    required_noise = noise_factor#68.0/db_factor
    required_takeoff_dist = 395.0/to_factor  #meters


    plyprops = plyproperties()

    matnames = plyprops.names

    name = "highmodulus_uni"
    idx = find(matnames -> matnames == name,matnames)
    uni_t = plyprops.plies[idx[1]].t

    name = "highmodulus_weave"
    idx = find(matnames -> matnames == name,matnames)
    weave_t = plyprops.plies[idx[1]].t

    name = "taylor_foam"
    idx = find(matnames -> matnames == name,matnames)
    foam_t = plyprops.plies[idx[1]].t

    # Prop inputs
    p_n_akima = 8
    Rtip = 0.5
    p_radii=[0.1,0.4,0.7,1.0]#[0.1,0.25,0.4,0.6,0.8,1.0]
    full_twist = ones(p_radii)*30.0
    full_chord = ones(p_radii)*0.07
    full_p_x_offset = ones(p_radii)*0.25
    full_p_aerocenter = ones(p_radii)*0.275

    p_lam_t = [
    ones(p_radii)*uni_t*10;#t1 highmodulus_uni
    ones(p_radii)*weave_t*85;#t2 highmodulus_weave
    ]
    p_usedmaterials = ["highmodulus_uni","highmodulus_weave"]
    p_orientation = [0.0,45.0]

    #Wing inputs
    # w_n_linear = 5
    # w_chord = [0.805,0.705,0.68]
    # w_halfspan = 10.5
    # w_nondim_halfspany = [0.0,.2,1.0]
    # full_w_x_offset = ones(w_nondim_halfspany)*0.25
    # full_w_aerocenter = ones(w_nondim_halfspany)*0.275
    # w_dihedral_d = zeros(w_chord)
    # w_sweep_d = ones(w_chord)*0.0

    #TODO: Technam inputs
    w_n_linear = 5
    w_halfspan = 11.4
    w_chord = [1.8/7,1.8/7,1.0/7]*w_halfspan
    w_nondim_halfspany = [0.0,4.0/7,1.0]
    full_w_x_offset = ones(w_nondim_halfspany)*0.25
    full_w_x_offset[end] = 1.0
    full_w_aerocenter = ones(w_nondim_halfspany)*0.275
    w_dihedral_d = zeros(w_chord)
    w_sweep_d = ones(w_chord)*0.0


    w_lam_t = [
    ones(w_nondim_halfspany)*uni_t*20;#t1 highmodulus_uni
    ones(w_nondim_halfspany)*weave_t*20;#t2 highmodulus_weave
    ones(w_nondim_halfspany)*foam_t*2;#t2 highmodulus_weave
    ones(w_nondim_halfspany)*weave_t*20;#t2 highmodulus_weave
    ones(w_nondim_halfspany)*foam_t*2;#t2 highmodulus_weave
    ]

    w_usedmaterials = ["highmodulus_uni","highmodulus_weave","taylor_foam","highmodulus_weave","taylor_foam"]
    w_orientation = [0.0,45.0,45.0] #uni, weave, webweave

    designparams = designparameters(;
    p_n_akima = p_n_akima,
    p_radii=p_radii,
    p_webloc=[0.999],
    altitude=0.0,
    plyprops=plyprops,
    p_usedmaterials = p_usedmaterials,
    p_orientation=p_orientation, #an p_orientation for each material at each section
    # structural design variables
    p_lam_t = p_lam_t,
    # geometric design variables
    p_chord = full_chord,
    p_sweep_d = fill(0.0,length(p_radii)-1),
    p_twist_d = full_twist,
    p_pitch_d = 1.0,
    p_dihedral_d = fill(0.0,length(p_radii)-1),
    p_airfoilthickness = fill(1.0,length(p_radii)), #not used
    p_aerocenter = full_p_aerocenter,
    p_x_offset = full_p_x_offset,
    # operating point design variables
    velocity = 25.0,
    # propulsion design variables
    p_Rtip = Rtip,
    blades = 2,
    numprops = numprops,
    machtip = 0.8,
    p_rpm = 1000.0,
    kv = 50.271513933940295,
    i0 = 2.0,
    batt_mass = 100.0,
    # Wing Design variables
    w_n_linear = w_n_linear, w_webloc=[0.5],
    w_lam_t = w_lam_t,
    w_chord = w_chord,
    w_aoa_d = 12.0,
    w_twist_d = ones(w_chord)*0.0,
    w_halfspan = w_halfspan,
    w_nondim_halfspany = w_nondim_halfspany,
    w_sweep_d = w_sweep_d,
    w_dihedral_d = w_dihedral_d,
    w_usedmaterials = w_usedmaterials,
    w_orientation = w_orientation,
    w_aerocenter = full_w_aerocenter,
    w_x_offset = full_w_x_offset,

    n_CP = n_CP, # 60.0,
    payload = payload, # 684.0,
    max_mass = max_mass, # 1500.0,
    gravity = gravity, # 9.81,
    h_set = h_set, # 20.0, #meters
    specif_energy = specif_energy, # 300.0, #wh/kg
    required_range = required_range, # 150000.0, #km
    required_takeoff_dist = required_takeoff_dist, # 395.0, #meters
    altitude_accel = altitude_accel, # 0.0,
    altitude_climb = altitude_climb, # 0.0,
    altitude_cruise = altitude_cruise, # 3000.0,
    required_noise = required_noise,
    required_cruise = required_cruise,

    )

    optxbounds = OptParams.initdict()
    getoptxbounds!(optxbounds,designparams)

    # get inputs to objective function
    x0,lb,ub = OptParams.assembleinput(designvars,optxbounds)
    rangevar = OptParams.getrangedict(designvars,optxbounds)



#             number_blades = 5.0
#
#
#             #Slow
#                 kv = 5.0
#                 i0 = 2.2037300540130143
#
#                 p_lam_t = [0.0007782, 0.000674115, 0.00050298, 0.000238061, 0.000218, 0.000218, 0.000218, 0.000218]
#                 p_chord = [0.217695, 0.179515, 0.0962492, 0.0129004]
#                 p_twist_d = [39.721, 18.0451, 10.7105, 4.20071]
#                 p_pitch_d_accel = 9.021458633187988
#                 p_pitch_d_climb = 2.6779561382353223
#                 p_pitch_d_cruise = 9.955675862211843
#                 velocity_accel = 21.907192616714354
#                 velocity_climb = 15.022172822158165
#                 velocity_cruise = 26.141124740636286
#                 p_Rtip = 0.1470747752685889
#                 p_rpm_accel = 4550.787583268113
#                 p_rpm_climb = 3617.25555613192
#                 p_rpm_cruise = 4206.390824864598
#                 w_aoa_d_accel = 11.805258303183681
#                 w_aoa_d_climb = 12.338223392213303
#                 w_aoa_d_cruise = 7.444036776377783
#                 batt_mass_accel = 0.01
#                 batt_mass_climb = 1.78807293443778
#                 batt_mass_cruise = 48.4944436895349
#
#                 number_blades = 2.0
#
#                 #SLOW
#                     kv = 5.0
#                     i0 = 5.999996165805586
#
#                 p_lam_t = [0.0106565, 0.00233204, 0.000521854, 0.000643721, 0.000235849, 0.000218, 0.000221644, 0.00023394]
#     p_chord = [0.152893, 0.189842, 0.0122381, 0.0126047]
#     p_twist_d = [43.6566, 25.3336, 14.8314, 2.15019]
#     p_pitch_d_accel = -1.8829344927066909
#     p_pitch_d_climb = -6.623043953201483
#     p_pitch_d_cruise = 3.1389505191969
#     velocity_accel = 10.516263781338967
#     velocity_climb = 22.110235464233497
#     velocity_cruise = 22.34991575727592
#     p_Rtip = 1.4607155806562915
#     p_rpm_accel = 1363.801616202054
#     p_rpm_climb = 1762.6245635737948
#     p_rpm_cruise = 785.088616955683
#     number_blades = 2.018354828584329
#     kv = 5.0000000099351745
#     i0 = 5.999999998436644
#     w_aoa_d_accel = 25.7671847059628
#     w_aoa_d_climb = 24.085518334858087
#     w_aoa_d_cruise = 7.6874985270624485
#     batt_mass_accel = 0.01
#     batt_mass_climb = 1.2565648978309747
#     batt_mass_cruise = 43.031131277642885
#
#     p_lam_t = [0.003252, 0.00142364, 0.000362052, 0.00099775, 0.000235849, 0.000218, 0.000221119, 0.000230127]
#     p_chord = [0.175241, 0.162084, 0.0107099, 0.0100064]
#     p_twist_d = [31.4471, 26.7358, 19.7654, 9.63553]
#     p_pitch_d_accel = -7.719866630174941
#     p_pitch_d_climb = -13.960686428934565
#     p_pitch_d_cruise = 1.173855411701484
#     velocity_accel = 10.178355866604997
#     velocity_climb = 10.271952373597921
#     p_Rtip = 1.3668277834692588
#     p_rpm_accel = 1145.630443055776
#     p_rpm_climb = 1665.6227960199028
#     p_rpm_cruise = 788.352329791602
#     number_blades = 2.445602851739818
#     kv = 6.105269032410999
#     i0 = 4.968452555691652
#     w_aoa_d_accel = 26.193250259194144
#     w_aoa_d_climb = 23.474670351125177
#     w_aoa_d_cruise = 7.79790809369224
#     batt_mass_accel = 0.01
#     batt_mass_climb = 0.6598943736785684
#     batt_mass_cruise = 99.85217871333208
#
#     p_lam_t = [0.00349585, 0.00154944, 0.000152, 0.00099775, 0.000235849, 0.000218, 0.000221119, 0.000230127]
# p_chord = [0.175239, 0.161664, 0.0110125, 0.0100079]
# p_twist_d = [32.2629, 27.3062, 21.6001, 10.2174]
# p_pitch_d_accel = -7.719866630174941
# p_pitch_d_climb = -13.924415137156194
# p_pitch_d_cruise = 1.173855411701484
# velocity_accel = 10.178355866604997
# velocity_climb = 10.197874857065022
# p_Rtip = 1.3668274103862978
# p_rpm_accel = 1145.630443055776
# p_rpm_climb = 1674.2859783769227
# p_rpm_cruise = 788.352329791602
# number_blades = 2.1438870567106965
# kv = 6.105269032410999
# i0 = 4.968452555691652
# w_aoa_d_accel = 26.145740407471234
# w_aoa_d_climb = 23.522365831369434
# w_aoa_d_cruise = 7.797908093691838
# batt_mass_accel = 0.01
# batt_mass_climb = 0.43827470796137974
# batt_mass_cruise = 98.31020185664771
#
# number_blades = 2.0
#
# p_lam_t = [0.00296032, 0.0015862, 0.000162432, 0.000559344, 0.000218, 0.000218, 0.000218, 0.000218]
# p_chord = [0.152775, 0.209865, 0.0599548, 0.0146912]
# p_twist_d = [43.7258, 36.3345, 24.8117, 15.5926]
# p_pitch_d_accel = 1.4774821755019851
# p_pitch_d_climb = 0.07394619362544558
# p_pitch_d_cruise = 1.2375953673499482
# velocity_accel = 19.472605774436627
# velocity_climb = 21.727538432606607
# velocity_cruise = 44.93261088033517
# p_Rtip = 0.7373054880998322
# p_rpm_accel = 1521.2024594724965
# p_rpm_climb = 2662.377310410668
# p_rpm_cruise = 2030.4950665640044
# number_blades = 2.0#2.3818948496117045
# kv = 9.01975436914855
# i0 = 5.942697466865178
# w_aoa_d_accel = 26.783254994297305
# w_aoa_d_climb = 15.682348122387333
# w_aoa_d_cruise = 1.9046025617909106
# batt_mass_accel = 0.01
# batt_mass_climb = 0.44404948094245933
# batt_mass_cruise = 119.39074611391598
#
# p_lam_t = [0.00491215, 0.00207895, 0.00109735, 0.000951505, 0.000218, 0.000218, 0.000218, 0.000218]
# p_chord = [0.151192, 0.325948, 0.192366, 0.0520607]
# p_twist_d = [45.0, 42.1689, 28.5539, 24.6288]
# p_pitch_d_accel = 0.5160959129483517
# p_pitch_d_climb = -6.379369415208218
# p_pitch_d_cruise = 1.2554747468011076
# velocity_accel = 17.504799071392068
# velocity_climb = 16.011054271640106
# velocity_cruise = 45.556861878613994
# p_Rtip = 0.6333333333333337
# p_rpm_accel = 1700.7399651777887
# p_rpm_climb = 2214.9896316661266
# p_rpm_cruise = 2033.7975016015014
# kv = 9.010532347513744
# i0 = 5.983991910449766
# w_aoa_d_accel = 26.4773214375485
# w_aoa_d_climb = 24.947976808858428
# w_aoa_d_cruise = 6.09880815849964
# batt_mass_accel = 0.0012680910284242287
# batt_mass_climb = 0.2174067264581435
# batt_mass_cruise = 122.52852597671871
#
# #With own af data
# p_lam_t = [0.00893571, 0.00234997, 0.00182375, 0.000152, 0.00083697, 0.000953682, 0.000876195, 0.000540092]
# p_chord = [0.106509, 0.150221, 0.0692545, 0.05]
# p_twist_d = [45.0, 39.3542, 30.3407, 28.5409]
# p_pitch_d_accel = -14.801534303269628
# p_pitch_d_climb = -16.663949276805756
# p_pitch_d_cruise = 1.116079551638493
# velocity_accel = 9.944097279374834
# velocity_climb = 7.73005569315936
# p_rpm_accel = 3785.4038793504687
# p_rpm_climb = 4097.5298786080875
# p_rpm_cruise = 1942.562723001659
# kv = 5.0
# i0 = 3.9951737452421248
# w_aoa_d_accel = 32.04188683229715
# w_aoa_d_climb = 31.337251317158273
# w_aoa_d_cruise = 4.691196096841078
# batt_mass_accel = 0.001
# batt_mass_climb = 1.1260003355278403
# batt_mass_cruise = 74.94512715737756
#
# p_lam_t = [0.00888127, 0.00239683, 0.00121549, 0.000152, 0.00083697, 0.000953682, 0.000876195, 0.000540092]
# p_chord = [0.116602, 0.143469, 0.0767974, 0.0262864]
# p_twist_d = [45.0, 38.8659, 30.5658, 29.6107]
# p_pitch_d_accel = -15.245276832060593
# p_pitch_d_climb = -16.990543523778836
# p_pitch_d_cruise = 1.1216558284591813
# velocity_accel = 9.951767417608675
# velocity_climb = 7.696453282205797
# p_rpm_accel = 3771.664591019246
# p_rpm_climb = 4097.529878608089
# p_rpm_cruise = 1943.4767749266518
# kv = 5.000001441672897
# i0 = 3.995010265667068
# w_aoa_d_accel = 32.07085586412711
# w_aoa_d_climb = 31.354963988641288
# w_aoa_d_cruise = 4.698885062498646
# batt_mass_accel = 0.001
# batt_mass_climb = 1.009896688608764
# batt_mass_cruise = 75.72755761609312
#
# p_lam_t = [0.00888127, 0.00239683, 0.00121549, 0.000152, 0.00083697, 0.000953682, 0.000876195, 0.000540092]
# p_chord = [0.116602, 0.143469, 0.0767974, 0.0262864]
# p_twist_d = [45.0, 38.8659, 30.5658, 29.6107]
# p_pitch_d_accel = -15.245276832060593
# p_pitch_d_climb = -16.990543523778836
# p_pitch_d_cruise = 1.1216558284591813
# velocity_accel = 9.951767417608675
# velocity_climb = 7.696453282205797
# p_rpm_accel = 3771.664591019246
# p_rpm_climb = 4097.529878608089
# p_rpm_cruise = 1943.4767749266518
# kv = 5.000001441672897
# i0 = 3.995010265667068
# w_aoa_d_accel = 32.07085586412711
# w_aoa_d_climb = 31.354963988641288
# w_aoa_d_cruise = 4.698885062498646
# batt_mass_accel = 0.001
# batt_mass_climb = 1.009896688608764
# batt_mass_cruise = 75.72755761609312
#
# #fixed aoa
# p_lam_t = [0.0090107, 0.00229688, 0.00114354, 0.000152023, 0.00083697, 0.000953682, 0.000218, 0.000218]
# p_chord = [0.102455, 0.150463, 0.0727664, 0.0242614]
# p_twist_d = [45.0, 39.5646, 30.6596, 28.9851]
# p_pitch_d_accel = -15.291890229352216
# p_pitch_d_climb = -16.82857992276188
# p_pitch_d_cruise = 1.1252752346571662
# velocity_accel = 7.719281971361207
# velocity_climb = 11.864801002989012
# p_rpm_accel = 3804.787843784844
# p_rpm_climb = 4097.529878608089
# p_rpm_cruise = 1943.5718130812047
# kv = 5.0
# i0 = 3.9791093851566774
# w_aoa_d_accel = 31.23708813422462
# w_aoa_d_climb = 32.507899772141734
# w_aoa_d_cruise = 4.649161517426887
# batt_mass_accel = 0.001
# batt_mass_climb = 0.7805441172811979
# batt_mass_cruise = 74.80296024119511
#
# #Fixed Structures
# p_lam_t = [0.007538, 0.00173269, 0.00088642, 0.000152, 0.00083697, 0.000953682, 0.000218, 0.000218]
# p_chord = [0.151763, 0.151763, 0.151763, 0.151763]
# p_twist_d = [45.0, 39.1948, 30.4775, 29.6842]
# p_pitch_d_accel = -15.296567262681904
# p_pitch_d_climb = -16.56620424726087
# p_pitch_d_cruise = 1.1254471622631674
# velocity_accel = 8.248197902866446
# velocity_climb = 12.52998048101345
# p_rpm_accel = 3811.247887669918
# p_rpm_climb = 4097.529878608089
# p_rpm_cruise = 1943.6393645230926
# kv = 5.008349580516639
# i0 = 3.966580647765332
# w_aoa_d_accel = 30.770162029306185
# w_aoa_d_climb = 32.18321792176881
# w_aoa_d_cruise = 4.655881487309745
# batt_mass_accel = 0.001
# batt_mass_climb = 0.8400014116756471
# batt_mass_cruise = 74.58121110817127
#
# w_halfspan = 5.7
#
# #Example GIF
# p_lam_t = [0.007538, 0.00173269, 0.00088642, 0.000152, 0.00083697, 0.000953682, 0.000218, 0.000218]
# p_chord = [0.151763, 0.151763, 0.151763, 0.151763]
# p_twist_d = [29.0, 29.0, 29.0, 29.0]
# p_pitch_d_accel = -15.296567262681904
# p_pitch_d_climb = -16.56620424726087
# p_pitch_d_cruise = 1.1254471622631674
# velocity_accel = 8.248197902866446
# velocity_climb = 12.52998048101345
# p_rpm_accel = 3811.247887669918
# p_rpm_climb = 4097.529878608089
# p_rpm_cruise = 1943.6393645230926
# kv = 5.008349580516639
# i0 = 3.966580647765332
# w_aoa_d_accel = 30.770162029306185
# w_aoa_d_climb = 32.18321792176881
# w_aoa_d_cruise = 4.655881487309745
# batt_mass_accel = 0.001
# batt_mass_climb = 0.8400014116756471
# batt_mass_cruise = 74.58121110817127
#
# w_halfspan = 5.7
#
# #Fixed Structures
# p_lam_t = [0.007538, 0.00173269, 0.00088642, 0.000152, 0.00083697, 0.000953682, 0.000218, 0.000218]*1.5
# p_chord = [0.101763, 0.153117, 0.0724272, 0.0237822]*1.5
# p_twist_d = [45.0, 39.1948, 30.4775, 29.6842]*0.9
# p_pitch_d_accel = -15.296567262681904*1.1
# p_pitch_d_climb = -16.56620424726087*1.1
# p_pitch_d_cruise = 1.1254471622631674*1.1
# velocity_accel = 8.248197902866446*1.1
# velocity_climb = 12.52998048101345*1.1
# p_rpm_accel = 3811.247887669918*1.1
# p_rpm_climb = 4097.529878608089*1.1
# p_rpm_cruise = 1943.6393645230926*1.1
# kv = 5.008349580516639
# i0 = 3.966580647765332
# w_aoa_d_accel = 30.770162029306185*1.1
# w_aoa_d_climb = 32.18321792176881*1.1
# w_aoa_d_cruise = 4.655881487309745*1.1
# batt_mass_accel = 0.001*1.1
# batt_mass_climb = 0.8400014116756471*1.1
# batt_mass_cruise = 74.58121110817127*1.1

w_halfspan = 5.7
p_chord = ones(4)*.2
p_twist_d = ones(4)*30.0
p_pitch_d_accel = -5.0
p_pitch_d_climb = -3.0
p_pitch_d_cruise = 0.0
velocity_accel = 13.0
velocity_climb = 21.0
velocity_cruise = 45.556861878613994
p_rpm_accel = 10000.0
p_rpm_climb = 10000.0
p_rpm_cruise = 5000.0
kv = 10.0
i0 = 5.0
w_aoa_d_accel = 15.0
w_aoa_d_climb = 30.0
w_aoa_d_cruise = 5.0
batt_mass_accel = 1.0
batt_mass_climb = 1.0
batt_mass_cruise = 20.0
p_Rtip = 0.6
number_blades = 2.0

p_lam_t = [0.0103505, 0.000285035, 0.000734952, 0.000348269, 0.000388793, 0.0017972, 0.000717596, 0.000407097]
p_chord = [0.131096, 0.151933, 0.0718815, 0.0189429]
p_twist_d = [44.8067, 38.1564, 30.961, 29.4853]
p_pitch_d_accel = -15.281648370110014
p_pitch_d_climb = -16.16372411483507
p_pitch_d_cruise = 1.3095227963163993
velocity_accel = 8.148039950692041
velocity_climb = 14.409762122371184
p_Rtip = 0.5814151924837737
p_rpm_accel = 4383.71038341194
p_rpm_climb = 4256.8755486191085
p_rpm_cruise = 2225.53942029128
kv = 5.000183261626615
i0 = 3.3656590406139397
w_aoa_d_accel = 31.226885842207395
w_aoa_d_climb = 28.68626884140169
w_aoa_d_cruise = 4.829112321652306
batt_mass_accel = 0.001
batt_mass_climb = 2.1438224245890556
batt_mass_cruise = 88.06496124285786

p_lam_t2 = [0.00317757, 0.000785375, 0.000970208, 0.000323186, 0.00101117, 0.00187457, 0.000479939, 0.000345597]
p_chord2 = [0.127766, 0.142388, 0.0683798, 0.0267748]
p_twist_d2 = [44.966, 39.5187, 30.8776, 28.7751]
p_pitch_d_accel2 = -16.298004244335992
p_pitch_d_climb2 = -16.26869739354685
p_pitch_d_cruise2 = 1.2379918784894843
p_Rtip2 = 0.6333333333333337
velocity_accel2 = 7.85487134395491
velocity_climb2 = 12.545434643740208
p_rpm_accel2 = 4095.136682148108
p_rpm_climb2 = 4097.501219539165
p_rpm_cruise2 = 1935.8454972489915
kv2 = 5.000390427715033
i02 = 3.9838408230865894
w_aoa_d_accel2 = 30.519865300480163
w_aoa_d_climb2 = 32.156697969032514
w_aoa_d_cruise2 = 4.655850131710671
batt_mass_accel2 = 0.0010000086985410138
batt_mass_climb2 = 0.8418729122777846
batt_mass_cruise2 = 74.2385048375922


blend = 0.2

p_lam_t = p_lam_t*blend + p_lam_t2*(1-blend)
p_chord = p_chord*blend + p_chord2*(1-blend)
p_twist_d = p_twist_d*blend + p_twist_d2*(1-blend)
p_pitch_d_accel = p_pitch_d_accel*blend + p_pitch_d_accel2*(1-blend)
p_pitch_d_climb = p_pitch_d_climb*blend + p_pitch_d_climb2*(1-blend)
p_pitch_d_cruise = p_pitch_d_cruise*blend + p_pitch_d_cruise2*(1-blend)
p_Rtip = p_Rtip*blend + p_Rtip2*(1-blend)
velocity_accel = velocity_accel*blend + velocity_accel2*(1-blend)
velocity_climb = velocity_climb*blend + velocity_climb2*(1-blend)
p_rpm_accel = p_rpm_accel*blend + p_rpm_accel2*(1-blend)
p_rpm_climb = p_rpm_climb*blend + p_rpm_climb2*(1-blend)
p_rpm_cruise = p_rpm_cruise*blend + p_rpm_cruise2*(1-blend)
kv = kv*blend + kv2*(1-blend)
i0 = i0*blend + i02*(1-blend)
w_aoa_d_accel = w_aoa_d_accel*blend + w_aoa_d_accel2*(1-blend)
w_aoa_d_climb = w_aoa_d_climb*blend + w_aoa_d_climb2*(1-blend)
w_aoa_d_cruise = w_aoa_d_cruise*blend + w_aoa_d_cruise2*(1-blend)
batt_mass_accel = batt_mass_accel*blend + batt_mass_accel2*(1-blend)
batt_mass_climb = batt_mass_climb*blend + batt_mass_climb2*(1-blend)
batt_mass_cruise = batt_mass_cruise*blend + batt_mass_cruise2*(1-blend)

    OptParams.setvar!(:p_lam_t,x0,optxbounds,rangevar,p_lam_t)
    OptParams.setvar!(:p_chord,x0,optxbounds,rangevar,p_chord)
    OptParams.setvar!(:p_twist_d,x0,optxbounds,rangevar,p_twist_d)
    OptParams.setvar!(:p_pitch_d_accel,x0,optxbounds,rangevar,p_pitch_d_accel)
    OptParams.setvar!(:p_pitch_d_climb,x0,optxbounds,rangevar,p_pitch_d_climb)
    OptParams.setvar!(:p_pitch_d_cruise,x0,optxbounds,rangevar,p_pitch_d_cruise)
    OptParams.setvar!(:velocity_accel,x0,optxbounds,rangevar,velocity_accel)
    OptParams.setvar!(:velocity_climb,x0,optxbounds,rangevar,velocity_climb)
    OptParams.setvar!(:velocity_cruise,x0,optxbounds,rangevar,velocity_cruise)
    OptParams.setvar!(:p_Rtip,x0,optxbounds,rangevar,p_Rtip)
    OptParams.setvar!(:p_rpm_accel,x0,optxbounds,rangevar,p_rpm_accel)
    OptParams.setvar!(:p_rpm_climb,x0,optxbounds,rangevar,p_rpm_climb)
    OptParams.setvar!(:p_rpm_cruise,x0,optxbounds,rangevar,p_rpm_cruise)
    OptParams.setvar!(:kv,x0,optxbounds,rangevar,kv)
    OptParams.setvar!(:i0,x0,optxbounds,rangevar,i0)
    OptParams.setvar!(:w_aoa_d_accel,x0,optxbounds,rangevar,w_aoa_d_accel)
    OptParams.setvar!(:w_aoa_d_climb,x0,optxbounds,rangevar,w_aoa_d_climb)
    OptParams.setvar!(:w_aoa_d_cruise,x0,optxbounds,rangevar,w_aoa_d_cruise)
    OptParams.setvar!(:batt_mass_accel,x0,optxbounds,rangevar,batt_mass_accel)
    OptParams.setvar!(:batt_mass_climb,x0,optxbounds,rangevar,batt_mass_climb)
    OptParams.setvar!(:batt_mass_cruise,x0,optxbounds,rangevar,batt_mass_cruise)
    OptParams.setvar!(:number_blades,x0,optxbounds,rangevar,number_blades)

    OptParams.setvar!(:w_halfspan,x0,optxbounds,rangevar,w_halfspan)

    # sparce vars 900kg and ~50m/s vinf all
    # xwarm = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 3881.27, 2112.42, 434.783, 3999.49, 4499.16, 2783.33, 152.197, 15.5839, -1814.98, 526.502, 11.2334, 5151.27, 2163.88, 1201.38, 1201.38, 4685.41, 3446.64, 20000.0, 1.48825, 1.85602, 1.0, 2.97122, 3.99704, 14.1136, 168.285, 169.996, 160.905, 114.286, 0.0440764, 0.0809586, 0.0139054, 0.01, 0.01, 0.01, 0.312975, 1.00904, 0.207683, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 63.476]

    # x0[1:length(xwarm)] = xwarm[1:end]

    table_af_wing_blown = JLD.load("$path/airfoils/af_prop_E212.jld")
    table_3D_wing_blown = table_af_wing_blown["NDtable"]
    spline_3D_blown = AirfoilPrep.SplineND_from_tableND(table_3D_wing_blown)


    table_af_wing_free = JLD.load("$path/airfoils/af_prop_E212.jld")
    table_3D_wing_free = table_af_wing_free["NDtable"]
    spline_3D_freestream = AirfoilPrep.SplineND_from_tableND(table_3D_wing_free)

    # # Load the prop airfoil ND table
    # ClarkY_non_extrap = JLD.load("$(path)/airfoils/af_prop_E212.jld")
    # ClarkY_non_extrap = ClarkY_non_extrap["NDtable"]
    # AR = 1.0
    # omega = p_rpm_climb*pi/30
    # TSR = omega * p_Rtip / velocity_cruise
    #
    # grid_alphas=[i for i in -180:1.0:180]
    # Rtip = OptParams.getvar(:p_Rtip,x0,optxbounds,rangevar)
    # chord_prop_pt = OptParams.getvar(:p_chord,x0,optxbounds,rangevar)
    # rpt = linspace(.2,.911,length(chord_prop_pt)) * Rtip
    # r = linspace(.2,.911,p_n_akima) * Rtip
    # chord_prop = Akima.interp_noderiv(rpt, chord_prop_pt, r)
    # chord_prop[chord_prop.<1E-3] = 1E-3
    # af = Array{Any}(p_n_akima)
    # for i = 1:p_n_akima
    #
    #    r_over_R = r[i]
    #    c_over_r = chord_prop[i]/r[i]
    #    NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(ClarkY_non_extrap,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)
    #    afspl = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)
    #    af[i] = afspl
    # end
    # JLD.save("$(path)/airfoils/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))af_rerun.jld", "af", af)
    af = JLD.load("$(path)/airfoils/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))af_rerun.jld")
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



    func(x) = objective(x,optxbounds,rangevar,designparams,designvars,constraints,
    ctlparams,af,p_xaf,p_yaf,p_xafstrain,p_yafstrain,w_xaf,w_yaf,w_xafstrain,
    w_yafstrain,spline_3D_freestream,spline_3D_blown;detailed=detailed,
    savefinalvtk=savefinalvtk, verification=verification, plots=plots)



    function fungrad(x)
        function func1(x)
            j,C,_ = func(x)
            return [j;C]
        end
        f, dfdx = Gradients.forwarddiff(func1,x)
        lenx = length(x)

        return f[1],f[2:end],dfdx[1,:],dfdx[2:end,:],false
    end

    function pfungrad(x)
        function func1(x)
            j,C,_ = func(x)
            return [j;C]
        end
        f, dfdx = Gradients.centraldiff(func1,x;h=1E-5,parallel = true)

        return f[1],f[2:end],dfdx[1,:],dfdx[2:end,:],false
    end


    # ----- Define Optimizer Options ----- #
    options = Dict{String, Any}()
    options["Derivative level"] = 0
    options["Verify level"] = 0
    # options["Function precision"] = 1.00E-4
    # options["Difference interval"] = 1e-4
    # options["Central difference interval"] = 1e-4
    options["Iterations limit"] = 1e8
    options["Major iterations limit"] = 1000
    # options["Minor iterations limit"]= 1e8
    options["Major optimality tolerance"] = 1e-8 #Should be scaled so it is optimal with a solid 2 digits
    # options["Minor optimality  tolerance"] = 1e-6
    options["Major feasibility tolerance"] = 5e-6
    # options["Minor feasibility tolerance"] = 1e-6
    # options["Minor print level"] = 1E8
    # options["Print frequency"] = 1E8
    options["Scale option"] = 1
    # options["Scale tolerance"] = .95
    # options["Verify level"] = 3 #only if specifying gradients

    printfile = "./SNOPT_out/SNOPT_printfile_p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise)).out"
    sumfile = "./SNOPT_out/SNOPT_sumfile_p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise)).out"
    time0 = 0.0
    #if true
        verification=false
        savefinalvtk=false
        detailed = true
        plots = false

        detailed_output_accel,detailed_output_climb,detailed_output_cruise,
        design_accel,design_climb,design_cruise, detailed_system, f, c = objective(x0,optxbounds,rangevar,
        designparams,designvars,constraints,ctlparams,af,p_xaf,p_yaf,p_xafstrain,
        p_yafstrain,w_xaf,w_yaf,w_xafstrain,w_yafstrain,spline_3D_freestream,spline_3D_blown;verification=verification,
        savefinalvtk=savefinalvtk,
        detailed =detailed,
        plots = plots)

        # tic()
        # finaloutput,p_geom,p_struct, f, c = objective(x0,optxbounds,rangevar,
        # designparams,designvars,constraints,ctlparams,af,p_xaf,p_yaf,p_xafstrain,
        # p_yafstrain,w_xaf,w_yaf,w_xafstrain,w_yafstrain,spline_3D_freestream,spline_3D_blown;verification=verification,
        # savefinalvtk=savefinalvtk,
        # detailed =detailed,
        # plots = plots)
        # time0 = toc()

        JLD.save("$(path)/OptSummaries/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))correct_span.jld",
        "detailed_output_accel", detailed_output_accel,
        "detailed_output_climb", detailed_output_climb,
        "detailed_output_cruise", detailed_output_cruise,
        "design_accel", design_accel,
        "design_climb", design_climb,
        "design_cruise", design_cruise,
        "detailed_system", detailed_system)
    #end

    dfdx1 = 0.0
    time1 = 0.0
    nameout = "./SNOPT_out/printout.txt"
    if false#true
        verification=false
        savefinalvtk=false
        detailed = false
        plots = false

        f1,c1,dfdx1,dcdx1,state = fungrad(x0)
        tic()
        f1,c1,dfdx1,dcdx1,state = fungrad(x0)
        time1 = toc()
    end

    if false#true
        verification=false
        savefinalvtk=false
        detailed = false
        plots = false

        # f2,c2,dfdx2,dcdx2,state = pfungrad(x0)
        tic()
        f2,c2,dfdx2,dcdx2,state = pfungrad(x0)
        time2 = toc()
        println("x0: $x0")
        println("f2: $(f2)")
        println("c2: $(c2)")
        println("dfdx2: $(dfdx2)")
        for i = 1:length(c2)
            if maximum(abs.(dcdx2[i,:]))>30.0 # || mean(abs.(dcdx2[i,:]))<1E-4
                println("dc$i:
                $(dcdx2[i,:])")
            end
        end
        println(dfdx1-dfdx2)

        open(nameout,"a")do x
            write(x,"x0: $x0
            f2: $f2
            dfdx2: $(dfdx2)
            $(dfdx1-dfdx2)
            nproc: $(nprocs())
            time0: $time0
            time1: $time1
            time2: $time2
            ")
        end


    end

# x0 = [0.207583, 0.16042, 0.0980674, 0.0341654, 0.0347087, 0.0203864, 0.0538978, 0.0121262, 3.15971, 2.49569, 1.83168, 0.66743, 3.73101, 3.05527, 2.26264, 1.56159, 406.152, -7.77481, -113.17, 50.0, 0.13628, 29.4291, 20.0, 88.1265, 8814.21, 2.40003, 26.5281, 2.87881, 0.4, 28.236, 3.47529, 0.001, 0.001, 0.5]

    optimize = false
    if optimize

        verification=false
        savefinalvtk=false
        detailed = false
        plots = false
        if nprocs()==1
            optfunc = func
        else
            optfunc = pfungrad
        end

        xopt, fopt, exitmsg = Snopt.snopt(optfunc, x0, lb, ub, options,
        printfile = printfile, sumfile = sumfile)


        println(exitmsg)
        println("fopt: $fopt")
        println("xopt: $xopt")

        verification=false
        savefinalvtk=true
        detailed = true
        plots = false

        detailed_output_accel,detailed_output_climb,detailed_output_cruise,
        design_accel,design_climb,design_cruise, detailed_system, f, c = objective(xopt,optxbounds,rangevar,
        designparams,designvars,constraints,ctlparams,af,p_xaf,p_yaf,p_xafstrain,
        p_yafstrain,w_xaf,w_yaf,w_xafstrain,w_yafstrain,spline_3D_freestream,spline_3D_blown;detailed=detailed, savefinalvtk=savefinalvtk, verification=verification, plots=plots)

        JLD.save("$(path)/OptSummaries/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))correct_span.jld",
        "detailed_output_accel", detailed_output_accel,
        "detailed_output_climb", detailed_output_climb,
        "detailed_output_cruise", detailed_output_cruise,
        "design_accel", design_accel,
        "design_climb", design_climb,
        "design_cruise", design_cruise,
        "detailed_system", detailed_system)

        println("----------Detailed Output----------")
        for field in fieldnames(detailed_system)
            # if field in keys(rangevar)
            if length(string(field)) < 7
                println(field,":\t\t\t",getfield(detailed_system,field))
            elseif length(string(field)) < 15
                println(field,":\t\t",getfield(detailed_system,field))
            else
                println(field,":\t",getfield(detailed_system,field))
            end
            # end
        end
    end

    C = []
    J = []
    Jin = []
    Cin = []
    variable = []
    # lb = x0*0.97
    # ub = x0*1.03

    if false #sweep # sweep through the variables

        verification=false
        savefinalvtk=false
        detailed = false
        plots = false


        N_points = 20 #discretization for sweeping
        for i = 1:length(x0)#[37] #[1,13,19,25,26,28,29,30,33,37,38] #[1,13,19,25,26,27,28,29,30,31,34,37,38]

            xsweep = copy(x0)
            J1,C1,fail1 = func(xsweep)
            variable = linspace(x0[i]*0.99,x0[i]*1.01,N_points)
            J = zeros(N_points,length(J1))
            C = zeros(N_points,length(C1))
            for j = 1:N_points
                xsweep[i] = variable[j]
                Jin,Cin,_ = func(xsweep)
                J[j,:] = Jin
                C[j,:] = Cin
                xsweep = copy(x0)
            end

            # figure("objective")
            # clf()
            # figure("contraints")
            # clf()

            figure("objective")
            PyPlot.plot(variable/maximum(abs.(variable)),J,"-",label = "J variable $i") #/maximum(abs(J))
            legend()

            # PyPlot.figure("contraints")
            # for k = 1:length(C[1,:])
            #
            #     # if !isapprox(norm(C[:,k]), norm(C[:,k+1]); atol = norm(C[:,k])*1E-1)
            #     # if (maximum(C[:,k]/maximum(abs.(C[:,k])))<=1.0 && minimum(C[:,k]/maximum(abs.(C[:,k])))>-.9) #&& minimum((C[:,k]./maximum(abs.(C[:,k]))))<-0.75
            #         PyPlot.plot(variable,C[:,k]/maximum(abs.(C[:,k])),label = "C$k var $i")
            #     # end
            #     # end
            # end
            # PyPlot.legend()
            pause(.100001)
        end
    end

    # end #optimize

    # end #module