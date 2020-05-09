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

printiter = 0
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
    printfreq = 50,
    savevtk = 5000,
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
    :velocity_cruise,
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
    cruise_factor = 50.0
    required_cruise = cruise_factor * 0.447 #MPH TO m/s


    # designvars, constraints, ctlparams = getoptparams()
    numprops = 16
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



            number_blades = 5.0


            #Slow
                kv = 5.0
                i0 = 2.2037300540130143

                p_lam_t = [0.0007782, 0.000674115, 0.00050298, 0.000238061, 0.000218, 0.000218, 0.000218, 0.000218]
                p_chord = [0.217695, 0.179515, 0.0962492, 0.0129004]
                p_twist_d = [39.721, 18.0451, 10.7105, 4.20071]
                p_pitch_d_accel = 9.021458633187988
                p_pitch_d_climb = 2.6779561382353223
                p_pitch_d_cruise = 9.955675862211843
                velocity_accel = 21.907192616714354
                velocity_climb = 15.022172822158165
                velocity_cruise = 26.141124740636286
                p_Rtip = 0.1470747752685889
                p_rpm_accel = 4550.787583268113
                p_rpm_climb = 3617.25555613192
                p_rpm_cruise = 4206.390824864598
                w_aoa_d_accel = 11.805258303183681
                w_aoa_d_climb = 12.338223392213303
                w_aoa_d_cruise = 7.444036776377783
                batt_mass_accel = 0.01
                batt_mass_climb = 1.78807293443778
                batt_mass_cruise = 48.4944436895349

                number_blades = 2.0

                #SLOW
                    kv = 5.0
                    i0 = 5.999996165805586

                p_lam_t = [0.0106565, 0.00233204, 0.000521854, 0.000643721, 0.000235849, 0.000218, 0.000221644, 0.00023394]
    p_chord = [0.152893, 0.189842, 0.0122381, 0.0126047]
    p_twist_d = [43.6566, 25.3336, 14.8314, 2.15019]
    p_pitch_d_accel = -1.8829344927066909
    p_pitch_d_climb = -6.623043953201483
    p_pitch_d_cruise = 3.1389505191969
    velocity_accel = 10.516263781338967
    velocity_climb = 22.110235464233497
    velocity_cruise = 22.34991575727592
    p_Rtip = 1.4607155806562915
    p_rpm_accel = 1363.801616202054
    p_rpm_climb = 1762.6245635737948
    p_rpm_cruise = 785.088616955683
    number_blades = 2.018354828584329
    kv = 5.0000000099351745
    i0 = 5.999999998436644
    w_aoa_d_accel = 25.7671847059628
    w_aoa_d_climb = 24.085518334858087
    w_aoa_d_cruise = 7.6874985270624485
    batt_mass_accel = 0.01
    batt_mass_climb = 1.2565648978309747
    batt_mass_cruise = 43.031131277642885

    p_lam_t = [0.00305852, 0.00137265, 0.000290901, 0.000996405, 0.000235849, 0.000218, 0.00022111, 0.000230092]
p_chord = [0.175304, 0.162308, 0.010835, 0.010002]
p_twist_d = [31.4144, 25.9963, 20.3517, 9.77923]
p_pitch_d_accel = -7.719228969357121
p_pitch_d_climb = -13.902685653464765
p_pitch_d_cruise = 1.1738434828282451
velocity_accel = 10.178092602026814
velocity_climb = 10.382760710177195
p_Rtip = 1.3667788139004986
p_rpm_accel = 1145.7358822436008
p_rpm_climb = 1666.0538848152646
p_rpm_cruise = 788.3481791457784
number_blades = 2.6434199441191963
kv = 6.103139971395562
i0 = 4.971025979554923
w_aoa_d_accel = 26.06240250515709
w_aoa_d_climb = 23.488736431800533
w_aoa_d_cruise = 7.818359355419493
batt_mass_accel = 0.01
batt_mass_climb = 0.5628561589513709
batt_mass_cruise = 100.64018900946813

p_lam_t = [0.0015526, 0.000152, 0.000471676, 0.000944611, 0.000235849, 0.000218, 0.000218, 0.000230092]
p_chord = [0.17515, 0.199292, 0.0229421, 0.01]
p_twist_d = [45.0, 42.4987, 36.3488, 11.1557]
p_pitch_d_accel = -7.724322717198662
p_pitch_d_climb = -23.676212247437675
p_pitch_d_cruise = 1.1959698743438472
velocity_accel = 10.777470512914551
velocity_climb = 13.317126291886758
p_Rtip = 0.7152893871652168
p_rpm_accel = 1193.2347473339564
p_rpm_climb = 2372.4328660875353
p_rpm_cruise = 789.8170769898898
number_blades = 4.94638522287553
kv = 6.103139971395562
i0 = 4.971025979554923
w_aoa_d_accel = 23.9786588926515
w_aoa_d_climb = 17.182660261944864
w_aoa_d_cruise = 7.735694908900511
batt_mass_accel = 0.01
batt_mass_climb = 0.06797890138891173
batt_mass_cruise = 54.462714647113415

number_blades = 5.0

p_lam_t = [0.00154862, 0.000152, 0.000471676, 0.000944611, 0.000235849, 0.000218, 0.000218, 0.000230092]
p_chord = [0.17515, 0.199226, 0.0230771, 0.01]
p_twist_d = [45.0, 42.5124, 35.8923, 11.2191]
p_pitch_d_accel = -7.724316126556464
p_pitch_d_climb = -23.698007219675056
p_pitch_d_cruise = 1.1959717286527558
velocity_accel = 10.77745929922272
velocity_climb = 12.543799230008656
p_Rtip = 0.7152904095770524
p_rpm_accel = 1193.3153054871473
p_rpm_climb = 2371.090894423183
p_rpm_cruise = 789.8170769898898
kv = 6.103139971395562
i0 = 4.971025979554923
w_aoa_d_accel = 23.972135997530387
w_aoa_d_climb = 17.18488697226424
w_aoa_d_cruise = 7.735829528475066
batt_mass_accel = 0.01
batt_mass_climb = 0.19926503808056195
batt_mass_cruise = 54.28437606222326


p_lam_t = [0.000152, 0.000152, 0.000152, 0.00016314, 0.000218, 0.000218, 0.000697237, 0.000218]
p_chord = [0.159116, 0.0881603, 0.0145218, 0.010409]
p_twist_d = [45.0, 41.3985, 28.554, 13.948]
p_pitch_d_accel = -6.702376950041459
p_pitch_d_climb = -0.9984657091922631
p_pitch_d_cruise = -2.9083478654690595
velocity_accel = 13.309245111241173
velocity_climb = 19.398578866562552
p_Rtip = 0.7152941176470594
p_rpm_accel = 1665.7033096670643
p_rpm_climb = 2761.0476758587279
p_rpm_cruise = 2030.5534516552439
kv = 8.1779358810975
i0 = 5.992269337186693
w_aoa_d_accel = 16.168181853486463
w_aoa_d_climb = 12.909173218251121
w_aoa_d_cruise = 1.9182052665194256
batt_mass_accel = 0.09288778086913847
batt_mass_climb = 0.7261166000961358
batt_mass_cruise = 124.77395766068845

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000163118, 0.000218, 0.000218, 0.000218, 0.000913548]
p_chord = [0.28597, 0.0885522, 0.05, 0.05]
p_twist_d = [45.0, 30.7086, 24.3697, 19.4367]
p_pitch_d_accel = -3.774924630424246
p_pitch_d_climb = -10.2393668560264
p_pitch_d_cruise = -1.3475341076979703
velocity_accel = 8.764174986241933
velocity_climb = 6.320916505351116
p_Rtip = 0.33529411764745026
p_rpm_accel = 5099.461390140229
p_rpm_climb = 7739.750204244843
p_rpm_cruise = 3028.0037195241225
kv = 10.841867649681772
i0 = 2.8460209616810417
w_aoa_d_accel = 31.03721205499132
w_aoa_d_climb = 26.311544050740377
w_aoa_d_cruise = 16.071648067942363
batt_mass_accel = 0.001
batt_mass_climb = 2.102751401285317
batt_mass_cruise = 96.21276517636045

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000163118, 0.000218, 0.000218, 0.000218, 0.000913548]
p_chord = [0.285924, 0.0885824, 0.05, 0.05]
p_twist_d = [45.0, 30.7137, 24.3516, 19.4907]
p_pitch_d_accel = -3.773989753177802
p_pitch_d_climb = -10.235233083425914
p_pitch_d_cruise = -1.3475341076979703
velocity_accel = 8.765340281260269
velocity_climb = 6.320092951780466
velocity_cruise = 22.35352481498706
p_Rtip = 0.33529411764705896
p_rpm_accel = 5099.8281374735725
p_rpm_climb = 7739.778659594037
p_rpm_cruise = 3028.0037195241225
number_blades = 4.997574103806742
kv = 10.84009771365462
i0 = 2.851496269925953
w_aoa_d_accel = 31.039082546571656
w_aoa_d_climb = 26.312602094893933
w_aoa_d_cruise = 16.07181091676642
batt_mass_accel = 0.001
batt_mass_climb = 2.1043603794491745
batt_mass_cruise = 96.19651969264403

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.00185717, 0.00191666, 0.00112342, 0.000314455]
p_chord = [0.152223, 0.109652, 0.0569723, 0.00359267]
p_twist_d = [35.3435, 30.8014, 26.2829, 33.1773]
p_pitch_d_accel = -5.774189754091159
p_pitch_d_climb = -8.43706605322609
p_pitch_d_cruise = 0.3927881340066803
velocity_accel = 7.86858511497682
velocity_climb = 7.494701447823569
velocity_cruise = 27.405722633936534
p_Rtip = 0.33447319165810563
p_rpm_accel = 7099.906676421568
p_rpm_climb = 7750.822325660634
p_rpm_cruise = 3061.887233770996
number_blades = 3.302512307166864
kv = 17.251245104090344
i0 = 4.864222401894478
w_aoa_d_accel = 29.026820159111328
w_aoa_d_climb = 21.90158858674213
w_aoa_d_cruise = 11.993998118569932
batt_mass_accel = 0.001
batt_mass_climb = 4.669035164679397
batt_mass_cruise = 74.6940061340688

#Fixed aoa

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.00162832, 0.00113173, 0.000218, 0.000991903]
p_chord = [0.130139, 0.0949486, 0.0493457, 0.0267967]
p_twist_d = [43.3092, 31.4552, 26.1137, 22.8278]
p_pitch_d_accel = -6.911559376171299
p_pitch_d_climb = -8.141674283707156
p_pitch_d_cruise = 0.408273469365135
velocity_accel = 5.759738766141724
velocity_climb = 8.967515560924905
velocity_cruise = 27.350858929169924
p_Rtip = 0.33529411764705885
p_rpm_accel = 7295.248491277464
p_rpm_climb = 7739.778659593088
p_rpm_cruise = 3063.381408240985
kv = 16.795228903846805
i0 = 4.938530749637611
w_aoa_d_accel = 27.92491191178969
w_aoa_d_climb = 30.38775292459795
w_aoa_d_cruise = 11.898797362085341
batt_mass_accel = 0.001
batt_mass_climb = 1.6426429748046207
batt_mass_cruise = 74.25527832166281

#Fixed Structure
p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.00129098, 0.000885536, 0.000218, 0.000991905]
p_chord = [0.107865, 0.0790412, 0.0466251, 0.0346808]
p_twist_d = [45.0, 34.7346, 28.405, 25.654]
p_pitch_d_accel = -9.065587554037227
p_pitch_d_climb = -10.141359701799603
p_pitch_d_cruise = 0.42165248419260576
velocity_accel = 5.7053095839546195
velocity_climb = 9.000641137406797
velocity_cruise = 31.634482005986875
p_Rtip = 0.33529411764705896
p_rpm_accel = 7352.179191083504
p_rpm_climb = 7739.778659593059
p_rpm_cruise = 3063.381408240985
kv = 16.795228903846805
i0 = 4.9389283421031855
w_aoa_d_accel = 27.857277084663412
w_aoa_d_climb = 30.390720953090746
w_aoa_d_cruise = 9.16183064450044
batt_mass_accel = 0.001
batt_mass_climb = 1.6450662302354349
batt_mass_cruise = 71.19518902540432

number_blades = 4.0
w_halfspan = 5.7






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
    # JLD.save("$(path)/airfoils/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))af.jld", "af", af)
    af = JLD.load("$(path)/airfoils/p8_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))af_rerun.jld")
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
    options["Major optimality tolerance"] = 1e-6 #Should be scaled so it is optimal with a solid 2 digits
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
    if true #true
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
    end

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
