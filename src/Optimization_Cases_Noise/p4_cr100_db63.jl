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

printiter = 1
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
    numprops = 4
    to_factor = 1
    db_factor = 1

    noise_factor = 63.0
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




            number_blades = 2.0

            #SLOW
                kv = 5.0
                i0 = 5.999996165805586

                p_lam_t = [0.0110025, 0.00233067, 0.000691286, 0.000641728, 0.000235849, 0.000218, 0.000221644, 0.00023394]
                p_chord = [0.152889, 0.189812, 0.0120205, 0.0125594]
                p_twist_d = [43.6009, 25.2406, 15.1695, 1.96521]
                p_pitch_d_accel = -1.8830676001098614
                p_pitch_d_climb = -6.651289051243336
                p_pitch_d_cruise = 3.1389545772338217
                velocity_accel = 10.516264138531302
                velocity_climb = 10.671170506667663
                velocity_cruise = 22.35123708929302
                p_Rtip = 1.4607960856634958
                p_rpm_accel = 1363.7927957326167
                p_rpm_climb = 1771.769821957893
                p_rpm_cruise = 785.090064714826
                number_blades = 2.0
                kv = 5.000126704265029
                i0 = 5.999753473854883
                w_aoa_d_accel = 26.5521851166395
                w_aoa_d_climb = 24.112851297141972
                w_aoa_d_cruise = 7.730874684275981
                batt_mass_accel = 0.01
                batt_mass_climb = 0.1549107266114388
                batt_mass_cruise = 42.91158400899825

                p_lam_t = [0.0116438, 0.00350416, 0.00404958, 0.000152, 0.000235849, 0.000218, 0.000218, 0.000233965]
                p_chord = [0.145164, 0.186175, 0.0685406, 0.010077]
                p_twist_d = [45.0, 26.8695, 21.4999, 22.2934]
                p_pitch_d_accel = -6.855551167943682
                p_pitch_d_climb = -12.667144982209622
                p_pitch_d_cruise = -0.8224628373467217
                velocity_accel = 10.67040554661978
                velocity_climb = 8.423250535771844
                p_Rtip = 1.3088151699357786
                p_rpm_accel = 1214.5450961829085
                p_rpm_climb = 1982.7872875183948
                p_rpm_cruise = 784.1909878114266
                kv = 5.0
                i0 = 6.0
                w_aoa_d_accel = 25.68969593406325
                w_aoa_d_climb = 26.322257801432617
                w_aoa_d_cruise = 7.805358279717483
                batt_mass_accel = 0.01
                batt_mass_climb = 0.706638078458127
                batt_mass_cruise = 44.12768373997392

                p_lam_t = [0.00961238, 0.00303307, 0.00430381, 0.00026858, 0.000235849, 0.000218, 0.000218, 0.000233965]
                p_chord = [0.145177, 0.186635, 0.0607687, 0.01]
                p_twist_d = [44.8803, 29.6608, 20.7856, 22.4135]
                p_pitch_d_accel = -6.701873561022049
                p_pitch_d_climb = -10.403335281957924
                p_pitch_d_cruise = 7.0036484354036705
                velocity_accel = 13.042183088287917
                velocity_climb = 21.78601550567518
                velocity_cruise = 44.69992450382849
                p_Rtip = 0.9794038430607419
                p_rpm_accel = 1258.7053098738095
                p_rpm_climb = 2464.3015238179573
                p_rpm_cruise = 1595.3065865367794
                w_aoa_d_accel = 25.68902157680819
                w_aoa_d_climb = 25.190975655524863
                w_aoa_d_cruise = 1.935788588268422
                batt_mass_climb = 12.281978088817898
                batt_mass_cruise = 110.34095803201556

                p_lam_t = [0.0100125, 0.00290661, 0.00490218, 0.000221752, 0.000235849, 0.000218, 0.000218, 0.000233965]
                p_chord = [0.145235, 0.186503, 0.0578565, 0.01]
                p_twist_d = [44.9288, 29.6574, 21.339, 22.4122]
                p_pitch_d_accel = -6.701390018772187
                p_pitch_d_climb = -9.756573609373271
                p_pitch_d_cruise = 7.003588042667996
                velocity_accel = 13.04175815965444
                velocity_climb = 12.420638496541775
                p_Rtip = 0.9782710111556879
                p_rpm_accel = 1258.9059632304622
                p_rpm_climb = 2518.2836370838063
                p_rpm_cruise = 1595.2947452820074
                kv = 5.000172553257532
                i0 = 6.0
                w_aoa_d_accel = 26.26830137134156
                w_aoa_d_climb = 25.702645930032226
                w_aoa_d_cruise = 1.933487483214442
                batt_mass_climb = 10.299573659674847
                batt_mass_cruise = 110.056240384561

                p_lam_t = [0.00981777, 0.00285966, 0.00456369, 0.00021496, 0.000235849, 0.000218, 0.000218, 0.000233965]
                p_chord = [0.145225, 0.186455, 0.0581635, 0.0101314]
                p_twist_d = [44.9398, 29.4376, 21.0365, 21.975]
                p_pitch_d_accel = -6.701359519499593
                p_pitch_d_climb = -9.778003389464745
                p_pitch_d_cruise = 7.003599347824445
                velocity_accel = 13.041731261532203
                velocity_climb = 9.552320864282073
                p_Rtip = 0.9782820488983894
                p_rpm_accel = 1258.9133261036197
                p_rpm_climb = 2517.7358944162934
                p_rpm_cruise = 1695.2973509664116
                number_blades = 2.082690127964101
                kv = 5.0
                i0 = 6.0
                w_aoa_d_accel = 27.12764594271704
                w_aoa_d_climb = 25.777505793785213
                w_aoa_d_cruise = 1.9364873496700061
                batt_mass_climb = 9.592819384903665
                batt_mass_cruise = 113.32237641452463

                number_blades = 2.0

                p_lam_t = [0.00497069, 0.00424201, 0.00527295, 0.000152, 0.000237324, 0.000218, 0.000218, 0.000752425]
                p_chord = [0.247114, 0.154571, 0.0575167, 0.05]
                p_twist_d = [45.0, 25.901, 18.6034, 16.7131]
                p_pitch_d_accel = -3.8533189527175358
                p_pitch_d_climb = -5.840289110324262
                p_pitch_d_cruise = -1.7864041045164314
                velocity_accel = 12.335741882949291
                velocity_climb = 11.1926238674717
                p_Rtip = 1.0277531687913324
                p_rpm_accel = 2377.592248196533
                p_rpm_climb = 2525.0248031305086
                p_rpm_cruise = 2129.2945710134904
                kv = 5.0
                i0 = 6.0
                w_aoa_d_accel = 33.83328207880271
                w_aoa_d_climb = 33.72722962344346
                w_aoa_d_cruise = 4.931345208684481
                batt_mass_accel = 0.001
                batt_mass_climb = 0.734068685217686
                batt_mass_cruise = 81.34301773985665

                p_lam_t = [0.00956401, 0.00198993, 0.004933, 0.000214831, 0.000292573, 0.00169026, 0.000954541, 0.000348048]
                p_chord = [0.157688, 0.162584, 0.0514416, 0.0503022]
                p_twist_d = [43.6212, 31.7104, 21.3084, 17.8097]
                p_pitch_d_accel = -7.837039738548938
                p_pitch_d_climb = -7.911860729864201
                p_pitch_d_cruise = -2.158856562329411
                velocity_accel = 13.786395569608965
                velocity_climb = 13.96677333240264
                p_Rtip = 0.9723762776008794
                p_rpm_accel = 2252.1998946840804
                p_rpm_climb = 2453.8872751359872
                p_rpm_cruise = 1992.28868303773
                kv = 5.039351052000106
                i0 = 5.964987710510177
                w_aoa_d_accel = 32.41155853831902
                w_aoa_d_climb = 28.5254701822914
                w_aoa_d_cruise = 5.12818786826107
                batt_mass_accel = 0.001
                batt_mass_climb = 0.5233667563276049
                batt_mass_cruise = 77.26709454794035

                p_lam_t = [0.00984977, 0.00308691, 0.00499031, 0.000152, 0.000986107, 0.000538009, 0.000652653, 0.000941326]
p_chord = [0.159518, 0.158507, 0.0510347, 0.0554556]
p_twist_d = [44.3889, 31.3959, 20.6764, 17.6741]
p_pitch_d_accel = -7.807252291719941
p_pitch_d_climb = -7.36771537532478
p_pitch_d_cruise = -2.159212862691299
velocity_accel = 13.712244632049552
velocity_climb = 13.15817582976719
p_Rtip = 0.9826920993377842
p_rpm_accel = 2256.732029683491
p_rpm_climb = 2424.240334716867
p_rpm_cruise = 1992.2412824198907
kv = 5.01442890998258
i0 = 5.979993488546374
w_aoa_d_accel = 32.313261651072175
w_aoa_d_climb = 29.608834394233764
w_aoa_d_cruise = 4.805903637871549
batt_mass_accel = 0.001
batt_mass_climb = 0.4249322546210109
batt_mass_cruise = 74.73979067638315

p_lam_t = [0.00924278, 0.00276194, 0.0045593, 0.000152, 0.000986107, 0.000538009, 0.000652653, 0.000941326]
p_chord = [0.161059, 0.147857, 0.05, 0.050001]
p_twist_d = [44.4805, 31.8139, 21.1869, 18.1616]
p_pitch_d_accel = -7.800340650460656
p_pitch_d_climb = -7.71267619456311
p_pitch_d_cruise = -2.1592483722596993
velocity_accel = 13.704885091769464
velocity_climb = 13.639890212631029
p_Rtip = 0.9735350164339817
p_rpm_accel = 2255.515518003605
p_rpm_climb = 2254.7324735133266
p_rpm_cruise = 1992.2345074810514
kv = 5.009912079852231
i0 = 5.986217387734277
w_aoa_d_accel = 32.055544637423544
w_aoa_d_climb = 29.630762928127933
w_aoa_d_cruise = 4.728930029756277
batt_mass_accel = 0.001
batt_mass_climb = 0.23783319950612603
batt_mass_cruise = 77.30902405630627

p_lam_t = [0.00423355, 0.00884946, 0.00364156, 0.00394408, 0.000986107, 0.000538009, 0.000652653, 0.000941326]
p_chord = [0.29722, 0.109427, 0.0364156, 0.0394408]
p_twist_d = [45.0, 27.5191, 17.9309, 15.565]
p_pitch_d_accel = -7.494284764030797
p_pitch_d_climb = -7.493047565452269
p_pitch_d_cruise = -2.162049426923569
velocity_accel = 12.431513034533541
velocity_climb = 12.20078748620543
p_Rtip = 1.1400000000000265
p_rpm_accel = 2276.4054881155234
p_rpm_climb = 2276.7236950785177
p_rpm_cruise = 1991.5633509596828
kv = 5.0
i0 = 5.905753652089577
w_aoa_d_accel = 33.98800763161594
w_aoa_d_climb = 32.51689588185419
w_aoa_d_cruise = 4.728006707891872
batt_mass_accel = 0.001
batt_mass_climb = 0.606607323441761
batt_mass_cruise = 75.47003413359957

p_lam_t = [0.00423286, 0.00885023, 0.00364046, 0.00394408, 0.000986107, 0.000538009, 0.000652653, 0.000941326]
p_chord = [0.29722, 0.109427, 0.0364046, 0.0394474]
p_twist_d = [45.0, 27.5201, 17.9346, 15.5517]
p_pitch_d_accel = -7.4943483586753885
p_pitch_d_climb = -7.49313034112774
p_pitch_d_cruise = -2.162049426923569
velocity_accel = 12.431920780034782
velocity_climb = 12.201715119635146
p_Rtip = 1.140000000000013
p_rpm_accel = 2276.4054881155234
p_rpm_climb = 2276.5604908668347
p_rpm_cruise = 1991.5633509596828
kv = 5.0
i0 = 5.905753652089577
w_aoa_d_accel = 33.98826120550194
w_aoa_d_climb = 32.51689588185419
w_aoa_d_cruise = 4.727756420463748
batt_mass_accel = 0.001
batt_mass_climb = 0.6065352965407658
batt_mass_cruise = 75.46805000652417

p_lam_t = [0.00423203, 0.00884696, 0.0036414, 0.00393914, 0.000986107, 0.000538009, 0.000652653, 0.000941326]
p_chord = [0.29722, 0.109408, 0.036414, 0.0393921]
p_twist_d = [45.0, 27.5174, 17.9359, 15.558]
p_pitch_d_accel = -7.495012742669088
p_pitch_d_climb = -7.494590124896626
p_pitch_d_cruise = -2.161758583776195
velocity_accel = 12.4323319064538
velocity_climb = 12.202048398616535
p_Rtip = 1.140000000000002
p_rpm_accel = 2276.350238795155
p_rpm_climb = 2276.423874641913
p_rpm_cruise = 1991.5633509596828
kv = 5.0
i0 = 5.905753652089577
w_aoa_d_accel = 33.988617220416636
w_aoa_d_climb = 32.519560344720595
w_aoa_d_cruise = 4.7274920717181566
batt_mass_accel = 0.001
batt_mass_climb = 0.606400687886306
batt_mass_cruise = 75.46741562989975

# fixed aoa

p_lam_t = [0.00432476, 0.00912701, 0.00366804, 0.00394804, 0.000986107, 0.000538009, 0.000652653, 0.000752223]
p_chord = [0.279417, 0.10735, 0.0366804, 0.0395359]
p_twist_d = [45.0, 27.4279, 17.9776, 15.4204]
p_pitch_d_accel = -7.761700785799464
p_pitch_d_climb = -7.217733978199294
p_pitch_d_cruise = -2.154636186710618
velocity_accel = 10.822550333849438
velocity_climb = 13.86311711292353
p_Rtip = 1.140000000000011
p_rpm_accel = 2276.405488114341
p_rpm_climb = 2276.4054881211737
p_rpm_cruise = 1991.984974491003
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 33.85890923013074
w_aoa_d_climb = 34.09214033964515
w_aoa_d_cruise = 4.701087227707984
batt_mass_accel = 0.001
batt_mass_climb = 0.5375184321032225
batt_mass_cruise = 75.43855086316978

#Fixed Structure
p_lam_t = [0.00347749, 0.00592931, 0.00354566, 0.000260414, 0.000986107, 0.000538009, 0.000652653, 0.000752223]
p_chord = [0.266859, 0.107527, 0.0354566, 0.0473888]
p_twist_d = [45.0, 27.4194, 17.5924, 15.7654]
p_pitch_d_accel = -7.195041126976776
p_pitch_d_climb = -6.53693397809769
p_pitch_d_cruise = -2.158526126696043
velocity_accel = 11.12211949591222
velocity_climb = 14.172919150981603
p_Rtip = 1.1400000000000003
p_rpm_accel = 2276.4054881156057
p_rpm_climb = 2276.405488115606
p_rpm_cruise = 1991.447591388328
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 33.51439470678863
w_aoa_d_climb = 33.83069866305622
w_aoa_d_cruise = 4.725901259614525
batt_mass_accel = 0.001
batt_mass_climb = 0.5650079840894046
batt_mass_cruise = 75.36033819643158

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
    af = JLD.load("$(path)/airfoils/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db65af.jld")
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
    options["Major feasibility tolerance"] = 5e-7
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
            time2: $time2c
            ")
        end


    end



    optimize = true
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