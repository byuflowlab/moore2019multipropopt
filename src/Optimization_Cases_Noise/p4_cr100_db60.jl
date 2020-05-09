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

    noise_factor = 60.0
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

                p_lam_t = [0.00386092, 0.00196405, 0.000176076, 0.0023944, 0.000237267, 0.000218, 0.000218, 0.000218]
                p_chord = [0.267913, 0.499966, 0.308411, 0.101915]
                p_twist_d = [14.9921, 14.9227, 20.5563, 18.1262]
                p_pitch_d_accel = -6.517916551276772
                p_pitch_d_climb = -7.974759470306106
                p_pitch_d_cruise = -2.3493018816611517
                velocity_accel = 14.873982337237894
                velocity_climb = 22.868299430199304
                p_Rtip = 0.7564430537650972
                p_rpm_accel = 1500.5612753892744
                p_rpm_climb = 2548.754972022308
                p_rpm_cruise = 1832.7804831698604
                kv = 5.0
                i0 = 6.0
                w_aoa_d_accel = 32.28462151540162
                w_aoa_d_climb = 20.06398120626823
                w_aoa_d_cruise = 4.695911014699417
                batt_mass_accel = 0.0011340365402889713
                batt_mass_climb = 0.62563497932973
                batt_mass_cruise = 83.83622510408681

                p_lam_t = [0.00261936, 0.00129239, 0.000318593, 0.00170887, 0.000230147, 0.000218001, 0.000218001, 0.000218001]
                p_chord = [0.288046, 0.425305, 0.28, 0.0898071]
                p_twist_d = [26.0897, 25.9517, 22.8351, 20.4363]
                p_pitch_d_accel = -6.647026285638015
                p_pitch_d_climb = -10.672558805274194
                p_pitch_d_cruise = -2.300927070925263
                velocity_accel = 14.51895130793713
                velocity_climb = 18.220358135886958
                p_Rtip = 0.6840352458199409
                p_rpm_accel = 1531.5879143108375
                p_rpm_climb = 2384.0128503597116
                p_rpm_cruise = 1832.2465482579948
                kv = 14.581424458602998
                i0 = 4.382129707990689
                w_aoa_d_accel = 28.33674283078558
                w_aoa_d_climb = 23.297484812179277
                w_aoa_d_cruise = 4.055398893068711
                batt_mass_accel = 0.0010842396020939447
                batt_mass_climb = 0.3969214237930061
batt_mass_cruise = 63.44944501191863

p_lam_t = [0.0026193, 0.00129236, 0.000318715, 0.00170883, 0.000230147, 0.000218001, 0.000218001, 0.000218001]
p_chord = [0.28804, 0.425295, 0.279994, 0.0898061]
p_twist_d = [26.0902, 25.9522, 22.8355, 20.4358]
p_pitch_d_accel = -6.647319782389796
p_pitch_d_climb = -10.672598280600084
p_pitch_d_cruise = 5.3009124014705638
velocity_accel = 14.518855091500484
velocity_climb = 18.22037699212219
p_Rtip = 0.6840272603827366
p_rpm_accel = 1531.574598032766
p_rpm_climb = 2384.014703171727
p_rpm_cruise = 1832.2471717828264
kv = 14.581758193276482
i0 = 4.3820180208445265
w_aoa_d_accel = 28.33605591192824
w_aoa_d_climb = 23.297394080617337
w_aoa_d_cruise = 4.055298335979353
batt_mass_accel = 0.0010842374049440699
batt_mass_climb = 0.39691133205068174
batt_mass_cruise = 63.4479205247701

p_lam_t = [0.00261924, 0.00129233, 0.000319103, 0.00170923, 0.000230147, 0.000218001, 0.000218001, 0.000218001]
p_chord = [0.288041, 0.425284, 0.28001, 0.0898167]
p_twist_d = [26.0907, 25.9527, 22.8356, 20.4354]
p_pitch_d_accel = -6.647349151428637
p_pitch_d_climb = -10.672590474284435
p_pitch_d_cruise = 5.300915909257823
velocity_accel = 14.518793855567504
velocity_climb = 18.22029018360142
p_Rtip = 0.6840255986089138
p_rpm_accel = 1531.565703917874
p_rpm_climb = 2384.0045495425175
p_rpm_cruise = 1832.247273538625
kv = 14.581509815608378
i0 = 4.382059962060628
w_aoa_d_accel = 28.336358263874
w_aoa_d_climb = 23.29739701783215
w_aoa_d_cruise = 4.0552906124872115
batt_mass_accel = 0.0010842352213479296
batt_mass_climb = 0.3969013025620974
batt_mass_cruise = 63.44640544184767

p_lam_t = [0.00261924, 0.00129233, 0.000319103, 0.00170923, 0.000233097, 0.000218001, 0.000218001, 0.000218001]
p_chord = [0.388041, 0.525284, 0.48001, 0.1898167]
p_twist_d = [26.0907, 25.9527, 22.8356, 20.4354]
p_pitch_d_accel = -6.647349151428637
p_pitch_d_climb = -10.672590474284435
p_pitch_d_cruise = 5.300915909257823
velocity_accel = 14.518793855567504
velocity_climb = 18.22029018360142
p_Rtip = 0.6840255986089138
p_rpm_accel = 2331.565703917874
p_rpm_climb = 2484.0045495425175
p_rpm_cruise = 1932.247273538625
kv = 14.581509815608378
i0 = 4.382059962060628
w_aoa_d_accel = 28.336358263874
w_aoa_d_climb = 23.29739701783215
w_aoa_d_cruise = 4.0552906124872115
batt_mass_accel = 0.0010842352213479296
batt_mass_climb = 0.3969013025620974
batt_mass_cruise = 63.44640544184767

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

p_lam_t = [0.000998319, 0.00607591, 0.00316333, 0.000787891, 0.00422525, 0.000741222, 0.00264572, 0.000364267]
p_chord = [0.519416, 0.116888, 0.0329347, 0.00918495]
p_twist_d = [44.2252, 29.4, 18.7897, 15.5361]
p_pitch_d_accel = -9.05236747554496
p_pitch_d_climb = -4.464432227744045
p_pitch_d_cruise = -2.169305055183001
velocity_accel = 13.946310159309812
velocity_climb = 22.972118069152756
p_Rtip = 1.1369896804464492
p_rpm_accel = 1851.1217311432706
p_rpm_climb = 1862.767079703493
p_rpm_cruise = 1993.4422608180894
kv = 5.012895298618748
i0 = 5.831467036499793
w_aoa_d_accel = 34.02958953207343
w_aoa_d_climb = 16.97913824557014
w_aoa_d_cruise = 4.605608183381458
batt_mass_accel = 0.010586193353069576
batt_mass_climb = 0.5482054063403234
batt_mass_cruise = 84.69366087101847

p_lam_t = [0.000152, 0.00821835, 0.00328633, 0.0015566, 0.00292193, 0.00290789, 0.00328633, 0.00241586]
p_chord = [0.55662, 0.0821835, 0.0328633, 0.0241262]
p_twist_d = [41.6739, 27.0823, 18.9412, 16.9996]
p_pitch_d_accel = -8.766515269869545
p_pitch_d_climb = -8.2786021673207
p_pitch_d_cruise = -2.1716476477527604
velocity_accel = 13.726644092060685
velocity_climb = 15.310710547534425
p_Rtip = 1.14
p_rpm_accel = 2060.817255348445
p_rpm_climb = 2076.718047753872
p_rpm_cruise = 1993.1660718678183
kv = 5.029702530837963
i0 = 5.916558236678992
w_aoa_d_accel = 32.24098574009515
w_aoa_d_climb = 26.728947243247433
w_aoa_d_cruise = 4.5492581301163675
batt_mass_accel = 0.001
batt_mass_climb = 0.01
batt_mass_cruise = 72.87035180738101

#Fixed aoa

p_lam_t = [0.000152, 0.00795315, 0.00314186, 0.00165438, 0.00257911, 0.00314105, 0.00314186, 0.00165437]
p_chord = [0.788397, 0.0795315, 0.0314186, 0.0165439]
p_twist_d = [45.0, 27.3359, 19.1007, 16.7498]
p_pitch_d_accel = -10.015406076262057
p_pitch_d_climb = -9.393896617570487
p_pitch_d_cruise = -2.1588512496592607
velocity_accel = 11.878192124860647
velocity_climb = 13.761110874058451
p_Rtip = 1.1399998385049939
p_rpm_accel = 2174.1528445820027
p_rpm_climb = 2168.6516335200963
p_rpm_cruise = 1995.041296665383
kv = 5.000007867955921
i0 = 5.449966668568358
w_aoa_d_accel = 34.150746625689685
w_aoa_d_climb = 33.24200762093198
w_aoa_d_cruise = 4.51127352104008
batt_mass_accel = 0.001
batt_mass_climb = 0.11740129813400288
batt_mass_cruise = 72.80029293805853

#Fixed Structure
p_lam_t = [0.000152, 0.00760407, 0.00281375, 0.00169663, 0.0024269, 0.000440446, 0.00281375, 0.00217664]
p_chord = [0.796986, 0.0785016, 0.0281375, 0.0276113]
p_twist_d = [45.0, 27.2802, 18.5533, 18.2285]
p_pitch_d_accel = -9.555860129609712
p_pitch_d_climb = -8.988769703695597
p_pitch_d_cruise = -2.1618306385408936
velocity_accel = 12.041891953629294
velocity_climb = 13.777939547958114
p_Rtip = 1.1400000000000001
p_rpm_accel = 2230.742455288636
p_rpm_climb = 2215.650594741242
p_rpm_cruise = 1994.6847506931902
kv = 5.0
i0 = 5.687202853797737
w_aoa_d_accel = 33.86443345276192
w_aoa_d_climb = 33.583113043962626
w_aoa_d_cruise = 4.524605600810268
batt_mass_accel = 0.001
batt_mass_climb = 0.5397311442102332
batt_mass_cruise = 72.24920430233269

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
    af = JLD.load("$(path)/airfoils/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))af.jld")
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
