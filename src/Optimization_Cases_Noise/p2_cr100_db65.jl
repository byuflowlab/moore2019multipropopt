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
    # :batt_mass_accel,
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
    numprops = 2
    to_factor = 1
    db_factor = 1

    noise_factor = 65.0
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

                p_lam_t = [0.0144977, 0.0041207, 0.00380252, 0.000152, 0.000235849, 0.000218, 0.000221649, 0.000267198]
                p_chord = [0.150925, 0.259938, 0.0547936, 0.01]
                p_twist_d = [40.2171, 24.9573, 15.4371, 12.2059]
                p_pitch_d_accel = -0.757024268372988
                p_pitch_d_climb = -4.227096682620433
                p_pitch_d_cruise = 3.131668618865019
                velocity_accel = 10.494871737917872
                velocity_climb = 10.038471203817535
                velocity_cruise = 22.350000000000406
                p_Rtip = 1.5016953607341559
                p_rpm_accel = 1451.8867462038745
                p_rpm_climb = 1696.18856249119
                p_rpm_cruise = 784.9979435271417
                w_aoa_d_accel = 27.117721582410013
                w_aoa_d_climb = 27.13183207420816
                w_aoa_d_cruise = 7.762172889745735
                batt_mass_accel = 0.01
                batt_mass_climb = 0.7502820899030893
                batt_mass_cruise = 44.26222220520083

                p_lam_t = [0.0144977, 0.00428491, 0.00413361, 0.000152, 0.000235849, 0.000218, 0.000221649, 0.000267198]
p_chord = [0.150926, 0.25996, 0.0561959, 0.0106229]
p_twist_d = [40.3493, 23.6495, 16.5499, 12.837]
p_pitch_d_accel = -0.7570366197630303
p_pitch_d_climb = -4.273533934261378
p_pitch_d_cruise = 3.1318767648056327
velocity_accel = 10.50115357735056
velocity_climb = 10.508385065847797
p_Rtip = 1.5017206693175817
p_rpm_accel = 1451.8853062368696
p_rpm_climb = 1702.004862847982
p_rpm_cruise = 784.997118577869
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 27.096542960779846
w_aoa_d_climb = 27.132800575422046
w_aoa_d_cruise = 7.770341903091186
batt_mass_climb = 0.3139155238745049
batt_mass_cruise = 44.210697202200514

p_lam_t = [0.0145278, 0.00416494, 0.00403764, 0.000152, 0.000235849, 0.000218, 0.000221649, 0.000267198]
p_chord = [0.150923, 0.259222, 0.0585371, 0.0110631]
p_twist_d = [40.3462, 24.4475, 15.3903, 12.8269]
p_pitch_d_accel = -0.7570513647195785
p_pitch_d_climb = -4.069606183292288
p_pitch_d_cruise = 3.1318846028190266
velocity_accel = 10.501260318853793
velocity_climb = 10.266257940503746
p_Rtip = 1.501729344875304
p_rpm_accel = 1451.8850895404096
p_rpm_climb = 1695.474202231303
p_rpm_cruise = 784.9972813630643
w_aoa_d_accel = 27.11320997863539
w_aoa_d_climb = 27.11447859601246
w_aoa_d_cruise = 7.766772071054681
batt_mass_climb = 0.4860451299610589
batt_mass_cruise = 44.63688514593232

p_lam_t = [0.0105731, 0.0262581, 0.00268972, 0.000152, 0.000218, 0.000218, 0.000218, 0.000282581]
p_chord = [0.105731, 0.498632, 0.0788773, 0.01]
p_twist_d = [45.0, 29.7836, 30.7799, 4.07252]
p_pitch_d_accel = -3.458984659697543
p_pitch_d_climb = 7.935396853796824
p_pitch_d_cruise = 14.26965595440019
velocity_accel = 21.301556472253356
velocity_climb = 27.095332167790886
velocity_cruise = 52.15133544597434
p_Rtip = 0.11467394626910472
p_rpm_accel = 2358.6342244212633
p_rpm_climb = 3619.3615537175065
p_rpm_cruise = 2487.4860292849053
w_aoa_d_accel = 12.26283930004526
w_aoa_d_climb = 15.106955892487754
w_aoa_d_cruise = 5.1488919489075
w_halfspan = 5.5287547694939585
batt_mass_climb = 14.367634347616098
batt_mass_cruise = 197.64118027592212

velocity_cruise = 45.0

p_lam_t = [0.00222326, 0.00281947, 0.00277319, 0.000929092, 0.00378516, 0.00195539, 0.00213457, 0.000225858]
p_chord = [0.271108, 0.212022, 0.182992, 0.0506781]
p_twist_d = [44.8915, 33.0477, 21.6509, 17.9445]
p_pitch_d_accel = -0.924740799739214
p_pitch_d_climb = -6.989465107034364
p_pitch_d_cruise = -2.2629817876582354
velocity_accel = 16.276966955425443
velocity_climb = 16.458264019189215
p_Rtip = 0.9215544938969867
p_rpm_accel = 1674.8838127274337
p_rpm_climb = 2118.227869571815
p_rpm_cruise = 2119.6049057252862
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 31.400850534737774
w_aoa_d_climb = 26.58610456443759
w_aoa_d_cruise = 4.8765509211359666
batt_mass_climb = 1.474223396970738
batt_mass_cruise = 90.3811395095503

p_lam_t = [0.000769917, 0.000346722, 0.000873985, 0.000351561, 0.0178791, 0.000287277, 0.000218, 0.0133924]
p_chord = [0.342426, 0.474738, 0.5, 0.472257]
p_twist_d = [45.0, 42.1658, 34.6507, 33.2096]
p_pitch_d_accel = 3.53313748255421
p_pitch_d_climb = -9.074088524447475
p_pitch_d_cruise = 3.5304954718815056
velocity_accel = 27.497867407288954
velocity_climb = 26.83350635091321
p_Rtip = 0.6812521094758671
p_rpm_accel = 2027.4800713877564
p_rpm_climb = 1973.5089027205986
p_rpm_cruise = 1692.1811295846044
kv = 5.0
i0 = 6.0000000000102975
w_aoa_d_accel = 11.605561986849098
w_aoa_d_climb = 11.993179304944821
w_aoa_d_cruise = 4.630038442857595
batt_mass_climb = 0.01
batt_mass_cruise = 16.934780966528994

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.0298305, 0.000552959, 0.000445033, 0.00944955]
p_chord = [0.304807, 0.79349, 0.342862, 0.36803]
p_twist_d = [45.0, 44.9926, 32.184, 34.9786]
p_pitch_d_accel = 2.490638226204294
p_pitch_d_climb = -11.190915299963157
p_pitch_d_cruise = 5.152147247145817
velocity_accel = 27.850111073746945
velocity_climb = 25.1777387031133
p_Rtip = 0.9174996029980375
p_rpm_accel = 1474.4458860420234
p_rpm_climb = 1389.9830668422587
p_rpm_cruise = 1027.6163176923305
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 11.878795614238367
w_aoa_d_climb = 15.267407568490794
w_aoa_d_cruise = 4.988402508795004
batt_mass_climb = 0.01
batt_mass_cruise = 69.77485442002452

p_lam_t = [0.012852, 0.000152, 0.00099182, 0.000152, 0.00724276, 0.0118719, 0.000496507, 0.010141]
p_chord = [0.12852, 0.395259, 0.0935259, 0.112919]
p_twist_d = [27.1948, 45.0, 26.8326, 14.0146]
p_pitch_d_accel = -8.888140333313912
p_pitch_d_climb = -12.97137745101056
p_pitch_d_cruise = -2.294140956952014
velocity_accel = 14.881030197399127
velocity_climb = 21.48714867233348
p_Rtip = 0.7725439532868138
p_rpm_accel = 1899.6655508008566
p_rpm_climb = 2668.3088402592334
p_rpm_cruise = 2121.4545554506954
kv = 5.136712906792837
i0 = 5.859045434420684
w_aoa_d_accel = 33.5154367732306
w_aoa_d_climb = 16.93408864189143
w_aoa_d_cruise = 4.3573376792476255
batt_mass_climb = 0.718467619858963
batt_mass_cruise = 77.21141337799982

p_lam_t = [0.0114561, 0.000881415, 0.000974668, 0.000297692, 0.000997415, 0.00686835, 0.00516608, 0.00165297]
p_chord = [0.210343, 0.348536, 0.184174, 0.159022]
p_twist_d = [44.9557, 40.0206, 25.7819, 33.485]
p_pitch_d_accel = 4.284567791567943
p_pitch_d_climb = -8.236394914491875
p_pitch_d_cruise = -2.6211610495310658
velocity_accel = 15.587162494247968
velocity_climb = 25.77395281148558
p_Rtip = 0.6754552859774366
p_rpm_accel = 1886.7123691467687
p_rpm_climb = 3300.2719800043424
p_rpm_cruise = 2408.519355731492
kv = 5.000529527786261
i0 = 5.999578762447254
w_aoa_d_accel = 32.16005283505898
w_aoa_d_climb = 23.096097835767015
w_aoa_d_cruise = 4.666212152975534
batt_mass_climb = 0.732671646235988
batt_mass_cruise = 82.1184464722631

p_lam_t = [0.0126445, 0.000152, 0.000681693, 0.000853767, 0.000532737, 0.00334368, 0.003328, 0.0033495]
p_chord = [0.132243, 0.40074, 0.0997776, 0.115934]
p_twist_d = [34.1531, 44.9999, 26.2506, 26.3513]
p_pitch_d_accel = -8.800520686355899
p_pitch_d_climb = -12.434239613486493
p_pitch_d_cruise = -2.294935240244655
velocity_accel = 14.800581975845887
velocity_climb = 21.908638756746264
p_Rtip = 0.776023153347384
p_rpm_accel = 2203.227371850445
p_rpm_climb = 2907.5740112278
p_rpm_cruise = 2121.3268313329736
kv = 5.073519434984816
i0 = 5.949795491223198
w_aoa_d_accel = 33.14987956278917
w_aoa_d_climb = 17.942707176478386
w_aoa_d_cruise = 4.6587273762243475
batt_mass_climb = 1.5930986310672346
batt_mass_cruise = 93.62371774344145

p_lam_t = [0.0132024, 0.000152, 0.000962997, 0.00165747, 0.00178382, 0.00336276, 0.00130714, 0.000218]
p_chord = [0.132024, 0.401106, 0.110586, 0.11794]
p_twist_d = [44.8114, 44.1184, 26.0984, 21.1484]
p_pitch_d_accel = -8.7817454406531
p_pitch_d_climb = -11.402049722002099
p_pitch_d_cruise = -2.2951290702117855
velocity_accel = 14.955080283053897
velocity_climb = 21.644623481806526
p_Rtip = 0.7866220275948844
p_rpm_accel = 2207.5260613079163
p_rpm_climb = 3143.1645843530455
p_rpm_cruise = 2121.2944429578906
kv = 5.0
i0 = 5.736989068716926
w_aoa_d_accel = 33.11721525763472
w_aoa_d_climb = 18.977963784469193
w_aoa_d_cruise = 4.685219774463085
batt_mass_climb = 0.741717018704567
batt_mass_cruise = 94.67662908384328

p_lam_t = [0.0135224, 0.000156558, 0.00268904, 0.000152, 0.000218, 0.00368144, 0.000279357, 0.000218]
p_chord = [0.139617, 0.3853, 0.105027, 0.129809]
p_twist_d = [45.0, 39.6983, 26.4305, 25.9365]
p_pitch_d_accel = -8.640823013758482
p_pitch_d_climb = -12.161900188152941
p_pitch_d_cruise = -2.2964760934622315
velocity_accel = 15.480722697563957
velocity_climb = 18.25354375051405
p_Rtip = 0.8120555017247091
p_rpm_accel = 2230.7834937256853
p_rpm_climb = 3195.717206018132
p_rpm_cruise = 2121.0800729835614
kv = 5.0
i0 = 5.912219495319556
w_aoa_d_accel = 33.14008538347226
w_aoa_d_climb = 23.914178390049628
w_aoa_d_cruise = 4.9141849210500785
batt_mass_climb = 0.7913576753939001
batt_mass_cruise = 101.18051355077026

p_lam_t = [0.01372, 0.000152, 0.00180032, 0.000152, 0.000218, 0.00375444, 0.000900508, 0.000218]
p_chord = [0.137201, 0.385153, 0.100022, 0.109661]
p_twist_d = [45.0, 31.9281, 20.2575, 15.8918]
p_pitch_d_accel = -7.812408706212089
p_pitch_d_climb = -7.811049362438349
p_pitch_d_cruise = -2.3049319409054947
velocity_accel = 13.823602198847928
velocity_climb = 13.709425718655707
p_Rtip = 0.9669676655260483
p_rpm_accel = 2683.8691799526364
p_rpm_climb = 2683.752879087068
p_rpm_cruise = 2116.0790054663253
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.400816381376096
w_aoa_d_climb = 31.62065353546675
w_aoa_d_cruise = 4.799913878461609
batt_mass_climb = 0.8998226399926497
batt_mass_cruise = 81.21165070617333

#fixed aoa
p_lam_t = [0.0142781, 0.000152, 0.00188924, 0.000152, 0.000218, 0.003802, 0.000889278, 0.000218]
p_chord = [0.142781, 0.377246, 0.103279, 0.105682]
p_twist_d = [45.0, 43.7044, 28.8957, 28.9261]
p_pitch_d_accel = -17.061087024445733
p_pitch_d_climb = -16.565033626266022
p_pitch_d_cruise = -1.855192384839103
velocity_accel = 12.32160737714068
velocity_climb = 15.503142775545795
p_Rtip = 0.922466733998593
p_rpm_accel = 2813.2204154552664
p_rpm_climb = 2813.220415455302
p_rpm_cruise = 1542.696118180986
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.15815515116242
w_aoa_d_climb = 32.66528362606995
w_aoa_d_cruise = 4.759955753408475
batt_mass_climb = 0.8274430880697909
batt_mass_cruise = 74.00662699365216

#Fixed Structure
p_lam_t = [0.0122773, 0.000152, 0.000673398, 0.00111581, 0.000218569, 0.00308722, 0.00126445, 0.000218]
p_chord = [0.12889, 0.376057, 0.103255, 0.104337]
p_twist_d = [45.0, 43.4668, 28.9718, 28.5676]
p_pitch_d_accel = -16.71549011901409
p_pitch_d_climb = -16.223776067761857
p_pitch_d_cruise = -1.8584256514623925
velocity_accel = 12.46421941125525
velocity_climb = 15.56169996976199
p_Rtip = 0.9208988562728082
p_rpm_accel = 2817.9949495175792
p_rpm_climb = 2818.010071379018
p_rpm_cruise = 1541.9368051275924
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.217186903840854
w_aoa_d_climb = 32.731626825976996
w_aoa_d_cruise = 4.799268992835256
batt_mass_climb = 0.8536579285163621
batt_mass_cruise = 73.54756854997733


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
    af = JLD.load("$(path)/airfoils/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db60af.jld")
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
