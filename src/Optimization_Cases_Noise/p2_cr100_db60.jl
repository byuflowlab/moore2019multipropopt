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

p_lam_t = [0.0138243, 0.00221821, 0.00292125, 0.000847897, 0.000218, 0.00276569, 0.000558416, 0.000991551]
p_chord = [0.148325, 0.251722, 0.104045, 0.0500003]
p_twist_d = [45.0, 28.5267, 20.0406, 16.8619]
p_pitch_d_accel = -3.710912237966906
p_pitch_d_climb = -5.366196278030933
p_pitch_d_cruise = -2.3131717130396776
velocity_accel = 13.994908934316992
velocity_climb = 13.31330255979673
p_Rtip = 1.015269927034646
p_rpm_accel = 2391.7734349536427
p_rpm_climb = 2556.426772122603
p_rpm_cruise = 2122.836317947841
w_aoa_d_accel = 32.55800045808513
w_aoa_d_climb = 32.47467481756906
w_aoa_d_cruise = 4.829469543115608
batt_mass_climb = 0.8583782045440271
batt_mass_cruise = 79.38145510209087

p_lam_t = [0.00930054, 0.000152063, 0.000154004, 0.000153933, 0.000993078, 0.0080352, 0.00334944, 0.00917818]
p_chord = [0.0539427, 0.469051, 0.0504748, 0.0918895]
p_twist_d = [43.0637, 31.4249, 25.87, 2.96723]
p_pitch_d_accel = -9.948838157281076
p_pitch_d_climb = -13.469131945806891
p_pitch_d_cruise = -2.2865989754292566
velocity_accel = 13.345888227910004
velocity_climb = 23.136658212682246
p_Rtip = 0.9273307083599485
p_rpm_accel = 1787.4709387014916
p_rpm_climb = 2401.732857902315
p_rpm_cruise = 2122.7803059009098
kv = 5.0000004052835445
i0 = 5.999999647440097
w_aoa_d_accel = 39.8258946787143
w_aoa_d_climb = 15.418741937168534
w_aoa_d_cruise = 4.197239148366411
batt_mass_climb = 0.01028227003485989
batt_mass_cruise = 73.66968499472041

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

p_lam_t = [0.00102089, 0.00240222, 0.000152, 0.00126711, 0.00490195, 0.000813647, 0.00117604, 0.000218]
p_chord = [0.345811, 0.388963, 0.0924136, 0.129313]
p_twist_d = [45.0, 42.4899, 23.3364, 11.1928]
p_pitch_d_accel = -5.52055187708344
p_pitch_d_climb = -9.795349396466072
p_pitch_d_cruise = -2.3125088348437948
velocity_accel = 14.54142275665122
velocity_climb = 22.61982867459151
p_Rtip = 0.829756326723427
p_rpm_accel = 1913.4307238218826
p_rpm_climb = 2452.1797935106297
p_rpm_cruise = 2118.1005204391095
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 33.07611838738169
w_aoa_d_climb = 15.899405694920365
w_aoa_d_cruise = 4.3548594797530535
batt_mass_climb = 0.6600771355898088
batt_mass_cruise = 77.97063194466823

p_lam_t = [0.000152, 0.000874961, 0.00102097, 0.00061527, 0.00418596, 0.00245375, 0.000703365, 0.000397152]
p_chord = [0.342623, 0.386695, 0.0900426, 0.133361]
p_twist_d = [45.0, 41.7869, 22.7627, 9.7627]
p_pitch_d_accel = -4.913514824569793
p_pitch_d_climb = -8.922976082655826
p_pitch_d_cruise = -2.3125088348437948
velocity_accel = 14.541625888231374
velocity_climb = 23.77762890187704
p_Rtip = 0.8597885553191547
p_rpm_accel = 1903.5538188952223
p_rpm_climb = 2336.2344450153546
p_rpm_cruise = 2118.1005204391095
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.923794457058655
w_aoa_d_climb = 14.954024217616144
w_aoa_d_cruise = 5.212370116883085
batt_mass_climb = 0.6482243174133087
batt_mass_cruise = 87.6928149665623

#fixed aoa
p_lam_t = [0.000156411, 0.00081511, 0.000874745, 0.000623624, 0.00418703, 0.00243309, 0.000767288, 0.000218]
p_chord = [0.351612, 0.388643, 0.0898165, 0.134077]
p_twist_d = [45.0, 41.5987, 22.6368, 9.989]
p_pitch_d_accel = -5.082594183477591
p_pitch_d_climb = -9.394208479113368
p_pitch_d_cruise = -2.310127729915182
velocity_accel = 13.864769969090068
velocity_climb = 22.292663856115674
p_Rtip = 0.8646036314533656
p_rpm_accel = 1884.5652403818103
p_rpm_climb = 2326.1781210541153
p_rpm_cruise = 2116.1044008464883
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.83875887467559
w_aoa_d_climb = 17.395463391087702
w_aoa_d_cruise = 5.652344847543516
batt_mass_climb = 0.6660735573785961
batt_mass_cruise = 89.15968986297709

#Fixed Structure
p_lam_t = [0.000156411, 0.00081511, 0.000874745, 0.000623624, 0.00339963, 0.00195669, 0.000218, 0.000249417]
p_chord = [0.35375, 0.388884, 0.0926167, 0.130399]
p_twist_d = [45.0, 42.0176, 23.0043, 9.8978]
p_pitch_d_accel = -4.58952220703476
p_pitch_d_climb = -9.372825175867224
p_pitch_d_cruise = -2.3148988340573773
velocity_accel = 13.924015620548477
velocity_climb = 22.103653360731137
p_Rtip = 0.8428968428452923
p_rpm_accel = 1931.5288841979334
p_rpm_climb = 2416.4192269583764
p_rpm_cruise = 2115.3384253036766
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.896240763350846
w_aoa_d_climb = 17.632150843057136
w_aoa_d_cruise = 4.806202876760801
batt_mass_climb = 0.40193561878888495
batt_mass_cruise = 81.32291033985385

#Try 2 63
p_lam_t = [0.010194, 0.000152, 0.00101874, 0.000179783, 0.000218, 0.00230649, 0.000218, 0.000218]
p_chord = [0.12866, 0.474082, 0.112861, 0.0694404]
p_twist_d = [44.9997, 44.8837, 23.2551, 18.4328]
p_pitch_d_accel = -3.5150527244375525
p_pitch_d_climb = -9.327895096375247
p_pitch_d_cruise = -2.3963543266946785
velocity_accel = 14.082893838693128
velocity_climb = 19.735983020258633
p_Rtip = 0.8220909295534947
p_rpm_accel = 1947.1680121564243
p_rpm_climb = 2405.081001812947
p_rpm_cruise = 2105.3013037438486
kv = 5.0
i0 = 6.0
w_aoa_d_accel = 32.63482477740545
w_aoa_d_climb = 21.9449591894635
w_aoa_d_cruise = 4.511418716869237
batt_mass_climb = 0.7584790106312591
batt_mass_cruise = 101.102614154679


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
