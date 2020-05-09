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
    numprops = 16
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





    p_lam_t = [0.00296032, 0.0015862, 0.000162432, 0.000559344, 0.000218, 0.000218, 0.000218, 0.000218]
    p_chord = [0.152775, 0.209865, 0.0599548, 0.0146912]
    p_twist_d = [43.7258, 36.3345, 24.8117, 15.5926]
    p_pitch_d_accel = 1.4774821755019851
    p_pitch_d_climb = 0.07394619362544558
    p_pitch_d_cruise = 1.2375953673499482
    velocity_accel = 19.472605774436627
    velocity_climb = 19.727538432606607
    velocity_cruise = 44.93261088033517
    p_Rtip = 0.7373054880998322
    p_rpm_accel = 1521.2024594724965
    p_rpm_climb = 1662.377310410668
    p_rpm_cruise = 2030.4950665640044
    number_blades = 2.0#2.3818948496117045
    kv = 9.01975436914855
    i0 = 5.942697466865178
    w_aoa_d_accel = 26.783254994297305
    w_aoa_d_climb = 10.682348122387333
    w_aoa_d_cruise = 1.9046025617909106
    batt_mass_accel = 0.01
    batt_mass_climb = 0.44404948094245933
    batt_mass_cruise = 119.39074611391598

    p_lam_t = [0.00336306, 0.0024627, 0.00221994, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
    p_chord = [0.155659, 0.201413, 0.0971402, 0.0100891]
    p_twist_d = [45.0, 33.5675, 22.6151, 14.0563]
    p_pitch_d_accel = 1.4829935332414101
    p_pitch_d_climb = -7.026136800337877
    p_pitch_d_cruise = 1.2383248802411353
    velocity_accel = 19.24363792774069
    velocity_climb = 12.153346213425108
    velocity_cruise = 44.81388239980551
    p_Rtip = 0.7782175243251072
    p_rpm_accel = 1543.53387273026
    p_rpm_climb = 2019.2972133656845
    p_rpm_cruise = 2030.550389292149
    kv = 7.219833930609759
    i0 = 6.0
    w_aoa_d_accel = 26.8000091392822
    w_aoa_d_climb = 18.92926786023156
    w_aoa_d_cruise = 1.9213751825354857
    batt_mass_accel = 0.01
    batt_mass_climb = 0.2641049737151147
    batt_mass_cruise = 120.95934410306408



    p_lam_t = [0.00265251, 0.000843842, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
    p_chord = [0.155659, 0.220257, 0.0117648, 0.0100891]
    p_twist_d = [45.0, 40.023, 28.6333, 1.88758]
    p_pitch_d_accel = 1.4829935332414101
    p_pitch_d_climb = -1.635119122794496
    p_pitch_d_cruise = -2.9277256334203314
    velocity_accel = 19.24363792774069
    velocity_climb = 13.643713798726896
    velocity_cruise = 44.73775558213799
    p_Rtip = 0.7118124546271428
    p_rpm_accel = 1976.2496570673752
    p_rpm_climb = 1943.6196478113399
    p_rpm_cruise = 2030.550389292149
    number_blades = 5.0#4.018472113116541
    kv = 8.155714770621731
    i0 = 5.997437563806129
    w_aoa_d_accel = 11.40529865732635
    w_aoa_d_climb = 20.23462263172518
    w_aoa_d_cruise = 1.931104840946755
    batt_mass_accel = 0.01
    batt_mass_climb = 0.16314428345511667
    batt_mass_cruise = 143.03066784675627

    p_lam_t = [0.0013721, 0.00095776, 0.00337691, 0.00235292, 0.0027877, 0.00075835, 0.00598907, 0.0462041]
    p_chord = [0.0969963, 0.44679, 0.175222, 0.0699541]
    p_twist_d = [33.2043, 45.0, 44.6104, 25.5642]
    p_pitch_d_accel = -24.34541997897509
    p_pitch_d_climb = -21.744777377644418
    p_pitch_d_cruise = -0.7531226039535881
    velocity_accel = 13.398350068887677
    velocity_climb = 16.704126149187054
    p_Rtip = 0.3493207679907946
    p_rpm_accel = 2740.550560953994
    p_rpm_climb = 2687.0237427401817
    p_rpm_cruise = 2553.5940186886087
    kv = 6.519522326064706
    i0 = 3.7616679058545
    w_aoa_d_accel = 37.67264350429828
    w_aoa_d_climb = 26.9287601118783
    w_aoa_d_cruise = 5.140434566691834
    batt_mass_accel = 0.001
    batt_mass_climb = 0.01
    batt_mass_cruise = 70.60902845083359

    p_lam_t = [0.00407943, 0.000338291, 0.000879098, 0.00467982, 0.000426051, 0.000218, 0.000504417, 0.0248299]
p_chord = [0.0649888, 0.463765, 0.170932, 0.166558]
p_twist_d = [45.0, 45.0, 44.2166, 34.369]
p_pitch_d_accel = -23.930194689275076
p_pitch_d_climb = -15.922557278863426
p_pitch_d_cruise = -0.6799300627541867
velocity_accel = 13.33078168102204
velocity_climb = 22.016510486356125
p_Rtip = 0.33530132186135064
p_rpm_accel = 2632.0094240780013
p_rpm_climb = 2529.2914565674614
p_rpm_cruise = 2565.677954227275
kv = 7.885840214630293
i0 = 0.40315507200272044
w_aoa_d_accel = 35.14598685378657
w_aoa_d_climb = 21.01061694925074
w_aoa_d_cruise = 4.648581714715574
batt_mass_accel = 0.001
batt_mass_climb = 0.12956306274306026
batt_mass_cruise = 9.765314584002086

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.00401272, 0.00245349, 0.00124912, 0.00674598]
p_chord = [0.0818858, 0.439626, 0.162802, 0.164482]
p_twist_d = [43.1119, 45.0, 40.9834, 28.0382]
p_pitch_d_accel = -14.35044805020054
p_pitch_d_climb = -9.255212493712818
p_pitch_d_cruise = 2.031738811143535
velocity_accel = 15.432111856296942
velocity_climb = 22.784818091923572
p_Rtip = 0.3078280145257298
p_rpm_accel = 2750.098427870975
p_rpm_climb = 2525.507729126974
p_rpm_cruise = 2798.7293653391066
kv = 9.00111664281726
i0 = 3.288202520284668
w_aoa_d_accel = 30.240597703491854
w_aoa_d_climb = 16.135812299419914
w_aoa_d_cruise = 4.643906823208933
batt_mass_accel = 0.001
batt_mass_climb = 3.3976635350712883
batt_mass_cruise = 121.5834491141962

#fixed aoa

p_lam_t = [0.00359077, 0.000370604, 0.000906367, 0.000691579, 0.00359077, 0.000218, 0.000690707, 0.000218]
p_chord = [0.0359077, 0.296071, 0.101652, 0.15762]
p_twist_d = [44.5486, 36.5624, 28.6058, 30.186]
p_pitch_d_accel = -11.328197156055754
p_pitch_d_climb = 4.510882190328494
p_pitch_d_cruise = 7.167826758086549
velocity_accel = 13.400420973021937
velocity_climb = 16.08917053705653
p_Rtip = 0.33529411764705896
p_rpm_accel = 3474.562027089696
p_rpm_climb = 2729.8834716163
p_rpm_cruise = 2935.6069223748345
kv = 13.384752042877295
i0 = 5.141624400462437
w_aoa_d_accel = 27.95529523309207
w_aoa_d_climb = 21.43192872111851
w_aoa_d_cruise = 4.666048465648842
batt_mass_accel = 0.004121895602014338
batt_mass_climb = 1.2005745202271156
batt_mass_cruise = 110.01534137452711

p_lam_t = [0.00375735, 0.000152, 0.000152682, 0.000183119, 0.00226924, 0.000218, 0.000443209, 0.000747398]
p_chord = [0.0399423, 0.312993, 0.101194, 0.277143]
p_twist_d = [44.9429, 42.696, 33.732, 20.4753]
p_pitch_d_accel = -15.631059277149353
p_pitch_d_climb = -2.0116425784519776
p_pitch_d_cruise = 7.195405538495134
velocity_accel = 12.570574777335706
velocity_climb = 14.117479035333085
p_Rtip = 0.3114304549748029
p_rpm_accel = 3532.538632035684
p_rpm_climb = 2898.756591852537
p_rpm_cruise = 2920.9993594464054
kv = 9.161093464688875
i0 = 4.163284474847829
w_aoa_d_accel = 32.54207169141614
w_aoa_d_climb = 24.493128055653564
w_aoa_d_cruise = 4.376143614693676
batt_mass_accel = 0.001
batt_mass_climb = 1.4873429916276262
batt_mass_cruise = 111.33064625794498


p_lam_t = [0.00341, 0.000152, 0.000589782, 0.000980885, 0.00340333, 0.000218, 0.000218, 0.000218]
p_chord = [0.0341, 0.323338, 0.0994817, 0.279841]
p_twist_d = [44.4665, 39.6343, 29.9125, 16.5639]
p_pitch_d_accel = -12.524522835042525
p_pitch_d_climb = 0.9685613637119833
p_pitch_d_cruise = 7.16686685338314
velocity_accel = 12.357786588283695
velocity_climb = 14.095539367627607
p_Rtip = 0.3352403923309577
p_rpm_accel = 3262.901280334454 +100.0
p_rpm_climb = 2671.765870925417 + 100.0
p_rpm_cruise = 2917.433833310762
kv = 7.607033706073792
i0 = 3.247269758179325
w_aoa_d_accel = 32.95463166244449
w_aoa_d_climb = 24.585469345540336
w_aoa_d_cruise = 4.473201307294152
batt_mass_accel = 0.001
batt_mass_climb = 1.3610062426857141
batt_mass_cruise = 104.24633528700593


#from supcomp 7/27
p_lam_t = [0.00371678, 0.000152, 0.000616299, 0.0010045, 0.00371645, 0.000218, 0.000641126, 0.000999863]
p_chord = [0.0371678, 0.340298, 0.10521, 0.29571]
p_twist_d = [42.3581, 39.8477, 30.282, 18.871]
p_pitch_d_accel = -13.705232863474155
p_pitch_d_climb = -1.3157966269654875
p_pitch_d_cruise = 7.153600857292216
velocity_accel = 12.167361678110272
velocity_climb = 13.149903343403887
p_Rtip = 0.33529411764705885
p_rpm_accel = 3442.557114540962
p_rpm_climb = 2876.346557900861
p_rpm_cruise = 2917.385269765737
kv = 7.564690720996396
i0 = 2.914988305549223
w_aoa_d_accel = 33.03769640188165
w_aoa_d_climb = 26.71777790852381
w_aoa_d_cruise = 4.402494361501817
batt_mass_accel = 0.001
batt_mass_climb = 1.2416120927206575
batt_mass_cruise = 98.70817710152097

#Fixed structure
p_lam_t = [0.00317119, 0.000152, 0.000602642, 0.000999455, 0.00305703, 0.000218, 0.00061662, 0.000985486]
p_chord = [0.0317119, 0.340654, 0.105173, 0.29573]
p_twist_d = [44.1595, 39.7716, 30.1964, 18.8546]
p_pitch_d_accel = -13.645922098333733
p_pitch_d_climb = -1.2065592868395605
p_pitch_d_cruise = 7.152977490363259
velocity_accel = 12.116310917896278
velocity_climb = 13.219746999857184
p_Rtip = 0.33528351477776414
p_rpm_accel = 3445.527963753064
p_rpm_climb = 2872.2248324549314
p_rpm_cruise = 2917.309073626354
kv = 7.472533298794895
i0 = 3.0717978416348415
w_aoa_d_accel = 33.016389441651214
w_aoa_d_climb = 26.46762602857269
w_aoa_d_cruise = 4.391175616531284
batt_mass_accel = 0.001
batt_mass_climb = 1.1797102485834798
batt_mass_cruise = 97.59207399490616



    number_blades = 3.0





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
#
# xwarm = [0.01000, 0.01000, 0.01000, 0.01000, 0.01000, 0.01000, 0.01000, 0.01000, 21.62153, 7.95808, 20.34969, 5.00338, 4.33754, 4.49540, 4.49814, 3.09121, -248.53136, -23.91318, -30.31537, 132.97473, 0.24305, 33.52939, 306.22862, 45.14119, 2476.55374, 4.99968, 7.26288, 4.03376, 0.33206, 14.33827, 4.75615, 0.00100, 0.07764, 13.59082]
#     x0[1:length(xwarm)] = xwarm[1:end]

    table_af_wing_blown = JLD.load("$path/airfoils/af_prop_E212.jld")
    table_3D_wing_blown = table_af_wing_blown["NDtable"]
    spline_3D_blown = AirfoilPrep.SplineND_from_tableND(table_3D_wing_blown)


    table_af_wing_free = JLD.load("$path/airfoils/af_prop_E212.jld")
    table_3D_wing_free = table_af_wing_free["NDtable"]
    spline_3D_freestream = AirfoilPrep.SplineND_from_tableND(table_3D_wing_free)

    # # Load the prop airfoil ND table
    # ClarkY_non_extrap = JLD.load("$(path)/airfoils/af_prop_E212.jld")
    # ClarkY_non_extrap = ClarkY_non_extrap["NDtable"]
    #
    # omega = p_rpm_cruise*pi/30
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
    af = JLD.load("$(path)/airfoils/p8_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))af.jld")
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

        for i = 1:3

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

            x0 = xopt

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
