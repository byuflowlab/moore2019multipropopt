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
    cruise_factor = 150.0
    required_cruise = cruise_factor * 0.447 #MPH TO m/s


    # designvars, constraints, ctlparams = getoptparams()
    numprops = 32
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

    p_lam_t = [0.000152, 0.000152, 0.000152, 0.00016314, 0.000218, 0.000218, 0.000697237, 0.000218]
    p_chord = [0.159116, 0.0881603, 0.0145218, 0.010409]
    p_twist_d = [45.0, 41.3985, 28.554, 13.948]
    p_pitch_d_accel = -6.702376950041459
    p_pitch_d_climb = -0.9984657091922631
    p_pitch_d_cruise = -2.9083478654690595
    velocity_accel = 13.309245111241173
    velocity_climb = 18.398578866562552
    p_Rtip = 0.7152941176470594
    p_rpm_accel = 2665.7033096670643
    p_rpm_climb = 2761.0476758587279
    p_rpm_cruise = 2030.5534516552439
    kv = 8.1779358810975
    i0 = 5.992269337186693
    w_aoa_d_accel = 16.168181853486463
    w_aoa_d_climb = 10.909173218251121
    w_aoa_d_cruise = 1.9182052665194256
    batt_mass_accel = 0.09288778086913847
    batt_mass_climb = 0.7261166000961358
    batt_mass_cruise = 124.77395766068845

    number_blades = 5.0



    p_lam_t = [0.000152, 0.000152, 0.000152, 0.000163126, 0.000218016, 0.000218, 0.000218, 0.000913548]
    p_chord = [0.300727, 0.12517, 0.0811, 0.0583803]
    p_twist_d = [41.6387, 39.1838, 34.5475, 13.5254]
    p_pitch_d_accel = -10.300854536492995
    p_pitch_d_climb = -13.655025374197056
    p_pitch_d_cruise = 6.936335111625699
    velocity_accel = 11.714582571100875
    velocity_climb = 13.207746915078056
    p_Rtip = 0.29979770350277596
    p_rpm_accel = 6047.892316807409
    p_rpm_climb = 4671.32644106729
    p_rpm_cruise = 3194.9027880993262
    kv = 34.94171095301614
    i0 = 4.546764300896736
    w_aoa_d_accel = 18.8984454233384
    w_aoa_d_climb = 24.1842091870094
    w_aoa_d_cruise = 5.242166025930194
    batt_mass_accel = 0.01673233853656761
    batt_mass_climb = 0.6940059656702224
    batt_mass_cruise = 114.55639556051761

    p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
    p_chord = [0.155186, 0.0945468, 0.05034, 0.0503816]
    p_twist_d = [33.2087, 35.7593, 33.1212, 24.2418]
    p_pitch_d_accel = -10.349044233803216
    p_pitch_d_climb = -16.667030190277142
    p_pitch_d_cruise = 9.154851737404712
    velocity_accel = 12.154931321487691
    velocity_climb = 14.864115457086813
    p_Rtip = 0.31333653365907904
    p_rpm_accel = 6629.347586057403
    p_rpm_climb = 7079.519370726989
    p_rpm_cruise = 5115.715593225067
    kv = 7.8378786280275206
    i0 = 3.3653501817076554
    w_aoa_d_accel = 31.63405185679079
    w_aoa_d_climb = 20.091121003196104
    w_aoa_d_cruise = 2.1968990061950286
    batt_mass_accel = 0.007854470776490605
    batt_mass_climb = 0.8945104698106916
    batt_mass_cruise = 161.03397177760934

    p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
p_chord = [0.155375, 0.0929051, 0.0484615, 0.0504848]
p_twist_d = [33.0738, 35.0526, 32.8843, 23.6201]
p_pitch_d_accel = -10.347428800333091
p_pitch_d_climb = -16.577056332686908
p_pitch_d_cruise = 9.154850549322514
velocity_accel = 12.15287894115304+1
velocity_climb = 13.36126090832502+1
p_Rtip = 0.3164524521940784
p_rpm_accel = 6629.382747865235+100
p_rpm_climb = 7086.861666833731+100
p_rpm_cruise = 5115.715568904886+100
kv = 7.8426259536951966
i0 = 3.3534461734149135
w_aoa_d_accel = 20.543391337772253
w_aoa_d_climb = 20.226825748936097
w_aoa_d_cruise = 2.126464674327615+1
batt_mass_accel = 0.0041794634507749655
batt_mass_climb = 0.7629857891516558
batt_mass_cruise = 158.6819866406862

#fixed aoa
p_lam_t = [0.00030417, 0.000310081, 0.000306877, 0.000252632, 0.000435382, 0.000787447, 0.000438089, 0.000361828]
p_chord = [0.01, 0.0427193, 0.0288075, 0.00118958]
p_twist_d = [44.3463, 33.3237, 22.5238, 15.5652]
p_pitch_d_accel = -1.6496417148005986
p_pitch_d_climb = -2.2526543724076076
p_pitch_d_cruise = 25.40229576521106
velocity_accel = 10.175362346224656
velocity_climb = 13.418219869940812
p_Rtip = 0.16503786262789163
p_rpm_accel = 8504.122456928331+ 900.0
p_rpm_climb = 7641.321253592031+ 900.0
p_rpm_cruise = 5577.194744033578+ 900.0
kv = 59.79044776110737
i0 = 2.8318004504437426
w_aoa_d_accel = 31.47205613527791
w_aoa_d_climb = 27.89394360279746
w_aoa_d_cruise = 4.0188076706202254
batt_mass_accel = 0.0020250043511871405
batt_mass_climb = 0.1274690349528058
batt_mass_cruise = 133.263665023343778

p_lam_t = [0.00030417, 0.000310081, 0.000306877, 0.000182507, 0.000279337, 0.000787447, 0.000438089, 0.000258195]
p_chord = [0.01, 0.0394062, 0.0217687, 0.00159402]
p_twist_d = [44.8178, 14.8528, 7.38811, 8.47636]
p_pitch_d_accel = 10.30848523036365
p_pitch_d_climb = 10.83881296264575
p_pitch_d_cruise = 32.5399632624233
velocity_accel = 16.419512842542964
velocity_climb = 15.249211108201559
p_Rtip = 0.18939383052897163
p_rpm_accel = 8138.723711847338
p_rpm_climb = 8147.2213991921135
p_rpm_cruise = 7977.238693273193
kv = 26.57102123412695
i0 = 2.6812740414389493
w_aoa_d_accel = 20.960545490355038
w_aoa_d_climb = 26.850339265710122
w_aoa_d_cruise = 1.3060820541181184
batt_mass_accel = 0.0012857487756668846
batt_mass_climb = 0.042747795536339436
batt_mass_cruise = 40.757102104512256

#Try 100 case as warmstart
p_lam_t = [0.000152008, 0.000152, 0.000152, 0.000152, 0.000218011, 0.000218, 0.000218159, 0.000218]
p_chord = [0.0709209, 0.0643777, 0.0209693, 0.065042]
p_twist_d = [38.3983, 38.024, 27.9631, 26.09]
p_pitch_d_accel = -9.90027456226873
p_pitch_d_climb = -11.289835122782023
p_pitch_d_cruise = 5.794831324824401
velocity_accel = 5.99658939977826
velocity_climb = 18.79501119133917
p_Rtip = 0.1726501662373227
p_rpm_accel = 15026.298396934739
p_rpm_climb = 14636.548518896087
p_rpm_cruise = 9999.999083479384
kv = 5.237707327236063
i0 = 2.7374477496284517
w_aoa_d_accel = 29.436833938528682
w_aoa_d_climb = 29.88453534239406
w_aoa_d_cruise = 2.387265208712961
batt_mass_accel = 0.001
batt_mass_climb = 0.08767773505016585
batt_mass_cruise = 195.43375439257395


p_lam_t = [0.000152008, 0.000152, 0.000152, 0.000152, 0.000218011, 0.000218, 0.000218159, 0.000218]
p_chord = [0.0709209, 0.0643777, 0.0209693, 0.065042]
p_twist_d = [38.3983, 38.024, 27.9631, 26.09]
p_pitch_d_accel = -9.90027456226873
p_pitch_d_climb = -11.289835122782023
p_pitch_d_cruise = 5.794831324824401
velocity_accel = 5.99658939977826
velocity_climb = 18.79501119133917
p_Rtip = 0.1726501662373227
p_rpm_accel = 15026.298396934739
p_rpm_climb = 14636.548518896087
p_rpm_cruise = 9999.999083479384
kv = 5.237707327236063
i0 = 2.7374477496284517
w_aoa_d_accel = 29.436833938528682
w_aoa_d_climb = 29.88453534239406
w_aoa_d_cruise = 2.387265208712961
batt_mass_accel = 0.001
batt_mass_climb = 0.08767773505016585
batt_mass_cruise = 195.43375439257395

p_lam_t = [0.000152, 0.000152001, 0.000152001, 0.000152, 0.000218, 0.000218002, 0.000218002, 0.000218]
p_chord = [0.0603394, 0.0447857, 0.0278678, 0.0389109]
p_twist_d = [41.5167, 34.7123, 29.1934, 30.0209]
p_pitch_d_accel = -9.913929511250501
p_pitch_d_climb = -12.711700172021517
p_pitch_d_cruise = 5.798604286799855
velocity_accel = 5.833141314579384
velocity_climb = 7.903091499318553
p_Rtip = 0.1727272727272728
p_rpm_accel = 15022.255518590646
p_rpm_climb = 14924.818699646814
p_rpm_cruise = 10000.0
kv = 6.7414388605388496
i0 = 1.1718516026158536
w_aoa_d_accel = 27.74311701357958
w_aoa_d_climb = 28.12344028041343
w_aoa_d_cruise = 1.977989281680936
batt_mass_accel = 0.001
batt_mass_climb = 2.33390543268507
batt_mass_cruise = 193.87414193417365




p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
p_chord = [0.0612406, 0.0464926, 0.03355, 0.0599469]
p_twist_d = [44.9741, 37.4358, 34.0795, 33.1921]
p_pitch_d_accel = -9.779273219389722
p_pitch_d_climb = -18.42213456899951
p_pitch_d_cruise = 5.798632648812015
velocity_accel = 6.535270095945317
velocity_climb = 7.385353771354927
p_Rtip = 0.17272681000651594
p_rpm_accel = 15024.316255304815
p_rpm_climb = 15024.316470339769
p_rpm_cruise = 9999.998271690989
kv = 6.698807617032466
i0 = 1.0728380295038447
w_aoa_d_accel = 30.717682757054376
w_aoa_d_climb = 30.586151857634714
w_aoa_d_cruise = 1.963963803933461
batt_mass_accel = 0.0010021059157188737
batt_mass_climb = 2.3309630878898395
batt_mass_cruise = 188.74505245322715

#Fixed Structure
p_lam_t = [0.000430107, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
p_chord = [0.102634, 0.0354151, 0.0310808, 0.0310067]
p_twist_d = [44.997, 45.0, 39.0873, 39.6742]
p_pitch_d_accel = -19.866102627280487
p_pitch_d_climb = -23.47643537594124
p_pitch_d_cruise = -4.734995186628237
velocity_accel = 5.881549223569948
velocity_climb = 7.186074148464279
p_Rtip = 0.17272727272727287
p_rpm_accel = 13250.149012431748
p_rpm_climb = 15024.276221562994
p_rpm_cruise = 9996.76972777564
kv = 9.278769941274271
i0 = 1.5181522183333138
w_aoa_d_accel = 30.10356655708142
w_aoa_d_climb = 31.165854917561827
w_aoa_d_cruise = 1.9795274634405933
batt_mass_accel = 0.001
batt_mass_climb = 2.247693285801534
batt_mass_cruise = 196.25587677749502



number_blades = 5.0
velocity_cruise = 150*0.447
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
    options["Major optimality tolerance"] = 1e-7 #Should be scaled so it is optimal with a solid 2 digits
    # options["Minor optimality  tolerance"] = 1e-6
    options["Major feasibility tolerance"] = 5e-6
    # options["Minor feasibility tolerance"] = 1e-6
    # options["Minor print level"] = 1E8
    # options["Print frequency"] = 1E8
    options["Scale option"] = 1
    options["Scale tolerance"] = .98
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

        for i = 1:1
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
