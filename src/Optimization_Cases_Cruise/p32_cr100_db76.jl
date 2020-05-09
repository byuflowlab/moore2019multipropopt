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

    p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
    p_chord = [0.326523, 0.0782953, 0.0869311, 0.0540222]
    p_twist_d = [45.0, 45.0, 44.9963, 43.0019]
    p_pitch_d_accel = -17.84083883893648
    p_pitch_d_climb = -25.0
    p_pitch_d_cruise = -0.1939052889102814
    velocity_accel = 10.070363699390153
    velocity_climb = 9.112159844132988
    p_Rtip = 0.33529411764705896
    p_rpm_accel = 3728.1188216845803
    p_rpm_climb = 4788.215359867585
    p_rpm_cruise = 2509.8050053693323
    kv = 6.1169974384665675
    i0 = 3.2352764130220324
    w_aoa_d_accel = 31.55695290142486
    w_aoa_d_climb = 29.920729687705446
    w_aoa_d_cruise = 4.7631567415964575
    batt_mass_accel = 0.001
    batt_mass_climb = 1.316145235752836
    batt_mass_cruise = 120.23016667996093

    p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
p_chord = [0.297983, 0.0795527, 0.084264, 0.0496087]
p_twist_d = [44.9891, 44.9957, 44.9887, 42.4323]
p_pitch_d_accel = -18.02215851665644
p_pitch_d_climb = -24.99901929675474
p_pitch_d_cruise = -0.13992840752574429
velocity_accel = 10.036014942261168
velocity_climb = 9.385040827888599
p_Rtip = 0.33529411764705896
p_rpm_accel = 3749.17886299392
p_rpm_climb = 4887.572254097
p_rpm_cruise = 2518.161446846273
kv = 6.1169974384665675
i0 = 3.2352764130220324
w_aoa_d_accel = 31.56448120393928
w_aoa_d_climb = 28.15845524966909
w_aoa_d_cruise = 4.759791288872362
batt_mass_accel = 0.001
batt_mass_climb = 1.6915365371493547
batt_mass_cruise = 119.07798147685504

p_lam_t = [0.000152, 0.000152, 0.000152, 0.000152, 0.000218, 0.000218, 0.000218, 0.000218]
p_chord = [0.297983, 0.0795527, 0.084264, 0.0496087]
p_twist_d = [44.9891, 44.9957, 44.9887, 42.4323]
p_pitch_d_accel = -18.02215851665644
p_pitch_d_climb = -24.99901929675474
p_pitch_d_cruise = -0.13992840752574429
velocity_accel = 10.036014942261168
velocity_climb = 9.385040827888599
p_Rtip = 0.33529411764705896
p_rpm_accel = 3749.17886299392
p_rpm_climb = 4887.572254097
p_rpm_cruise = 2518.161446846273
kv = 6.1169974384665675
i0 = 3.2352764130220324
w_aoa_d_accel = 31.56448120393928
w_aoa_d_climb = 28.15845524966909
w_aoa_d_cruise = 4.759791288872362
batt_mass_accel = 0.001
batt_mass_climb = 1.6915365371493547
batt_mass_cruise = 119.07798147685504

#fixed aoa
p_lam_t = [0.000152039, 0.000152, 0.000152709, 0.000152, 0.000218056, 0.000219017, 0.000219017, 0.000218012]
p_chord = [0.0746714, 0.0544777, 0.0832914, 0.138318]
p_twist_d = [42.1725, 42.3992, 40.0619, 42.6916]
p_pitch_d_accel = -6.1983587030869485
p_pitch_d_climb = -1.370391032611906
p_pitch_d_cruise = 13.480031833536383
velocity_accel = 15.842742296443285
velocity_climb = 32.13056847988888
p_Rtip = 0.1727272727272728
p_rpm_accel = 4345.681730822606
p_rpm_climb = 5479.251904169189
p_rpm_cruise = 4436.589982528049
kv = 7.385362299063935
i0 = 1.2372387909197842
w_aoa_d_accel = 22.380290133243466
w_aoa_d_climb = 10.874001379029613
w_aoa_d_cruise = 4.1290775383082154
batt_mass_accel = 0.001
batt_mass_climb = 1.7863030982999202
batt_mass_cruise = 161.6834420201758

p_lam_t = [0.000152039, 0.000152, 0.000152709, 0.000152, 0.000218056, 0.000219017, 0.000219017, 0.000218012]
p_chord = [0.0746714, 0.0544777, 0.0832914, 0.138318]
p_twist_d = [42.1725, 42.3992, 40.0619, 42.6916]
p_pitch_d_accel = -6.1983587030869485
p_pitch_d_climb = -1.370391032611906
p_pitch_d_cruise = 13.480031833536383
velocity_accel = 15.842742296443285
velocity_climb = 32.13056847988888
p_Rtip = 0.1727272727272728
p_rpm_accel = 4345.681730822606
p_rpm_climb = 5479.251904169189
p_rpm_cruise = 4436.589982528049
kv = 7.385362299063935
i0 = 1.2372387909197842
w_aoa_d_accel = 22.380290133243466
w_aoa_d_climb = 10.874001379029613
w_aoa_d_cruise = 4.1290775383082154
batt_mass_accel = 0.001
batt_mass_climb = 1.7863030982999202
batt_mass_cruise = 161.6834420201758

#Fixed AOA
p_lam_t = [0.000153553, 0.000152, 0.000154222, 0.000152, 0.000220227, 0.000218, 0.000221188, 0.000218]
p_chord = [0.0716154, 0.0570286, 0.0339698, 0.0604225]
p_twist_d = [44.9963, 27.6546, 30.6431, 30.0685]
p_pitch_d_accel = -9.958294771395655
p_pitch_d_climb = -9.875033717965351
p_pitch_d_cruise = 13.795670884292571
velocity_accel = 5.0
velocity_climb = 7.04829957046534
p_Rtip = 0.17272727272728414
p_rpm_accel = 15024.276221562472
p_rpm_climb = 15024.276221562008
p_rpm_cruise = 5917.325926589708
kv = 5.0
i0 = 0.8324888812466171
w_aoa_d_accel = 28.58870344149924
w_aoa_d_climb = 30.29724182384655
w_aoa_d_cruise = 4.3410025055206765
batt_mass_accel = 0.001
batt_mass_climb = 2.992518823651641
batt_mass_cruise = 137.0187593015944

#Fixed Structure
p_lam_t = [0.000153551, 0.000152, 0.00015422, 0.000152, 0.000220225, 0.000218, 0.000221186, 0.000218]
p_chord = [0.064936, 0.0468867, 0.0331557, 0.0604225]
p_twist_d = [45.0, 27.4481, 30.5663, 30.3685]
p_pitch_d_accel = -9.958294771395655
p_pitch_d_climb = -9.605670006720334
p_pitch_d_cruise = 13.795337021203192
velocity_accel = 5.0
velocity_climb = 7.1278035023102255
p_Rtip = 0.17272727272727284
p_rpm_accel = 15024.276221562472
p_rpm_climb = 15024.276221562985
p_rpm_cruise = 5917.325926589708
kv = 5.0
i0 = 0.8324888812466171
w_aoa_d_accel = 28.783632735342362
w_aoa_d_climb = 30.419796125540014
w_aoa_d_cruise = 4.343209920189757
batt_mass_accel = 0.001
batt_mass_climb = 2.913937555510022
batt_mass_cruise = 136.30187570347059

    number_blades = 5.0





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
