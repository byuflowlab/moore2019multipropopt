#------- SET UP DESIGN PARAMETERS --------#

struct plyproperties
    names#::Array{String,1}
    plies::Array{Composites.material,1}
end

struct geometry
    t::Array{Float64,1}
    s::Array{Float64,1}
    x::Array{Float64,1} #x-location of airfoil section (positive x = downstream)
    y::Array{Float64,1} #y-location of airfoil section (positive y = spanwise)
    z::Array{Float64,1} #z-location of airfoil section (positive z = up)
    chord::Array{Float64,1}
    normalchord::Array{Float64,1} #normal chord length
    twist::Array{Float64,1} #airfoil twist in degrees
    sweep::Array{Float64,1}
    dihedral::Array{Float64,1}
    airfoilthickness::Array{Float64,1}
    xaf::Array{Float64,2}
    yaf::Array{Float64,2}
end

struct compstructure
    material::Array{Array{Composites.material,1},1}
    laminate::Array{Composites.laminate,1}
    A::Array{Array{Float64,2},1}
    B::Array{Array{Float64,2},1}
    D::Array{Array{Float64,2},1}
    bucklingstrain::Array{Float64,2}
    precompinput::Array{PreComp.input,1}
    precompoutput::Array{PreComp.output,1}
end

struct designparameters
    # configuration assumptions
    p_n_akima::Int64
    p_radii::Array{Float64,1}
    p_webloc::Array{Float64,1}

    # operating point assumptions
    altitude::Float64 #m

    # material properties
    plyprops::plyproperties
    p_usedmaterials
    p_orientation::Array{Float64,1}

    # Propeller Design
    p_lam_t::Array{Float64,1}

    p_chord::Array{Float64,1}
    p_sweep_d::Array{Float64,1}
    p_twist_d::Array{Float64,1}
    p_pitch_d::Float64
    p_dihedral_d::Array{Float64,1}
    p_airfoilthickness::Array{Float64,1}
    p_aerocenter::Array{Float64,1}
    p_x_offset::Array{Float64,1}
    velocity::Float64
    p_Rtip::Float64
    blades::Float64
    numprops::Int64
    machtip::Float64
    p_rpm::Float64
    kv::Float64
    i0::Float64
    batt_mass::Float64

    # Wing Design
    w_n_linear::Int64
    w_webloc::Array{Float64,1}
    w_lam_t::Array{Float64,1}
    w_chord::Array{Float64,1}
    w_aoa_d::Float64
    w_twist_d::Array{Float64,1}
    w_halfspan::Float64
    w_nondim_halfspany::Array{Float64,1}
    w_sweep_d::Array{Float64,1}
    w_dihedral_d::Array{Float64,1}
    w_usedmaterials #array of strings???
    w_orientation::Array{Float64,1}
    w_x_offset::Array{Float64,1}
    w_aerocenter::Array{Float64,1}

    n_CP::Float64
    payload::Float64
    max_mass::Float64
    gravity::Float64
    h_set::Float64
    specif_energy::Float64
    required_range::Float64
    required_takeoff_dist::Float64
    altitude_accel::Float64
    altitude_climb::Float64
    altitude_cruise::Float64
    required_noise::Float64
    required_cruise::Float64

end

function designparameters(;
    p_n_akima = 20,
    p_radii::Array{Float64,1}=Float64[0.0,0.2,.4,.6,.8,1.0],
    p_webloc::Array{Float64,1}=[1.0],
    altitude::Float64=18288.0,
    plyprops::plyproperties=plyproperties(),
    p_usedmaterials = ["highmodulus_uni","highmodulus_weave","taylor_foam"],
    p_orientation=[0.0,-45],
    # structural design variables
    p_lam_t = ones(length(p_radii)*length(p_usedmaterials)),
    # geometric design variables
    p_chord = Float64[0.1, 0.1, 0.1, 0.1, 0.1],
    p_sweep_d = Float64[0.0,0.0,0.0,0.0,0.0],
    p_twist_d = Float64[30.0, 30.0, 30.0, 30.0],
    p_pitch_d = 0.0,
    p_dihedral_d = fill(0.0,5),
    p_airfoilthickness = Float64[0.147538, 0.132307, 0.127404, 0.130604, 0.245836, 0.0629818],
    p_aerocenter = ones(p_radii)*0.25,
    p_x_offset = zeros(p_radii),
    # operating point design variables
    velocity = 15.40662113981876,
    # propulsion design variables
    p_Rtip = 2.25,
    blades = 2,
    numprops = 4,
    machtip = 0.8,
    p_rpm = 1020.45,
    kv = 4.271513933940295,
    i0 = 0.1,
    batt_mass = 1E8,

    w_n_linear = 20,
    w_halfspan = 10.0,
    w_webloc::Array{Float64,1}=[0.40],
    w_usedmaterials = ["highmodulus_uni","highmodulus_weave","taylor_foam"],
    w_nondim_halfspany = linspace(0,w_halfspan,length(w_chord)),
    w_lam_t = ones(length(w_nondim_halfspany)*length(w_usedmaterials)),
    w_chord = [0.5,0.5],
    w_aoa_d = 0.0,
    w_twist_d = [0.0,0.0],
    w_sweep_d = [0.0,0.0],
    w_dihedral_d = [0.0,0.0],
    w_orientation = [0.0,45.0,45.0], #uni, weave, webweave
    w_x_offset = ones(w_nondim_halfspany)*0.25,
    w_aerocenter = ones(w_nondim_halfspany)*0.275,

    n_CP = 60.0,
    payload = 684.0,
    max_mass = 1500.0,
    gravity = 9.81,
    h_set = 20.0, #meters
    specif_energy = 300.0, #wh/kg
    required_range = 150000.0, #km
    required_takeoff_dist = 395.0, #meters
    altitude_accel = 0.0,
    altitude_climb = 0.0,
    altitude_cruise = 3000.0,
    required_noise = 68.0,
    required_cruise = 28.0,

    )

    parameters = designparameters(p_n_akima,p_radii,p_webloc,altitude,plyprops,p_usedmaterials,p_orientation,p_lam_t,
    p_chord,p_sweep_d,p_twist_d,p_pitch_d,p_dihedral_d,p_airfoilthickness,p_aerocenter,p_x_offset,velocity,p_Rtip,
    blades,numprops,machtip,p_rpm,kv,i0,batt_mass,w_n_linear,w_webloc,w_lam_t,w_chord,w_aoa_d,w_twist_d,w_halfspan,w_nondim_halfspany,w_sweep_d,w_dihedral_d,w_usedmaterials,
    w_orientation,w_x_offset,w_aerocenter,n_CP,payload,max_mass,gravity,h_set,
    specif_energy,required_range,required_takeoff_dist,altitude_accel,altitude_climb,
    altitude_cruise,required_noise,required_cruise)
    return parameters
end


struct propwingdesign
    #Design Variables at this scope
    batt_mass::Float64
    w_aoa::Float64
    p_rpm::Float64
    p_pitch::Float64
    velocity::Float64

    #Design Parameters at this scope
    constraints
    af
    p_xaf
    p_yaf
    p_xafstrain
    p_yafstrain
    w_xaf
    w_yaf
    w_xafstrain
    w_yafstrain
    spline_3D_freestream
    spline_3D_blown

    # Design variables unique to the system function
    p_lam_tin
    p_orientation::Array{Float64,1}
    p_chord_pts::Array{Float64,1}
    p_twist_pts::Array{Float64,1}
    p_Rtip::Float64
    kv::Float64
    i0::Float64
    w_chord_pts::Array{Float64,1}
    w_twist_pts::Array{Float64,1}
    w_halfspan::Float64
    w_sweep_pts::Array{Float64,1}
    w_dihedral_pts::Array{Float64,1}
    w_lam_tin
    w_orientation::Array{Float64,1}
    w_airfoilthickness::Array{Float64,1}

    #Design variables unique to the system function
    etaprop_center::Array{Float64,1}
    p_n_akima::Int64
    p_radii::Array{Float64,1}
    altitude::Float64
    numprops::Int64
    p_x_offset::Array{Float64,1}
    p_aerocenter::Array{Float64,1}
    p_usedmaterials
    p_webloc::Array{Float64,1}
    plyprops
    blades::Float64
    w_n_linear::Int64
    w_x_offset::Array{Float64,1}
    w_aerocenter::Array{Float64,1}
    w_nondim_halfspany::Array{Float64,1}
    w_usedmaterials
    w_webloc::Array{Float64,1}
    prop_tilt::Float64

    #Control parameters
    detailed::Bool
    VTKfilename
    savefinalvtk::Bool
    verification::Bool
    plots::Bool
    n_CP::Int64

end

function propwingdesign(;
    #Design Variables at this scope
    batt_mass = 0.0,
    w_aoa = 0.0,
    p_rpm = 0.0,
    p_pitch = 0.0,
    velocity = 0.0,

    #Design Parameters at this scope
    constraints = [],
    af = [],
    p_xaf = [],
    p_yaf = [],
    p_xafstrain = [],
    p_yafstrain = [],
    w_xaf = [],
    w_yaf = [],
    w_xafstrain = [],
    w_yafstrain = [],
    spline_3D_freestream = [],
    spline_3D_blown = [],

    # Design variables unique to the system function
    p_lam_tin= [],
    p_orientation = zeros(2),
    p_chord_pts = zeros(2),
    p_twist_pts = zeros(2),
    p_Rtip =0.0,
    kv =0.0,
    i0 =0.0,
    w_chord_pts = zeros(2),
    w_twist_pts = zeros(2),
    w_halfspan =0.0,
    w_sweep_pts = zeros(2),
    w_dihedral_pts = zeros(2),
    w_lam_tin= [],
    w_orientation = zeros(2),
    w_airfoilthickness = zeros(2),

    #Design variables unique to the system function
    etaprop_center = zeros(2),
    p_n_akima = 10 ,
    p_radii = zeros(2),
    altitude =0.0,
    numprops = 10 ,
    p_x_offset = zeros(2),
    p_aerocenter = zeros(2),
    p_usedmaterials = [],
    p_webloc = zeros(2),
    plyprops = [],
    blades = 10 ,
    w_n_linear = 10 ,
    w_x_offset = zeros(2),
    w_aerocenter = zeros(2),
    w_nondim_halfspany = zeros(2),
    w_usedmaterials = [],
    w_webloc = zeros(2),
    prop_tilt = 0.0,

    #Control parameters
    detailed = false,
    VTKfilename = "notdefined.vtk",
    savefinalvtk = false,
    verification = false,
    plots = false,
    n_CP = 60,
    )

    return propwingdesign(batt_mass, w_aoa, p_rpm, p_pitch, velocity,
        constraints, af, p_xaf, p_yaf, p_xafstrain, p_yafstrain, w_xaf, w_yaf,
        w_xafstrain, w_yafstrain, spline_3D_freestream, spline_3D_blown,
        p_lam_tin, p_orientation, p_chord_pts, p_twist_pts, p_Rtip, kv, i0,
        w_chord_pts, w_twist_pts, w_halfspan, w_sweep_pts, w_dihedral_pts,
         w_lam_tin, w_orientation, w_airfoilthickness, etaprop_center,
        p_n_akima, p_radii, altitude, numprops, p_x_offset, p_aerocenter,
        p_usedmaterials, p_webloc, plyprops, blades, w_n_linear, w_x_offset,
        w_aerocenter, w_nondim_halfspany, w_usedmaterials, w_webloc, prop_tilt,
        detailed, VTKfilename, savefinalvtk, verification, plots, n_CP)
end

struct propwing_out
    #Propulsion
    ETA_total::Float64
    ETA_m::Float64
    ETA_p::Float64
    # Propeller
    CT::Float64
    CQ::Float64
    Np::Array{Float64,1}
    Tp::Array{Float64,1}
    p_alphas::Array{Float64,1}
    uprop::Array{Float64,1}
    vprop::Array{Float64,1}
    total_thrust::Float64
    advance_ratio::Float64
    torque::Float64
    db_test::Array{Float64,1}

    # Motor
    V_m::Float64
    I_m::Float64
    # Masses
    propwing_mass::Float64
    mass_prop::Float64
    mass_ESC::Float64
    mass_motor::Float64
    mass_wing::Float64

    # Wing
    panels#::Array{Float64,1}
    Lift::Float64
    Drag::Float64
    viscousDrag::Float64
    a::Float64
    Fp#::Array{Float64,1}
    cllocal::Array{Float64,1}
    cl::Array{Float64,1}
    Vinfeff::Array{Float64,1}
    alphaeff::Array{Float64,1}

    #Structures
    p_c_stress #::Array{Float64,2}
    w_c_stress #::Array{Float64,2}
    p_c_buckling #::Array{Float64,2}
    w_c_buckling #::Array{Float64,2}

    fail::Bool
end

function propwing_out(;
    #Propulsion
    ETA_total = 0.0,
    ETA_m = 0.0,
    ETA_p = 0.0,
    # Propeller
    CT = 0.0,
    CQ = 0.0,
    Np = zeros(2),
    Tp = zeros(2),
    p_alphas = zeros(2),
    uprop = zeros(2),
    vprop = zeros(2),
    total_thrust = 0.0,
    advance_ratio = 0.0,
    torque = 0.0,
    db_test = zeros(2),

    # Motor
    V_m = 0.0,
    I_m = 0.0,
    # Masses
    propwing_mass = 0.0,
    mass_prop = 0.0,
    mass_ESC = 0.0,
    mass_motor = 0.0,
    mass_wing = 0.0,

    # Wing
    panels = zeros(2),
    Lift = 0.0,
    Drag = 0.0,
    viscousDrag = 0.0,
    a = 0.0,
    Fp = [],# = zeros(2)
    cllocal = zeros(2),
    cl = zeros(2),
    Vinfeff = zeros(2),
    alphaeff = zeros(2),

    #Structures
    p_c_stress = zeros(2),
    w_c_stress = zeros(2),
    p_c_buckling = zeros(2),
    w_c_buckling = zeros(2),
    fail = false,
)

    return propwing_out(ETA_total,ETA_m,ETA_p,CT,CQ,Np,Tp,p_alphas,uprop,vprop,total_thrust,
    advance_ratio,torque,db_test,V_m,I_m,propwing_mass,mass_prop,mass_ESC,
    mass_motor,mass_wing,panels,Lift,Drag,viscousDrag,a,Fp,cllocal,cl,Vinfeff,alphaeff,
    p_c_stress,w_c_stress,p_c_buckling,w_c_buckling, fail)


end

struct system_output
    range_true::Float64
    dist_accel::Float64
    dist_climb::Float64

    energy_accel::Float64
    energy_climb::Float64
    energy_cruise::Float64

    t_accel::Float64
    t_transition::Float64
    t_climb::Float64
    t_cruise::Float64

    total_mass::Float64
    total_weight::Float64

    gama::Float64
    radius::Float64
    xtr::Float64
    ytr::Float64
    H_left::Float64
    x_left::Float64

    optconvals::Dict{Symbol,Array{Float64,1}}

    payload::Float64
    max_mass::Float64
    gravity::Float64
    h_set::Float64
    specif_energy::Float64
    required_takeoff_dist::Float64
    altitude_accel::Float64
    altitude_climb::Float64
    altitude_cruise::Float64
    required_noise::Float64
    etaprop_center::Array{Float64,1}
    batt_mass_accel::Float64
    batt_mass_climb::Float64
    batt_mass_cruise::Float64
    Fl::Float64
end

function system_output(;
    range_true = 0.0,
    dist_accel = 0.0,
    dist_climb = 0.0,

    energy_accel = 0.0,
    energy_climb = 0.0,
    energy_cruise = 0.0,

    t_accel = 0.0,
    t_transition = 0.0,
    t_climb = 0.0,
    t_cruise = 0.0,

    total_mass = 0.0,
    total_weight = 0.0,

    gama = 0.0,
    radius = 0.0,
    xtr = 0.0,
    ytr = 0.0,
    H_left = 0.0,
    x_left = 0.0,

    optconvals = Dict{Symbol,Array{Float64,1}}(),

    payload = 0.0,
    max_mass = 0.0,
    gravity = 0.0,
    h_set = 0.0,
    specif_energy = 0.0,
    required_takeoff_dist = 0.0,
    altitude_accel = 0.0,
    altitude_climb = 0.0,
    altitude_cruise = 0.0,
    required_noise = 0.0,
    etaprop_center = [0.0,0.0],
    batt_mass_accel = 0.0,
    batt_mass_climb = 0.0,
    batt_mass_cruise = 0.0,
    Fl = 0.0,
    )

    return system_output(range_true,dist_accel,dist_climb, energy_accel,
    energy_climb,energy_climb, t_accel,t_transition,t_climb,t_cruise,
    total_mass,total_weight, gama,radius,xtr,ytr,H_left,x_left, optconvals,
    payload,max_mass,gravity,h_set,specif_energy,required_takeoff_dist,
    altitude_accel,altitude_climb,altitude_cruise,required_noise,
    etaprop_center,batt_mass_accel,batt_mass_climb,batt_mass_cruise,Fl)

end

# This file sets up default values and bounds for design variables
function getoptxbounds!(paramdict,designparams)

    ##############################################################################
    #----------------------- Prop Design Variables ------------------------#
    ##############################################################################

    lb_scale = 1.0 #1E-6
    ub_scale = 1E3
    # Skin Fabric Ply Thickness Design Variables
    p_nstations = length(designparams.p_radii)
    matnames = designparams.plyprops.names
    p_nplies = length(designparams.p_usedmaterials)

    name = :p_lam_t
    default = zeros(designparams.p_lam_t)
    lowerbound = zeros(designparams.p_lam_t)
    upperbound = zeros(designparams.p_lam_t)
    scaling = zeros(designparams.p_lam_t)
    for i = 1:p_nplies
        idx = find(matnames -> matnames == designparams.p_usedmaterials[i],matnames)
        material = designparams.plyprops.plies[idx[1]]

        default1 = designparams.p_lam_t[(i-1)*p_nstations+1:(i)*p_nstations]
        lowerbound1 = fill(designparams.plyprops.plies[idx[1]].t*lb_scale,p_nstations)
        upperbound1 = fill(designparams.plyprops.plies[idx[1]].t*ub_scale,p_nstations)
        scaling1 = fill(0.01/designparams.plyprops.plies[idx[1]].t,p_nstations)

        default[(i-1)*p_nstations+1:(i)*p_nstations] = default1
        lowerbound[(i-1)*p_nstations+1:(i)*p_nstations] = lowerbound1
        upperbound[(i-1)*p_nstations+1:(i)*p_nstations] = upperbound1
        scaling[(i-1)*p_nstations+1:(i)*p_nstations] = scaling1

    end
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)


    name = :p_orientation
    default = designparams.p_orientation
    lowerbound = fill(-90.0,length(default))
    upperbound = fill(90.0,length(default))
    scaling = fill(1.0./90.0,length(default))
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :p_chord
    default = designparams.p_chord
    lowerbound = fill(0.01,p_nstations)
    lowerbound[end] = 0.001 #tip, works well for splines
    upperbound = fill(0.8,p_nstations)
    scaling = fill(1E2,p_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :p_sweep_d
    default = designparams.p_sweep_d
    lowerbound = fill(0.0,p_nstations-1)
    upperbound = fill(45.0,p_nstations-1)
    scaling = fill(pi/180.0,p_nstations-1)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    OptParams.addparam!(paramdict,:p_pitch_d_accel, designparams.p_pitch_d*10.0, -25.0, 45.0, 1E1)
    OptParams.addparam!(paramdict,:p_pitch_d_climb, designparams.p_pitch_d*10.0, -25.0, 45.0, 1E0)
    OptParams.addparam!(paramdict,:p_pitch_d_cruise, designparams.p_pitch_d*-2, -5.0, 80.0, 1E2)

    name = :p_twist_d
    default = designparams.p_twist_d
    lowerbound = fill(-0.0,p_nstations)
    upperbound = fill(45.0,p_nstations)
    scaling = fill(0.1,p_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :p_dihedral_d
    default = designparams.p_dihedral_d
    lowerbound = fill(0.0,p_nstations-1)
    upperbound = fill(10.0,p_nstations-1)
    scaling = fill(pi/180.0,p_nstations-1)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :p_airfoilthickness
    default = designparams.p_airfoilthickness
    lowerbound = fill(0.050,p_nstations)
    upperbound = fill(0.250,p_nstations)
    scaling = fill(10.0,p_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)


    OptParams.addparam!(paramdict,:p_Rtip,designparams.p_Rtip,0.1,5.0,1E-1)

    ##############################################################################
    #----------------------- Wing Design Variables ------------------------#
    ##############################################################################

    # Skin Fabric Ply Thickness Design Variables
    w_nstations = length(designparams.w_chord)

    name = :w_lam_t
    default = zeros(designparams.w_lam_t)
    lowerbound = zeros(designparams.w_lam_t)
    upperbound = zeros(designparams.w_lam_t)
    scaling = zeros(designparams.w_lam_t)
    for i = 1:length(designparams.w_usedmaterials)
        idx = find(matnames -> matnames == designparams.w_usedmaterials[i],matnames)
        material = designparams.plyprops.plies[idx[1]]

        default1 = designparams.w_lam_t[(i-1)*w_nstations+1:(i)*w_nstations]
        lowerbound1 = fill(designparams.plyprops.plies[idx[1]].t*lb_scale,w_nstations)
        upperbound1 = fill(designparams.plyprops.plies[idx[1]].t*ub_scale,w_nstations)
        scaling1 = fill(0.01/designparams.plyprops.plies[idx[1]].t,w_nstations)

        default[(i-1)*w_nstations+1:(i)*w_nstations] = default1
        lowerbound[(i-1)*w_nstations+1:(i)*w_nstations] = lowerbound1
        upperbound[(i-1)*w_nstations+1:(i)*w_nstations] = upperbound1
        scaling[(i-1)*w_nstations+1:(i)*w_nstations] = scaling1

    end
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :w_orientation
    default = designparams.w_orientation
    lowerbound = fill(-90.0,length(default))
    upperbound = fill(90.0,length(default))
    scaling = fill(1.0./90.0,length(default))
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :w_chord
    default = designparams.w_chord
    lowerbound = fill(0.1,w_nstations)
    upperbound = fill(5.0,w_nstations)
    scaling = fill(100.0/designparams.w_halfspan,w_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    OptParams.addparam!(paramdict,:w_aoa_d_accel,designparams.w_aoa_d*1.1,2.0,40.0,1E-2)
    OptParams.addparam!(paramdict,:w_aoa_d_climb,designparams.w_aoa_d*0.5,5.0,40.0,1E0)
    OptParams.addparam!(paramdict,:w_aoa_d_cruise,designparams.w_aoa_d*0.5,0.20,20.0,1E0)

    name = :w_twist_d
    default = designparams.w_twist_d
    lowerbound = fill(-10.0,w_nstations)
    upperbound = fill(20.0,w_nstations)
    scaling = fill(10.0,w_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :w_halfspan
    default = designparams.w_halfspan
    lowerbound = 5.699
    upperbound = 12.0
    scaling = 1.0/designparams.w_halfspan
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :w_sweep_d
    default = designparams.w_sweep_d
    lowerbound = fill(0.0,w_nstations)
    upperbound = fill(20.0,w_nstations)
    scaling = fill(10.0,w_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :w_dihedral_d
    default = designparams.w_dihedral_d
    lowerbound = fill(0.0,w_nstations)
    upperbound = fill(20.0,w_nstations)
    scaling = fill(10.0,w_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)

    name = :w_airfoilthickness
    default = designparams.p_airfoilthickness
    lowerbound = fill(0.050,p_nstations)
    upperbound = fill(0.250,p_nstations)
    scaling = fill(10.0,p_nstations)
    OptParams.addparam!(paramdict,name,default,lowerbound,upperbound,scaling)




    ##############################################################################
    #-----------------  Operating Point Design Variables ------------------------#
    ##############################################################################

    OptParams.addparam!(paramdict,:velocity_accel,designparams.velocity/1.4,5.0,28.0,1E1)
    OptParams.addparam!(paramdict,:velocity_climb,designparams.velocity*1.0,5.0,60.0,1E-2)
    OptParams.addparam!(paramdict,:velocity_cruise,designparams.velocity*1.0,20.0,100.0,1E1)

    ##############################################################################
    #-------------------------  Propulsion Design Variables ---------------------#
    ##############################################################################

    OptParams.addparam!(paramdict,:p_rpm_accel,designparams.p_rpm*4.5,200.0,40000.0,1E-1)
    OptParams.addparam!(paramdict,:p_rpm_climb,designparams.p_rpm*4.6,200.0,40000.0,1E-2)
    OptParams.addparam!(paramdict,:p_rpm_cruise,designparams.p_rpm*2.5,100.0,10000.0,1E0)
    OptParams.addparam!(paramdict,:kv,designparams.kv,5.0,5000.0,1E-0)
    OptParams.addparam!(paramdict,:i0,designparams.i0,0.1,6.0,1E0)
    OptParams.addparam!(paramdict,:batt_mass_accel,designparams.batt_mass,0.001,100.0,1E0)
    OptParams.addparam!(paramdict,:batt_mass_climb,designparams.batt_mass,0.01,100.0,1E-1)
    OptParams.addparam!(paramdict,:batt_mass_cruise,designparams.batt_mass,5.0,1000.0,1E-1)
    OptParams.addparam!(paramdict,:number_blades,2.0,2.0,5.0,1E-0)
    OptParams.addparam!(paramdict,:radius,100.0,2.0,2000.0,1E-1)
    OptParams.addparam!(paramdict,:gama_d,30.0,0.1,89.9,1E-1)
end

#------ SET UP OPTIMIZATION PARAMETERS -------#

mutable struct controlparams
    objective::Symbol
    printfreq::Int
    savevtk::Int
    printobjective::Bool
    printdetailedoutput::Bool
    printconstraintviolations::Bool
end

function controlparams(;
    objective::Symbol = :mass,
    printfreq::Int = 1,
    savevtk::Int = 3,
    printobjective::Bool = true,
    printdetailedoutput::Bool = true,
    printconstraintviolations::Bool = true,
    )
    controlparams(objective,printfreq,savevtk,printobjective,printdetailedoutput,
    printconstraintviolations)
end

# Define design variables
function getoptparams()
    designvar =  [
    :p_lam_t, #1-12
    # :p_orientation,
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
    # :kv, #28
    # :i0, #29


    # :w_orientation,
    :w_chord, #30-32
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
    # :w_materialfailure,
    # :w_localbuckling,

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
    :p_alphas,
    ]

    # define run controls
    ctlparams = controlparams(
    objective= :range, #:mass,#:range #:ETA_total_cruise,
    printfreq = 50,
    savevtk = 500,
    printobjective=true,
    printdetailedoutput=true,
    printconstraintviolations=true,
    )

    return designvar, constraints, ctlparams
end


#------- MATERIAL PROPERTIES -------#
function plyproperties()
    names = ["highmodulus_uni",
    "standard_uni",
    "MR60H",
    "T3900_uni",
    "T700_uni",
    "ELT5500",
    "UDCarbon",
    "highmodulus_weave",
    "standard_weave",
    "T3900_weave",
    "T700_weave",
    "Gelcoat",
    "Triax",
    "Saertex",
    "taylor_foam",
    "SNL_foam"]

    e1  =zeros(16)
    e2  =zeros(16)
    g12 =zeros(16)
    anu =zeros(16)
    rho =zeros(16)
    xt  =zeros(16)
    xc  =zeros(16)
    yt  =zeros(16)
    yc  =zeros(16)
    s   =zeros(16)
    plythickness =zeros(16)

    # "highmodulus_uni"
    e1[1]  = 175.0e9
    e2[1]  = 8.0e9
    g12[1] = 5.0e9
    anu[1] = 0.30
    rho[1] = 1600.0
    xt[1]  = min(1.0,1355.656/1682.011)*1000e6 #mean -> A-basis
    xc[1]  = min(1.0,1103.943/1396.504)*850e6 #mean -> A-basis
    yt[1]  = min(1.0,39.226/52.975)*40e6 #mean -> A-basis
    yc[1]  = min(1.0,235.434/282.439)*200e6 #mean -> A-basis
    s[1]   = min(1.0,142.411/159.516)*60e6 #mean -> A-basis
    plythickness[1] = 0.152e-3
    # "standard_uni"
    e1[2]  = 135.0e9
    e2[2]  = 10.0e9
    g12[2] = 5.0e9
    anu[2] = 0.30
    rho[2] = 1600.0
    xt[2]  = min(1.0,1355.656/1682.011)*1500e6 #mean -> A-basis
    xc[2]  = min(1.0,1103.943/1396.504)*1200e6 #mean -> A-basis
    yt[2]  = min(1.0,39.226/52.975)*50e6 #mean -> A-basis
    yc[2]  = min(1.0,235.434/282.439)*250e6 #mean -> A-basis
    s[2]   = min(1.0,142.411/159.516)*70e6 #mean -> A-basis
    plythickness[2] = 0.152e-3
    # "MR60H"
    e1[3]  = (165e9+150e9)/2.0
    e2[3]  = 8.56e9
    g12[3] = 4.39e9
    anu[3] = 0.326
    rho[3] = 1810.0
    xt[3]  = min(1.0,1355.656/1682.011)*3190e6 #mean -> A-basis
    xc[3]  = min(1.0,1103.943/1396.504)*1440e6 #mean -> A-basis
    yt[3]  = min(1.0,39.226/52.975)*82.0e6 #mean -> A-basis
    yc[3]  = min(1.0,235.434/282.439)*200.0e6 #mean -> A-basis
    s[3]   = min(1.0,142.411/159.516)*141e6 #mean -> A-basis
    plythickness[3] = 0.152e-3
    # "T3900_uni"
    e1[4]  = (148e9+131e9)/2.0
    e2[4]  = (9.7e9+9.7e9)/2.0
    g12[4] = 4.83e9
    anu[4] = 0.33
    rho[4] = 1573.0
    xt[4]  = min(1.0,1355.656/1682.011)*2830e6 #CTD mean -> A-basis
    xc[4]  = min(1.0,1103.943/1396.504)*1772e6 #CTD mean -> A-basis
    yt[4]  = min(1.0,39.226/52.975)*56.9e6 #CTD mean -> A-basis
    yc[4]  = min(1.0,235.434/282.439)*303e6 #CTD mean -> A-basis
    s[4]   = min(1.0,142.411/159.516)*89.6e6 #CTD mean -> A-basis
    plythickness[4] = 0.191e-3
    # "T700_uni"
    e1[5]=120.8e9
    e2[5]=11.57e9
    g12[5]=5.219e9
    anu[5]=0.350
    rho[5]=1525.0
    xt[5]=1356e6
    xc[5]=1104e6
    yt[5]=39.23e6
    yc[5]=235.4e6
    s[5]=142.4e6
    plythickness[5] = 0.152e-3
    # "ELT5500"
    e1[6]=41.8e9
    e2[6]=14.00e9
    g12[6]=2.630e9
    anu[6]=0.280
    rho[6]=1920.0
    xt[6]=972.0e6
    xc[6]=702.0e6
    yt[6]=100.0e6 #made up
    yc[6]=100.0e6 #made up
    s[6]=100.0e6 #made up
    plythickness[6] = 0.47e-3
    # "UDCarbon"
    e1[7]=114.5e9
    e2[7]=8.39e9
    g12[7]=5.990e9
    anu[7]=0.270
    rho[7]=1220.0
    xt[7]=1546.0e6
    xc[7]=1047.0e6
    yt[7]=100.0e6 #made up
    yc[7]=100.0e6 #made up
    s[7]=100.0e6 #made up
    plythickness[7] = 0.47e-3

    # FABRICS
    # "highmodulus_weave"
    e1[8]  = 85.0e9
    e2[8]  = 85.0e9
    g12[8] = 5.0e9
    anu[8] = 0.10
    rho[8] = 1600.0
    xt[8]  = min(1.0,1355.656/1682.011)*350e6 #mean -> A-basis
    xc[8]  = min(1.0,1103.943/1396.504)*150e6 #mean -> A-basis
    yt[8]  = min(1.0,39.226/52.975)*350e6 #mean -> A-basis
    yc[8]  = min(1.0,235.434/282.439)*150e6 #mean -> A-basis
    s[8]   = min(1.0,142.411/159.516)*35e6 #mean -> A-basis
    plythickness[8] = 0.218e-3
    # "standard_weave"
    e1[9]  = 70.0e9
    e2[9]  = 70.0e9
    g12[9] = 5.0e9
    anu[9] = 0.10
    rho[9] = 1600.0
    xt[9]  = min(1.0,1355.656/1682.011)*600e6 #mean -> A-basis
    xc[9]  = min(1.0,1103.943/1396.504)*570e6 #mean -> A-basis
    yt[9]  = min(1.0,39.226/52.975)*600e6 #mean -> A-basis
    yc[9]  = min(1.0,235.434/282.439)*570e6 #mean -> A-basis
    s[9]   = min(1.0,142.411/159.516)*90e6 #mean -> A-basis
    plythickness[9] = 0.218e-3
    # "T3900_weave"
    e1[10]  = (70.3e9+71.0e9)/2.0
    e2[10]  = (68.9e9+67.6e9)/2.0
    g12[10] = 4.6e9
    anu[10] = 0.032
    rho[10] = 1551.0
    xt[10]  = min(1.0,701.302/803.236)*1055e6 #CTD mean -> A-basis
    xc[10]  = min(1.0,549.748/749.955)*676e6 #CTD mean -> A-basis
    yt[10]  = min(1.0,557.575/722.602)*945e6 #CTD mean -> A-basis
    yc[10]  = min(1.0,604.067/741.866)*614e6 #CTD mean -> A-basis
    s[10]   = min(1.0,138.440/154.888)*79.3e6 #CTD mean -> A-basis
    plythickness[10] = 0.218e-3
    # "T700_weave"
    e1[11]=55.82e9
    e2[11]=52.10e9
    g12[11]=4.295e9
    anu[11]=0.085
    rho[11]=1501.0
    xt[11]=701.32e6
    yt[11]=557.59e6
    xc[11]=549.77e6
    yc[11]=604.08e6
    s[11]=138.44e6
    plythickness[11] = 0.218e-3
    # "Gelcoat"
    e1[12]=3.44e9
    e2[12]=3.44e9
    g12[12]=1.38e9
    anu[12]=0.3
    rho[12]=1235.0
    xt[12]=100.0e6 #made up
    xc[12]=100.0e6 #made up
    yt[12]=100.0e6 #made up
    yc[12]=100.0e6 #made up
    s[12]=100.0e6 #made up
    plythickness[12] = 0.05e-3
    # "Triax"
    e1[13]=27.7e9
    e2[13]=13.65e9
    g12[13]=7.2e9
    anu[13]=0.39
    rho[13]=1850.0
    xt[13]=700.0e6
    xc[13]=100.0e6 #made up
    yt[13]=700.0e6
    yc[13]=100.0e6 #made up
    s[13]=100.0e6 #made up
    plythickness[13] = 0.94e-3
    # "Saertex"
    e1[14]=13.6e9
    e2[14]=13.3e9
    g12[14]=11.8e9
    anu[14]=0.49
    rho[14]=1780.0
    xt[14]=144.0e6
    xc[14]=213.0e6
    yt[14]=144.0e6
    yc[14]=213.0e6
    s[14]=100.0e6 #made up
    plythickness[14] = 1.0e-3

    #FOAMS
    # "taylor_foam"
    e1[15]=48.0e6
    e2[15]=48.0e6
    g12[15]=28.0e6
    anu[15]=0.3
    rho[15]=75.0
    xt[15]=100.0e6 #made up
    yt[15]=100.0e6 #made up
    xc[15]=100.0e6 #made up
    yc[15]=100.0e6 #made up
    s[15]=100.0e6 #made up
    plythickness[15]=1.0E-3
    # "SNL_foam"
    e1[16]=256.0e6
    e2[16]=256.0e6
    g12[16]=22.0e6
    anu[16]=0.3
    rho[16]=200.0
    xt[16]=100.0e6 #made up
    yt[16]=100.0e6 #made up
    xc[16]=100.0e6 #made up
    yc[16]=100.0e6 #made up
    s[16]=100.0e6 #made up
    plythickness[16]=1.0E-3

    return plyproperties(names,Composites.material(e1,e2,g12,anu,rho,xt,xc,yt,yc,s,plythickness))
end
