 # && (nprocs()==1 || (myid()==2))
function objective(x,optxbounds,rangevar,designparams,designvars,constraints,
    ctlparams,af,p_xaf,p_yaf,p_xafstrain,p_yafstrain,w_xaf,w_yaf,w_xafstrain,
    w_yafstrain,spline_3D_freestream,spline_3D_blown; detailed = false,
    savefinalvtk = false, verification = false, plots = false)

    optconvals = Dict{Symbol,Array{Float64,1}}() # dictionary for constraints
    fail = false #hasn't failed yet

    #Unpack multipoint variables
    batt_mass_accel = OptParams.getvar(:batt_mass_accel,x,optxbounds,rangevar)
    batt_mass_climb = OptParams.getvar(:batt_mass_climb,x,optxbounds,rangevar)
    batt_mass_cruise = OptParams.getvar(:batt_mass_cruise,x,optxbounds,rangevar)
    w_aoa_accel = OptParams.getvar(:w_aoa_d_accel,x,optxbounds,rangevar)*pi/180
    w_aoa_climb = OptParams.getvar(:w_aoa_d_climb,x,optxbounds,rangevar)*pi/180
    w_aoa_cruise = OptParams.getvar(:w_aoa_d_cruise,x,optxbounds,rangevar)*pi/180
    p_rpm_accel = OptParams.getvar(:p_rpm_accel,x,optxbounds,rangevar)
    p_rpm_climb = OptParams.getvar(:p_rpm_climb,x,optxbounds,rangevar)
    p_rpm_cruise = OptParams.getvar(:p_rpm_cruise,x,optxbounds,rangevar)
    p_pitch_accel = OptParams.getvar(:p_pitch_d_accel,x,optxbounds,rangevar)*pi/180
    p_pitch_climb = OptParams.getvar(:p_pitch_d_climb,x,optxbounds,rangevar)*pi/180
    p_pitch_cruise = OptParams.getvar(:p_pitch_d_cruise,x,optxbounds,rangevar)*pi/180
    velocity_accel = OptParams.getvar(:velocity_accel,x,optxbounds,rangevar)
    velocity_climb = OptParams.getvar(:velocity_climb,x,optxbounds,rangevar)
    velocity_cruise = OptParams.getvar(:velocity_cruise,x,optxbounds,rangevar)
    number_blades = OptParams.getvar(:number_blades,x,optxbounds,rangevar)

    p_Rtip = OptParams.getvar(:p_Rtip,x,optxbounds,rangevar)

    p_lam_t = OptParams.getvar(:p_lam_t,x,optxbounds,rangevar)
    p_chord = OptParams.getvar(:p_chord,x,optxbounds,rangevar)
    w_halfspan = OptParams.getvar(:w_halfspan,x,optxbounds,rangevar)
    w_lam_t = OptParams.getvar(:w_lam_t,x,optxbounds,rangevar)
    w_chord = OptParams.getvar(:w_chord,x,optxbounds,rangevar)
    # radius = OptParams.getvar(:radius,x,optxbounds,rangevar)
    # gama = OptParams.getvar(:gama_d,x,optxbounds,rangevar)*pi/180

    global printiter
    printiter+=1

    n_CP = designparams.n_CP # 60
    payload = designparams.payload # 684.0
    max_mass = designparams.max_mass # 1500.0
    gravity = designparams.gravity # 9.81
    h_set = designparams.h_set # 20.0 #meters
    specif_energy = designparams.specif_energy # 300.0 #wh/kg
    required_range = designparams.required_range # 150000.0 #km
    required_takeoff_dist = designparams.required_takeoff_dist # 395.0 #meters
    altitude_accel = designparams.altitude_accel # 0.0
    altitude_climb = designparams.altitude_climb # 0.0
    altitude_cruise = designparams.altitude_cruise # 3000.0
    required_noise = designparams.required_noise
    required_cruise = designparams.required_cruise
    flap_zero_aoa = -14.0*pi/180

    p_spacing = w_halfspan*2/(designparams.numprops+1)
    etaprop_center = collect(linspace(-w_halfspan+p_spacing,w_halfspan-p_spacing,numprops))/w_halfspan
    p_spacing = etaprop_center[2]*w_halfspan - etaprop_center[1]*w_halfspan
    ##############################################################################
    #----------------------------- ACCELERATION --------------------------#
    ##############################################################################

    design_accel = propwingdesign(;
        #Design Variables at this scope
        batt_mass = batt_mass_accel,
        w_aoa = w_aoa_accel,
        p_rpm = p_rpm_accel,
        p_pitch = p_pitch_accel,
        velocity = velocity_accel,

        #Design Parameters at this scope
        constraints = constraints,
        af = af,
        p_xaf = p_xaf,
        p_yaf = p_yaf,
        p_xafstrain = p_xafstrain,
        p_yafstrain = p_yafstrain,
        w_xaf = w_xaf,
        w_yaf = w_yaf,
        w_xafstrain = w_xafstrain,
        w_yafstrain = w_yafstrain,
        spline_3D_freestream = spline_3D_freestream,
        spline_3D_blown = spline_3D_blown,

        # Design variables unique to the system function
        p_lam_tin = p_lam_t,
        p_orientation = OptParams.getvar(:p_orientation,x,optxbounds,rangevar),
        p_chord_pts = p_chord,
        p_twist_pts = OptParams.getvar(:p_twist_d,x,optxbounds,rangevar)*pi/180 + p_pitch_accel,
        p_Rtip = p_Rtip,
        kv = OptParams.getvar(:kv,x,optxbounds,rangevar),
        i0 = OptParams.getvar(:i0,x,optxbounds,rangevar),
        w_chord_pts = OptParams.getvar(:w_chord,x,optxbounds,rangevar),
        w_twist_pts = OptParams.getvar(:w_twist_d,x,optxbounds,rangevar)*pi/180+w_aoa_accel,
        w_halfspan = OptParams.getvar(:w_halfspan,x,optxbounds,rangevar),
        w_sweep_pts = OptParams.getvar(:w_sweep_d,x,optxbounds,rangevar)*pi/180,
        w_dihedral_pts = OptParams.getvar(:w_dihedral_d,x,optxbounds,rangevar)*pi/180,
        w_lam_tin = OptParams.getvar(:w_lam_t,x,optxbounds,rangevar),
        w_orientation = OptParams.getvar(:w_orientation,x,optxbounds,rangevar),
        w_airfoilthickness = OptParams.getvar(:w_airfoilthickness,x,optxbounds,rangevar),

        #Design parameters unique to the system function
        etaprop_center = etaprop_center,
        p_n_akima = designparams.p_n_akima,
        p_radii = designparams.p_radii,
        altitude = altitude_accel,
        numprops = designparams.numprops,
        p_x_offset = designparams.p_x_offset,
        p_aerocenter = designparams.p_aerocenter,
        p_usedmaterials = designparams.p_usedmaterials,
        p_webloc = designparams.p_webloc,
        plyprops = designparams.plyprops,
        blades = number_blades,
        w_n_linear = designparams.w_n_linear,
        w_x_offset = designparams.w_x_offset,
        w_aerocenter = designparams.w_aerocenter,
        w_nondim_halfspany = designparams.w_nondim_halfspany,
        w_usedmaterials = designparams.w_usedmaterials,
        w_webloc = designparams.w_webloc,
        prop_tilt = w_aoa_accel+flap_zero_aoa,

        #Control parameters
        detailed = detailed,
        VTKfilename = "./VTK_output/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))_testaccel",
        savefinalvtk = savefinalvtk,
        verification = verification,
        plots = plots,
        n_CP = n_CP,

    )
    if ((printiter%(ctlparams.printfreq)==0.0) || (printiter == 1)) && (nprocs()==1 || (myid()==2))  ; println("ACCELERATION"); end

    #Call the PropWing System
    detailed_output_accel = system(ctlparams,design_accel;run_noise=true)

    ETA_total_accel = detailed_output_accel.ETA_total
    thrust_accel = detailed_output_accel.total_thrust
    Lift_accel = detailed_output_accel.Lift
    Drag_accel = detailed_output_accel.Drag
    propwing_mass_accel = detailed_output_accel.propwing_mass

    a = detailed_output_accel.a
    db_test_accel = detailed_output_accel.db_test
    p_c_stress_accel = detailed_output_accel.p_c_stress
    w_c_stress_accel = detailed_output_accel.w_c_stress

    p_c_buckling_accel = detailed_output_accel.p_c_buckling
    w_c_buckling_accel = detailed_output_accel.w_c_buckling
    cllocal_accel = detailed_output_accel.cllocal


    ##############################################################################
    #----------------------------- TRANSITION/CLIMB --------------------------#
    ##############################################################################

    design_climb = propwingdesign(;
        #Design Variables at this scope
        batt_mass = batt_mass_climb,
        w_aoa = w_aoa_climb,
        p_rpm = p_rpm_climb,
        p_pitch = p_pitch_climb,
        velocity = velocity_climb,

        #Design Parameters at this scope
        constraints = constraints,
        af = af,
        p_xaf = p_xaf,
        p_yaf = p_yaf,
        p_xafstrain = p_xafstrain,
        p_yafstrain = p_yafstrain,
        w_xaf = w_xaf,
        w_yaf = w_yaf,
        w_xafstrain = w_xafstrain,
        w_yafstrain = w_yafstrain,
        spline_3D_freestream = spline_3D_freestream,
        spline_3D_blown = spline_3D_blown,

        # Design variables unique to the system function
        p_lam_tin = OptParams.getvar(:p_lam_t,x,optxbounds,rangevar),
        p_orientation = OptParams.getvar(:p_orientation,x,optxbounds,rangevar),
        p_chord_pts = OptParams.getvar(:p_chord,x,optxbounds,rangevar),
        p_twist_pts = OptParams.getvar(:p_twist_d,x,optxbounds,rangevar)*pi/180 + p_pitch_climb,
        p_Rtip = OptParams.getvar(:p_Rtip,x,optxbounds,rangevar),
        kv = OptParams.getvar(:kv,x,optxbounds,rangevar),
        i0 = OptParams.getvar(:i0,x,optxbounds,rangevar),
        w_chord_pts = OptParams.getvar(:w_chord,x,optxbounds,rangevar),
        w_twist_pts = OptParams.getvar(:w_twist_d,x,optxbounds,rangevar)*pi/180+w_aoa_climb,
        w_halfspan = OptParams.getvar(:w_halfspan,x,optxbounds,rangevar),
        w_sweep_pts = OptParams.getvar(:w_sweep_d,x,optxbounds,rangevar)*pi/180,
        w_dihedral_pts = OptParams.getvar(:w_dihedral_d,x,optxbounds,rangevar)*pi/180,
        w_lam_tin = OptParams.getvar(:w_lam_t,x,optxbounds,rangevar),
        w_orientation = OptParams.getvar(:w_orientation,x,optxbounds,rangevar),
        w_airfoilthickness = OptParams.getvar(:w_airfoilthickness,x,optxbounds,rangevar),

        #Design variables unique to the system function
        etaprop_center = etaprop_center,
        p_n_akima = designparams.p_n_akima,
        p_radii = designparams.p_radii,
        altitude = altitude_climb,
        numprops = designparams.numprops,
        p_x_offset = designparams.p_x_offset,
        p_aerocenter = designparams.p_aerocenter,
        p_usedmaterials = designparams.p_usedmaterials,
        p_webloc = designparams.p_webloc,
        plyprops = designparams.plyprops,
        blades = number_blades,
        w_n_linear = designparams.w_n_linear,
        w_x_offset = designparams.w_x_offset,
        w_aerocenter = designparams.w_aerocenter,
        w_nondim_halfspany = designparams.w_nondim_halfspany,
        w_usedmaterials = designparams.w_usedmaterials,
        w_webloc = designparams.w_webloc,
        prop_tilt = w_aoa_climb+flap_zero_aoa,

        #Control parameters
        detailed = detailed,
        VTKfilename = "./VTK_output/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))_testclimb",
        savefinalvtk = savefinalvtk,
        verification = verification,
        plots = plots,
        n_CP = n_CP,

    )
    if ((printiter%(ctlparams.printfreq)==0.0) || (printiter == 1)) && (nprocs()==1 || (myid()==2))  ; println("TRANSITION/CLIMB"); end

    #Call the PropWing System
    detailed_output_climb = system(ctlparams,design_climb;run_noise=true)

    ETA_total_climb = detailed_output_climb.ETA_total
    thrust_climb = detailed_output_climb.total_thrust
    Lift_climb = detailed_output_climb.Lift
    Drag_climb = detailed_output_climb.Drag
    propwing_mass_climb = detailed_output_climb.propwing_mass

    a = detailed_output_climb.a
    db_test_climb = detailed_output_climb.db_test
    p_c_stress_climb = detailed_output_climb.p_c_stress
    w_c_stress_climb = detailed_output_climb.w_c_stress

    p_c_buckling_climb = detailed_output_climb.p_c_buckling
    w_c_buckling_climb = detailed_output_climb.w_c_buckling
    cllocal_climb = detailed_output_climb.cllocal


    ##############################################################################
    #----------------------------- CRUISE --------------------------#
    ##############################################################################
    design_cruise = propwingdesign(;
        #Design Variables at this scope
        batt_mass = batt_mass_cruise,
        w_aoa = w_aoa_cruise,
        p_rpm = p_rpm_cruise,
        p_pitch = p_pitch_cruise,
        velocity = velocity_cruise,

        #Design Parameters at this scope
        constraints = constraints,
        af = af,
        p_xaf = p_xaf,
        p_yaf = p_yaf,
        p_xafstrain = p_xafstrain,
        p_yafstrain = p_yafstrain,
        w_xaf = w_xaf,
        w_yaf = w_yaf,
        w_xafstrain = w_xafstrain,
        w_yafstrain = w_yafstrain,
        spline_3D_freestream = spline_3D_freestream,
        spline_3D_blown = spline_3D_blown,

        # Design variables unique to the system function
        p_lam_tin = OptParams.getvar(:p_lam_t,x,optxbounds,rangevar),
        p_orientation = OptParams.getvar(:p_orientation,x,optxbounds,rangevar),
        p_chord_pts = OptParams.getvar(:p_chord,x,optxbounds,rangevar),
        p_twist_pts = OptParams.getvar(:p_twist_d,x,optxbounds,rangevar)*pi/180 + p_pitch_cruise,
        p_Rtip = OptParams.getvar(:p_Rtip,x,optxbounds,rangevar),
        kv = OptParams.getvar(:kv,x,optxbounds,rangevar),
        i0 = OptParams.getvar(:i0,x,optxbounds,rangevar),
        w_chord_pts = OptParams.getvar(:w_chord,x,optxbounds,rangevar),
        w_twist_pts = OptParams.getvar(:w_twist_d,x,optxbounds,rangevar)*pi/180+w_aoa_cruise,
        w_halfspan = OptParams.getvar(:w_halfspan,x,optxbounds,rangevar),
        w_sweep_pts = OptParams.getvar(:w_sweep_d,x,optxbounds,rangevar)*pi/180,
        w_dihedral_pts = OptParams.getvar(:w_dihedral_d,x,optxbounds,rangevar)*pi/180,
        w_lam_tin = OptParams.getvar(:w_lam_t,x,optxbounds,rangevar),
        w_orientation = OptParams.getvar(:w_orientation,x,optxbounds,rangevar),
        w_airfoilthickness = OptParams.getvar(:w_airfoilthickness,x,optxbounds,rangevar),

        #Design variables unique to the system function
        etaprop_center = etaprop_center,
        p_n_akima = designparams.p_n_akima,
        p_radii = designparams.p_radii,
        altitude = altitude_cruise,
        numprops = designparams.numprops,
        p_x_offset = designparams.p_x_offset,
        p_aerocenter = designparams.p_aerocenter,
        p_usedmaterials = designparams.p_usedmaterials,
        p_webloc = designparams.p_webloc,
        plyprops = designparams.plyprops,
        blades = number_blades,
        w_n_linear = designparams.w_n_linear,
        w_x_offset = designparams.w_x_offset,
        w_aerocenter = designparams.w_aerocenter,
        w_nondim_halfspany = designparams.w_nondim_halfspany,
        w_usedmaterials = designparams.w_usedmaterials,
        w_webloc = designparams.w_webloc,
        prop_tilt = w_aoa_cruise+flap_zero_aoa,

        #Control parameters
        detailed = detailed,
        VTKfilename = "./VTK_output/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))_testcruise",
        savefinalvtk = savefinalvtk,
        verification = verification,
        plots = plots,
        n_CP = n_CP,

    )


    if ((printiter%(ctlparams.printfreq)==0.0) || (printiter == 1)) && (nprocs()==1 || (myid()==2))  ; println("CRUISE"); end

    #Call the PropWing System
    detailed_output_cruise = system(ctlparams,design_cruise;run_noise=false)


    ETA_total_cruise = detailed_output_cruise.ETA_total
    thrust_cruise = detailed_output_cruise.total_thrust
    Lift_cruise = detailed_output_cruise.Lift
    Drag_cruise = detailed_output_cruise.Drag
    propwing_mass_cruise = detailed_output_cruise.propwing_mass

    a = detailed_output_cruise.a
    db_test_cruise = detailed_output_cruise.db_test
    p_c_stress_cruise = detailed_output_cruise.p_c_stress
    w_c_stress_cruise = detailed_output_cruise.w_c_stress

    p_c_buckling_cruise = detailed_output_cruise.p_c_buckling
    w_c_buckling_cruise = detailed_output_cruise.w_c_buckling
    cllocal_cruise = detailed_output_cruise.cllocal

    ##############################################################################
    #----------------------------- CHECK FOR FAILURE --------------------------#
    ##############################################################################

    # if detailed_output_accel.fail || detailed_output_climb.fail || detailed_output_cruise.fail; fail = true; end
    if detailed_output_cruise.fail; fail = true; end

    ##############################################################################
    #----------------------------- Flight Performance --------------------------#
    ##############################################################################

    total_batt_mass = (batt_mass_accel+batt_mass_climb+batt_mass_cruise)
    total_mass =total_batt_mass+payload + OptParams.softmax([propwing_mass_accel;propwing_mass_climb;propwing_mass_cruise])
    total_weight = total_mass*gravity

    #ACCELERATION
    p_accel = (thrust_accel-Drag_accel)*velocity_accel
    t_accel = total_weight*velocity_accel^2/(2*gravity*(p_accel/ETA_total_accel)*ETA_total_accel)
    dist_accel = total_weight*velocity_accel^3/(3*gravity*(p_accel/ETA_total_accel)*ETA_total_accel)
    energy_accel = t_accel * (p_accel)/ETA_total_accel #joules

    # energy_accel = batt_mass_accel*specif_energy*3600 #W*s
    # t_accel = energy_accel / ((p_accel)/ETA_total_accel)  # (W*s)/W
    # dist_accel = t_accel*velocity_accel*2/3

    # OLD TRANSITION AND CLIMB
    gama = pi/2
    if (thrust_climb*cos(w_aoa_climb+flap_zero_aoa)-Drag_climb)/total_weight > .9999
        gama = pi/2 #cap out at straight up to avoid numerical error
    elseif (thrust_climb*cos(w_aoa_climb+flap_zero_aoa)-Drag_climb)/total_weight < -.9999
        gama = -pi/2 #cap out at straight down to avoid numerical error
        fail = true
    else
        try
            gama = asin((thrust_climb*cos(w_aoa_climb+flap_zero_aoa)-Drag_climb)/total_weight)
        catch
            fail = true
            warn("flight path angle bad, thrust to weight: $(thrust_climb*cos(w_aoa_climb+flap_zero_aoa)-Drag_climb)")
        end
    end
    Fl = Lift_climb - total_weight + thrust_climb*sin(w_aoa_climb+flap_zero_aoa)
    radius = total_mass*velocity_climb^2/Fl
    arclength = radius*gama
    t_transition = arclength/velocity_climb
    xtr = radius*cos(gama-pi/2)
    ytr = radius*sin(gama-pi/2)+radius

    t_climb = 0.0
    if ytr > h_set #if we get to altitude during the climb
        H_left = h_set-ytr
        x_left = H_left/tan(gama)
        dist_climb = xtr
        energy_climb = t_transition * (thrust_climb*velocity_climb)/ETA_total_climb #joules
    else #include the steady state flight path angle to get to altitude
        H_left = h_set-ytr
        x_left = H_left/tan(gama)
        t_climb = sqrt(x_left^2+H_left^2)/velocity_climb
        t_climb += t_transition
        dist_climb = xtr+x_left
        energy_climb = t_climb * (thrust_climb*velocity_climb)/ETA_total_climb #joules
    end

    if gama<0.0 || radius<0.0; fail = true; end

    # # NEW TRANSITION AND CLIMB
    #
    # Thrust_climb_target = (total_weight*sin(gama)+Drag_climb)/cos(w_aoa_climb+flap_zero_aoa)
    #
    # Lift_transition_target =  total_mass*velocity_climb^2/radius + total_weight*cos(gama) - thrust_climb*sin(w_aoa_climb+flap_zero_aoa)
    # Fl = Lift_transition_target
    # arclength = radius*gama
    # t_transition = arclength/velocity_climb
    # xtr = radius*cos(gama-pi/2)
    # ytr = radius*sin(gama-pi/2)+radius
    #
    # t_climb = 0.0
    # if ytr > h_set #if we get to altitude during the climb
    #     H_left = h_set-ytr
    #     x_left = H_left/tan(gama)
    #     dist_climb = xtr
    #     energy_climb = t_transition * (thrust_climb*velocity_climb)/ETA_total_climb #joules
    # else #include the steady state flight path angle to get to altitude
    #     H_left = h_set-ytr
    #     x_left = H_left/tan(gama)
    #     t_climb = sqrt(x_left^2+H_left^2)/velocity_climb
    #     t_climb += t_transition
    #     dist_climb = xtr+x_left
    #     energy_climb = t_climb * (thrust_climb*velocity_climb)/ETA_total_climb #joules
    # end



    # if energy_accel<0.0; energy_accel = 0.0; end #remove incentive to push these in the negative
    # if energy_climb<0.0; energy_climb = 0.0; end
    #RANGE
    # capacity_left = (batt_mass*specif_energy*3600) - max(0.0,energy_accel) - max(0.0,energy_climb)
    # R = capacity_left/(thrust_cruise * velocity_cruise / ETA_total_cruise)*velocity_cruise #meters
    # range = batt_mass_cruise*specif_energy*3600 * ETA_total_cruise / thrust_cruise #equivalent to above

    range = batt_mass_cruise/10*ETA_total_cruise/(thrust_cruise/1000)
    range_true = batt_mass_cruise*specif_energy*3600*ETA_total_cruise/thrust_cruise
    rangeplot = batt_mass_cruise*specif_energy*3600*ETA_total_cruise*Lift_cruise/(Drag_cruise*total_weight)
    t_cruise = range_true/velocity_cruise #seconds
    energy_cruise = t_cruise * (thrust_cruise*velocity_cruise)/ETA_total_cruise #joules

    #system constraints
    if :Lift in constraints #doesn't matter for accelerations since it requires Lift
        optconvals[:Lift_accel] = [((total_weight-thrust_accel*sin(w_aoa_accel+flap_zero_aoa)) - Lift_accel)/10000.0]
        optconvals[:Lift_climb] = [((total_weight*cos(gama)-thrust_climb*sin(w_aoa_climb+flap_zero_aoa)) - Lift_climb)/10000.0] #Activate constraint before it dips
        optconvals[:Lift_cruise] = [((total_weight-thrust_cruise*sin(w_aoa_cruise)) - Lift_cruise)/10000.0]
    end

    if :range in constraints #doesn't matter for accelerations since it requires Lift
        range_constraint = required_range-range_true
        optconvals[:range] = [range_constraint]*1E-6#[log(abs(range_constraint+1.0))*sign(range_constraint)]
    end

    if :cruisespeed in constraints
        optconvals[:cruisespeed] = [required_cruise - velocity_cruise]*1E-3
    end

    if :takeoff_dist in constraints #doesn't matter for accelerations since it requires Lift
        takeoff_dist =  dist_climb+dist_accel #max(0.0,dist_accel)+max(0.0,dist_climb)
        takeoff_constraint = takeoff_dist-required_takeoff_dist

        optconvals[:takeoff_dist] = [takeoff_constraint]*1E-4 #[log(abs(takeoff_constraint+1.0))*sign(takeoff_constraint)]#[(max(0.0,dist_accel)+max(0.0,dist_climb) - required_takeoff_dist)/1000.0]

        # optconvals[:y_transition] = [(0 - ytr)/10.0]

        ####optconvals[:x_transition] = [(0 - xtr)/100.0]

        # optconvals[:x_left] = [(0.0 - 1/x_left)*100000.0]
    end

    if :radius in constraints
        optconvals[:positiveradius] = [(0 - 1/radius)*10000.0]
        optconvals[:radiussize] = [(radius-10000)/10000.0]
        # optconvals[:Lift_transition_target] = [Lift_transition_target - Lift_climb]*1E-3
    end

    if :flightpathangle in constraints
        optconvals[:flightpathangle] = [(0.1 - gama*180/pi)]*10.0 #this should ensure positive climb
        optconvals[:positive_dist_accel] = [(0 - log(abs(dist_accel+1.0))*sign(dist_accel))]
        # optconvals[:Thrust_climb_target] = [Thrust_climb_target-thrust_climb]*1E-3
    end


    if :capacity in constraints
        optconvals[:capacity_accel] = [energy_accel/(specif_energy*3600)-(batt_mass_accel)]/1E6
        optconvals[:capacity_climb] = [energy_climb/(specif_energy*3600)-(batt_mass_climb)]/1E6
    end


    #internal constraints
    if :w_materialfailure in constraints
        optconvals[:w_stress_accel] = w_c_stress_accel/100.0
        optconvals[:w_stress_climb] = w_c_stress_climb/100.0
        optconvals[:w_stress_cruise] = w_c_stress_cruise/100.0
    end

    if :p_materialfailure in constraints
        optconvals[:p_stress_accel] = p_c_stress_accel
        optconvals[:p_stress_climb] = p_c_stress_climb
        optconvals[:p_stress_cruise] = p_c_stress_cruise
    end

    if :w_localbuckling in constraints
        optconvals[:w_buckling_accel] = w_c_buckling_accel
        optconvals[:w_buckling_climb] = w_c_buckling_climb
        optconvals[:w_buckling_cruise] = w_c_buckling_cruise
    end

    if :p_localbuckling in constraints
        optconvals[:p_buckling_accel] = p_c_buckling_accel
        optconvals[:p_buckling_climb] = p_c_buckling_climb
        optconvals[:p_buckling_cruise] = p_c_buckling_cruise
    end

    if :machtip in constraints
        optconvals[:machtip_accel] = [(p_rpm_accel*pi/30.0*p_Rtip)/a - designparams.machtip]
        optconvals[:machtip_climb] = [(p_rpm_climb*pi/30.0*p_Rtip)/a - designparams.machtip]
        optconvals[:machtip_cruise] = [(p_rpm_cruise*pi/30.0*p_Rtip)/a - designparams.machtip]
    end

    if :Drag in constraints
        optconvals[:Drag_accel] = [(Drag_accel - thrust_accel*cos(w_aoa_accel+flap_zero_aoa))/1000.0] #numprops already included
        optconvals[:Drag_climb] = [(Drag_climb - thrust_climb*cos(w_aoa_climb+flap_zero_aoa))/1000.0] #numprops already included
        optconvals[:Drag_cruise] = [(Drag_cruise - thrust_cruise*cos(w_aoa_cruise))/100.0] #numprops already included
    end

    if :span in constraints
        optconvals[:span] = [(w_halfspan - 11.4/2)/100.0]
    end

    if :stall in constraints

        stallcon_accel = (cllocal_accel - 2.4)./1000
        stallcon_climb = (cllocal_climb - 2.4)./1000
        stallcon_cruise = (cllocal_cruise - 1.2)./1000

        #Ensure tips don't stall first

        stallcon_accel[end-6:end] = (cllocal_accel[end-6:end] - 2.0)./1000
        stallcon_climb[end-6:end] = (cllocal_climb[end-6:end] - 2.0)./1000
        stallcon_cruise[end-6:end] = (cllocal_cruise[end-6:end] - 1.0)./1000

        optconvals[:cl_max_accel] = stallcon_accel
        optconvals[:cl_max_climb] = stallcon_climb
        optconvals[:cl_max_cruise] = stallcon_cruise
    end

    if :noise in constraints
        optconvals[:noise_accel] = (db_test_accel - required_noise)*1E-2
        optconvals[:noise_climb] = (db_test_climb - required_noise)*1E-2
        # optconvals[:noise_cruise] = (db_test_cruise - required_noise)*1E-2
    end

    p_lam_t_pts = reshape(p_lam_t,length(designparams.p_radii),length(designparams.p_usedmaterials))

    if :p_chordlamthick in constraints
        c_lamt = zeros(p_lam_t_pts)
        for i = 1:length(p_lam_t_pts[:,1])
            c_lamt[i,:] = p_lam_t_pts[i,:]-p_chord[i]*0.1  # laminate thickness less than total airfoil thickness
        end
        optconvals[:p_chordlamthick] = reshape(c_lamt,prod(size(c_lamt)))
    end

    w_lam_t_pts = reshape(w_lam_t,length(designparams.w_nondim_halfspany),length(designparams.w_usedmaterials))

    if :w_chordlamthick in constraints
        c_lamt = zeros(w_lam_t_pts)
        for i = 1:length(w_lam_t_pts[:,1])
            c_lamt[i,:] = w_lam_t_pts[i,:]-w_chord[i]*0.1  # laminate thickness less than total airfoil thickness
        end
        optconvals[:w_chordlamthick] = reshape(c_lamt,prod(size(c_lamt)))
    end

    if :totalmass in constraints
        optconvals[:totalmass] = [(total_mass - max_mass)/1000.0]
    end

    if :p_alphas in constraints #Only constrain tip alpha
        optconvals[:p_max_alphas_accel] =  (detailed_output_accel.p_alphas - 19.0*pi/180)*1E1
        optconvals[:p_max_alphas_climb] = (detailed_output_climb.p_alphas-19.0*pi/180)*1E1
        optconvals[:p_max_alphas_cruise] = (detailed_output_cruise.p_alphas-19.0*pi/180)*1E1

        optconvals[:p_tipmin_alphas_accel] = [(0.1*pi/180 - detailed_output_accel.p_alphas[end-1])*1E0]
        optconvals[:p_tipmin_alphas_climb] = [(0.1*pi/180 - detailed_output_climb.p_alphas[end-1])*1E0]
        optconvals[:p_tipmin_alphas_cruise] = [(0.1*pi/180 - detailed_output_cruise.p_alphas[end-1])*1E0]

        optconvals[:p_min_alphas_accel] = (-5.0*pi/180 - detailed_output_accel.p_alphas)*1E1
        optconvals[:p_min_alphas_climb] = (-5.0*pi/180 - detailed_output_climb.p_alphas)*1E1
        optconvals[:p_min_alphas_cruise] = (-10.0*pi/180 - detailed_output_cruise.p_alphas)*1E1
    end

    if :p_separation in constraints #Only constrain tip alpha
        optconvals[:p_separation] =  [(p_Rtip - p_spacing/2.0)*10.0]
    end

    if :alphaeff in constraints #Only constrain tip alpha
        optconvals[:alphaeff_accel] = (detailed_output_accel.alphaeff*180/pi - 16.0)
        optconvals[:alphaeff_climb] = (detailed_output_climb.alphaeff*180/pi - 16.0)
        optconvals[:alphaeff_cruise] = (detailed_output_cruise.alphaeff*180/pi - 16.0)
    end


    ##############################################################################
    #------------------------------- ASSEMBLE OUTPUT ----------------------------#
    ##############################################################################

    c = Float64[]
    for con in keys(optconvals)
        c = append!(c,optconvals[con])
        # println("$con $(length(optconvals[con]))")
    end

    if ctlparams.objective == :totalmass
        f = (total_mass)/1000
    elseif ctlparams.objective == :range
        # f = -1/(1+e^(-range/1000000))

        f = -range #rangeplot/1E4#-range#-log(abs(range/1000.0+1.0))*sign(range) #
    elseif ctlparams.objective ==:ETA_total_cruise
        f = -ETA_total_cruise*10
    elseif ctlparams.objective ==:mass
        f = total_mass*1E-3
    elseif ctlparams.objective ==:constraints
        f = sum(c)*1E-2
    elseif ctlparams.objective ==:endurance
        f = -t_cruise/(24*60*60) # scale to be in days endurance
    elseif ctlparams.objective ==:takeoff
        f = (dist_climb+dist_accel)*1E-2
    elseif ctlparams.objective ==:work
        f = -(gama+1/radius)
    else
        error("objective key is not recognized, please define it or use a recognized one")
    end



    if fail
        warn("Objective/Constraint Function Failed")
    end

    if ((printiter%(35)==0.0) || (printiter == 1))  && (nprocs()==1 || (myid()==2))



        if ctlparams.printobjective
            println("------------OBJCON INPUT/OUTPUT------------")
            println("f: ",f)
            println("x: ",x)
            for var in designvars
                println(string(var)," = ",OptParams.getvar(var,x,optxbounds,rangevar))
            end

        end
        # print constraint violations
        if ctlparams.printconstraintviolations
            println("----------ACTIVE CONSTRAINTS----------")
            for con in keys(optconvals)
                printmaxviol(optconvals[con],con)
            end
        end

        # print some useful information
        if ctlparams.printdetailedoutput
            println("----------PARAMETERS OF INTEREST----------")
            println("
            range: $(round.(range,5))
            range_true: $(round.(range_true/1000,2))km
            rangeplot: $(round.(rangeplot/1000,2))km
            t_cruise: $(round(t_cruise/60,3)) min

            total_mass: $(round.(total_mass,2)) kg
            total_weight: $(round.(total_weight,2)) N

            dist_accel: $(round.(dist_accel,2)) m
            t_accel: $(round(t_accel,3))  s

            dist_climb: $(round.(dist_climb,2)) m
            t_transition: $(round(t_transition,3)) s
            t_climb: $(round(t_climb,3)) s


            flight path angle: $(gama*180/pi)
            radius: $(round.(radius,3))
            xtr: $(round.(xtr,3))
            ytr: $(round.(ytr,3))
            H_left: $(round.(H_left,3))
            x_left: $(round.(x_left,3))

            energy_accel%: $(round.(energy_accel/(total_batt_mass*specif_energy*3600)*100,2))%
            energy_climb%: $(round.(energy_climb/(total_batt_mass*specif_energy*3600)*100,2))%
            energy_cruise%: $(round.(energy_cruise/(total_batt_mass*specif_energy*3600)*100,2))%

            Propeller Angles of Attack
            accel: $(detailed_output_accel.p_alphas*180/pi)
            climb: $(detailed_output_climb.p_alphas*180/pi)
            cruise: $(detailed_output_cruise.p_alphas*180/pi)

            ")
        end
    end

    if ((printiter%(ctlparams.printfreq)==0.0) || (printiter == 1))  && (nprocs()==1 || (myid()==2))

        # ctlparams.printfreq = Int(round(ctlparams.printfreq*1.71))

        detailed_system = system_output(;
            range_true = range_true,
            dist_accel = dist_accel,
            dist_climb = dist_climb,

            energy_accel = energy_accel,
            energy_climb = energy_climb,
            energy_cruise = energy_cruise,

            t_accel = t_accel,
            t_transition = t_transition,
            t_climb = t_climb,
            t_cruise = t_cruise,

            total_mass = total_mass,
            total_weight = total_weight,

            gama = gama,
            radius = radius,
            xtr = xtr,
            ytr = ytr,
            H_left = H_left,
            x_left = x_left,

            optconvals = optconvals,

            payload = payload,
            max_mass = max_mass,
            gravity = gravity,
            h_set = h_set,
            specif_energy = specif_energy,
            required_takeoff_dist = required_takeoff_dist,
            altitude_accel = altitude_accel,
            altitude_climb = altitude_climb,
            altitude_cruise = altitude_cruise,
            required_noise = required_noise,
            etaprop_center = etaprop_center,
            batt_mass_accel = batt_mass_accel,
            batt_mass_climb = batt_mass_climb,
            batt_mass_cruise = batt_mass_cruise,
            Fl = Fl,
            )

        JLD.save("$(path)/OptSummaries/Journal_GIF/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))correct_span$(printiter).jld",
        "detailed_output_accel", detailed_output_accel,
        "detailed_output_climb", detailed_output_climb,
        "detailed_output_cruise", detailed_output_cruise,
        "design_accel", design_accel,
        "design_climb", design_climb,
        "design_cruise", design_cruise,
        "detailed_system", detailed_system,
        "objective", f,
        "constraints",c)
    end

    if savefinalvtk
        for var in designvars
            varnums = OptParams.getvar(var,x,optxbounds,rangevar)
            open("./SNOPT_out/p$(numprops)_cr$(round(Int,required_cruise/0.447))_db$(round(Int,required_noise))_design.txt","a")do x
                write(x,"$(string(var)) = $(varnums)
")
            end
        end
        println("Design File Written")
    end

    if detailed
        detailed_system = system_output(;
            range_true = range_true,
            dist_accel = dist_accel,
            dist_climb = dist_climb,

            energy_accel = energy_accel,
            energy_climb = energy_climb,
            energy_cruise = energy_cruise,

            t_accel = t_accel,
            t_transition = t_transition,
            t_climb = t_climb,
            t_cruise = t_cruise,

            total_mass = total_mass,
            total_weight = total_weight,

            gama = gama,
            radius = radius,
            xtr = xtr,
            ytr = ytr,
            H_left = H_left,
            x_left = x_left,

            optconvals = optconvals,

            payload = payload,
            max_mass = max_mass,
            gravity = gravity,
            h_set = h_set,
            specif_energy = specif_energy,
            required_takeoff_dist = required_takeoff_dist,
            altitude_accel = altitude_accel,
            altitude_climb = altitude_climb,
            altitude_cruise = altitude_cruise,
            required_noise = required_noise,
            etaprop_center = etaprop_center,
            batt_mass_accel = batt_mass_accel,
            batt_mass_climb = batt_mass_climb,
            batt_mass_cruise = batt_mass_cruise,
            Fl = Fl,

            )
        return detailed_output_accel,detailed_output_climb,detailed_output_cruise,
        design_accel,design_climb,design_cruise, detailed_system, f, c
    end

    return f,c,fail

end

function system(ctlparams,design;run_noise=true)

    fail = false #hasn't failed yet

    batt_mass = design.batt_mass
    w_aoa = design.w_aoa
    p_rpm = design.p_rpm
    p_pitch = design.p_pitch
    velocity = design.velocity
    constraints = design.constraints
    af = design.af
    p_xaf = design.p_xaf
    p_yaf = design.p_yaf
    p_xafstrain = design.p_xafstrain
    p_yafstrain = design.p_yafstrain
    w_xaf = design.w_xaf
    w_yaf = design.w_yaf
    w_xafstrain = design.w_xafstrain
    w_yafstrain = design.w_yafstrain
    spline_3D_freestream = design.spline_3D_freestream
    spline_3D_blown = design.spline_3D_blown
    altitude = design.altitude
    detailed = design.detailed
    VTKfilename = design.VTKfilename
    savefinalvtk = design.savefinalvtk
    verification = design.verification
    plots = design.plots


    #note p_ is for propeller, so it makes it easier for Taylor to incorporate this into his code

    rho,mu,a = Atmosphere.atmospherefit(altitude) #Use average altitude for each evaluation

    ##############################################################################
    #---------------------------------- PROPELLER ------------------------------#
    ##############################################################################

    # Unpack geometric design variables
    p_n_akima = design.p_n_akima
    p_radii_pts = design.p_radii

    p_radii = collect(linspace(p_radii_pts[1],p_radii_pts[end],p_n_akima))

    # operating point assumptions
    altitude = design.altitude
    numprops = design.numprops

    p_x_offset_pts = design.p_x_offset
    p_aerocenter_pts = design.p_aerocenter

    p_x_offset = zeros(p_n_akima)
    p_aerocenter = zeros(p_n_akima)

    for i = 1:p_n_akima
        p_x_offset[i] = linear_interp(p_radii_pts,p_x_offset_pts,p_radii[i])
        p_aerocenter[i] = linear_interp(p_radii_pts,p_aerocenter_pts,p_radii[i])
    end

    # Propeller Design
    p_lam_tin = design.p_lam_tin
    p_lam_t_pts = reshape(p_lam_tin,length(p_radii_pts),length(design.p_usedmaterials)) #map p_lam_t into matrix each row is a blade section, each column is a thickness
    p_orientation = design.p_orientation
    p_chord_pts = design.p_chord_pts
    p_sweep = zeros(p_n_akima)
    p_twist_pts = design.p_twist_pts
    p_dihedral = zeros(p_n_akima)
    p_airfoilthickness = ones(p_n_akima)
    p_Rtip = design.p_Rtip
    kv = design.kv
    i0 = design.i0

    # Apply akima splines to propeller
    p_chord = Akima.interp_noderiv(p_radii_pts, p_chord_pts, p_radii)
    p_chord[p_chord.<1E-3] = 1E-3
    p_twist = Akima.interp_noderiv(p_radii_pts, p_twist_pts, p_radii)
    p_twist[p_twist.<1E-3] = 1E-3
    p_lam_t = zeros(length(p_radii),length(design.p_usedmaterials))
    for i = 1:length(design.p_usedmaterials)
        p_lam_t[:,i] = Akima.interp_noderiv(p_radii_pts, p_lam_t_pts[:,i], p_radii)
    end
    p_lam_t[p_lam_t.<1E-3] = 1E-3

    #Translate prop sweep, dihedral into xyz locations
    p_sloc,p_xloc,p_yloc,p_zloc = spancoord(p_radii,p_Rtip,p_dihedral,p_sweep,p_twist,p_x_offset,p_chord)
    p_nchord = chordlengths(p_chord,p_sweep) #Normal chord needed for structures

    # Package into geometry object
    p_geom = geometry(p_radii,p_sloc,p_xloc,p_yloc,p_zloc,p_chord,p_nchord,
    p_twist,p_sweep,p_dihedral,p_airfoilthickness,p_xaf,p_yaf)



    ##############################################################################
    #---------------------------------- Propeller STRUCTURES ------------------------------#
    ##############################################################################

    p_struct = assemble_structures(p_layup,p_n_akima,p_twist,p_sloc,p_xaf,p_yaf,p_nchord,p_lam_t,design,p_x_offset,p_orientation,p_xloc,p_zloc,p_geom,design.p_webloc[1]) #only 1 we[1]) #only 1 web

    ##############################################################################
    #----------------------------- PROPELLER AERODYNAMICS -------------------------#
    ##############################################################################

    omega = p_rpm*pi/30
    R_station = p_radii*p_Rtip

    B = design.blades
    precone = 0.0
    xtowing = 0.2

    function propaero(omega;minimize=true)
        inflow = CCBlade.generalinflow(velocity, omega, R_station[2:end-1], precone, rho, mu, a;tilt=design.prop_tilt)
        rotor = CCBlade.Rotor(R_station[2:end-1], p_chord[2:end-1], p_twist[2:end-1], af, R_station[1], R_station[end], B, precone)

        # loads and thrust
        turbine = false
        Np, Tp, uprop, vprop, p_alphas = CCBlade.distributedloads(rotor, inflow, turbine)
        Tsub, Qsub = CCBlade.thrusttorqueintegrate(rotor, Np, Tp)
        thrust = B * Tsub
        torque = B * Qsub
        # thrust, torque = CCBlade.thrusttorque(rotor, [inflow], turbine)
        ETA_p, CT, CQ = CCBlade.nondim(thrust, torque, velocity, omega, rho, R_station[end], precone, turbine)

        abar = abs.(mean(uprop/velocity))  # an average induction used for contraction rate #TODO: if uprop/velocity is -1, the sqrt breaks, thrust constraint should keep uprop positive
        contraction_ratio = sqrt( (1 + abar)/ (1 + abar*(1 + xtowing/(sqrt(p_Rtip^2 + xtowing^2)))))
        reff = R_station*contraction_ratio #[Rhub; r; Rtip] * contraction_ratio
        uprop = [0; uprop; 0]
        vprop = [0; vprop; 0]

        if plots

            figname = "alpha_prop"
            figure(figname)
            plot(R_station,[0.0;p_alphas;0.0]*180/pi,".-",label = "Angle of Attack")
            plot(R_station,p_twist*180/pi,".-",label = "Geometric Twist (includes pitch)")
            xlabel("radius (m)")
            ylabel("AOA (deg)")
            legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # savefig("./figures/$figname.pdf",transparent = true)


        end

        # CP = torque*omega/(0.5*rho*pi*(p_Rtip)^2*velocity^3)
        TSR = omega * p_Rtip / velocity
        if minimize
            figure("output")
            plot(TSR,CP,".b")
            xlabel("TSR")
            ylabel("CP")
            return -CP
        else
            return Np, Tp, uprop, vprop, thrust, torque, ETA_p, CT, CQ, p_alphas
        end
    end

    if verification
        omega,CP = BrentMin.brent_min(0.0,30.0*2*pi/60,propaero;singleoutput=true)
        Np, Tp, uprop, vprop, thrust, torque, ETA_p, CT, CQ, p_alphas = propaero(omega;minimize=false)
        Np = -Np #Switch becuase turbine blade flies upside down
    else
        Np, Tp, uprop, vprop, thrust, torque, ETA_p, CT, CQ, p_alphas = propaero(omega;minimize=false)
    end

    Tp = -Tp #Switch because airfoil defined from leading to trailing by precomp and structures
    TSR = omega * p_Rtip / velocity
    RPMk = omega/(2*pi)*60
    advance_ratio = velocity/(omega/(2*pi)*p_Rtip*2)
    ETA_p = (thrust*velocity)/(omega*torque)

    ##############################################################################
    #----------------------------------- POWER ---------------------------------#
    ##############################################################################

    V_m,I_m,Po_m,ETA_m,mass_motor,Ro_m,Vmax_m = MotorPower.motorVI(torque,omega,kv,i0,numMotors=1,
    Case=1,mtype = "astroflight")

    mass_ESC, _ = MotorPower.ESC(V_m, I_m;esctype="AstroFlightHighVoltage")
    ETA_total = ETA_p*ETA_m#(thrust*velocity)/(V_m*I_m)
    ##############################################################################
    #----------------------------- PROP STRUCTURAL STRAINS -------------------------#
    ##############################################################################

    #TODO calculate centripetal force

    #-------- CALCULATE AERODYNAMIC CENTERS, MOMENT ARMS (AeroCenter to ShearCenter etc) -------#
    if maximum(p_xaf)>1.0 || maximum(w_xaf)>1.0
        warn("airfoil chord not of unit length, must be for proper structural calculations")
    end

    Fp = zeros(3,length(Np)+2)

    Fp[3,:] = [0.0;Np;0.0] #Lift is z in aero, y in beam
    Fp[1,:] = [0.0;Tp;0.0] #Drag is x in aero, z in beam
    # Fp[2,:] = #TODO: centripetal force #compression is y in aero, x in beam]

    p_strain, p_delta, p_freq, mass_prop,p_xafstrain,p_yafstrain = beam_strain_wrapper(Fp,p_n_akima,ctlparams,printiter,plots,p_xaf,p_yaf,p_xafstrain,p_yafstrain,p_twist,p_x_offset,p_aerocenter,p_nchord,p_struct,R_station,p_xloc,p_yloc,p_zloc)

    # calculate ply stress in prop
    p_shear = zeros(length(p_strain[:,1])) #TODO: add shear
    p_c_stress, p_stress, p_stresslayer = stress_wrapper(design.p_usedmaterials,design.plyprops,R_station,p_orientation,p_strain,p_shear)
    p_c_buckling = bucklingcon(p_struct,-p_strain,design.p_webloc[1],p_geom.normalchord)



    ##############################################################################
    #----------------------------- WING ANALYSIS -------------------------#
    ##############################################################################

    ##############################################################################
    #---------------------------------- Wing Design ------------------------------#
    ##############################################################################

    w_n_linear = design.w_n_linear
    w_x_offset_pts = design.w_x_offset
    w_aerocenter_pts = design.w_aerocenter

    w_nondim_halfspany_pts = design.w_nondim_halfspany
    w_nondim_halfspany = collect(linspace(w_nondim_halfspany_pts[1],w_nondim_halfspany_pts[end],w_n_linear))

    w_x_offset = zeros(w_n_linear)
    w_aerocenter = zeros(w_n_linear)

    for i = 1:w_n_linear
        w_x_offset[i] = linear_interp(w_nondim_halfspany_pts,w_x_offset_pts,w_nondim_halfspany[i])
        w_aerocenter[i] = linear_interp(w_nondim_halfspany_pts,w_aerocenter_pts,w_nondim_halfspany[i])
    end

    w_halfspan = design.w_halfspan
    w_chord_pts = design.w_chord_pts
    w_twist_pts = design.w_twist_pts
    w_sweep_pts = design.w_sweep_pts
    w_dihedral_pts = design.w_dihedral_pts
    w_lam_tin = design.w_lam_tin
    w_lam_t_pts = reshape(w_lam_tin,length(w_nondim_halfspany_pts),length(design.w_usedmaterials)) #map w_lam_t into matrix each row is a blade section, each column is a thickness
    w_orientation = design.w_orientation
    w_airfoilthickness = design.w_airfoilthickness

    # Apply linear interp to wing
    w_chord = linear_interp(w_nondim_halfspany_pts, w_chord_pts, w_nondim_halfspany)
    w_chord[w_chord.<1E-6] = 1E-6

    w_twist = linear_interp(w_nondim_halfspany_pts, w_twist_pts, w_nondim_halfspany)
    w_twist[w_twist.<1E-6] = 1E-6

    w_sweep = linear_interp(w_nondim_halfspany_pts, w_sweep_pts, w_nondim_halfspany)
    w_sweep[w_sweep.<1E-6] = 1E-6

    w_dihedral = linear_interp(w_nondim_halfspany_pts, w_dihedral_pts, w_nondim_halfspany)
    w_dihedral[w_dihedral.<1E-6] = 1E-6

    w_lam_t = zeros(length(w_nondim_halfspany),length(design.w_usedmaterials))
    for i = 1:length(design.w_usedmaterials)
        w_lam_t[:,i] = linear_interp(w_nondim_halfspany_pts, w_lam_t_pts[:,i], w_nondim_halfspany)
    end
    w_lam_t[w_lam_t.<1E-6] = 1E-6

    #Translate wing sweep, dihedral into xyz locations
    w_sloc,w_xloc,w_yloc,w_zloc = spancoord(w_nondim_halfspany,w_halfspan,w_dihedral,w_sweep,w_twist-w_aoa,0.0,w_chord)
    w_nchord = chordlengths(w_chord,w_sweep) #Normal chord needed for structures

    # Package into geometry object
    w_geom = geometry(w_nondim_halfspany,w_sloc,w_xloc,w_yloc,w_zloc,w_chord,w_nchord,
    w_twist,w_sweep,w_dihedral,w_airfoilthickness,w_xaf,w_yaf)

    ##############################################################################
    #---------------------------------- Wing Aero ------------------------------#
    ##############################################################################


    prop_center = design.etaprop_center*w_halfspan
    proploc = zeros(length(prop_center),3)
    for i = 1:length(prop_center)
        if prop_center[i]<0 sign=-1; else sign=1; end
        proploc[i,1] = linear_interp(w_yloc*sign,w_xloc,prop_center[i])
        proploc[i,2] = prop_center[i]
        proploc[i,3] = linear_interp(w_yloc*sign,w_zloc,prop_center[i])
    end

    #NOTE: in the structures frame, the forces are normal to the freestream, so the wing geometry includes aoa, but aerodynamically not the case including the prop on wing effect (unless prop is always normal to freestream).  Right now, the w_x y and z locations include the aoa, so slightly incorrect
    Lift, Drag, Fp, cllocal, cl, Vinfeff, alphaeff, viscousDrag, panels, fail = winganalysis(w_halfspan,w_chord, w_aoa, w_twist-w_aoa,w_sloc,w_xloc, w_yloc, w_zloc, rho, mu, a, velocity, R_station, uprop, vprop, proploc,spline_3D_freestream,spline_3D_blown;n_CP = design.n_CP, cw=1.0, plots=plots)

    ##############################################################################
    #---------------------------------- Wing STRUCTURES ------------------------------#
    ##############################################################################

    if  (:w_materialfailure in constraints) || (:w_localbuckling in constraints)
        w_struct= assemble_structures(w_layup,w_n_linear,w_twist,w_sloc,w_xaf,w_yaf,w_nchord,w_lam_t,design,w_x_offset,w_orientation,w_xloc,w_zloc,w_geom,design.w_webloc[1]) #only 1 web


        w_Fp = zeros(3,length(Fp[1,:])+1)
        w_Fp[:,1:length(Fp[1,:])] = Fp #TODO: this is taking the in-between load/length and applying it to the left and putting a 0 at the end
        w_Fp[1,:] = -w_Fp[1,:] #TODO: include props pulling on wings not here though since they aren't a distributed load, also any payload and motor/battery distributed load
        # w_Fp[2,:] = 0.0
        w_Fp[3,:] = -w_Fp[3,:] #TODO: VLM Frame Of Reference must be reverse?

        w_strain, w_delta, w_freq, mass_wing, w_xafstrain,w_yafstrain = beam_strain_wrapper(w_Fp,w_n_linear,ctlparams,printiter,plots,w_xaf,w_yaf,w_xafstrain,w_yafstrain,w_twist,w_x_offset,w_aerocenter,w_nchord,w_struct,w_nondim_halfspany*w_halfspan,w_xloc,w_yloc,w_zloc)

        # calculate ply stress in wing
        w_shear = zeros(length(w_strain[:,1])) #TODO: add shear
        w_c_stress, w_stress, w_stresslayer = stress_wrapper(design.w_usedmaterials,design.plyprops,w_nondim_halfspany*w_halfspan,w_orientation,w_strain,w_shear)

        # interpolate buckling strains to node locations

        w_c_buckling = bucklingcon(w_struct,-w_strain,design.w_webloc[1],w_geom.normalchord)
    else

            w_strainlocx = zeros(w_xafstrain)
            w_strainlocy = zeros(w_yafstrain)
            w_c_stress = zeros(w_xafstrain)
            w_stress = zeros(w_xafstrain)
            w_strain = zeros(w_xafstrain)
            w_c_buckling = zeros(w_xafstrain)
            mass_wing = 0.0
            w_stresslayer = zeros(2,length(w_xafstrain[:,1]),length(w_xafstrain[1,:]))
            rot = zeros(2,2)
            for i = 1:w_n_linear
                # find leading edge
                lei = indmin(abs.(w_xaf[i,:]))
                # shift so leading edge is first
                xpc = circshift(w_xaf[i,:],-(lei-1))
                ypc = circshift(w_yaf[i,:],-(lei-1))

                rot[1,:] = [cos(-w_twist[i]),-sin(-w_twist[i])]
                rot[2,:] = [sin(-w_twist[i]),cos(-w_twist[i])]

                xy = [w_xafstrain[i,:]*w_nchord[i],w_yafstrain[i,:]*w_nchord[i]]
                xyrot = rot*xy
                # determine locations at which to calculate strain
                w_strainlocx[i,:] = xyrot[1] +w_xloc[i]
                w_strainlocy[i,:] = xyrot[2] +w_zloc[i]
            end

    end

    ##############################################################################
    #----------------------------- PROP NOISE ANALYSIS -------------------------#
    ##############################################################################
    # Calculate propeller locations
    if numprops < 16
        n_db = 10
        db_test = zeros(n_db)
        yloc = linspace(10*p_Rtip,60*p_Rtip,n_db) # all props
    elseif numprops == 16
        n_db = 1
        db_test = zeros(n_db)
        yloc = 26.0*p_Rtip# linspace(24*p_Rtip,26*p_Rtip,n_db) # 16 props
    elseif numprops == 32
        n_db = 1
        db_test = zeros(n_db)
        yloc = 55.0*p_Rtip #linspace(50*p_Rtip,55*p_Rtip,n_db) #32 props
    end

    blades = design.blades
    blades = round(Int, B)

    if run_noise
        for i = 1:n_db
            x_test = proploc[:,2] # x-locations of turbines (m)
            y_test = ones(numprops)*0.0 # y-locations of turbines (m)
            obs_test = [0., yloc[i], 1.5] # x-, y-, and z-location of the observer (m)
            winddir_test = 0. # wind direction (deg)
            rpm_test = ones(numprops)*p_rpm # rotation rate of the tubrines (rpm)
            windvel_test = ones(numprops)*velocity # wind velocity (m/s)
            h_test = 15.24 # height of the turbine hub (m)
            noise_corr = 0.8697933840957954 # correction factor for noise

            alpha =  abs.([p_alphas[1];p_alphas;p_alphas[end]]*180/pi)
            AR = p_Rtip/mean(p_chord) #TODO: actual AR
            psi = 14.0 # solid angle (deg)

            db_test[i] = 0.0

            try

                db_test[i] = BPM.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, blades, h_test, R_station, p_chord, p_chord*0.25, alpha, mu/rho, a, psi, AR, noise_corr)
            catch
                # fail = true
                # warn("BPM failed prop angle of attack: $alpha")
                # println("prop tip mach number: $(maximum(rpm_test*(2*pi)/60*p_Rtip/a))")
            end
        end
    end

    ##############################################################################
    #---------------------------------- MASS ------------------------------#
    ##############################################################################


    propwing_mass = (mass_prop*design.blades + mass_ESC + mass_motor)*numprops + mass_wing*2.0


    ##############################################################################
    #----------------------------- PRINT RELEVANT OUTPUT ------------------------#
    ##############################################################################
    global printiter

    if ((printiter%(35)==0.0) || (printiter == 1))  && (nprocs()==1 || (myid()==2))
        # print objective and design variable values

        # print some useful information
        if ctlparams.printdetailedoutput
            println("----------PARAMETERS OF INTEREST----------")
            println("
    V_m: $(round.(V_m,3))
    I_m: $(round.(I_m,3))

    mass_prop: $(round.(mass_prop,3))
    mass_ESC: $(round.(mass_ESC,3))
    mass_motor: $(round.(mass_motor,3))
    mass_wing $(round.(mass_wing,3))

    ETA_m: $(round.(ETA_m,3))
    ETA_p: $(round.(ETA_p,3))
    ETA_total: $(round.(ETA_total,3))

    Np: $(round.(Np,3))
    Tp: $(round.(Tp,3))
    uprop: $(round.(uprop,3))
    vprop: $(round.(vprop,3))
    total thrust: $(round.(thrust*numprops, 3))
    advance_ratio: $(round.(advance_ratio,5))
    torque: $(round.(torque,3))
    Velocity: $(round.(velocity,3))

    db_test $(round.(db_test,3))

    Lift $(round.(Lift,3))
    Drag $(round.(Drag,3))
    viscousDrag $(round.(viscousDrag,3))

            ")

        end

        if plots
            figure("chords")
            plot(p_radii_pts, p_chord_pts,".",label = "original")
            plot(p_radii, p_chord, "-",label = "spline")
            xlabel("Radial Position (r/R)")
            ylabel("Chord (m)")
            legend(loc = "best")

            figure("twists")
            plot(p_radii_pts, p_twist_pts*180/pi,".",label = "original")
            plot(p_radii, p_twist*180/pi, "-",label = "spline")
            xlabel("Radial Position (r/R)")
            ylabel("Twist (deg)")
            legend(loc = "best")
        end

        for i = 1:p_n_akima
            # Calculate aerodynamic center from precalculated x-aero center and chord line
            idxle = indmin(abs.(p_xaf[i,:]))
            idxte = indmax(abs.(p_xaf[i,:]))

            yle = p_yaf[i,idxle]
            yte = p_yaf[i,idxte]

            p_acx = p_aerocenter[i]
            p_acy = (yte-yle)/(1)*p_acx + yle #y/mx+b

        end

    end

    p_c_bucklingplot = reshape(p_c_buckling,size(p_strain))
    w_c_bucklingplot = reshape(w_c_buckling,size(w_strain))

    if ((printiter%(ctlparams.savevtk)==0.0) || (printiter == 1) || (savefinalvtk == true))  && (nprocs()==1 || (myid()==2))

        # ctlparams.savevtk = Int(round(ctlparams.savevtk*1.71))

        save_vtk_opt(p_stress,p_strain,p_stresslayer,p_c_bucklingplot,p_xloc,p_yloc,p_zloc,numprops,w_stress,w_strain,w_stresslayer,w_c_bucklingplot,w_xloc,w_yloc,w_zloc;
        file_name=VTKfilename,
        p_radii = p_radii,
        chord = p_chord,
        p_twist_d = p_twist*180/pi,
        p_dihedral_d = [p_dihedral;p_dihedral[end]]*180/pi,
        p_sweep_d = [p_sweep;p_sweep[end]]*180/pi,
        Rtip = p_Rtip,
        p_airfoilx = p_xaf,
        p_airfoily = p_yaf,
        p_strainx = p_xafstrain,
        p_strainy = p_yafstrain,
        w_twist_d = (w_twist-w_aoa)*180/pi, #Plot level and have freestream change with aoa
        w_dihedral_d = w_dihedral*180/pi,
        w_sweep_d = w_sweep*180/pi,
        w_strainx = w_xafstrain,
        w_strainy = w_yafstrain,
        spanloc = w_sloc/w_sloc[end],
        halfspan = w_sloc[end],
        printiter = printiter,
        proploc = proploc, #xyz offset from port to starboard
        numblades = round(Int,blades),
        plots = plots )


        # println("
        # p_stress = $p_stress
        # p_strain = $p_strain
        # p_stresslayer = $p_stresslayer
        # p_c_bucklingplot = $p_c_bucklingplot
        # p_xloc = $p_xloc
        # p_yloc = $p_yloc
        # p_zloc = $p_zloc
        # numprops = $numprops
        # w_stress = $w_stress
        # w_strain = $w_strain
        # w_stresslayer = $w_stresslayer
        # w_c_bucklingplot = $w_c_bucklingplot
        # w_xloc = $w_xloc
        # w_yloc = $w_yloc
        # w_zloc = $w_zloc
        # file_name=$(VTKfilename)
        # p_radii = $(p_radii)
        # chord = $(p_chord)
        # p_twist_d = $(p_twist*180/pi)
        # p_dihedral_d = $([p_dihedral;p_dihedral[end]]*180/pi)
        # p_sweep_d = $([p_sweep;p_sweep[end]]*180/pi)
        # Rtip = $(p_Rtip)
        # p_airfoilx = $(p_xaf)
        # p_airfoily = $(p_yaf)
        # p_strainx = $(p_strainlocx)
        # p_strainy = $(p_strainlocy)
        # w_twist_d = $((w_twist-w_aoa)*180/pi)
        # w_dihedral_d = $(w_dihedral*180/pi)
        # w_sweep_d = $(w_sweep*180/pi)
        # w_strainx = $(w_strainlocx)
        # w_strainy = $(w_strainlocy)
        # spanloc = $(w_sloc/w_sloc[end])
        # halfspan = $(w_sloc[end])
        # printiter = $(printiter)
        # proploc = $(proploc)
        # plots = $(plots)
        # ")

    end



    p_lam_tsave = reshape(p_lam_t,prod(size(p_lam_t)))
    w_lam_tsave = reshape(w_lam_t,prod(size(w_lam_t)))

    finaloutput = designparameters(;

    p_orientation=p_orientation, #an p_orientation for each material at each section
    # structural design variables
    p_lam_t = p_lam_tsave,
    # geometric design variables
    p_chord = p_chord,
    p_sweep_d = p_sweep*180/pi,
    p_twist_d = p_twist*180/pi,
    p_pitch_d = p_pitch*180/pi,
    p_dihedral_d = p_dihedral*180/pi,
    p_airfoilthickness = p_airfoilthickness, #not used

    # operating point design variables
    velocity = velocity,
    # propulsion design variables
    p_Rtip = p_Rtip,

    p_rpm = p_rpm,
    kv = kv,
    i0 = i0,
    # Wing Design variables
    w_n_linear = w_n_linear,
    w_lam_t = w_lam_tsave,
    w_chord = w_chord,
    w_aoa_d = w_aoa*180/pi,
    w_twist_d = w_twist*180/pi,
    w_halfspan = w_halfspan,
    w_nondim_halfspany = w_nondim_halfspany,
    w_sweep_d = w_sweep*180/pi,
    w_dihedral_d = w_dihedral*180/pi,
    w_orientation = w_orientation,

    )

    detailed_output = propwing_out(;
    # Propulsion
    ETA_total = ETA_total,
    ETA_m = ETA_m,
    ETA_p = ETA_p,
    # Propeller
    CT = CT,
    CQ = CQ,
    Np = Np,
    Tp = Tp,
    p_alphas = p_alphas,
    uprop = uprop,
    vprop = vprop,
    total_thrust = thrust*numprops,
    advance_ratio = advance_ratio,
    torque = torque,
    db_test = db_test,

    # Motor
    V_m = V_m,
    I_m = I_m,
    # Masses
    propwing_mass = propwing_mass,
    mass_prop = mass_prop,
    mass_ESC = mass_ESC,
    mass_motor = mass_motor,
    mass_wing = mass_wing,

    # Wing
    panels = panels,
    Lift = Lift,
    Drag = Drag,
    viscousDrag = viscousDrag,
    a = a,
    Fp = Fp,
    cllocal = cllocal,
    cl = cl,
    Vinfeff = Vinfeff,
    alphaeff = alphaeff,

    #Structures
    p_c_stress = p_c_stress,
    w_c_stress = w_c_stress,
    p_c_buckling = p_c_buckling,
    w_c_buckling = w_c_buckling,

    fail = fail,

    )


    return detailed_output
end #objcon
