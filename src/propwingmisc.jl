
"""
viscous Drag: right now assumes a single aifoil across the wing.

Arguments:

Returns:

"""
function viscousdrag(spline_3D_freestream,spline_3D_blown,alpha_array,Re_array,
    M_array,chord_array,q_array,ds_array,blown_array)

    cdp = zeros(alpha_array)
    for i = 1:length(alpha_array)
        if blown_array[i]==0
            vars = [alpha_array[i],Re_array[i],M_array[i]]
            cdp[i] = AirfoilPrep.interpND(spline_3D_freestream[2],vars)
        elseif blown_array[i]==1
            vars = [alpha_array[i],Re_array[i],M_array[i]]
            cdp[i] = AirfoilPrep.interpND(spline_3D_blown[2],vars)
        else
            error("blown_array does not contain values of 0 for not blown and 1 for blown")
        end
    end

    Dp = sum(cdp.*q_array.*chord_array.*ds_array)*2 #VLM is in half span, so x2 for other half

    halfspan_dp = cdp.*q_array.*chord_array
    return Dp, halfspan_dp
end

"""helper function"""
function interp1(xpt, ypt, x)
    intf = Interpolations.interpolate((xpt,), ypt, Interpolations.Gridded(Interpolations.Linear()))
    y = zeros(x)
    idx = (x .> xpt[1]) .& (x.< xpt[end])
    y[idx] = intf[x[idx]]
    return y
end


"""
interpolate propwash onto wing

Dnew: desired diameter (rounding up or down to match panel spacing)
ypropc: y-location for center of propeller
cw: rotating clockwise
"""
function project_wake_onto_wing(Dnew, ypropc, R_station, uprop, vprop, cw, CPy)

    # prop velocities as function of y (get both sides of prop)
    Dprop = R_station[end]*2.0
    yprop = [ypropc - reverse(R_station)*Dnew/Dprop; ypropc + R_station*Dnew/Dprop]   # adjust width to panel size
    Vxprop = [2*reverse(uprop); 2*uprop]  # factor of two from far field
    Vzprop = cw * [reverse(vprop); -vprop]  # factor of 2 from far-field and factor of 0.5 from Veldhuis cancel out, times by negative 1 if counter-clockwise

    fail = false
    Vxext = zeros(CPy)
    Vzext = zeros(CPy)
    # interpolate onto wing
    try
        Vxext = interp1(yprop, Vxprop, CPy)
        Vzext = interp1(yprop, Vzprop, CPy)
    catch
        fail = true
        warn("prop on wing failed: yprop: $yprop
        Vxprop: $Vxprop
        CPy: $CPy")
    end
    # Vyext = zeros(CP.y)

    return Vxext, Vzext, fail
end

# """
# wing analysis # with newer VLM
# """
# function winganalysis(w_halfspan,w_chord, w_aoa, w_twist,w_sloc,w_xloc, w_yloc, w_zloc, rho, mu, a, velocity, R_station, uprop, vprop, proploc,spline_3D_freestream,spline_3D_blown;n_CP = 60, cw=1, plots = false)
#     #Setup and run VLM
#     center_chord_array = (w_chord[2:end]+w_chord[1:end-1])/2.0
#     center_w_yloc = cumsum(w_yloc[2:end]-w_yloc[1:end-1])
#
#     ds_array = w_sloc[2:end]-w_sloc[1:end-1]
#     beta = 0.0
#     Omega = [0.0; 0.0; 0.0]
#
#     npanels = round.(Int,ones(length(w_chord)-1)*n_CP/length(w_chord))
#     duplicate = false
#     spacing = "uniform"
#     panels = VLM.linearsections(w_xloc, w_yloc, w_zloc, w_chord, w_twist, npanels, duplicate, spacing)
#     ds_array_VLM = ones(sum(npanels))*ds_array[1]*length(center_w_yloc)/sum(npanels)  #ASSUMES LINEAR SPACING
#
#     Sref = sum(center_chord_array.*ds_array)*2.0
#     cref = 2.0/Sref*sum(center_chord_array.^2.0.*ds_array)
#     bref = w_halfspan*2
#     rcg = [0.50, 0.0, 0.0]
#     ref = VLM.Reference(Sref, cref, bref, rcg)
#
#     # setup
#
#     # use only starboard props since symmetric
#     yprop_center = proploc[div(end,2)+1:end,2]
#     # # -------- leave commented out: only used to show how the default method of wake control point interaction works poorly -----------
#     # Vx, Vz = project_wake_onto_wing(2*R_station[end], yprop_center[1], R_station, uprop, vprop, cw, CP.y)
#     # Vext = VLM.velocity(Vx, zeros(CP.y), Vz)
#     # L, Di, cl, cllocal, Vinfeff = VLM.run(wing, fs, Vext)
#     # return L, Di, CP, cl, cllocal, Vinfeff
#     # # ----------------------------------------------
#
#     # ---------- interpolate prop wakes onto wing ------------
#
#     # arbitrary decision: could center on control points (diameters are odd number of panel width) or on QC lines (diameters are even number of panel width).
#     # I chose center of control points.
#
#     #Get control point center locations
#     N = length(panels)
#     CPy = zeros(N)
#     for i = 1:N  # CP
#         rcp = VLM.control_point(panels[i])
#         CPy[i] = rcp[2]
#     end
#
#     # panel spacing  (NOTE: this assumes panel spacing is equal, which you must ensure yourself when VLM.wingsection is created.  If you only have one wing section, then it's automatic.)
#     panel_width = CPy[2] - CPy[1]
#
#     # nearest diameter that fits over complete panels (must be an odd number)
#     Deff = 2*R_station[end]
#
#     npanels2 = Deff/panel_width
#     fewerpanels = floor(npanels2)
#     if mod(fewerpanels, 2) == 0 # even number
#         fewerpanels -= 1
#     end
#     Dsmaller = fewerpanels * panel_width
#     Dbigger = Dsmaller + 2*panel_width  # needs to be 2, one on each side, to remain centered
#
#     # iterate over all the props and add up velocities (for bigger and smaller diameter separately)
#     Vxsmall = zeros(CPy)
#     Vzsmall = zeros(CPy)
#     Vxbig = zeros(CPy)
#     Vzbig = zeros(CPy)
#     Vy = zeros(CPy)
#     fail = false
#
#     nprops = length(yprop_center)
#     for i = 1:nprops
#
#         # adjust prop center to center on nearest control point  (this shouldn't change during optimization because we don't change number of props or wing span during optimization)
#         idx = indmin(abs.(CPy - yprop_center[i]))
#         ypc = CPy[idx]
#
#         Vxs, Vzs, fail1 = project_wake_onto_wing(Dsmaller, ypc, R_station, uprop, vprop, cw, CPy)
#         Vxb, Vzb, fail2 = project_wake_onto_wing(Dbigger, ypc, R_station, uprop, vprop, cw, CPy)
#
#         if fail1 || fail2; fail=true; end
#
#         Vxsmall += Vxs
#         Vzsmall += Vzs
#         Vxbig += Vxb
#         Vzbig += Vzb
#     end
#
#     #TODO: verify coordinate systems and interpolation function
#
#     # Wrapper functions for the propeller full wakes
#     intf_bigz = Interpolations.interpolate((CPy,), Vzbig, Interpolations.Gridded(Interpolations.Linear())) #don't re-make the spline each time VLM calls it
#     intf_bigx = Interpolations.interpolate((CPy,), Vxbig, Interpolations.Gridded(Interpolations.Linear()))
#     function bigwake(r)
#         vz = intf_bigz[r[2]]
#         vx = intf_bigx[r[2]]
#         vy = 0.0
#
#         # vx = vx*cos(w_aoa)+vz*sin(w_aoa)
#         # vz = vz*cos(w_aoa)+vx*sin(w_aoa)
#         return [vx,vy,vz]
#     end
#
#     intf_smallz = Interpolations.interpolate((CPy,), Vzsmall, Interpolations.Gridded(Interpolations.Linear()))
#     intf_smallx = Interpolations.interpolate((CPy,), Vxsmall, Interpolations.Gridded(Interpolations.Linear()))
#     function smallwake(r)
#         vz = intf_smallz[r[2]]
#         vx = intf_smallx[r[2]]
#         vy = 0.0
#
#         # vx = vx*cos(-w_aoa)+vz*sin(-w_aoa)
#         # vz = vz*cos(-w_aoa)+vx*sin(-w_aoa)
#         return [vx,vy,vz]
#     end
#
#     fs_small = VLM.Freestream(velocity, w_aoa, beta, Omega, smallwake) #NOTE: prop is assumed to be normal to alpha, not alpha+twist
#     fs_big = VLM.Freestream(velocity, w_aoa, beta, Omega, bigwake)
#
#     # ------ wing analysis ------
#     symmetric = true
#     CF1, CM1, ymid1, zmid1, cl1, dCF1, dCM1, Fp1, alphaeff1, Vinfeff1, cllocal1 = VLM.run(panels, ref, fs_small, symmetric) #TODO: include extra swirl Drag
#     CD1, CY1, CL1 = CF1
#     # Cl, Cm, Cn = CM
#
#     CF2, CM2, ymid2, zmid2, cl2, dCF2, dCM2, Fp2, alphaeff2, Vinfeff2, cllocal2 = VLM.run(panels, ref, fs_big, symmetric) #TODO: include extra swirl Drag
#     CD2, CY2, CL2 = CF2
#     # Cl, Cm, Cn = CM
#
#     # linearly interpolate
#     fraction = (Deff - Dsmaller)/(Dbigger - Dsmaller)
#     CL = CL1 + fraction * (CL2 - CL1)
#     CD = CD1 + fraction * (CD2 - CD1)
#     Fp = Fp1 + fraction * (Fp2 - Fp1)
#     cl = cl1 + fraction * (cl2 - cl1)
#     cllocal = cllocal1 + fraction * (cllocal2 - cllocal1)
#     Vinfeff = Vinfeff1 + fraction * (Vinfeff2 - Vinfeff1)
#     alphaeff = alphaeff1 + fraction * (alphaeff2 - alphaeff1)
#
#     ####################################
#     ########### Viscous Drag ###########
#     ####################################
#     center_chord_array_CPy = linear_interp(center_w_yloc,center_chord_array, CPy) # translate wing geom into VLM frame
#     Re_array = rho/mu*Vinfeff.*center_chord_array_CPy
#     M_array = ones(Re_array).*Vinfeff/a
#     q_array = 0.5*rho*Vinfeff.^2
#
#     blown_array = zeros(Vinfeff)
#     for i = 1:length(Vinfeff)
#         if Vinfeff[i]>velocity
#             blown_array[i] = 1 #turn on if blown
#         end
#     end
#
#     #TODO: include prop induced Drag from swirl
#     viscousDrag, halfspan_dp = viscousdrag(spline_3D_freestream,spline_3D_blown,alphaeff,Re_array,
#     M_array,center_chord_array_CPy,q_array,ds_array_VLM,blown_array)
#
#     # ########### Outputs ###########
#     dyn_press = 0.5*rho*velocity^2
#     Lift = CL * dyn_press*Sref
#     Drag = CD * dyn_press*Sref + viscousDrag
#
#     Fp[1,:] += halfspan_dp
#
#
#
#     #Convert back to structures dimensions
#     Fp_struct = zeros(3,length(center_w_yloc))
#     for i = 1:3
#         Fnew = linear_interp(CPy,Fp[i,:].*ds_array_VLM, center_w_yloc) #
#         Fp_struct[i,:] = Fnew./ds_array
#     end
#
#     if plots
#         figname = "Distrubuted_Lift"
#         figure(figname)
#         plot(CPy,Fp[3,:].*ds_array_VLM,label="VLM F.O.R")
#         plot(center_w_yloc,Fp_struct[3,:].*ds_array,label="Structures F.O.R")
#         xlabel("semispan (m)")
#         ylabel("Lift/span (N/m)")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#         figname = "Distrubuted_Drag"
#         figure(figname)
#         plot(CPy,Fp[1,:].*ds_array_VLM,label="VLM F.O.R")
#         plot(center_w_yloc,Fp_struct[1,:].*ds_array,label="Structures F.O.R")
#         xlabel("semispan (m)")
#         ylabel("total Drag/span (N/m)")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#         figname = "cl"
#         figure(figname)
#         plot(CPy,cl,".-",label = "cl")
#         plot(CPy,cllocal,".-",label = "cllocal")
#         xlabel("semispan (m)")
#         ylabel("cl")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#         figname = "AxialVelocity"
#         figure(figname)
#         plot(CPy,Vxbig,".-",label = "Axial Velocity")
#         # plot(CPy,velocity*ones(CPy),".-",label = "Freestream")
#         xlabel("semispan (m)")
#         ylabel("Velocity (m/s)")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#         figname = "TangentialVelocity"
#         figure(figname)
#         plot(CPy,Vzbig,".-",label = "Tangential Velocity")
#         # plot(CPy,velocity*ones(CPy),".-",label = "Freestream")
#         xlabel("semispan (m)")
#         ylabel("Velocity (m/s)")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#         figname = "Vinf_Eff"
#         figure(figname)
#         plot(CPy,Vinfeff,".-",label = "Effective Vinf")
#         plot(CPy,velocity*ones(CPy),".-",label = "Freestream")
#         xlabel("semispan (m)")
#         ylabel("Vinf (m/s)")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#         figname = "alphaeff"
#         figure(figname)
#         plot(CPy,alphaeff*180/pi,".-",label = "Effective AOA")
#         plot(CPy,StatsBase.mode(alphaeff)*ones(CPy)*180/pi,".-",label = "Wing AOA")
#         xlabel("semispan (m)")
#         ylabel("AOA (deg)")
#         legend(loc="center left",bbox_to_anchor=(1, 0.5))
#         # savefig("./figures/$figname.pdf",transparent = true)
#
#
#
#     end
#
#     return Lift, Drag, Fp, cllocal, cl, Vinfeff, alphaeff, viscousDrag, panels, fail
# end


function oldwinganalysis(sectionspan, chord, twist, tc, sweep, dihedral, N, alpha, rho, Vinf, rprop, uprop, vprop, etaprop_center, cw; plots = false)

        # setup
        beta = 0.0
        wing = VLM.wingsection(sectionspan, chord, twist, tc, sweep, dihedral, N)
        fs = VLM.fs_def(alpha, beta, rho, Vinf)

        QC, TE, CP, LE = VLM.geometry(wing)

        # dimensionalize prop locations
        semispan = sum(sectionspan)
        yprop_center = etaprop_center * semispan

        # # -------- leave commented out: only used to show how the default method of wake control point interaction works poorly -----------
        # Vx, Vz = project_wake_onto_wing(2*rprop[end], yprop_center[1], rprop, uprop, vprop, cw, CP.y)
        # Vext = VLM.velocity(Vx, zeros(CP.y), Vz)
        # L, Di, cl, cllocal, Vinfeff = VLM.run(wing, fs, Vext)
        # return L, Di, CP, cl, cllocal, Vinfeff
        # # ----------------------------------------------

        # ---------- interpolate prop wakes onto wing ------------

        # arbitrary decision: could center on control points (diameters are odd number of panel width) or on QC lines (diameters are even number of panel width).
        # I chose center of control points.

        # panel spacing  (NOTE: this assumes panel spacing is equal, which you must ensure yourself when VLM.wingsection is created.  If you only have one wing section, then it's automatic.)
        panel_width = CP.y[2] - CP.y[1]

        # nearest diameter that fits over complete panels (must be an odd number)
        Deff = 2*rprop[end]
        npanels = Deff/panel_width
        fewerpanels = floor(npanels)
        if mod(fewerpanels, 2) == 0 # even number
            fewerpanels -= 1
        end
        Dsmaller = fewerpanels * panel_width
        Dbigger = Dsmaller + 2*panel_width  # needs to be 2, one on each side, to remain centered

        # iterate over all the props and add up velocities (for bigger and smaller diameter separately)
        Vxsmall = zeros(CP.y)
        Vzsmall = zeros(CP.y)
        Vxbig = zeros(CP.y)
        Vzbig = zeros(CP.y)
        Vy = zeros(CP.y)

        nprops = length(yprop_center)
        for i = 1:nprops

            # adjust prop center to center on nearest control point  (this shouldn't change during optimization because we don't change number of props or wing span during optimization)
            idx = indmin(abs.(CP.y - yprop_center[i]))
            ypc = CP.y[idx]

            Vxs, Vzs = project_wake_onto_wing(Dsmaller, ypc, rprop, uprop, vprop, cw, CP.y)
            Vxb, Vzb = project_wake_onto_wing(Dbigger, ypc, rprop, uprop, vprop, cw, CP.y)

            Vxsmall += Vxs
            Vzsmall += Vzs
            Vxbig += Vxb
            Vzbig += Vzb
        end
        Vextsmall = VLM.velocity(Vxsmall, Vy, Vzsmall)
        Vextbig = VLM.velocity(Vxbig, Vy, Vzbig)

        # ------ wing analysis ------
        L1, Di1, cl1, cllocal1, Vinfeff1, alphaeff1, Vextn1 = VLM.run(wing, fs, Vextsmall)
        L2, Di2, cl2, cllocal2, Vinfeff2, alphaeff2, Vextn2 = VLM.run(wing, fs, Vextbig)
        # linearly interpolate
        fraction = (Deff - Dsmaller)/(Dbigger - Dsmaller)
        L = L1 + fraction * (L2 - L1)
        Di = Di1 + fraction * (Di2 - Di1)
        cl = cl1 + fraction * (cl2 - cl1)
        cllocal = cllocal1 + fraction * (cllocal2 - cllocal1)
        Vinfeff = Vinfeff1 + fraction * (Vinfeff2 - Vinfeff1)
        alphaeff = alphaeff1 + fraction * (alphaeff2 - alphaeff1)
        Vextn = Vextn1 + fraction * (Vextn2 - Vextn1)



        if plots
            # figname = "Distrubuted_Lift"
            # figure(figname)
            # plot(CP.y,Fp[3,:].*ds_array_VLM,label="VLM F.O.R")
            # plot(center_w_yloc,Fp_struct[3,:].*ds_array,label="Structures F.O.R")
            # xlabel("semispan (m)")
            # ylabel("Lift/span (N/m)")
            # legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # # savefig("./figures/$figname.pdf",transparent = true)
            #
            # figname = "Distrubuted_Drag"
            # figure(figname)
            # plot(CP.y,Fp[1,:].*ds_array_VLM,label="VLM F.O.R")
            # plot(center_w_yloc,Fp_struct[1,:].*ds_array,label="Structures F.O.R")
            # xlabel("semispan (m)")
            # ylabel("total Drag/span (N/m)")
            # legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # # savefig("./figures/$figname.pdf",transparent = true)

            figname = "cl"
            figure(figname)
            plot(CP.y,cl,".-",label = "cl")
            plot(CP.y,cllocal,".-",label = "cllocal")
            xlabel("semispan (m)")
            ylabel("cl")
            legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # savefig("./figures/$figname.pdf",transparent = true)

            figname = "AxialVelocity"
            figure(figname)
            plot(CP.y,Vxbig,".-",label = "Axial Velocity")
            # plot(CP.y,velocity*ones(CP.y),".-",label = "Freestream")
            xlabel("semispan (m)")
            ylabel("Velocity (m/s)")
            legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # savefig("./figures/$figname.pdf",transparent = true)

            figname = "TangentialVelocity"
            figure(figname)
            plot(CP.y,Vzbig,".-",label = "Tangential Velocity")
            # plot(CP.y,velocity*ones(CP.y),".-",label = "Freestream")
            xlabel("semispan (m)")
            ylabel("Velocity (m/s)")
            legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # savefig("./figures/$figname.pdf",transparent = true)

            figname = "Vinf_Eff"
            figure(figname)
            plot(CP.y,Vinfeff,".-",label = "Effective Vinf")
            plot(CP.y,velocity*ones(CP.y),".-",label = "Freestream")
            xlabel("semispan (m)")
            ylabel("Vinf (m/s)")
            legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # savefig("./figures/$figname.pdf",transparent = true)

            figname = "alphaeff"
            figure(figname)
            plot(CP.y,alphaeff*180/pi,".-",label = "Effective AOA")
            plot(CP.y,StatsBase.mode(alphaeff)*ones(CP.y)*180/pi,".-",label = "Wing AOA")
            xlabel("semispan (m)")
            ylabel("AOA (deg)")
            legend(loc="center left",bbox_to_anchor=(1, 0.5))
            # savefig("./figures/$figname.pdf",transparent = true)



        end

        return L, Di, CP, cl, cllocal, Vinfeff, alphaeff, Vextn
    end


function winganalysis(w_halfspan,w_chord, w_aoa, w_twist,w_sloc,w_xloc, w_yloc, w_zloc, rho, mu, a, velocity, R_station, uprop, vprop, proploc,spline_3D_freestream,spline_3D_blown;n_CP = 60, cw=1, plots = false)

    # uprop = uprop-velocity
    # vprop = zeros(vprop)

    center_chord_array = (w_chord[2:end]+w_chord[1:end-1])/2.0
    center_w_yloc = cumsum(w_yloc[2:end]-w_yloc[1:end-1])

    npanels = round.(Int,ones(length(w_chord)-1)*n_CP/length(w_chord))
    ds_array = w_sloc[2:end]-w_sloc[1:end-1]
    ds_array_VLM = ones(sum(npanels))*ds_array[1]*length(center_w_yloc)/sum(npanels)  #ASSUMES LINEAR SPACING


    sectionspan = w_yloc[2:end] - w_yloc[1:end-1]
    chord = w_chord
    twist = w_twist
    tc = zeros(w_twist)
    sweep = zeros(w_twist)
    dihedral = zeros(w_twist)
    N = npanels
    alpha = w_aoa
    rho = rho
    Vinf = velocity
    rprop = R_station
    uprop = uprop
    vprop = vprop
    etaprop_center = proploc[div(end,2)+1:end,2]/w_halfspan
    cw = cw


    L, Di, CP, cl, cllocal, Vinfeff, alphaeff, Vextn = oldwinganalysis(sectionspan, chord, twist, tc, sweep, dihedral, N, alpha, rho, Vinf, rprop, uprop, vprop, etaprop_center, cw; plots = plots)

    ####################################
    ########### Viscous Drag ###########
    ####################################
    center_chord_array_CPy = linear_interp(center_w_yloc,center_chord_array, CP.y) # translate wing geom into VLM frame
    Re_array = rho/mu*Vinfeff.*center_chord_array_CPy
    M_array = ones(Re_array).*Vinfeff/a
    q_array = 0.5*rho*Vinfeff.^2

    blown_array = zeros(Vinfeff)
    for i = 1:length(Vinfeff)
        if Vinfeff[i]>velocity
            blown_array[i] = 1 #turn on if blown
        end
    end

    #TODO: include prop induced Drag from swirl
    viscousDrag, halfspan_dp = viscousdrag(spline_3D_freestream,spline_3D_blown,alphaeff,Re_array,
    M_array,center_chord_array_CPy,q_array,ds_array_VLM,blown_array)

    swirl_loss = rho/2*2*pi*R_station[end]*mean(abs.(Vextn))^2*length(proploc) #TODO: full integral?

    # chordVLM = [p.chord for p in detailed_output0.panels]
    chordVLM = CP.chord

    Lift = L
    Drag = Di + viscousDrag + swirl_loss
    Fp = zeros(3,length(cl))
    Fp[1,:] = halfspan_dp+(Di+swirl_loss)*ds_array_VLM/w_halfspan #TODO before wing structures are included
    Fp[2,:] = cl.*chordVLM
    cllocal = cllocal
    cl = cl
    Vinfeff = Vinfeff
    alphaeff = alphaeff
    viscousDrag = viscousDrag
    panels = CP
    fail = false


    return Lift, Drag, Fp, cllocal, cl, Vinfeff, alphaeff, viscousDrag, panels, fail
end
