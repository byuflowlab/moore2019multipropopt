function assemble_structures(layup,n_pt,twist,sloc,xaf,yaf,nchord,lam_t,designparams,x_offset,orientation,xloc,zloc,geom,webloc)
    # assemble structural properties
    mat = Array{Array{Composites.material,1}}(n_pt)
    lam = Array{Composites.laminate,1}(n_pt)
    precompinput = Array{PreComp.input,1}(n_pt)
    precompoutput = Array{PreComp.output,1}(n_pt)
    twist_d = twist*180.0/pi
    twistrate_d = PreComp.tw_rate(n_pt,sloc[1:n_pt],twist_d)
    # strainx = zeros(xafstrain)
    # strainy = zeros(yafstrain)
    # rot = zeros(2,2)
    for i = 1:n_pt
        # find leading edge
        lei = indmin(abs.(xaf[i,:]))
        # shift so leading edge is first
        xpc = circshift(xaf[i,:],-(lei-1))
        ypc = circshift(yaf[i,:],-(lei-1))

        # assemble input
        precompinput[i],mat[i],lam[i] = layup(nchord[i],twist_d[i],
        twistrate_d[i],xpc,ypc,lam_t[i,:],designparams,x_offset[i],orientation)

        # calculate composite properties: stiffness, mass, etc
        precompoutput[i] = PreComp.properties(precompinput[i])

        # rot[1,:] = [cos(-twist[i]),-sin(-twist[i])]
        # rot[2,:] = [sin(-twist[i]),cos(-twist[i])]
        #
        # xy = [xafstrain[i,:]*nchord[i],yafstrain[i,:]*nchord[i]]
        # xyrot = rot*xy
        # # determine locations at which to calculate strain
        # strainx[i,:] = xyrot[1] +xloc[i]
        # strainy[i,:] = xyrot[2] +zloc[i]
    end

    # Get ABD Matrices
    A = Array{Array{Float64,2},1}(n_pt)
    B = Array{Array{Float64,2},1}(n_pt)
    D = Array{Array{Float64,2},1}(n_pt)
    for i = 1:n_pt
        Q = Composites.getQ(mat[i])
        A[i],B[i],D[i] = Composites.getABD(lam[i],Q)
    end

    # Get buckling strains
    bucklingstrain = zeros(Float64,n_pt,2) #currently hardcoded for a single web
    _,bucklingstrain[:,1] = Composites.localbuckling(A,D,webloc*geom.normalchord)
    _,bucklingstrain[:,2] = Composites.localbuckling(A,D,(1-webloc)*geom.normalchord)

    # Assemble structural properties
     return compstructure(mat,lam,A,B,D,bucklingstrain,precompinput,precompoutput)

 end


function p_layup(normalchord,twist_d,twistrate_d,xaf,yaf,
    p_lam_t,designparams,p_x_offset,p_orientation)

    t1 = p_lam_t[1]#uni
    t2 = p_lam_t[2]#weave

    mat = []

    # Materials Input
    e1 = zeros(length(designparams.p_usedmaterials))
    e2 = zeros(length(designparams.p_usedmaterials))
    g12 = zeros(length(designparams.p_usedmaterials))
    anu12 = zeros(length(designparams.p_usedmaterials))
    density = zeros(length(designparams.p_usedmaterials))

    for i = 1:length(designparams.p_usedmaterials)
        matnames = designparams.plyprops.names
        idx = find(matnames -> matnames == designparams.p_usedmaterials[i],matnames)
        material = designparams.plyprops.plies[idx[1]]
        push!(mat,material)
        e1[i] = material.e1
        e2[i] = material.e2
        g12[i] = material.g12
        anu12[i] = material.nu12
        density[i] = material.rho
    end

    web1 = designparams.p_webloc[1]
    xsec_nodeU=[0.0;web1;1.0]

    n_laminaU=[2,2]
    n_pliesU=[1,1,
        1,1]

    t_lamU=[t1,t2,
        t1,t2]

    tht_lamU=[p_orientation[1],p_orientation[2],
        p_orientation[1],p_orientation[2]]
    mat_lamU=[1,2,
        1,2]

    # Lower surface
    xsec_nodeL = xsec_nodeU
    n_laminaL = n_laminaU
    n_pliesL = n_pliesU
    t_lamL = t_lamU
    tht_lamL = tht_lamU
    mat_lamL = mat_lamU

    # Web
    loc_web = [web1]
    n_laminaW = [2]
    n_pliesW = [1,1]
    t_lamW = [0.0,0.0]
    tht_lamW = [p_orientation[1],p_orientation[2]]
    mat_lamW = [1,2]

    leloc = p_x_offset

    # assemble output
    lam = Composites.laminate{Int64,Float64}(mat_lamL[1:div(end,2)],n_pliesL[1:div(end,2)],t_lamL[1:div(end,2)],tht_lamL[1:div(end,2)])

    precompinput = PreComp.input(
    normalchord,
    -twist_d,-twistrate_d,
    leloc,xaf,yaf,
    e1,e2,g12,anu12,density,
    xsec_nodeU,n_laminaU,n_pliesU,t_lamU,tht_lamU,mat_lamU,
    xsec_nodeL,n_laminaL,n_pliesL,t_lamL,tht_lamL,mat_lamL,
    loc_web,n_laminaW,n_pliesW,t_lamW,tht_lamW,mat_lamW)

    return precompinput,mat,lam
end

function w_layup(normalchord,twist_d,twistrate_d,xaf,yaf,
    w_lam_t,designparams,w_x_offset,w_orientation)

    t1 = w_lam_t[1]#uni
    t2 = w_lam_t[2]#weave
    t3 = w_lam_t[3]#wing foam
    t4 = w_lam_t[4]#web fabric
    t5 = w_lam_t[5]#web foam

    mat = []

    # Materials Input
    e1 = zeros(length(designparams.w_usedmaterials))
    e2 = zeros(length(designparams.w_usedmaterials))
    g12 = zeros(length(designparams.w_usedmaterials))
    anu12 = zeros(length(designparams.w_usedmaterials))
    density = zeros(length(designparams.w_usedmaterials))

    for i = 1:length(designparams.w_usedmaterials)
        matnames = designparams.plyprops.names
        idx = find(matnames -> matnames == designparams.w_usedmaterials[i],matnames)
        material = designparams.plyprops.plies[idx[1]]
        push!(mat,material)
        e1[i] = material.e1
        e2[i] = material.e2
        g12[i] = material.g12
        anu12[i] = material.nu12
        density[i] = material.rho
    end

    web1 = designparams.w_webloc[1]
    xsec_nodeU=[0.0;web1;1.0]

    n_laminaU=[5,5]
    n_pliesU=[1,1,1,1,1,
        1,1,1,1,1]

    t_lamU=[t1,t2,t3,t2,t1,
        t1,t2,t3,t2,t1]

    tht_lamU=[w_orientation[1],w_orientation[2],0.0,w_orientation[2],w_orientation[1],
        w_orientation[1],w_orientation[2],0.0,w_orientation[2],w_orientation[1]]
    mat_lamU=[1,2,3,2,1,
        1,2,3,2,1]

    # Lower surface
    xsec_nodeL = xsec_nodeU
    n_laminaL = n_laminaU
    n_pliesL = n_pliesU
    t_lamL = t_lamU
    tht_lamL = tht_lamU
    mat_lamL = mat_lamU

    # Web
    loc_web = [web1]
    n_laminaW = [3]
    n_pliesW = [1,1,1]
    t_lamW = [t4,t5,t4]
    tht_lamW = [w_orientation[3],0.0,w_orientation[3]]
    mat_lamW = [2,3,2]

    leloc = w_x_offset

    # assemble output
    lam = Composites.laminate{Int64,Float64}(mat_lamL[1:div(end,2)],n_pliesL[1:div(end,2)],t_lamL[1:div(end,2)],tht_lamL[1:div(end,2)])

    precompinput = PreComp.input(
    normalchord,
    -twist_d,-twistrate_d,
    leloc,xaf,yaf,
    e1,e2,g12,anu12,density,
    xsec_nodeU,n_laminaU,n_pliesU,t_lamU,tht_lamU,mat_lamU,
    xsec_nodeL,n_laminaL,n_pliesL,t_lamL,tht_lamL,mat_lamL,
    loc_web,n_laminaW,n_pliesW,t_lamW,tht_lamW,mat_lamW)

    return precompinput,mat,lam
end

function beam_strain_wrapper(Fp,n_sec,ctlparams,printiter,plots,xaf,yaf,xafstrain,yafstrain,twist,x_offset,aerocenter,nchord,structure,spanyloc,xloc,yloc,zloc)

    momarmx_scac = zeros(n_sec) #shear center (sc) to aero center (ac) moment arm
    momarmy_scac = zeros(n_sec)
    momarmx_cmtc = zeros(n_sec)
    momarmy_cmtc = zeros(n_sec)
    momarmx_sccm = zeros(n_sec)
    momarmy_sccm = zeros(n_sec)

    if (printiter%(ctlparams.printfreq)==0.0 && plots)
        figure("af")
        clf()
    end


    strainx = zeros(xafstrain) # for strain locations about the shear center
    strainy = zeros(yafstrain)
    strainx_vis = zeros(xafstrain)
    strainy_vis = zeros(yafstrain)
    rot = zeros(2,2)

    for i = 1:n_sec
        # Calculate aerodynamic center from precalculated x-aero center and chord line
        idxle = indmin(abs.(xaf[i,:]))
        idxte = indmax(abs.(xaf[i,:]))

        yle = yaf[i,idxle]
        yte = yaf[i,idxte]

        acx = aerocenter[i]
        acy = (yte-yle)/(1)*acx + yle #y/mx+b

        y_ac = (acx+x_offset[i])*nchord[i] #x and y are flipped in precomp
        x_ac = acy*nchord[i]
        x_sc = x_ac + structure.precompoutput[i].x_sc
        y_sc = y_ac + structure.precompoutput[i].y_sc
        x_tc = x_ac + structure.precompoutput[i].x_tc
        y_tc = y_ac + structure.precompoutput[i].y_tc
        x_cm = x_ac + structure.precompoutput[i].x_cm
        y_cm = y_ac + structure.precompoutput[i].y_cm


        #moment arms defined positive as in Fig. 13 in precomp manual with tc above cm, x and y are flipped as well in precomp
        momarmx_scac[i] = x_ac-y_sc
        momarmy_scac[i] = y_ac-x_sc

        momarmx_cmtc[i] = y_tc-y_cm
        momarmy_cmtc[i] = x_tc-x_cm

        momarmx_sccm[i] = y_cm-y_sc
        momarmy_sccm[i] = x_cm-x_sc

        #-------- Calculate Rotated Strain Locations about the Shear Center --------#
        rot[1,:] = [cos(-twist[i]),-sin(-twist[i])]
        rot[2,:] = [sin(-twist[i]),cos(-twist[i])]

        xy = [xafstrain[i,:]*nchord[i]+x_offset[i]*nchord[i]-y_sc,yafstrain[i,:]*nchord[i]-x_sc] #Offest by the shear center, precomp FOR is swapped
        xyrot = rot*xy
        # determine locations at which to calculate strain
        strainx[i,:] = xyrot[1]
        strainy[i,:] = xyrot[2]

        # add offset for visualization
        xy = [xafstrain[i,:]*nchord[i]-0.125*nchord[i],yafstrain[i,:]*nchord[i]] #Offest by the shear center, precomp FOR is swapped
        xyrot = rot*xy
        strainx_vis[i,:] = xyrot[1] + xloc[i]
        strainy_vis[i,:] = xyrot[2] + zloc[i]

        # Plots
        if (i%5==0 || i==n_sec) && (printiter%(ctlparams.printfreq)==0.0 && plots)

            figure("af$i")
            # clf()
            plot((xaf[i,:]+x_offset[i])*nchord[i],yaf[i,:]*nchord[i])
            if i==n_sec
                plot(y_ac,x_ac,"x",label = "Aero Center")
                plot(y_sc,x_sc,"D",label = "Shear Center")
                plot(y_tc,x_tc,".",label = "Mass Center")
                plot(y_cm,x_cm,"+",label = "Tension Center")
            else
                plot(y_ac,x_ac,"x")
                plot(y_sc,x_sc,"D")
                plot(y_tc,x_tc,".")
                plot(y_cm,x_cm,"+")
            end
            axis("equal")
            legend(loc="center left", bbox_to_anchor=(1, 0.5))
            pause(0.001)

            figure("strain_af$i")
            # clf()
            plot(strainx[i,:],strainy[i,:])
            plot(xafstrain[i,:]*nchord[i]+x_offset[i]*nchord[i],yafstrain[i,:]*nchord[i])
            if i==n_sec
                plot(y_ac-y_sc,x_ac-x_sc,"x",label = "Aero Center")
                plot(y_sc-y_sc,x_sc-x_sc,"D",label = "Shear Center")
                plot(y_tc-y_sc,x_tc-x_sc,".",label = "Mass Center")
                plot(y_cm-y_sc,x_cm-x_sc,"+",label = "Tension Center")

                plot(y_ac,x_ac,"x",label = "Aero Center")
                plot(y_sc,x_sc,"D",label = "Shear Center")
                plot(y_tc,x_tc,".",label = "Mass Center")
                plot(y_cm,x_cm,"+",label = "Tension Center")
            else
                plot(y_ac-y_sc,x_ac-x_sc,"x")
                plot(y_sc-y_sc,x_sc-x_sc,"D")
                plot(y_tc-y_sc,x_tc-x_sc,".")
                plot(y_cm-y_sc,x_cm-x_sc,"+")

                plot(y_ac,x_ac,"x")
                plot(y_sc,x_sc,"D")
                plot(y_tc,x_tc,".")
                plot(y_cm,x_cm,"+")
            end
            axis("equal")
            legend(loc="center left", bbox_to_anchor=(1, 0.5))
            pause(0.001)

        end

    end



    #-------- TRANSLATE CCBLADE LOADS TO PRECOMP/BEAM LOADS -------#
    nodes = length(spanyloc)
    elements = nodes - 1

    #Beam goes from root to tip
    Py = Fp[3,:] #Lift is z in aero, y in beam
    Pz = Fp[1,:] #Drag is x in aero, z in beam
    Px = Fp[2,:] #compression is y in aero, x in beam

    # --- extract the point forces from distributed in order to apply moments from ac to sc etc -----
    DOF = 6
    F = zeros(DOF*nodes)
    for i = 1:elements
        start = (i-1)*DOF  # (0, 0) start of matrix

        _, _, Fsub = BeamFEA.beam_matrix(spanyloc[i+1] - spanyloc[i], [0.0,0.0], [0.0,0.0], [0.0,0.0],
        [0.0,0.0], [0.0,0.0], [0.0,0.0], Px[i:i+1], Py[i:i+1], Pz[i:i+1])

        idx = start+1:start+2*DOF
        F[idx] += Fsub
    end

    # These are in the beam frame of reference
    EIy = zeros(nodes)
    EIz = zeros(nodes)
    EA = zeros(nodes)
    GJ = zeros(nodes)
    rhoA = zeros(nodes)
    rhoJ = zeros(nodes)
    # Px = zeros(nodes) # Already Calculated
    # Py = zeros(nodes)
    # Pz = zeros(nodes)
    Fx = zeros(nodes)
    Fy = zeros(nodes)
    Fz = zeros(nodes)
    Mx = zeros(nodes)
    My = zeros(nodes)
    Mz = zeros(nodes)
    kx = zeros(nodes)
    ky = zeros(nodes)
    kz = zeros(nodes)
    kthetax = zeros(nodes)
    kthetay = zeros(nodes)
    kthetaz = zeros(nodes)

    idxf = 1
    for i = 1:nodes
        # these are all distributed properties
        mass = structure.precompoutput[i].mass
        iyy = structure.precompoutput[i].flap_iner
        ixx = structure.precompoutput[i].lag_iner
        d_sccm2 =momarmx_sccm[i]^2 + momarmy_sccm[i]^2
        izz = ixx+iyy + mass*d_sccm2 #perpendicular axis theorem and parallel axis theorm

        EIy[i] = structure.precompoutput[i].ei_flap
        EIz[i] = structure.precompoutput[i].ei_lag
        EA[i] = structure.precompoutput[i].ea
        GJ[i] = structure.precompoutput[i].gj
        rhoA[i] = mass
        rhoJ[i] = izz
        # Px[i] = 0.0 # Already Calculated
        # Py[i] = 0.0
        # Pz[i] = 0.0
        Fx[i] = 0.0 #TODO:Extra Point Loads centripetal force here
        Fy[i] = 0.0
        Fz[i] = 0.0
        Mx[i] = F[idxf+3]*momarmy_scac[i] + F[idxf+5]*momarmx_scac[i]
        My[i] = 0.0 #TODO:centripetal here, but tension center not aligned with shear center so...?
        Mz[i] = 0.0 #TODO:centripetal here, but tension center not aligned with shear center so...?

        idxf += DOF
    end

    #fixed at hub (future work may allow axial twisting via a passive pitching device)
    kx[1] = Inf
    ky[1] = Inf
    kz[1] = Inf
    kthetax[1] = Inf
    kthetay[1] = Inf
    kthetaz[1] = Inf

    #-------- FEA ANALYSIS -------#
    delta, freq, V, K, M, F = BeamFEA.fea_analysis(spanyloc, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
    Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)
    # println(minimum(delta))
    #TODO: use non-rigid strain calculation, i.e. composite curvature
    #TODO: include shear strain
    # get strain for each x strip along the blade
    strains = zeros(length(spanyloc),length(strainy[1,:]))
    for i = 1:length(strainy[1,:])
        strains[:,i], Nx, Vy, Vz, Tx, Myout, Mzout = BeamFEA.strain(spanyloc,
        strainy[:,i], strainx[:,i], EIy, EIz, EA,
        Px, Py, Pz, Fx, Fy, Fz, Mx, My, Mz, false)
    end

    mass_structure = trapz(spanyloc,rhoA)


    return strains, delta, freq, mass_structure,strainx_vis,strainy_vis
end


function bucklingcon(structure,strain,webloc,normalchord)
  # Used for constraining buckling
  sf = 1.5 #safety factor

  # buckling strain is negative, strain is always positive

  # NOT SMOOTH - will cause problems if webloc is design variable
  c_buckling = zeros(size(strain,1),size(strain,2))
  for i = 1:size(c_buckling,1)
    for iloc = 1:size(c_buckling,2)
      # if structure.strainlocx[i,iloc] < webloc*normalchord[i] #|| webloc == 0.0
          c_buckling[i,iloc] = (structure.bucklingstrain[i,1]+sf*strain[i,iloc])
      # else
      #     c_buckling[i,iloc] = (structure.bucklingstrain[i,2]+sf*strain[i,iloc])
      # end
    end
  end


  # # SMOOTH, but unnecessarily conservative
  # c_buckling = zeros(size(strain,1),size(strain,2),size(def.bucklingstrain,2)*2)
  # for i = 1:size(c_buckling,1)
  #   for iloc = 1:size(c_buckling,2)
  #     for ibuckle = 1:size(def.bucklingstrain,2)
  #       c_buckling[i,iloc,2*ibuckle-1] = (-def.bucklingstrain[i,ibuckle]-sf*strain[i,iloc])/def.bucklingstrain[i,ibuckle]
  #       c_buckling[i,iloc,2*ibuckle] = (-def.bucklingstrain[i,ibuckle]+sf*strain[i,iloc])/def.bucklingstrain[i,ibuckle]
  #     end
  #   end
  # end

  return c_buckling[:]/1E2
end

function stresscalc(strain::Array{Float64,2},shear::Array{Float64,1},
    mat::Composites.material,theta_d::Float64)
    # Get rotated material stiffness matrix
    qbar = Composites.getQ(mat,theta_d)
    # Determine Stresses from Strains in individual plys
    nnodes = size(strain,1)
    nloc = size(strain,2)
    plystress = zeros(nnodes,nloc,3)
    for i = 1:nnodes #Run entire length of wing #
        for iloc = 1:nloc # Run at every specified strain location
            # Convert shear stress to shear strain
            gam12 = shear[i]/qbar[3,3]-qbar[1,3]/qbar[3,3]*strain[i,iloc]
            # Calculate laminate stresses
            stress = qbar*[strain[i,iloc],0,gam12]
            # Transform stresses to ply axes
            plystress[i,iloc,:] = Composites.rotstress(stress,theta_d)
        end
    end
    return plystress
end #stress_calc

function stresscon(stress,mat::Composites.material,criteria = "tsaiwu";BVID=0.65,sf=1.5)
    # Used to constrain stress
    sigma1 = sf*stress[:,:,1]
    sigma2 = sf*stress[:,:,2]
    tau12 = sf*stress[:,:,3]
    sigma1tu = BVID*mat.xt
    sigma2tu = BVID*mat.yt
    sigma1cu = BVID*mat.xc
    sigma2cu = BVID*mat.yc
    tau12u = BVID*mat.s
    if criteria=="maxstress"
        c_stress = zeros(Float64,6,size(sigma1)...)
        for i = 1:size(sigma1,1)
            for j = 1:size(sigma1,2)
                c_stress[:,i,j],_ = Composites.maxstress(sigma1[i,j],sigma2[i,j],tau12[i,j],
                sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
            end
        end
    elseif criteria=="tsaiwu"
        c_stress = zeros(Float64,size(sigma1)...)
        for i = 1:size(sigma1,1)
            for j = 1:size(sigma1,2)
                c_stress[i,j],_ = Composites.tsaiwu(sigma1[i,j],sigma2[i,j],tau12[i,j],
                sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
                # maximum(c_stress)
            end
        end
    elseif criteria=="hashinrotem"
        c_stress = zeros(Float64,4,size(sigma1)...)
        for i = 1:size(sigma1,1)
            for j = 1:size(sigma1,2)
                c_stress[:,i,j],_ = Composites.hashinrotem(sigma1[i,j],sigma2[i,j],tau12[i,j],
                sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
            end
        end
    else
        error("Specified failure criteria has no implementation")
    end
    return c_stress[:]-1.0
end #stresscon

function stress_wrapper(usedmaterials,plyprops,spanlocy,orientation,strain,shear)

    c_stress = []
    stress = [] #used in vtk output
    stresslayer = zeros(length(usedmaterials),length(spanlocy),length(strain[1,:]))
    j=1
    for i = 1:length(usedmaterials)

        if !contains(usedmaterials[i],"foam") #don't check stress on foam
            matnames = plyprops.names
            idx = find(matnames -> matnames == usedmaterials[i],matnames)

            material = plyprops.plies[idx[1]]
            orien = orientation[j]
            j+=1 #so we don't have to have extra design variables in orientation

            stressi = stresscalc(strain,shear,material,orien)

            if i==1 #only calc once since it's the same for the whole structure (assuming thin layups)
                stress = sqrt.(stressi[:,:,1].^2+stressi[:,:,2].^2+stressi[:,:,3].^2) #TODO break up or calculate von-mises etc
            end

            # determine contraint values
            c_stressi = stresscon(stressi,material)
            stresslayer[i,:,:] = c_stressi

            c_stress = [c_stress;c_stressi]
        end
    end

    return c_stress, stress, stresslayer
end #stress_wrapper

function spancoord(spanyloc,Rtip,dihedral,sweep,twist,x_offset,chord)
    # Define number of sections
    nsections = length(spanyloc)
    # Main  sections
    xloc = zeros(nsections)
    yloc = zeros(nsections)
    zloc = zeros(nsections)
    xloc[2:end] = cumsum((spanyloc[2:end]-spanyloc[1:(end-1)]).*cos.(dihedral[1:nsections-1]).*sin.(sweep[1:nsections-1]))
    yloc[2:end] = cumsum((spanyloc[2:end]-spanyloc[1:(end-1)]).*cos.(dihedral[1:nsections-1]).*cos.(sweep[1:nsections-1]))
    zloc[2:end] = cumsum((spanyloc[2:end]-spanyloc[1:(end-1)]).*sin.(dihedral[1:nsections-1]).*cos.(sweep[1:nsections-1]))
    # Find y length of  by stepping back spanylochrough let
    yend = Rtip
    # Scale appropriately to match span
    scaling = yend/yloc[nsections]
    xloc = xloc*scaling - cos.(twist).*x_offset.*chord
    yloc = yloc*scaling
    zloc = zloc*scaling + sin.(twist).*x_offset.*chord
    # Get structural spanwise parameter
    sloc = zeros(nsections)
    for i = 2:nsections
        sloc[i] = sloc[i-1]+sqrt((xloc[i]-xloc[i-1])^2.0+
        (yloc[i]-yloc[i-1])^2.0+(zloc[i]-zloc[i-1])^2.0)
    end

    # Assemble full  parameters
    sloc = sloc
    xloc = xloc
    yloc = yloc
    zloc = zloc

    return sloc,xloc,yloc,zloc
end

function chordlengths(chord,sweep)
    # Determine number of sections
    nsections = length(sweep) + 1
    # Normal Chord lengths
    normalchord = zeros(nsections)
    normalchord[1] = chord[1]*cos(sweep[1])
    for i = 2:(nsections-1)
        normalchord[i] = chord[i]*cos((sweep[i-1]+sweep[i])/2.0)
    end
    normalchord[end] = chord[end]*cos(sweep[end])
    return normalchord
end

function printmaxviol(c,name)
    if !isempty(find(c.>-1e-4))
        println(string(name),": ",maximum(c))
    end
end

#
# (private)
# trapezoidal integration
#
function trapz(x::Array{Float64,1}, y::Array{Float64,1})  # integrate y w.r.t. x
    integral = 0.0
    for i = 1:length(x)-1
        integral += (x[i+1]-x[i])*0.5*(y[i] + y[i+1])
    end
    return integral
end

function linear_interp(x_array::Array{Float64,1},y_array::Array{Float64,1},xnew::Float64)


    i = 1
    if xnew<minimum(x_array) # cap at max and min
        return y_array[indmin(x_array)]
    elseif xnew>maximum(x_array)
        return y_array[indmax(x_array)]
    else
        for i = 1:length(x_array)-1
            if x_array[i]<=xnew && x_array[i+1]>=xnew
                break
            end
        end

        fraction = (xnew - x_array[i])/(x_array[i+1] - x_array[i])
        return y_array[i] + fraction * (y_array[i+1] - y_array[i])
    end

end

function linear_interp(x_array::Array{Float64,1},y_array::Array{Float64,1},xnew::Array{Float64,1})
    ynew = zeros(xnew)
    for i = 1:length(xnew)
        ynew[i] = linear_interp(x_array,y_array,xnew[i])
    end

    return ynew
end
