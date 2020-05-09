#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
#   * Modified  : Kevin Moore moorekevin711@gmail.com
=###############################################################################

function save_vtk_opt(p_stress,p_strain,p_stresslayer,p_c_bucklingplot,p_xloc,p_yloc,p_zloc,numprops,w_stress,w_strain,w_stresslayer,w_c_bucklingplot,w_xloc,w_yloc,w_zloc;
    file_name="./VTK_output/temp_kevin00",
    p_radii = [0.0,0.2,0.4,0.6,0.8,1.0],
    chord = [1.619, 1.473, 1.456, 0.789, 0.789, 0.789],
    p_twist_d = [1.619, 1.473, 1.456, 0.789, 0.789, 0.789], #twist includes pitch
    p_dihedral_d = pi/180*vcat([0 for i in 1:5], [23.20]),
    p_sweep_d = pi/180*vcat([23.20 for i in 1:6], []),
    Rtip = 1.25,
    p_airfoilx = 0.0, #TODO, input default
    p_airfoily = 0.0,
    p_strainx = 0.0,
    p_strainy = 0.0,
    w_twist_d = [1.619, 1.473, 1.456, 0.789, 0.789, 0.789],
    w_dihedral_d = pi/180*vcat([0 for i in 1:5], [23.20]),
    w_sweep_d = pi/180*vcat([23.20 for i in 1:6], []),
    w_strainx = 0.0,
    w_strainy = 0.0,
    spanloc = [0.0,1.0],
    halfspan = 10.0,
    printiter = 1,
    proploc = 0.0,
    numblades = 2,
    plots = false )


    numblades = Int(round(numblades))

    # GEOMETRY DEFINITION
    # dihedral = p_dihedral_d * pi/180
    # sweep = p_sweep_d * pi/180

    forward_offset = Rtip/10.0

    y_pos = copy(p_yloc)#Rtip*p_radii # span position of each section
    x_pos = zeros(p_yloc)-forward_offset#[y*tan(sweep[i]) for (i,y) in enumerate(y_pos)] #Already translated
    z_pos = zeros(p_yloc)#[y*tan(dihedral[i]) for (i,y) in enumerate(y_pos)]
    # PARAMETERS
    r = fill(1.0,length(p_radii))          # Expansion ratio in both surfaces of each airfoil

    azimuth_deg_spacing = round(Int, 360/numblades)
    azimuth_deg = 0.0
    for i = 1#:numprops
        for j = 1:numblades

            generateprop!(p_strainx,p_strainy,p_strain,p_stresslayer,p_stress,p_c_bucklingplot,p_airfoilx,p_airfoily,x_pos,y_pos,z_pos,proploc[i,:],p_sweep_d,p_twist_d,i,j,p_radii,azimuth_deg,file_name;plots=plots)
            azimuth_deg = azimuth_deg_spacing + azimuth_deg

        end
    end


    # ############################################
    # ########    Generate Wing  #################
    # ############################################
    #
    # y_pos = copy(w_yloc)#Rtip*w_radii # span position of each section
    # x_pos = zeros(w_yloc)#[y*tan(sweep[i]) for (i,y) in enumerate(y_pos)] #Already translated
    # z_pos = zeros(w_yloc)#[y*tan(dihedral[i]) for (i,y) in enumerate(y_pos)]
    # # PARAMETERS
    # r = fill(1.0,length(spanloc))          # Expansion ratio in both surfaces of each airfoil
    #
    # sections = fill(Tuple{Float64,Int64,Float64,Bool}[(1.0, 10, 1.0, false)],(length(spanloc)-1))
    #
    # # Leading edge position of each airfoil
    # Os = [ [x_pos[i], y_pos[i], z_pos[i]] for i in 1:size(w_strainx)[1]]
    # # Orientation of chord of each airfoil (yaw, pitch, roll)
    # orien = zeros(length(w_twist_d),3) #twist already included
    # for i = 1:length(orien[:,1])
    #     orien[i,:] = [0.0,w_sweep_d[i],270.0]
    #     if i == 1
    #         orien[i,:] = [0.0,0.0,270.0]
    #     end
    # end
    #
    # crosssections = []        # It will store here the cross sections for lofting
    # strains = []
    # stress1 = []
    # stress2 = []
    # stressplot = []
    # buckle = []
    # point_datas = []
    #
    # # Processes each airfoil geometry - STRESS/STRAIN CALCS
    # for i = 1:length(w_strainx[:,1])
    #
    #     # Read airfoil file
    #     x = w_strainx[i,:]
    #     y = w_strainy[i,:]
    #
    #     # Scales the airfoil acording to its chord length
    #     new_x = x#new_x #already scaled
    #     new_y = y#new_y
    #
    #     # plot_airfoil(new_x, new_y; style=styles[i], label="airfoil")
    #     if plots
    #         figure("wing_airfoil")
    #         plot(new_x,new_y)
    #     end
    #
    #     # Reformats into points
    #     npoints = size(new_x)[1]
    #     airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]
    #
    #     # Positions the airfoil along the blade in the right p_orientation
    #     Oaxis = VTKtools.rotation_matrix(orien[i,1], orien[i,2], orien[i,3])
    #     invOaxis = inv(Oaxis)
    #     airfoil = VTKtools.countertransform(airfoil, invOaxis, Os[i])
    #
    #     push!(crosssections, airfoil)
    #
    #     push!(strains, w_strain[i,:])
    #     push!(stress1, w_stresslayer[1,i,:])
    #     push!(stress2, w_stresslayer[2,i,:])
    #     push!(stressplot, w_stress[i,:])
    #     push!(buckle, w_c_bucklingplot[i,:])
    #     push!(point_datas, [j for j in npoints*(i-1)+1:npoints*i])
    # end
    #
    # # Generates cells in VTK Legacy format
    # points, vtk_cells, strain = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=strains)
    # _, _, stresses1 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stress1)
    # _, _, stresses2 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stress2)
    # _, _, stresses6 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stressplot)
    # _, _, bucklecon = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=buckle)
    # points, vtk_cells, point_datas = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=point_datas)
    #
    # # Formats the point data for generateVTK
    # data = []
    # push!(data, Dict(
    # "field_name" => "Strain",
    # "field_type" => "scalar",
    # "field_data" => strain))
    #
    # push!(data, Dict(
    # "field_name" => "uni_stress",
    # "field_type" => "scalar",
    # "field_data" => stresses1))
    #
    # push!(data, Dict(
    # "field_name" => "weave_stress",
    # "field_type" => "scalar",
    # "field_data" => stresses2))
    #
    # push!(data, Dict(
    # "field_name" => "stress",
    # "field_type" => "scalar",
    # "field_data" => stresses6))
    #
    # push!(data, Dict(
    # "field_name" => "Buckling_Constraint",
    # "field_type" => "scalar",
    # "field_data" => bucklecon))
    #
    # push!(data, Dict(
    # "field_name" => "PointIndex",
    # "field_type" => "scalar",
    # "field_data" => point_datas))
    #
    # # Generates the vtk file
    # VTKtools.generateVTK("$(file_name)wingR.$printiter", points; cells=vtk_cells, point_data=data, keep_points=false)
    #
    #
    # #generate other wing
    # Os = [ [x_pos[i], -y_pos[i], z_pos[i]] for i in 1:size(w_strainx)[1]]
    # # Orientation of chord of each airfoil (yaw, pitch, roll)
    # orien = zeros(length(w_twist_d),3) #twist already included
    # for i = 1:length(orien[:,1])
    #     orien[i,:] = [0.0,-w_sweep_d[i],270.0]
    #     if i == 1
    #         orien[i,:] = [0.0,0.0,270.0]
    #     end
    # end
    #
    # crosssections = []        # It will store here the cross sections for lofting
    # strains = []
    # stress1 = []
    # stress2 = []
    # stressplot = []
    # buckle = []
    # point_datas = []
    # for i = length(w_strainx[:,1]):-1:1
    #
    #     # Read airfoil file
    #     x = w_strainx[i,:]
    #     y = w_strainy[i,:]
    #
    #     # Scales the airfoil acording to its chord length
    #     new_x = x#new_x #already scaled
    #     new_y = y#new_y
    #
    #     # plot_airfoil(new_x, new_y; style=styles[i], label="airfoil")
    #     if plots
    #         figure("wing_airfoil")
    #         plot(new_x,new_y)
    #     end
    #
    #     # Reformats into points
    #     npoints = size(new_x)[1]
    #     airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]
    #
    #     # Positions the airfoil along the blade in the right p_orientation
    #     Oaxis = VTKtools.rotation_matrix(orien[i,1], orien[i,2], orien[i,3])
    #     invOaxis = inv(Oaxis)
    #     airfoil = VTKtools.countertransform(airfoil, invOaxis, Os[i])
    #
    #     push!(crosssections, airfoil)
    #
    #     push!(strains, w_strain[i,:])
    #     push!(stress1, w_stresslayer[1,i,:])
    #     push!(stress2, w_stresslayer[2,i,:])
    #     push!(stressplot, w_stress[i,:])
    #     push!(buckle, w_c_bucklingplot[i,:])
    #     push!(point_datas, [j for j in npoints*(i-1)+1:npoints*i])
    # end
    #
    # # Generates cells in VTK Legacy format
    # points, vtk_cells, strain = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=strains)
    # _, _, stresses1 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stress1)
    # _, _, stresses2 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stress2)
    # _, _, stresses6 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stressplot)
    # _, _, bucklecon = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=buckle)
    # points, vtk_cells, point_datas = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=point_datas)
    #
    # # Formats the point data for generateVTK
    # data = []
    # push!(data, Dict(
    # "field_name" => "Strain",
    # "field_type" => "scalar",
    # "field_data" => strain))
    #
    # push!(data, Dict(
    # "field_name" => "uni_stress",
    # "field_type" => "scalar",
    # "field_data" => stresses1))
    #
    # push!(data, Dict(
    # "field_name" => "weave_stress",
    # "field_type" => "scalar",
    # "field_data" => stresses2))
    #
    # push!(data, Dict(
    # "field_name" => "stress",
    # "field_type" => "scalar",
    # "field_data" => stresses6))
    #
    # push!(data, Dict(
    # "field_name" => "Buckling_Constraint",
    # "field_type" => "scalar",
    # "field_data" => bucklecon))
    #
    # push!(data, Dict(
    # "field_name" => "PointIndex",
    # "field_type" => "scalar",
    # "field_data" => point_datas))
    #
    # # Generates the vtk file
    # VTKtools.generateVTK("$(file_name)wingL.$printiter", points; cells=vtk_cells, point_data=data, keep_points=false)
    #
    #
    # # #TEST vector output
    # # p0 = [0,0,0.0]
    # # p1 = [0,1,0.0]
    # # c0 = [p0,p1,p1,p0]
    # # d0 = [1.0,1.0,1.0,1.0]
    # # points, vtk_cells, lift = VTKtools.lines2vtk([c0];
    # #                                        values=[d0])
    # # data = []
    # # push!(data, Dict(
    # # "field_name" => "Strain",
    # # "field_type" => "vector",
    # # "field_data" => lift))
    # #
    # # VTKtools.generateVTK("$(file_name)lift.$printiter", points;
    # #          cells=vtk_cells, point_data=data)

    println("VTK Saved")

end

function generateprop!(p_strainx,p_strainy,p_strain,stresslayer,stress,p_c_bucklingplot,p_airfoilx,p_airfoily,x_pos,y_pos,z_pos,proploc,p_sweep_d,p_twist_d,propnum,bladenum,p_radii,azimuth_deg,file_name;plots = false)

    sections = fill(Tuple{Float64,Int64,Float64,Bool}[(1.0, 10, 1.0, false)],((length(p_radii)-1)))#*2+1))

    crosssections = []        # It will store here the cross sections for lofting
    strains = []
    stress1 = []
    stress2 = []
    stressplot = []
    buckle = []
    point_datas = []

    # Leading edge position of each airfoil
    Os = [ [x_pos[i]+proploc[1], y_pos[i]+proploc[2], z_pos[i]]+proploc[3] for i in 1:size(p_airfoilx)[1]]
    # Orientation of chord of each airfoil (yaw, pitch, roll)
    orien = zeros(length(p_twist_d),3) #twist already included
    for i = 1:length(orien[:,1])
        orien[i,:] = [-90.0,0.0,270.0]
    end

    # Processes each airfoil geometry - STRESS/STRAIN CALCS
    for i = 1:length(p_strainx[:,1])

        # Read airfoil file
        x = p_strainx[i,:]
        y = p_strainy[i,:]

        # Scales the airfoil acording to its chord length
        new_x = x#new_x #already scaled
        new_y = y#new_y

        # plot_airfoil(new_x, new_y; style=styles[i], label="airfoil")
        if plots
            figure("prop_airfoil")
            plot(new_x,new_y)
        end

        # Reformats into points
        npoints = size(new_x)[1]
        airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]

        # Positions the airfoil along the blade in the right p_orientation
        Oaxis = VTKtools.rotation_matrix(orien[i,1], orien[i,2], orien[i,3])
        invOaxis = inv(Oaxis)
        airfoil = VTKtools.countertransform(airfoil, invOaxis, Os[i])

        spinaxis = VTKtools.rotation_matrix(0.0, 0.0, azimuth_deg)
        # invspinaxis = inv(spinaxis)
        airfoil = transform(airfoil, spinaxis, [proploc[1],proploc[2],proploc[3]])

        push!(crosssections, airfoil)

        push!(strains, p_strain[i,:])
        push!(stress1, stresslayer[1,i,:])
        push!(stress2, stresslayer[2,i,:])
        push!(stressplot, stress[i,:])
        push!(buckle, p_c_bucklingplot[i,:])
        push!(point_datas, [j for j in npoints*(i-1)+1:npoints*i])
    end

    # #generate other blade
    # Os = [ [x_pos[i]+proploc[1], -y_pos[i]+proploc[2], z_pos[i]]+proploc[3] for i in 1:size(p_airfoilx)[1]]
    # for i = length(p_strainx[:,1]):-1:1
    #
    #     # Read airfoil file
    #     x = p_strainx[i,:]
    #     y = p_strainy[i,:]
    #
    #     # Scales the airfoil acording to its chord length
    #     new_x = x#new_x #already scaled
    #     new_y = y#new_y
    #
    #     # plot_airfoil(new_x, new_y; style=styles[i], label="airfoil")
    #     if plots
    #         figure("prop_airfoil")
    #         plot(new_x,new_y)
    #     end
    #
    #     # Reformats into points
    #     npoints = size(new_x)[1]
    #     airfoil = Array{Float64, 1}[[-new_x[j], new_y[j], 0] for j in 1:npoints]
    #
    #     # Positions the airfoil along the blade in the right p_orientation
    #     Oaxis = VTKtools.rotation_matrix(orien[i,1], orien[i,2], orien[i,3])
    #     invOaxis = inv(Oaxis)
    #     airfoil = VTKtools.countertransform(airfoil, invOaxis, Os[i])
    #
    #     spinaxis = VTKtools.rotation_matrix(0.0, 0.0, 90.0)
    #     # invspinaxis = inv(spinaxis)
    #     airfoil = transform(airfoil, spinaxis, [proploc[1],proploc[2],proploc[3]])
    #
    #
    #     push!(crosssections, airfoil)
    #
    #     push!(strains, p_strain[i,:])
    #     push!(stress1, stresslayer[1,i,:])
    #     push!(stress2, stresslayer[2,i,:])
    #     push!(stressplot, stress[i,:])
    #     push!(buckle, p_c_bucklingplot[i,:])
    #     push!(point_datas, [j for j in npoints*(i-1)+1:npoints*i])
    # end


    # Generates cells in VTK Legacy format
    points, vtk_cells, strain = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=strains)
    _, _, stresses1 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stress1)
    _, _, stresses2 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stress2)
    _, _, stresses6 = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=stressplot)
    _, _, bucklecon = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=buckle)
    _, _, point_datas = VTKtools.multilines2vtkmulticells(crosssections, sections;point_datas=point_datas)
    # Formats the point data for generateVTK
    data = []
    push!(data, Dict(
    "field_name" => "Strain",
    "field_type" => "scalar",
    "field_data" => strain))

    push!(data, Dict(
    "field_name" => "uni_stress",
    "field_type" => "scalar",
    "field_data" => stresses1))

    push!(data, Dict(
    "field_name" => "weave_stress",
    "field_type" => "scalar",
    "field_data" => stresses2))

    push!(data, Dict(
    "field_name" => "stress",
    "field_type" => "scalar",
    "field_data" => stresses6))


    push!(data, Dict(
    "field_name" => "Buckling_Constraint",
    "field_type" => "scalar",
    "field_data" => bucklecon))

    push!(data, Dict(
    "field_name" => "PointIndex",
    "field_type" => "scalar",
    "field_data" => point_datas))

    # Generates the vtk file
    VTKtools.generateVTK("$(file_name)prop$(propnum)bladenum$(bladenum).$printiter", points; cells=vtk_cells, point_data=data, keep_points=false)


end


function transform(V::Array{Float64,1},
    M::Array{Float64,2}, T::Array{Float64,1})
    return M*(V-T)+T
end

function transform(Vs::Array{Array{Float64,1},1},
    M::Array{Float64,2}, T::Array{Float64,1})
    out = Array{Float64,1}[]
    for V in Vs
        push!(out, transform(V, M, T))
    end
    return out
end


#Generate prop with all airfoil points
# y_pos = p_yloc#Rtip*p_radii # span position of each section
# x_pos = p_xloc#[y*tan(sweep[i]) for (i,y) in enumerate(y_pos)]
# z_pos = p_zloc#[y*tan(dihedral[i]) for (i,y) in enumerate(y_pos)]
# # PARAMETERS
# r = fill(1.0,length(p_radii))          # Expansion ratio in both surfaces of each airfoil
#
# sections = fill(Tuple{Float64,Int64,Float64,Bool}[(1.0, 10, 1.0, false)],length(p_radii)-1)
#
# # Leading edge position of each airfoil
# Os = [ [x_pos[i], y_pos[i], z_pos[i]] for i in 1:size(p_airfoilx)[1]]
# # Orientation of chord of each airfoil (yaw, pitch, roll)
# orien = zeros(length(p_twist_d),3) #twist already included
# for i = 1:length(orien[:,1])
#     orien[i,:] = [p_twist_d[i],p_sweep_d[i],270.0]
# end
#
#
# crosssections = []        # It will store here the cross sections for lofting
# point_datas = []
#
# # Processes each airfoil geometry
# for i = 1:length(p_airfoilx[:,1])
#
#     # Read airfoil file
#     x = p_airfoilx[i,:]
#     y = p_airfoily[i,:]
#
#     # Scales the airfoil acording to its chord length
#     new_x = chord[i]*x#new_x
#     new_y = chord[i]*y#new_y
#
#     # plot_airfoil(new_x, new_y; style=styles[i], label="airfoil")
#     if plots
#         figure("wing_airfoil")
#         plot(new_x,new_y)
#     end
#
#     # Reformats into points
#     npoints = size(new_x)[1]
#     airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]
#
#     # Positions the airfoil along the blade in the right p_orientation
#     Oaxis = VTKtools.rotation_matrix(orien[i,1], orien[i,2], orien[i,3])
#     invOaxis = inv(Oaxis)
#     airfoil = VTKtools.countertransform(airfoil, invOaxis, Os[i])
#
#     push!(crosssections, airfoil)
#     push!(point_datas, [j for j in npoints*(i-1)+1:npoints*i])
#
# end
#
# # Generates cells in VTK Legacy format
# points, vtk_cells, pt_idx = VTKtools.multilines2vtkmulticells(crosssections, sections;
# point_datas=point_datas)
#
# data = []
# push!(data, Dict(
# "field_name" => "Point_index",
# "field_type" => "scalar",
# "field_data" => pt_idx
# )
# )
#
# # Generates the vtk file for geometry
# VTKtools.generateVTK("$(file_name)_geom$printiter", points; cells=vtk_cells, point_data=data, keep_points=false)
