path,_ = splitdir(@__FILE__)
using PyPlot
import JLD


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
include("$path/../oldvlm.jl")
import BPM
import AirfoilPrep

printiter = 4600
#TODO verify this works with parallel


include("$path/../objective.jl")
include("$path/../propstructmisc.jl")
include("$path/../propwingmisc.jl")
include("$path/../setup.jl")
include("$path/../save_vtk.jl")



detailed_output_accel = []
detailed_output_climb = []
detailed_output_cruise = []
design_accel = []
design_climb = []
design_cruise = []
detailed_system = []

rc("figure", figsize=(60*1.8, 60))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
# rc("figure", dpi=300)
rc("axes.spines", right=false, top=false)
rc("axes", linewidth=15)
# rc("figure.subplot", left=0.27, bottom=0.2, top=0.8, right=.95)
rc("figure.subplot", left=0.01, bottom=0.01, top=0.99, right=.99)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

ioff()

index = 70
failnum = 0
imagenum = 0
iterations = zeros(1)
objective1 = ones(1)*100
objectiveColor = ones(1)
indexcorrection = 0
background_transp = ones(1)
objcol = 0.0

while index<51521

    # try

close("all")

println("imagenum: $imagenum index: $index")
# success = true
# while success


if index == 2170
    index = 21070
    indexcorrection = 18900
end

filename = "$(path)/Journal_GIF/p8_cr100_db76correct_span$(index).jld"
summary = JLD.load(filename)
image = imread("$(path)/Journal_PNG/blade.$(lpad(imagenum,4,0)).png")
# success = false
# failnum = 0
# catch
#     failnum += 1
#     index += 1
#
# end

detailed_output_accel = summary["detailed_output_accel"]
detailed_output_climb = summary["detailed_output_climb"]
detailed_output_cruise = summary["detailed_output_cruise"]
design_accel = summary["design_accel"]
design_climb = summary["design_climb"]
design_cruise = summary["design_cruise"]
detailed_system = summary["detailed_system"]
obj = summary["objective"]

G = linspace(0,detailed_system.gama,100)-pi/2

xtr = detailed_system.dist_accel+detailed_system.radius*cos.(G)
ytr = detailed_system.radius*sin.(G)+detailed_system.radius

pathx = [0;detailed_system.dist_accel;xtr;detailed_system.dist_accel+detailed_system.dist_climb] # dist_climb includes detailed_system.xtr
plotx = [0,detailed_system.dist_accel,detailed_system.dist_accel+detailed_system.xtr,detailed_system.dist_accel+detailed_system.dist_climb]
pathy = [0;0;ytr; 15.24]

iterations = append!(iterations, iterations[end]+index-indexcorrection)

if pathx[end]<1.0

    objective1 = append!(objective1, objective1[end])

else

    objective1 = append!(objective1, pathx[end])

end

w_height = 0

imagenum += 1
index += 70

figname = "Journal_GIF$(lpad(imagenum,4,0))"
figure(figname)

# SETUP
axes()
axis("off")

#TOP TEXTS

# text(0.04, 0.95, "Props (#): 2 4    16 32     Noise (dBa): 60 65          Cruise (MPH): 50        150",
# verticalalignment="bottom", horizontalalignment="left",
# color=(0.0, 0.0, 0.0, 1.0), fontsize=180)
#
# text(0.212, 0.95, "8",
# verticalalignment="bottom", horizontalalignment="left",
# color=(0.0, 0.0, 0.0, 1.0), fontsize=180,
# bbox=Dict("facecolor"=>"none", "edgecolor"=>"black", "boxstyle"=>"round,pad=0.1","linewidth"=>"15"))
#
# text(0.57, 0.95, "70",
# verticalalignment="bottom", horizontalalignment="left",
# color=(0.0, 0.0, 0.0, 1.0), fontsize=180,
# bbox=Dict("facecolor"=>"none", "edgecolor"=>"black", "boxstyle"=>"round,pad=0.1","linewidth"=>"15"))
#
# text(0.843, 0.95, "100",
# verticalalignment="bottom", horizontalalignment="left",
# color=(0.0, 0.0, 0.0, 1.0), fontsize=180,
# bbox=Dict("facecolor"=>"none", "edgecolor"=>"black", "boxstyle"=>"round,pad=0.1","linewidth"=>"15"))
#


con_names = ["Weight",
"Noise",
"Range",
"Lift",
"Thrust",
"Tip Separation",
"Tip Speed",
"Stall",
"Cruise Speed",
"Composite Failure",
"Buckling",
"Battery Capacity"]

con_dic_names = [:totalmass,
:noise_accel,
:range,
:Lift_cruise,
:Drag_cruise,
:radiussize,
:machtip_accel,
:cl_max_accel,
:cruisespeed,
:p_stress_climb,
:p_buckling_climb,
:capacity_climb]

# for names in keys(detailed_system.optconvals)
#     println(names)
# end

constraint_y = collect(range(.82,-.07,length(con_names)))

objcol = 0

#CONSTRAINTS
for i = 1:length(con_names)
    try #catch if one of the constraints listed above is not in the list here
        background_transp = maximum(detailed_system.optconvals[con_dic_names[i]])*100
    catch
        background_transp = 0.0
    end

    # if background_transp< 0.0
    #     background_transp = 0.0
    # elseif background_transp > 1.0
    #     background_transp = 1.0
    # end

    if objcol < background_transp
        objcol = background_transp
    end

    coolwarm_constr_color = PyPlot.cm[:coolwarm](background_transp)

    # println(background_transp)
    #Background
    text(0.915, constraint_y[i], con_names[i],
    verticalalignment="bottom", horizontalalignment="center",
    color=Tuple(coolwarm_constr_color), fontsize=130,fontweight = 1000)
    # Actual text
    # text(0.915, constraint_y[i], con_names[i],
    # verticalalignment="bottom", horizontalalignment="center",
    # color=(0.0, 0.0, 0.0, 1.0), fontsize=140,fontweight = 540.0,)
end

objectiveColor = append!(objectiveColor, objcol)

#LABEL TEXTS
text(0.1, 0.91, "Props On Wing",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=160)


text(0.55, 0.90, "Takeoff Profile",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=160)

text(0.44, 0.468, "Optimization Progression",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=160)


text(0.915, 0.91, "Constraints",
verticalalignment="bottom", horizontalalignment="center",
color=(0.0, 0.0, 0.0, 1.0), fontsize=170)




# VTK Image
# Scaled by 1/2 when outputting in paraview
axes([-0.050, 0.04, 0.35, 0.68])
imshow(image)

text(210, 0, "Blade Design",
verticalalignment="bottom", horizontalalignment="center",
color=(0.0, 0.0, 0.0, 1.0), fontsize=160)

xoff = 40

text(285-xoff, 50, "Composite Failure",
verticalalignment="bottom", horizontalalignment="center",
color=(0.0, 0.0, 0.0, 1.0), fontsize=130)#,rotation = "vertical")

text(275-xoff, 321, "Buckling Failure",
verticalalignment="bottom", horizontalalignment="center",
color=(0.0, 0.0, 0.0, 1.0), fontsize=130)#,rotation = "vertical")

# TOP Bar Annotation
center = 178
offset = 95

text(205-xoff, center-offset, "Violated",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=95)

text(205-xoff, center, "Constrained",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=95)

text(205-xoff, center+offset, "Satisfied",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=95)

# BOTTOM Bar Annotation
center = 445
offset = 95

text(205-xoff, center-offset, "Violated",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=95)

text(205-xoff, center, "Constrained",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=95)

text(205-xoff, center+offset, "Satisfied",
verticalalignment="bottom", horizontalalignment="left",
color=(0.0, 0.0, 0.0, 1.0), fontsize=95)



axis("off")



# Wing and Props
axes([0.0, 0.56, 0.4, 0.6])
x_tips = [0.0,1/8,2/8,3/8,4/8,5/8,6/8,7/8,1.0]
x = ((x_tips[2:end]-x_tips[1:end-1])/2+x_tips[1:end-1])*11.4
y = ones(length(x))*w_height
areas = ones(length(x))*pi*design_accel.p_Rtip^2*85000
# areas = ones(length(x))*pi*.7125^2*85000 #tip to tip for tuning scaling factor
scatter(x,y,s=areas,alpha=0.5)
plot([0.1,11.3],ones(2)*w_height,linewidth=75.0)
axis("off")



# Takeoff
axes([0.47, 0.655, .34, 0.238])
plot(pathx[1:2],pathy[1:2],"-",linewidth=25.0,color = "#1a455e")
plot(pathx[2:102],pathy[2:102],"-",linewidth=25.0,color = "#28698f")
plot(pathx[102:103],pathy[102:103],"-",linewidth=25.0,color = "#348ABD")
tick_params(which="major",length=10*7,width=2*7,labelsize=16*7)
xlabel("Distance (m)",fontsize=140)
xticks(plotx, rotation = 35)
ylabel("Height (m)",fontsize=140)
yticks([0,15.24],["0","15.24\n (50 ft)"],rotation = 90)
ylim([0.0,20.0])
xlim([0.0,100.0])
axis("equal")


coolwarm_obj_color = PyPlot.cm[:coolwarm](objectiveColor)

# Objective
axes([0.38, 0.119, 0.40, 0.35])
for i = 1:length(iterations)-1
    semilogy(iterations[i:i+1]/1000/35, objective1[i:i+1],"-",c = Tuple(coolwarm_obj_color[i,:]),linewidth=25.0)
end
xlabel("Iterations (1k)",fontsize=140)
ylabel("Takeoff Distance (m)",fontsize=140)
yticks(rotation = 30)
# ylim([0,100])
# xlim([0,100])
tick_params(which="major",length=10*7,width=2*7,labelsize=16*7)
tick_params(which="minor",length=10*7,width=2*7,labelsize=16*7)






# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/Journal_GIF/Output/$figname.png")#,transparent = false)

    # catch
    #     index += 70
    #     failnum+=1
    #     println("$failnum")
    # end
end


# run(`cd Journal_GIF/Output ; pwd ; convert -delay 2 -loop 0 *.pdf Moore2018_3.gif`)
