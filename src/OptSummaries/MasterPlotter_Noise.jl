
path,_ = splitdir(@__FILE__)
import Interpolations
using PyPlot
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

include("$path/../setup.jl") # for types
# rc("figure", figsize=(4.0, 2.5))
rc("figure", figsize=(3.5, 3.2))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
# rc("figure.subplot", left=0.17, bottom=0.21, top=0.97, right=.65)
rc("figure.subplot", left=0.22, bottom=0.32, top=0.99, right=.95)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
close("all")

numprops = [2,4,8,16,32]

cr_factor = 100#reverse([50,100, 150])
db_factor = reverse([76,65,63,60])

db_factor_plot = reverse(["<76","65","63","60"])

DETAILED = false

detailed_output_accel = []
detailed_output_climb = []
detailed_output_cruise = []
design_accel = []
design_climb = []
design_cruise = []
detailed_system = []

ranges = zeros(length(numprops),length(db_factor))
CLs_accel = zeros(length(numprops),length(db_factor))
CLs_climb = zeros(length(numprops),length(db_factor))
CLs_cruise = zeros(length(numprops),length(db_factor))

Drags_accel = zeros(length(numprops),length(db_factor))
Drags_climb = zeros(length(numprops),length(db_factor))
Drags_cruise = zeros(length(numprops),length(db_factor))

CDs_accel = zeros(length(numprops),length(db_factor))
CDs_climb = zeros(length(numprops),length(db_factor))
CDs_cruise = zeros(length(numprops),length(db_factor))

total_weight = zeros(length(numprops),length(db_factor))
climb_thrust = zeros(length(numprops),length(db_factor))
gama = zeros(length(numprops),length(db_factor))

ETAs_accel = zeros(length(numprops),length(db_factor))
ETAs_climb = zeros(length(numprops),length(db_factor))
ETAs_cruise = zeros(length(numprops),length(db_factor))

BattMasses_accel = zeros(length(numprops),length(db_factor))
BattMasses_climb = zeros(length(numprops),length(db_factor))
BattMasses_cruise = zeros(length(numprops),length(db_factor))

Power_accel = zeros(length(numprops),length(db_factor))
Power_climb = zeros(length(numprops),length(db_factor))
Power_cruise = zeros(length(numprops),length(db_factor))

DBs_accel = zeros(length(numprops),length(db_factor))
DBs_climb = zeros(length(numprops),length(db_factor))
DBs_cruise = zeros(length(numprops),length(db_factor))

Thrust_accel = zeros(length(numprops),length(db_factor))
Thrust_climb = zeros(length(numprops),length(db_factor))
Thrust_cruise = zeros(length(numprops),length(db_factor))

Velocity_accel = zeros(length(numprops),length(db_factor))
Velocity_climb = zeros(length(numprops),length(db_factor))
Velocity_cruise = zeros(length(numprops),length(db_factor))

Rtips = zeros(length(numprops),length(db_factor))

Propulsion_Mass = zeros(length(numprops),length(db_factor))

takeoff_dist = zeros(length(numprops),length(db_factor))
dist_climb = zeros(length(numprops),length(db_factor))
dist_trans = zeros(length(numprops),length(db_factor))
dist_accel = zeros(length(numprops),length(db_factor))

tip_mach_accel = zeros(length(numprops),length(db_factor))
tip_mach_climb = zeros(length(numprops),length(db_factor))
tip_mach_cruise = zeros(length(numprops),length(db_factor))

aoa_accel = zeros(length(numprops),length(db_factor))
aoa_climb = zeros(length(numprops),length(db_factor))
aoa_cruise = zeros(length(numprops),length(db_factor))

mass_prop = zeros(length(numprops),length(db_factor))
blades_prop = zeros(length(numprops),length(db_factor))

span = zeros(length(numprops),length(db_factor))


mean_chord = zeros(length(numprops),length(db_factor))

for i = 1:length(numprops)
    for j = 1:length(db_factor)
    # try
        filename = "$(path)/p$(numprops[i])_cr$(cr_factor)_db$(db_factor[j])correct_span.jld"
        summary = JLD.load(filename)

        detailed_output_accel = summary["detailed_output_accel"]
        detailed_output_climb = summary["detailed_output_climb"]
        detailed_output_cruise = summary["detailed_output_cruise"]
        design_accel = summary["design_accel"]
        design_climb = summary["design_climb"]
        design_cruise = summary["design_cruise"]
        detailed_system = summary["detailed_system"]

        ranges[i,j] = detailed_system.range_true
        CLs_accel[i,j] = detailed_output_accel.Lift/(0.5*1.225*design_accel.velocity^2*design_accel.w_halfspan*2*design_accel.w_chord_pts[1]) #TODO: rho, Sref
        CLs_climb[i,j] = detailed_output_climb.Lift/(0.5*1.225*design_climb.velocity^2*design_climb.w_halfspan*2*design_climb.w_chord_pts[1]) #TODO: rho, Sref
        CLs_cruise[i,j] = detailed_output_cruise.Lift/(0.5*1.225*design_cruise.velocity^2*design_cruise.w_halfspan*2*design_cruise.w_chord_pts[1]) #TODO: rho, Sref

        CDs_accel[i,j] = detailed_output_accel.Drag/(0.5*1.225*design_accel.velocity^2*design_accel.w_halfspan*2*design_accel.w_chord_pts[1]) #TODO: rho, Sref
        CDs_climb[i,j] = detailed_output_climb.Drag/(0.5*1.225*design_climb.velocity^2*design_climb.w_halfspan*2*design_climb.w_chord_pts[1]) #TODO: rho, Sref
        CDs_cruise[i,j] = detailed_output_cruise.Drag/(0.5*1.225*design_cruise.velocity^2*design_cruise.w_halfspan*2*design_cruise.w_chord_pts[1]) #TODO: rho, Sref

        total_weight[i,j] = detailed_system.total_weight

        climb_thrust[i,j] = detailed_output_climb.total_thrust
        gama[i,j] = detailed_system.gama

        ETAs_accel[i,j] = detailed_output_accel.ETA_total
        ETAs_climb[i,j] = detailed_output_climb.ETA_total
        ETAs_cruise[i,j] = detailed_output_cruise.ETA_total

        BattMasses_accel[i,j] = detailed_system.batt_mass_accel
        BattMasses_climb[i,j] = detailed_system.batt_mass_climb
        BattMasses_cruise[i,j] = detailed_system.batt_mass_cruise

        Rtips[i,j] = design_cruise.p_Rtip

        span[i,j] = design_cruise.w_halfspan*2

        Propulsion_Mass[i,j] = max(detailed_output_accel.propwing_mass,detailed_output_climb.propwing_mass,detailed_output_cruise.propwing_mass)

        Power_accel[i,j] = detailed_output_accel.total_thrust*design_accel.velocity
        Power_climb[i,j] = detailed_output_climb.total_thrust*design_climb.velocity
        Power_cruise[i,j] = detailed_output_cruise.total_thrust*design_cruise.velocity

        Thrust_accel[i,j] = detailed_output_accel.total_thrust
        Thrust_climb[i,j] = detailed_output_climb.total_thrust
        Thrust_cruise[i,j] = detailed_output_cruise.total_thrust

        Drags_accel[i,j] = detailed_output_accel.Drag
        Drags_climb[i,j] = detailed_output_climb.Drag
        Drags_cruise[i,j] = detailed_output_cruise.Drag

        Velocity_accel[i,j] = design_accel.velocity
        Velocity_climb[i,j] = design_climb.velocity
        Velocity_cruise[i,j] = design_cruise.velocity

        DBs_accel[i,j] = maximum(detailed_output_accel.db_test)
        DBs_climb[i,j] = maximum(detailed_output_climb.db_test)
        DBs_cruise[i,j] = maximum(detailed_output_cruise.db_test)

        takeoff_dist[i,j] = detailed_system.dist_climb+detailed_system.dist_accel
        dist_accel[i,j] = detailed_system.dist_accel
        dist_trans[i,j] = detailed_system.xtr
        dist_climb[i,j] = detailed_system.x_left

        tip_mach_accel[i,j] = design_accel.p_Rtip*design_accel.p_rpm/60*2*pi/detailed_output_accel.a
        tip_mach_climb[i,j] = design_climb.p_Rtip*design_climb.p_rpm/60*2*pi/detailed_output_climb.a
        tip_mach_cruise[i,j] = design_cruise.p_Rtip*design_cruise.p_rpm/60*2*pi/detailed_output_cruise.a

        aoa_accel[i,j] = design_accel.w_aoa
        aoa_climb[i,j] = design_climb.w_aoa
        aoa_cruise[i,j] = design_cruise.w_aoa

        mean_chord[i,j] = mean(design_cruise.p_chord_pts)

        mass_prop[i,j] = detailed_output_cruise.mass_prop*design_cruise.numprops*design_cruise.blades #same for all three
        blades_prop[i,j] = design_cruise.blades #same for all three


        if DETAILED
            figname = "cl_p_db$(numprops[i])_cr$(cr_factor)_db$(db_factor[j])"
            figure(figname)
            # maxy = zeros(length(db_factor))
            plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),detailed_output_climb.cllocal,".-",label = "local cl")
            plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),detailed_output_climb.cl,".-",label = "freestream cl")
            xlabel("wing halfspan")
            ylabel("cl")
            # ylim([0,1.05*ceil(Int,maximum(maxy))])
            # legend(loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
            # savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

        end

    # end
    end
end

# High level results
figname = "Takeoff_Speed_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],takeoff_dist[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(takeoff_dist[:,i])
        plot(numprops,takeoff_dist[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Takeoff Distance (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Range_Props_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],ranges[1:end-1,i]/1000.0,"o-",label="$(db_factor_plot[i])")#"1/$(cr_factor[i])")
    else
        maxy[i] = maximum(ranges[:,i]/1000.0)
        plot(numprops,ranges[:,i]/1000.0,"o-",label="$(db_factor_plot[i])")#"1/$(cr_factor[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Range (km)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Accel_Speed_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],dist_accel[1:end-1,i],"o-",label="$(db_factor_plot[i])")#"1/$(cr_factor[i])")
    else
        maxy[i] = maximum(dist_accel[:,i])
        plot(numprops,dist_accel[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Ground Roll Distance (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "span_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],span[1:end-1,i],"o-",label="$(db_factor_plot[i])")#"1/$(cr_factor[i])")
    else
        maxy[i] = maximum(span[:,i])
        plot(numprops,span[:,i],"o-",label="$(db_factor_plot[i])")#"1/$(cr_factor[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Span (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "CL_climb_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],CLs_climb[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(CLs_climb[:,i])
        plot(numprops,CLs_climb[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Lift Coefficient (CL)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "CD_cruise_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],CDs_cruise[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(CDs_cruise[:,i])
        plot(numprops,CDs_cruise[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Drag Coefficient (CD)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Noise_climb_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],DBs_climb[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(DBs_climb[:,i])
        plot(numprops,DBs_climb[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Sound Pressure Level (dBa)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "FlightPathAngle_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],gama[1:end-1,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(gama[:,i]*180/pi)
        plot(numprops,gama[:,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Flight Path Angle (degrees)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Total_Mass_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],total_weight[1:end-1,i]/9.81,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(total_weight[:,i]/9.81)
        plot(numprops,total_weight[:,i]/9.81,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Total Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Propulsion_Mass_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Propulsion_Mass[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Propulsion_Mass[:,i])
        plot(numprops,Propulsion_Mass[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


Propulsion_Mass

figname = "Rtip_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Rtips[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Rtips[:,i])
        plot(numprops,Rtips[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Radius (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Battery_Cruise_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],BattMasses_cruise[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(BattMasses_cruise[:,i])
        plot(numprops,BattMasses_cruise[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Battery Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Battery_Total_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],BattMasses_cruise[1:end-1,i]+BattMasses_accel[1:end-1,i]+BattMasses_climb[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(BattMasses_cruise[:,i]+BattMasses_accel[:,i]+BattMasses_climb[:,i])
        plot(numprops,BattMasses_cruise[:,i]+BattMasses_accel[:,i]+BattMasses_climb[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Battery Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Propulsion_ETA_Cruise_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],ETAs_cruise[1:end-1,i]*100.0,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(ETAs_cruise[:,i]*100.0)
        plot(numprops,ETAs_cruise[:,i]*100.0,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Propulsion_ETA_Accel_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],ETAs_accel[1:end-1,i]*100.0,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(ETAs_accel[:,i]*100.0)
        plot(numprops,ETAs_accel[:,i]*100.0,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Propulsion_ETA_Climb_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],ETAs_climb[1:end-1,i]*100.0,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(ETAs_climb[:,i]*100.0)
        plot(numprops,ETAs_climb[:,i]*100.0,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Climb_Power_In_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Power_climb[1:end-1,i]./ETAs_climb[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Power_climb[:,i]./ETAs_climb[:,i])
        plot(numprops,Power_climb[:,i]./ETAs_climb[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Climb Power In (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Climb_Power_Out_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Power_climb[1:end-1,i]/1000.0,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Power_climb[:,i]/1000.0)
        plot(numprops,Power_climb[:,i]/1000.0,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Climb Power Out (kW)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Climb_Thrust_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Thrust_climb[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Thrust_climb[:,i])
        plot(numprops,Thrust_climb[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Total Thrust (N)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Climb_Thrust2Weight_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Thrust_climb[1:end-1,i]./(total_weight[1:end-1,i]),"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Thrust_climb[:,i]./(total_weight[:,i]))
        plot(numprops,Thrust_climb[:,i]./(total_weight[:,i]),"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Thrust/Weight")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Climb_NetThrust2Weight_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],(Thrust_climb[1:end-1,i]-Drags_climb[1:end-1,i])./(total_weight[1:end-1,i]),"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum((Thrust_climb[:,i]-Drags_climb[:,i])./(total_weight[:,i]))
        plot(numprops,(Thrust_climb[:,i]-Drags_climb[:,i])./(total_weight[:,i]),"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Net-Thrust/Weight")
ylim([0,0.9])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Takeoff_Energy_Percent_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    # if i = || i==2= 1 #Don't plot 16 props for tight db constraint, very bad
    #     to_percent = 100[1:end-1].0*1(end-1BattMasses_accel[:,i]+BattMasses_climb[:,i])./(BattMasses_accel[:,i]+BattMasses_climb[:,i]+BattMasses_cruise[:,i])
    # else
        to_percent = 100.0*(BattMasses_accel[:,i]+BattMasses_climb[:,i])./(BattMasses_accel[:,i]+BattMasses_climb[:,i]+BattMasses_cruise[:,i])
    # end
    semilogy(numprops,to_percent,"o-",label="$(db_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Takeoff Energy (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Velocity_cruise_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Velocity_cruise[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Velocity_cruise[:,i])
        plot(numprops,Velocity_cruise[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Velocity (m/s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Velocity_cruiseMPH_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Velocity_cruise[1:end-1,i]*2.237,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Velocity_cruise[:,i]*2.237)
        plot(numprops,Velocity_cruise[:,i]*2.237,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Velocity (mph)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "tip_mach_climb_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],tip_mach_climb[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(tip_mach_climb[:,i])
        plot(numprops,tip_mach_climb[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Tip Speed (mach)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "aoa_cruise_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],aoa_cruise[1:end-1,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(aoa_cruise[:,i]*180/pi)
        plot(numprops,aoa_cruise[:,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Angle of Attack (deg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "aoa_accel_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],aoa_accel[1:end-1,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(aoa_accel[:,i]*180/pi)
        plot(numprops,aoa_accel[:,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Angle of Attack (deg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "aoa_climb_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],aoa_climb[1:end-1,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(aoa_climb[:,i]*180/pi)
        plot(numprops,aoa_climb[:,i]*180/pi,"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Angle of Attack (deg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "prop_mean_chord_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],mean_chord[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(mean_chord[:,i])
        plot(numprops,mean_chord[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Mean Propeller Chord (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "mass_prop_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],mass_prop[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(mass_prop[:,i])
        plot(numprops,mass_prop[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "blades_prop_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],blades_prop[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(blades_prop[:,i])
        plot(numprops,blades_prop[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Blades (#)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Velocity_accel_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=1:length(db_factor)
    if false #i == 1 || i==2 #Don't plot 16 props for tight db constraint, very bad
        plot(numprops[1:end-1],Velocity_accel[1:end-1,i],"o-",label="$(db_factor_plot[i])")
    else
        maxy[i] = maximum(Velocity_accel[:,i])
        plot(numprops,Velocity_accel[:,i],"o-",label="$(db_factor_plot[i])")
    end
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Liftoff Speed (m/s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Relative_Distance_cr100_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=2#1:length(db_factor)

    plot(numprops,dist_accel[:,i],"o-",label="Accel")

    plot(numprops,dist_trans[:,i],"o-",label="Transition")

    plot(numprops,dist_climb[:,i],"o-",label="Climb")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Distance (m)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Percent_Relative_Distance_db65_db"
figure(figname)
maxy = zeros(length(db_factor))
for i=2#1:length(db_factor)

    plot(numprops,dist_accel[:,i]./takeoff_dist[:,i]*100,"o-",label="Accel")

    plot(numprops,dist_trans[:,i]./takeoff_dist[:,i]*100,"o-",label="Transition")

    plot(numprops,dist_climb[:,i]./takeoff_dist[:,i]*100,"o-",label="Climb")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Relative Distance (%)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "distance_noise_8"
figure(figname)
maxy = zeros(length(db_factor))
distance_noise_8 = zeros(length(db_factor))
for i=1:length(db_factor)
    distance_noise_8[i] = takeoff_dist[3,i]
end
db_factor[end] = round(Int,DBs_climb[3,4])
plot(db_factor,distance_noise_8,"ko-")
xticks(db_factor)
yticks([0.0;round(distance_noise_8[[1,2,4]],-1)])
xlabel("Sound Pressure Level (dBa)")
ylabel("Takeoff Distance (m)")
ylim([0,1.05*ceil(Int,1.1*maximum(distance_noise_8))])
# legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "distance_noise_16"
figure(figname)
maxy = zeros(length(db_factor))
distance_noise_16 = zeros(length(db_factor))
for i=1:length(db_factor)
    distance_noise_16[i] = takeoff_dist[4,i]
end
db_factor[end] = round(Int,DBs_climb[4,4])
plot(db_factor,distance_noise_16,"ko-")
xticks(db_factor)
yticks([0.0;round(distance_noise_16[[1,2,4]],-1)])
xlabel("Sound Pressure Level (dBa)")
ylabel("Takeoff Distance (m)")
ylim([0,1.05*ceil(Int,1.1*maximum(distance_noise_16))])
# legend(title = "Noise Constraint (dBa)",loc="upper center", ncol = length(db_factor), bbox_to_anchor=(0.4, -.20))#loc="center left", bbox_to_anchor=(1, 0.5))
# savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)
