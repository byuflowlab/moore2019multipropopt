
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
rc("figure", figsize=(3.5, 3.20))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
# rc("figure.subplot", left=0.27, bottom=0.2, top=0.8, right=.95)
rc("figure.subplot", left=0.22, bottom=0.32, top=0.99, right=.95)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
close("all")

numprops = [2,4,8,16,32]

cr_factor = ([50,100, 150])

cr_factor_plot = ([">22","45", "67"])

db_factor = 76 #[1,2,8,32,128]

DETAILED = false

detailed_output_accel = []
detailed_output_climb = []
detailed_output_cruise = []
design_accel = []
design_climb = []
design_cruise = []
detailed_system = []

ranges = zeros(length(numprops),length(cr_factor))
CLs_accel = zeros(length(numprops),length(cr_factor))
CLs_climb = zeros(length(numprops),length(cr_factor))
CLs_cruise = zeros(length(numprops),length(cr_factor))

CDs_accel = zeros(length(numprops),length(cr_factor))
CDs_climb = zeros(length(numprops),length(cr_factor))
CDs_cruise = zeros(length(numprops),length(cr_factor))

Drags_accel = zeros(length(numprops),length(cr_factor))
Drags_climb = zeros(length(numprops),length(cr_factor))
Drags_cruise = zeros(length(numprops),length(cr_factor))

total_weight = zeros(length(numprops),length(cr_factor))
climb_thrust = zeros(length(numprops),length(cr_factor))
gama = zeros(length(numprops),length(cr_factor))

ETAs_accel = zeros(length(numprops),length(cr_factor))
ETAs_climb = zeros(length(numprops),length(cr_factor))
ETAs_cruise = zeros(length(numprops),length(cr_factor))

ETAs_prop_accel = zeros(length(numprops),length(cr_factor))
ETAs_prop_climb = zeros(length(numprops),length(cr_factor))
ETAs_prop_cruise = zeros(length(numprops),length(cr_factor))

BattMasses_accel = zeros(length(numprops),length(cr_factor))
BattMasses_climb = zeros(length(numprops),length(cr_factor))
BattMasses_cruise = zeros(length(numprops),length(cr_factor))

MotorMasses_accel = zeros(length(numprops),length(cr_factor))
MotorMasses_climb = zeros(length(numprops),length(cr_factor))
MotorMasses_cruise = zeros(length(numprops),length(cr_factor))

Power_accel = zeros(length(numprops),length(cr_factor))
Power_climb = zeros(length(numprops),length(cr_factor))
Power_cruise = zeros(length(numprops),length(cr_factor))

Power_IN_accel = zeros(length(numprops),length(cr_factor))
Power_IN_climb = zeros(length(numprops),length(cr_factor))
Power_IN_cruise = zeros(length(numprops),length(cr_factor))

Power_Shaft_accel = zeros(length(numprops),length(cr_factor))
Power_Shaft_climb = zeros(length(numprops),length(cr_factor))
Power_Shaft_cruise = zeros(length(numprops),length(cr_factor))

DBs_accel = zeros(length(numprops),length(cr_factor))
DBs_climb = zeros(length(numprops),length(cr_factor))
DBs_cruise = zeros(length(numprops),length(cr_factor))

Thrust_accel = zeros(length(numprops),length(cr_factor))
Thrust_climb = zeros(length(numprops),length(cr_factor))
Thrust_cruise = zeros(length(numprops),length(cr_factor))

Velocity_accel = zeros(length(numprops),length(cr_factor))
Velocity_climb = zeros(length(numprops),length(cr_factor))
Velocity_cruise = zeros(length(numprops),length(cr_factor))

Rtips = zeros(length(numprops),length(cr_factor))

R_transition = zeros(length(numprops),length(cr_factor))

Propulsion_Mass = zeros(length(numprops),length(cr_factor))

takeoff_dist = zeros(length(numprops),length(cr_factor))
dist_climb = zeros(length(numprops),length(cr_factor))
dist_trans = zeros(length(numprops),length(cr_factor))
dist_accel = zeros(length(numprops),length(cr_factor))

tip_mach_accel = zeros(length(numprops),length(cr_factor))
tip_mach_climb = zeros(length(numprops),length(cr_factor))
tip_mach_cruise = zeros(length(numprops),length(cr_factor))

aoa_accel = zeros(length(numprops),length(cr_factor))
aoa_climb = zeros(length(numprops),length(cr_factor))
aoa_cruise = zeros(length(numprops),length(cr_factor))

t_accel = zeros(length(numprops),length(cr_factor))
t_transition = zeros(length(numprops),length(cr_factor))
t_climb = zeros(length(numprops),length(cr_factor))
t_cruise = zeros(length(numprops),length(cr_factor))

mass_prop = zeros(length(numprops),length(cr_factor))

span = zeros(length(numprops),length(cr_factor))


mean_chord = zeros(length(numprops),length(cr_factor))

blades_prop = zeros(length(numprops),length(cr_factor))

for i = 1:length(numprops)
    for j = 1:length(cr_factor)
    # try
        filename = "$(path)/p$(numprops[i])_cr$(cr_factor[j])_db$(db_factor)correct_span.jld"
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

        Drags_accel[i,j] = detailed_output_accel.Drag
        Drags_climb[i,j] = detailed_output_climb.Drag
        Drags_cruise[i,j] = detailed_output_cruise.Drag

        total_weight[i,j] = detailed_system.total_weight

        R_transition[i,j] = detailed_system.radius

        climb_thrust[i,j] = detailed_output_climb.total_thrust
        gama[i,j] = detailed_system.gama

        ETAs_accel[i,j] = detailed_output_accel.ETA_total
        ETAs_climb[i,j] = detailed_output_climb.ETA_total
        ETAs_cruise[i,j] = detailed_output_cruise.ETA_total

        ETAs_prop_accel[i,j] = detailed_output_accel.ETA_p
        ETAs_prop_climb[i,j] = detailed_output_climb.ETA_p
        ETAs_prop_cruise[i,j] = detailed_output_cruise.ETA_p

        BattMasses_accel[i,j] = detailed_system.batt_mass_accel
        BattMasses_climb[i,j] = detailed_system.batt_mass_climb
        BattMasses_cruise[i,j] = detailed_system.batt_mass_cruise

        MotorMasses_accel[i,j] = detailed_output_accel.mass_motor
        MotorMasses_climb[i,j] = detailed_output_climb.mass_motor
        MotorMasses_cruise[i,j] = detailed_output_cruise.mass_motor

        Rtips[i,j] = design_cruise.p_Rtip

        span[i,j] = design_cruise.w_halfspan*2

        Propulsion_Mass[i,j] = max(detailed_output_accel.propwing_mass,detailed_output_climb.propwing_mass,detailed_output_cruise.propwing_mass)

        Power_accel[i,j] = detailed_output_accel.total_thrust*design_accel.velocity
        Power_climb[i,j] = detailed_output_climb.total_thrust*design_climb.velocity
        Power_cruise[i,j] = detailed_output_cruise.total_thrust*design_cruise.velocity

        Power_IN_accel[i,j] = detailed_output_accel.V_m*detailed_output_accel.I_m*numprops[i]
        Power_IN_climb[i,j] = detailed_output_climb.V_m*detailed_output_climb.I_m*numprops[i]
        Power_IN_cruise[i,j] = detailed_output_cruise.V_m*detailed_output_cruise.I_m*numprops[i]

        Power_Shaft_accel[i,j] = Power_IN_accel[i,j] * detailed_output_accel.ETA_m
        Power_Shaft_climb[i,j] = Power_IN_climb[i,j] * detailed_output_climb.ETA_m
        Power_Shaft_cruise[i,j] = Power_IN_cruise[i,j] * detailed_output_cruise.ETA_m

        Thrust_accel[i,j] = detailed_output_accel.total_thrust
        Thrust_climb[i,j] = detailed_output_climb.total_thrust
        Thrust_cruise[i,j] = detailed_output_cruise.total_thrust

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

        t_accel[i,j] = detailed_system.t_accel
        t_transition[i,j] = detailed_system.t_transition
        t_climb[i,j] = detailed_system.t_climb
        t_cruise[i,j] = detailed_system.t_cruise

        mean_chord[i,j] = mean(design_cruise.p_chord_pts)

        mass_prop[i,j] = detailed_output_cruise.mass_prop*design_cruise.numprops*design_cruise.blades #same for all three

        blades_prop[i,j] = design_cruise.blades #same for all three


        if DETAILED
            splchord = Dierckx.Spline1D(design_climb.w_nondim_halfspany,design_climb.w_chord_pts,k = 1)
            splchordpts = splchord(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y))
            figname = "cl_p$(numprops[i])_cr$(cr_factor[j])"
            figure(figname)
            plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),detailed_output_climb.cllocal,".-",label = "local cl")
            plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),detailed_output_climb.cl,".-",label = "freestream cl")
            xlabel("wing halfspan")
            ylabel("cl")
            # legend(loc="center left", bbox_to_anchor=(1, 0.5))
            savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

            figname = "lift_p$(numprops[i])_cr$(cr_factor[j])"
            figure(figname)
            # plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),detailed_output_climb.cllocal,".-",label = "local cl")
            plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),detailed_output_climb.cl*0.5*1.225*design_climb.velocity^2.0.*splchordpts,".-",label = "freestream cl")
            xlabel("wing halfspan")
            ylabel("local lift")
            # legend(loc="center left", bbox_to_anchor=(1, 0.5))
            savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

            figname = "localq_p$(numprops[i])_cr$(cr_factor[j])"
            figure(figname)
            plot(detailed_output_climb.panels.y/maximum(detailed_output_climb.panels.y),0.5*1.225*detailed_output_climb.Vinfeff.^2,".-",label = "local q")
            xlabel("wing halfspan")
            ylabel("q")
            # legend(loc="center left", bbox_to_anchor=(1, 0.5))
            savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

            figname = "chord_p$(numprops[i])_cr$(cr_factor[j])"
            figure(figname)
            plot(design_climb.w_nondim_halfspany*design_climb.w_halfspan,design_climb.w_chord_pts,".-",label = "local cl")
            xlabel("wing halfspan (meters)")
            ylabel("chord (meters)")

            # legend(loc="center left", bbox_to_anchor=(1, 0.5))
            # ylim([0.1.2*0,maximum(design_climb.maxy1])
            savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



            figname = "geometricaoa_p$(numprops[i])_cr$(cr_factor[j])"
            figure(figname)
            plot(linspace(0,1,length(detailed_output_climb.alphaeff)),design_climb.w_aoa*ones(detailed_output_climb.alphaeff)*180/pi,".-",label = "geometric aoa")
            plot(linspace(0,1,length(detailed_output_climb.alphaeff)),detailed_output_climb.alphaeff*180/pi,".-",label = "effective aoa")
            plot(linspace(0,1,length(detailed_output_climb.alphaeff)),(detailed_output_climb.alphaeff-design_climb.w_aoa)*180/pi,".-",label = "effective - geometric aoa")
            xlabel("wing halfspan")
            ylabel("aoa (deg)")
            legend(loc="center left", bbox_to_anchor=(1, 0.5))
            savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


        end

    # end
    end
end

# High level results
figname = "Takeoff_Speed"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(takeoff_dist[:,i])
    plot(numprops,takeoff_dist[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Takeoff Distance (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Accel_Speed"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(dist_accel[:,i])
    plot(numprops,dist_accel[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Ground Roll Distance (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Range_Props"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(ranges[:,i]/1000.0)
    plot(numprops,ranges[:,i]/1000.0,"o-",label="$(cr_factor_plot[i])")#"1/$(cr_factor[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Range (km)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "span"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(span[:,i])
    plot(numprops,span[:,i],"o-",label="$(cr_factor_plot[i])")#"1/$(cr_factor[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Span (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "CL_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(CLs_climb[:,i])
    plot(numprops,CLs_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Lift Coefficient (CL)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "CL_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(CLs_cruise[:,i])
    plot(numprops,CLs_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Lift Coefficient (CL)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "CD_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(CDs_cruise[:,i])
    plot(numprops,CDs_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Drag Coefficient (CD)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Noise_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(DBs_climb[:,i])
    plot(numprops,DBs_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Sound Pressure Level (dBa)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "FlightPathAngle"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(gama[:,i]*180/pi)
    plot(numprops,gama[:,i]*180/pi,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Flight Path Angle (degrees)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Total_Mass"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(total_weight[:,i])
    plot(numprops,total_weight[:,i]/9.81,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Total Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Propulsion_Mass"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Propulsion_Mass[:,i])
    plot(numprops,Propulsion_Mass[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Motor_Mass"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)


    maxy[i] = maximum(MotorMasses_climb[:,i].*numprops)
    plot(numprops,MotorMasses_climb[:,i].*numprops,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Motor Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "Rtip"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Rtips[:,i])
    plot(numprops,Rtips[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Radius (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Battery_Cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(BattMasses_cruise[:,i])
    plot(numprops,BattMasses_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Battery Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Battery_Total"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(BattMasses_cruise[:,i])
    plot(numprops,BattMasses_cruise[:,i]+BattMasses_accel[:,i]+BattMasses_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Battery Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Battery_Mass_Fraction"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum((BattMasses_cruise[:,i]+BattMasses_accel[:,i]+BattMasses_climb[:,i])./(total_weight[:,i]/9.81)*100)
    plot(numprops,(BattMasses_cruise[:,i]+BattMasses_accel[:,i]+BattMasses_climb[:,i])./(total_weight[:,i]/9.81)*100,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Battery Mass Fraction (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "Propulsion_ETA_Cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(ETAs_cruise[:,i]*100.0)
    plot(numprops,ETAs_cruise[:,i]*100.0,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Propulsion_ETA_Accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(ETAs_accel[:,i]*100.0)
    plot(numprops,ETAs_accel[:,i]*100.0,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Propulsion_ETA_Climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(ETAs_climb[:,i]*100.0)
    plot(numprops,ETAs_climb[:,i]*100.0,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propulsion Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Propeller_ETA_Climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(ETAs_prop_climb[:,i]*100.0)
    plot(numprops,ETAs_prop_climb[:,i]*100.0,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Propeller_ETA_Cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(ETAs_prop_cruise[:,i]*100.0)
    plot(numprops,ETAs_prop_cruise[:,i]*100.0,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Efficiency (%)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "Climb_Thrust"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Thrust_climb[:,i])
    plot(numprops,Thrust_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Total Thrust (N)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Climb_Thrust2Weight"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Thrust_climb[:,i])
    plot(numprops,Thrust_climb[:,i]./(total_weight[:,i]),"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Thrust/Weight")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Climb_NetThrust2Weight"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum((Thrust_climb[:,i]-Drags_climb[:,i])./(total_weight[:,i]))
    plot(numprops,(Thrust_climb[:,i]-Drags_climb[:,i])./(total_weight[:,i]),"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Net-Thrust/Weight")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Transition_Radius"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(R_transition[:,i])
    plot(numprops,R_transition[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Transition Radius (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "Takeoff_Energy_Percent"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    to_percent = 100.0*(BattMasses_accel[:,i]+BattMasses_climb[:,i])./(BattMasses_accel[:,i]+BattMasses_climb[:,i]+BattMasses_cruise[:,i])
    semilogy(numprops,to_percent,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Takeoff Energy (%)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Velocity_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Velocity_cruise[:,i])
    plot(numprops,Velocity_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Velocity (m/s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Velocity_cruiseMPH"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Velocity_cruise[:,i])
    plot(numprops,Velocity_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Velocity (mph)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "time_accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(t_accel[:,i])
    plot(numprops,t_accel[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Acceleration Time (s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "time_transition"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(t_transition[:,i])
    plot(numprops,t_transition[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Transition Time (s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "time_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(t_climb[:,i])
    plot(numprops,t_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Climb Time (s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "time_tclimb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(t_climb[:,i]+t_transition[:,i])
    plot(numprops,t_climb[:,i]+t_transition[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Transition \& Climb Time (s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "time_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(t_cruise[:,i])
    plot(numprops,t_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Time (s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "tip_mach_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(tip_mach_climb[:,i])
    plot(numprops,tip_mach_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Tip Speed (mach)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "aoa_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(aoa_cruise[:,i])
    plot(numprops,aoa_cruise[:,i]*180/pi,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Angle of Attack (deg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "aoa_accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(aoa_accel[:,i])
    plot(numprops,aoa_accel[:,i]*180/pi,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Angle of Attack (deg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "aoa_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(aoa_climb[:,i])
    plot(numprops,aoa_climb[:,i]*180/pi,"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Wing Angle of Attack (deg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "prop_mean_chord"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(mean_chord[:,i])
    plot(numprops,mean_chord[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Mean Propeller Chord (m)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "mass_prop"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(mass_prop[:,i])
    plot(numprops,mass_prop[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Mass (kg)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "blades_prop"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(blades_prop[:,i])
    plot(numprops,blades_prop[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Propeller Blades (#)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Power_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_climb[:,i])
    plot(numprops,Power_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Climb Power Out (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Power_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_cruise[:,i])
    plot(numprops,Power_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Power Out (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Power_accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_accel[:,i])
    plot(numprops,Power_accel[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Acceleration Power Out (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Power_IN_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_IN_climb[:,i])
    plot(numprops,Power_IN_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Climb Power IN (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Power_IN_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_IN_cruise[:,i])
    plot(numprops,Power_IN_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Power IN (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Power_IN_accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_IN_accel[:,i])
    plot(numprops,Power_IN_accel[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Acceleration Power IN (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Power_Shaft_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_Shaft_climb[:,i])
    plot(numprops,Power_Shaft_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Climb Power Shaft (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Power_Shaft_cruise"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_Shaft_cruise[:,i])
    plot(numprops,Power_Shaft_cruise[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Cruise Power Shaft (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Power_Shaft_accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Power_Shaft_accel[:,i])
    plot(numprops,Power_Shaft_accel[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Acceleration Power Shaft (watts)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)



figname = "Velocity_accel"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Velocity_accel[:,i])
    plot(numprops,Velocity_accel[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Liftoff Speed (m/s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Velocity_climb"
figure(figname)
maxy = zeros(length(cr_factor))
for i=1:length(cr_factor)
    maxy[i] = maximum(Velocity_climb[:,i])
    plot(numprops,Velocity_climb[:,i],"o-",label="$(cr_factor_plot[i])")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Liftoff Speed (m/s)")
ylim([0,1.05*ceil(Int,maximum(maxy))])
# legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, 1.30))#loc="center left", bbox_to_anchor=(1, 0.5))
legend(title = "Cruise Constraint (m/s)",loc="upper center", ncol = length(cr_factor), bbox_to_anchor=(0.5, -.20))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)

figname = "Relative_Distance_cr100"
figure(figname)
for i=2#1:length(cr_factor)
    plot(numprops,dist_accel[:,i],"o-",label="Accel")
    plot(numprops,dist_trans[:,i],"o-",label="Transition")
    plot(numprops,dist_climb[:,i],"o-",label="Climb")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Distance (m)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(loc="center left", bbox_to_anchor=(1, 0.5))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)


figname = "Percent_Relative_Distance_cr100"
figure(figname)
for i=2#1:length(cr_factor)
    plot(numprops,dist_accel[:,i]./takeoff_dist[:,i]*100,"o-",label="Accel")
    plot(numprops,dist_trans[:,i]./takeoff_dist[:,i]*100,"o-",label="Transition")
    plot(numprops,dist_climb[:,i]./takeoff_dist[:,i]*100,"o-",label="Climb")
end
xticks(numprops)
xlabel("Propellers (#)")
ylabel("Relative Distance (%)")
# ylim([0,1.05*ceil(Int,maximum(maxy))])
legend(loc="center left", bbox_to_anchor=(1, 0.5))
savefig("$path/../../MDAO_Paper/figures/$figname.pdf",transparent = true)
