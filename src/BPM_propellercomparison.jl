using BPM

using PyPlot
rc("figure", figsize=(3.5, 3.5))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.15, bottom=0.05, top=0.97, right=.93)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
import Atmosphere
# function runme()
meshgrid(x,y) = (repmat(x',length(y),1),repmat(y,1,length(x)))

# INITIALIZE
nturbines = 8
x_test = ones(nturbines)*0.0 # x-locations of turbines (m)
y_test = ones(nturbines)*0.0 # y-locations of turbines (m)
obs_test = [0., 20.0, 1.5] # x-, y-, and z-location of the observer (m)
winddir_test = 0. # wind direction (deg)
rpm_test = ones(nturbines)*5000.0 # rotation rate of the tubrines (rpm)
windvel_test = ones(nturbines)*1.0 # wind velocity (m/s)
h_test = 15.24 # height of the turbine hub (m)
noise_corr = 0.8697933840957954 # correction factor for noise

alpha =  [5.10514, 5.10514, 6.8703, 6.96041, 6.07495, 4.58821, 3.20019, 3.20019]#[10.0,5.0,1.0]
R_station = [0.0778218, 0.177878, 0.277935, 0.377991, 0.478048, 0.578104, 0.678161, 0.778218]#[0.0,0.125,0.25]*2
p_chord = [0.155659, 0.19989, 0.210314, 0.172972, 0.125581, 0.0843676, 0.0460056, 0.0100891]#[.01,.05,.01]
AR = R_station[end]/mean(p_chord) #TODO: actual AR #10.0 #TODO: actual AR
psi = 14.0 # solid angle (deg)

area = pi*R_station[end]^2 * nturbines

Rtip_orig = R_station[end]

rho,mu,a = Atmosphere.atmospherefit(0.0) #Use average altitude for each evaluation


blades = 2

tipspeed = collect(linspace(0.1,0.8,10))
omega_correct_speed = tipspeed*a/R_station[end]
RPM_correct_speed = omega_correct_speed * 60 / (2*pi)
numprops = [2,4,8,16,32]#collect(linspace(2,32,5))
nBlades = collect(linspace(2,5,length(numprops)))
numprops2 = [2,4,8,16,32]
db_test1 = zeros(length(numprops),length(tipspeed))
checkrpm = zeros(length(numprops),length(tipspeed))

levels = [60,65,70,72.5]

X,Y=meshgrid(tipspeed,numprops)

for j = 1:length(tipspeed)
    for i = 1:length(numprops)
        nturbines = 8#round(Int,numprops[i])
        Rtip_scale = sqrt(area/nturbines/pi)/Rtip_orig
        blades = 4#round(Int,nBlades[i])
        R_station = [0.0778218, 0.177878, 0.277935, 0.377991, 0.478048, 0.578104, 0.678161, 0.778218]*Rtip_scale
        p_chord = [0.155659, 0.19989, 0.210314, 0.172972, 0.125581, 0.0843676, 0.0460056, 0.0100891]*Rtip_scale

        omega_correct_speed = tipspeed[j]*a/R_station[end]
        RPM_correct_speed = omega_correct_speed * 60 / (2*pi)

        x_test = ones(nturbines)*0.0 # x-locations of turbines (m)
        y_test = ones(nturbines)*0.0 # y-locations of turbines (m)
        windvel_test = ones(nturbines)*1.0 # wind velocity (m/s)
        AR = R_station[end]/mean(p_chord) #TODO: actual AR #10.0 #TODO: actual AR


        rpm_test = ones(nturbines)*RPM_correct_speed # rotation rate of the tubrines (rpm)
        checkrpm[i,j] = RPM_correct_speed
        db_test1[i,j] = BPM.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, blades, h_test, R_station, p_chord, p_chord*0.25, alpha, mu/rho, a, psi, AR, noise_corr)
        println(mean(p_chord))
    end
    # println(tipspeed[j])

end

println(maximum(db_test1))

close("all")
figname = "Noise_Tipspeed_Contour_scaledchord"
fig = figure(figname)
CS = PyPlot.contourf(Y,X,db_test1,interpolation = "spline16",origin = "lower",cmap = PyPlot.cm[:Blues],levels = linspace(28,maximum(db_test1),1000))#linspace(28,77,1000))
PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "SPL (dBa)",pad = 0.2,aspect = 20,ticks = [30,40,50,60,70])
CS2 = PyPlot.contour(Y,X,db_test1,colors = "white",levels = levels)#linspace(28,77,1000))
clabel( CS2,fmt = "%1.1f")
xticks(numprops2)
PyPlot.xlabel("Propellers (#)")
PyPlot.ylabel("Tipspeed (mach)")
# xlim([Inlet_pos,Outlet_pos])
# ylim([Bot_pos,Top_pos])
PyPlot.savefig("../MDAO_paper/figures/$figname.png",transparent = true,dpi = 600)

db_test2 = zeros(length(numprops),length(tipspeed))

for j = 1:length(tipspeed)
    for i = 1:length(numprops)
        nturbines = round(Int,numprops[i])
        Rtip_scale = sqrt(area/nturbines/pi)/Rtip_orig
        blades = round(Int,nBlades[i])
        R_station = [0.0778218, 0.177878, 0.277935, 0.377991, 0.478048, 0.578104, 0.678161, 0.778218]*Rtip_scale
        p_chord = [0.155659, 0.19989, 0.210314, 0.172972, 0.125581, 0.0843676, 0.0460056, 0.0100891]#*Rtip_scale

        omega_correct_speed = tipspeed[j]*a/R_station[end]
        RPM_correct_speed = omega_correct_speed * 60 / (2*pi)

        x_test = ones(nturbines)*0.0 # x-locations of turbines (m)
        y_test = ones(nturbines)*0.0 # y-locations of turbines (m)
        windvel_test = ones(nturbines)*1.0 # wind velocity (m/s)
        AR = R_station[end]/mean(p_chord) #TODO: actual AR #10.0 #TODO: actual AR


        rpm_test = ones(nturbines)*RPM_correct_speed # rotation rate of the tubrines (rpm)
        checkrpm[i,j] = RPM_correct_speed
        db_test2[i,j] = BPM.turbinepos(x_test, y_test, obs_test, winddir_test, windvel_test, rpm_test, blades, h_test, R_station, p_chord, p_chord*0.25, alpha, mu/rho, a, psi, AR, noise_corr)

    end
    println(tipspeed[j])
end


figname = "Noise_Tipspeed_Contour_constchord"
fig = figure(figname)
CS = PyPlot.contourf(Y,X,db_test2,interpolation = "spline16",origin = "lower",cmap = PyPlot.cm[:Blues],levels = linspace(28,maximum(db_test1),1000))#linspace(28,77,1000))
PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "SPL (dBa)",pad = 0.2,aspect = 20,ticks = [30,40,50,60,70])
CS2 = PyPlot.contour(Y,X,db_test2,colors = "white",levels = levels)#linspace(28,77,1000))
clabel( CS2,fmt = "%1.1f")
xticks(numprops2)
PyPlot.xlabel("Propellers (#)")
PyPlot.ylabel("Tipspeed (mach)")
# xlim([Inlet_pos,Outlet_pos])
# ylim([Bot_pos,Top_pos])
PyPlot.savefig("../MDAO_paper/figures/$figname.png",transparent = true,dpi = 600)


# figname = "CheckRPM"
# figure(figname)
# CS = PyPlot.contourf(Y,X,checkrpm,interpolation = "spline16",origin = "lower",cmap = PyPlot.cm[:Blues],levels = linspace(minimum(checkrpm),maximum(checkrpm),1000))
# PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "SPL (dba)",pad = 0.2,aspect = 20,ticks = round.(linspace(minimum(checkrpm),maximum(checkrpm),8),1))
# xticks(numprops2)
# PyPlot.xlabel("Propellers (#)")
# PyPlot.ylabel("Tipspeed (mach)")
# # xlim([Inlet_pos,Outlet_pos])
# # ylim([Bot_pos,Top_pos])
# # PyPlot.savefig("../MDAO_paper/figures/$figname.pdf",transparent = true)
