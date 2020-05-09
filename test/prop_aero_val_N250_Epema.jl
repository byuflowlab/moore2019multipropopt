# Propeller Aerodynamic Validation using UiUC Propeller data,


include("propwing.jl")

module SysOptInit_M
# export sysoptinit, getvar, printmaxviol

# ----- Define function for getting values within objective function ----- #
function getvar(key,x,optvars,optrange,optval,optscale)
    # Determine whether variable is a design variable, assign default if not
    if key in optvars
        value = x[optrange[key]]./optscale[key]'
    else
        value = optval[key][:,3]
    end
    # Make sure Float64 is returned instead of Array{Float64,1} where applicable
    if length(value) == 1
        value = value[1]
    end
    return value
end

# ----- Define function for setting values within objective function ----- #
function setvar(key,x,optvars,optrange,optval,optscale,value)
    # Determine whether variable is a design variable, assign default if not
    if key in optvars
        x[optrange[key]] = value.*optscale[key]'
    end

end

# ----- Define function for printing maximum constraint violations ----- #
function printmaxviol(c,name)
    if !isempty(find(c.>0.0))
        println("$name: ",maximum(c))
    end
end

function sysoptinit(optvars::Array{Symbol,1})
    # ----- Define lower, upper, and starting values ----- #
    optval = Dict{Symbol,Array{Float64,2}}()
    #---General Geometry Variables

    #---Power System Variables
    # Propeller Defaults are for the Epema N250

    #Propulsion
    optval[:velocity]  = [1.0 50.0 19.0]  #m/s
    optval[:proptwist] = [-15.0 80.0 54.9242;
    -15.0 75.0 49.6058;
    -15.0 60.0 44.6923;
    -15.0 45.0 40.1836;
    -15.0 45.0 36.0798;
    -15.0 45.0 32.3809;
    -15.0 45.0 29.0869;
    -15.0 45.0 26.1977;
    -15.0 45.0 23.7134]  #degrees # Propeller Defaults are for the Epema N250

    optval[:propchord] = [0.03 0.3 0.158788;
    0.03 0.3 0.15083;
    0.03 0.3 0.146158;
    0.03 0.3 0.145793;
    0.03 0.3 0.148978;
    0.03 0.3 0.153175;
    0.03 0.3 0.154064;
    0.03 0.3 0.145548;
    0.03 0.3 0.119747]  #chord/Rtip # Propeller Defaults are for the Epema N250


    optval[:RPMprop]   = [10.0 20000.0 4040.0] #RPM
    optval[:Rtip]      = [0.01 2.0 .203] # Radius meters, # Propeller Defaults are for the Epema N250

    # ----- Define scaling for each variable ----- #
    optscale = Dict()

    #Propulsion
    optscale[:velocity]  = [1e-3]
    optscale[:proptwist] = [1e-1]
    optscale[:propchord] = [1e3]
    optscale[:RPMprop]   = [1e-4] #RPM
    optscale[:Rtip]      = [1e0] # Radius meters
    optscale[:T_prop]      = [1e-2] # Radius meters

    # ----- Apply scaling ----- #
    for var in optvars
        optval[var] = optval[var].*optscale[var]'
    end

    # ----- Define Ranges ----- #
    optrange = Dict()
    cx = 0
    for var in optvars
        optrange[var] = cx+(1:size(optval[var],1))
        cx += size(optval[var],1)
    end
    vartotal = cx

    # ----- Assemble lower, upper, and starting optval ----- #
    lb = zeros(vartotal)
    ub = zeros(vartotal)
    x0 = zeros(vartotal)
    for var in optvars
        lb[optrange[var]] = optval[var][:,1]
        ub[optrange[var]] = optval[var][:,2]
        x0[optrange[var]] = optval[var][:,3]
    end

    # ----- Define Optimizer Options ----- #
    options = Dict{String, Any}()
    options["Derivative level"] = 0
    # options["Function precision"] = 3.00E-4
    options["Difference interval"] = 1e-5
    options["Central difference interval"] = 1e-6
    options["Iterations limit"] = 1e8
    options["Major iterations limit"] = 500
    options["Minor iterations limit"]= 1e8
    options["Major optimality tolerance"] = 1e-4
    options["Minor optimality  tolerance"] = 1e-4
    options["Major feasibility tolerance"] = 1e-4
    options["Minor feasibility tolerance"] = 1e-4
    options["Minor print level"] = 1
    options["Print frequency"] = 100
    options["Scale option"] = 2
    options["Scale tolerance"] = .95
    # options["Verify level"] = 3 #only if specifying gradients

    return x0,lb,ub,options,optrange,optval,optscale
end #function sysoptinit

end #module


# Do all calling relative to file's location
fileLoc,_ = splitdir(@__FILE__)


include("$(fileLoc)/../../hale-module-miscellany/params.jl")
push!(LOAD_PATH,"$(fileLoc)/../../hale-module-miscellany/")
using StandardAtmosphere

using CCBlade


using JLD
using PyPlot
close("all")

function prop_aero(x0,optvars,optrange,optval,optscale,T_req,af,rho,unitareas,pitch=0.0,optimize=false,numblades = 2)
    #Unpack the variables, sysoptinit takes care of which are being optimized
    theta = (SysOptInit_M.getvar(:proptwist,x0,optvars,optrange,optval,optscale)+pitch)*pi/180 #in degrees, convert to rad
    RPMp = SysOptInit_M.getvar(:RPMprop,x0,optvars,optrange,optval,optscale)
    Rtip = SysOptInit_M.getvar(:Rtip,x0,optvars,optrange,optval,optscale)
    chords = SysOptInit_M.getvar(:propchord,x0,optvars,optrange,optval,optscale)*Rtip
    velocity = SysOptInit_M.getvar(:velocity,x0,optvars,optrange,optval,optscale)

    r = linspace(.2,.911,length(chords)) * Rtip
    B = numblades
    Rhub = r[1]/2
    chord_prop = chords
    twist_prop = theta
    Vinf = velocity
    Omega = RPMp/60*2*pi
    xtowing = 0.5
    uwake, vwake, T, Q, reff = propwing.propanalysis(Rhub, Rtip, r, chord_prop, twist_prop, B, af, Vinf, Omega, rho, xtowing)
    # Tprop,Qprop,Uprop,Utprop,propvolume,Np,Tp,J_p,eff,CT,CQ = CCcall.cccall(RPMp,theta,chords,CCBlade,Rtip,r,rho,velocity,af,unitareas,numblades)
    Vhub = Vinf
    precone = 0.0
    turbine = false
    J_p = Vhub/(RPMp/60*2*Rtip)
    eff, CT, CQ = CCBlade.nondim(T, Q, Vhub, Omega, rho, Rtip, precone, turbine)
    return T,Q,uwake, vwake, 0.0, 0.0, 0.0,J_p,eff, CT, CQ
end



# push!(LOAD_PATH,"$(fileLoc)/../../Snopt.jl/src/")
# using Snopt
#--------- Initialization ---------#
optvars = [
# :proptwist;
:RPMprop;
:velocity
# :Rtip;
]

x0,lb,ub,options,optrange,optval,optscale = SysOptInit_M.sysoptinit(optvars)


#Initialize Airfoil Data
twists = SysOptInit_M.getvar(:proptwist,x0,optvars,optrange,optval,optscale)
Rtip = SysOptInit_M.getvar(:Rtip,x0,optvars,optrange,optval,optscale)
chords = SysOptInit_M.getvar(:propchord,x0,optvars,optrange,optval,optscale)*Rtip
velocity = SysOptInit_M.getvar(:velocity,x0,optvars,optrange,optval,optscale)
RPMprop = SysOptInit_M.getvar(:RPMprop,x0,optvars,optrange,optval,optscale)
r = collect(linspace(.2,.911,length(chords))) .* Rtip
J = velocity/(RPMprop/60*2*Rtip) #Vinf/(Omega/(2*pi)*Rtip*2)
alt = 222.0 #meters
numblades = 6
rho,mu,q,a = StandardAtmosphere.atmosphere(alt)
re = rho/mu.*chords.*sqrt.((RPMprop*pi/30.0*r).^2+velocity^2)

afgeomnames = ["$(fileLoc)/EpemaExperimental/R023.dat";"$(fileLoc)/EpemaExperimental/R023.dat";
"$(fileLoc)/EpemaExperimental/R0405.dat";"$(fileLoc)/EpemaExperimental/R0405.dat";
"$(fileLoc)/EpemaExperimental/R0405.dat";"$(fileLoc)/EpemaExperimental/R079.dat";
"$(fileLoc)/EpemaExperimental/R079.dat";"$(fileLoc)/EpemaExperimental/R079.dat";"$(fileLoc)/EpemaExperimental/R1.dat"]


N=100
M = 3
new_RPM=4040 #linspace(4000,20000,M)
pitch_array = [-5.0,0.0,5.0]
new_velocity = linspace(1,60,N)
J_p = zeros(N,M)
eta_p = zeros(N,M)
CT_p = zeros(N,M)
CQ_p = zeros(N,M)

SysOptInit_M.setvar(:RPMprop,x0,optvars,optrange,optval,optscale,new_RPM)

velocity = SysOptInit_M.getvar(:velocity,x0,optvars,optrange,optval,optscale)
RPMprop = SysOptInit_M.getvar(:RPMprop,x0,optvars,optrange,optval,optscale)

J = velocity/(RPMprop/60*2*Rtip)

re = rho/mu.*chords.*sqrt.((RPMprop*pi/30.0*r).^2+velocity^2)

# af,unitareas, Al, Au = CCBladeAirfoilTools.af_performance_from_geom(afgeomnames,chords,Rtip,r,J,re)
# JLD.save("$(fileLoc)/af_prop$(new_RPM).jld", "af", af, "unitareas",unitareas,"Al", Al, "Au", Au)
JLDin = JLD.load("$(fileLoc)/af_prop$(new_RPM).jld")
af = JLDin["af"]
unitareas = JLDin["unitareas"]

for j = 1:M

    T_req = 0.0
    pitch = pitch_array[j]+3 #degrees
    optimize=false
    numblades = 6
    args = (optvars,optrange,optval,optscale,T_req,af,rho,unitareas,pitch,optimize,numblades)

    for i = 1:N
        SysOptInit_M.setvar(:velocity,x0,optvars,optrange,optval,optscale,new_velocity[i])

        Tprop,Qprop,Uprop, Utprop, propvolume, Np, Tp,J_p[i,j],eta_p[i,j], CT_p[i,j], CQ_p[i,j] = prop_aero(x0,args...)

    end
end

using CSV
CT_N250_Epema_25deg = CSV.read("$(fileLoc)/EpemaExperimental/CT_N250_25deg.txt", delim = '\t',header = false)
CP_N250_Epema_25deg = CSV.read("$(fileLoc)/EpemaExperimental/CP_N250_25deg.txt", delim = '\t',header = false)
CQ_N250_Epema_25deg = copy(CP_N250_Epema_25deg)
CQ_N250_Epema_25deg[:,2] = CQ_N250_Epema_25deg[:,2]/(2*pi)
ETA_N250_Epema_25deg = CT_N250_Epema_25deg[:,2].*CT_N250_Epema_25deg[:,1]./CP_N250_Epema_25deg[:,2]

CT_N250_Epema_30deg = CSV.read("$(fileLoc)/EpemaExperimental/CT_N250_30deg.txt", delim = '\t',header = false)
CP_N250_Epema_30deg = CSV.read("$(fileLoc)/EpemaExperimental/CP_N250_30deg.txt", delim = '\t',header = false)
CQ_N250_Epema_30deg = copy(CP_N250_Epema_30deg)
CQ_N250_Epema_30deg[:,2] = CQ_N250_Epema_30deg[:,2]/(2*pi)
ETA_N250_Epema_30deg = CT_N250_Epema_30deg[:,2].*CT_N250_Epema_30deg[:,1]./CP_N250_Epema_30deg[:,2]

CT_N250_Epema_35deg = CSV.read("$(fileLoc)/EpemaExperimental/CT_N250_35deg.txt", delim = '\t',header = false)
CP_N250_Epema_35deg = CSV.read("$(fileLoc)/EpemaExperimental/CP_N250_35deg.txt", delim = '\t',header = false)
CQ_N250_Epema_35deg = copy(CP_N250_Epema_35deg)
CQ_N250_Epema_35deg[:,2] = CQ_N250_Epema_35deg[:,2]/(2*pi)
ETA_N250_Epema_35deg = CT_N250_Epema_35deg[:,2].*CT_N250_Epema_35deg[:,1]./CP_N250_Epema_35deg[:,2]


rc("figure", figsize=(6.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.15, bottom=0.18, top=0.97, right=0.72)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])


mygreen = "#009e73"
myblue = "#348abd"
myred = "#a60628"
plotColors = [myblue,myred,mygreen]

#Thrust Coefficient
PyPlot.figure()
for i = 1:M
    PyPlot.plot(J_p[:,i],CT_p[:,i], marker="o",color = plotColors[i],markeredgecolor = "k",label = "BEM $(30+pitch_array[i])deg RPM=$(new_RPM)")

end
PyPlot.plot(CT_N250_Epema_25deg[:,1],CT_N250_Epema_25deg[:,2],linestyle="-", color = plotColors[1],label = "Epema 25deg")
PyPlot.plot(CT_N250_Epema_30deg[:,1],CT_N250_Epema_30deg[:,2],linestyle="-", color = plotColors[2],label = "Epema 30deg")
PyPlot.plot(CT_N250_Epema_35deg[:,1],CT_N250_Epema_35deg[:,2],linestyle="-", color = plotColors[3],label = "Epema 35deg")
PyPlot.xlabel(L"J\, (V/(n*D_{prop})")
PyPlot.ylabel(L"CT\, (T/(\rho*n^2*{D_{prop}}^4))")
PyPlot.xlim([0,2.0])
PyPlot.ylim([0,.4])
PyPlot.legend(loc = 1)


#Torque Coefficient

PyPlot.figure()
for i = 1:M
    PyPlot.plot(J_p[:,i],CQ_p[:,i], marker="o",color = plotColors[i],markeredgecolor = "k",label = "BEM $(30+pitch_array[i])deg RPM=$(new_RPM)")

end
PyPlot.plot(CQ_N250_Epema_25deg[:,1],CQ_N250_Epema_25deg[:,2],linestyle="-", color = plotColors[1],label = "Epema 25deg")
PyPlot.plot(CQ_N250_Epema_30deg[:,1],CQ_N250_Epema_30deg[:,2],linestyle="-", color = plotColors[2],label = "Epema 30deg")
PyPlot.plot(CQ_N250_Epema_35deg[:,1],CQ_N250_Epema_35deg[:,2],linestyle="-", color = plotColors[3],label = "Epema 35deg")
PyPlot.xlabel(L"J\, (V/(n*D_{prop})")
PyPlot.ylabel(L"Cq\, (T/(\rho*n^2*{D_{prop}}^5)")
PyPlot.xlim([0,2.0])
PyPlot.ylim([0,.1])
PyPlot.legend(loc = 1)

#Efficiency
PyPlot.figure()
for i = 1:M
    PyPlot.plot(J_p[:,i],eta_p[:,i], linestyle="-",color = plotColors[i])#,markeredgecolor = "k",label = "BEM $(30+pitch_array[i])deg RPM=$(new_RPM)")

end
PyPlot.plot(CT_N250_Epema_25deg[:,1],ETA_N250_Epema_25deg,marker="o",linestyle="", color = plotColors[1])#,label = "Epema 25deg")
PyPlot.plot(CT_N250_Epema_30deg[:,1],ETA_N250_Epema_30deg,marker="o",linestyle="", color = plotColors[2])#,label = "Epema 30deg")
PyPlot.plot(CT_N250_Epema_35deg[:,1],ETA_N250_Epema_35deg,marker="o",linestyle="", color = plotColors[3])#,label = "Epema 35deg")

PyPlot.plot([],[],"ko",label = "Experimental")
PyPlot.plot([],[],"k-",label = "Blade Element Momentum")
PyPlot.xlabel("Advance Ratio")
PyPlot.ylabel("Efficiency")
# PyPlot.xlim([0,2.0])
# PyPlot.ylim([0,1.5])
PyPlot.legend(loc="center left", bbox_to_anchor=(.8, 0.8))
