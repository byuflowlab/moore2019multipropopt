using PyPlot

N = 100
chord = 1.0
ds = .26
halfspan = 26
Vinf = 25.85 #(props are removed, so constant Vinfeff = Vinf)
Vinfeff = Vinf
Sref = 26*2

y = linspace(0,halfspan,N)
ds_ell = y[2]-y[1]
cl0 = 1.2
cl = cl0*sqrt.(1-(y/halfspan).^2)

clave_ell = mean(cl)


clopt = [1.18335, 1.18348, 1.18367, 1.1839, 1.18416, 1.18444, 1.18473, 1.18503, 1.18534, 1.18566, 1.18597, 1.18629, 1.18661, 1.18693, 1.18725, 1.18756, 1.18787, 1.18818, 1.18849, 1.18879, 1.18909, 1.18938, 1.18966, 1.18994, 1.19022, 1.19048, 1.19074, 1.191, 1.19124, 1.19148, 1.19171, 1.19193, 1.19214, 1.19234, 1.19253, 1.19271, 1.19287, 1.19303, 1.19317, 1.1933, 1.19342, 1.19352, 1.19361, 1.19368, 1.19373, 1.19376, 1.19378, 1.19377, 1.19375, 1.1937, 1.19362, 1.19352, 1.1934, 1.19324, 1.19305, 1.19283, 1.19257, 1.19228, 1.19194, 1.19156, 1.19113, 1.19065, 1.19011, 1.18951, 1.18885, 1.18811, 1.18729, 1.18639, 1.18539, 1.18429, 1.18307, 1.18172, 1.18023, 1.17857, 1.17674, 1.17471, 1.17246, 1.16994, 1.16713, 1.16399, 1.16046, 1.15648, 1.15198, 1.14687, 1.14104, 1.13435, 1.12662, 1.11765, 1.10713, 1.09472, 1.07991, 1.06206, 1.04025, 1.01321, 0.979061, 0.934965, 0.876372, 0.795405, 0.676626, 0.481338]
y2 = linspace(0,halfspan,length(clopt))
clopt_ave = mean(clopt)

CL_ell = 2*sum(cl.*Vinfeff.*chord.*ds_ell)/(Vinf.*Sref)
CL_opt = 2*sum(clopt.*Vinfeff.*chord.*ds)/(Vinf.*Sref)

PyPlot.figure()
PyPlot.plot(y,cl,"b-",label = "Elliptical cl: CL = $(round.(CL_ell,3))")
PyPlot.plot(y,ones(cl)*clave_ell,"b--",label = "Elliptical cl Average")
PyPlot.plot(y2,clopt,"r-",label = "FB results cl: CL = $(round.(CL_opt,3))")
PyPlot.plot(y2,ones(clopt)*clopt_ave,"r--",label = "FB results cl Average")
PyPlot.xlabel("Halfspan")
PyPlot.ylabel("cl")
PyPlot.legend(loc = "best")
