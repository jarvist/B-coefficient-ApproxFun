using ApproxFun

#pts=points(c,100) #This is a built in function; Approxfun would then do eg. vals=cos(pts) to get the values
#show(pts) # Density looks like a 'U' with greater sampling density near the extremes of the range

# From: https://github.com/ApproxFun/ApproxFun.jl/issues/275 , courtesy of private communication with Sheehan Olver
# Least squares approximation of data on an evenly spaced grid with Chebyshev series
function vandermonde(S,n,x::AbstractVector)
    V=Array(Float64,length(x),n)
    for k=1:n
        V[:,k]=Fun([zeros(k-1);1],S)(x)
    end
    V
end

function loadB(filename)
    c=Chebyshev([0,0.08]) #Define Chebyshev domain in this range (to match data imported)

    # Standard two column data form
    df=readdlm(filename)
    
    pts=df[:,1] # Points
    vals=df[:,2] #Values at these points

    # For ...(this)... case, make sure `length(pts) >> n`.
    n=13 # This is a magic number, found to give a good fit to Pooya's data
    println("Attempting Vandermonde / Chebyshev fit with: Range: ",c," Points in fit: ",n," Data points:",length(pts))
    V=vandermonde(c,n,pts)
    # Are you ready for the magic?
    af=Fun(V\vals,c) # Approximate Function (af)
    # me is now an ApproxFun representation of the tabulated data. 
    # As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.
    return af,df
end

#B,df=loadB("jagged.dat") # Only 8 points, very sawtoothy
B,df=loadB("B.cgs.kT-.0259.nh3ch3")
# Pooya's data runs from densities between 0 .. 0.08
# Unit is in electrons / unit cell, so multiply by ~4.03e21 to get cm^-3

# Print out the fit vs. the raw data, so you can see the residual
function graphB(af,df)
    # Logscale version...
    plot(af,label="Chebyshev (fit) B")
    plot!(df[:,1],df[:,2],label="Tabulated (input) B")
    yaxis!("B coeff")
    xaxis!("Density")
    png("Bcoeff_linear.png")
    
    xaxis!(:log10)
    png("Bcoeff_log.png")
end

using Plots
default(show=true)
default(size=(1024,768)) # of rendered PNGs

graphB(B,df)

# OK; we have an ApproxFun function (B) fitted to the tabulated data

#Herz values from DOI: 10.1021/acs.accounts.5b00411
A=5e6
Bconst=0.9e-10
# B is our Approxfun fit; internally it's a polynomial, but you can differentiate etc. as if it were analytic

# We are now going to build our Ordinary Differential Equation model for n(t)
using ODE

type Model
    label
    ODE::Function # ODE of model
    t # holds calculated data
    xv
end

Models = [
    Model("Bcoeff",   (t,n) -> [ - B(n[1]/4.03e21) * n[1]*n[1] ; n[1]],[],[]), # Pure bimolecular; Pooya B(n)
    Model("Bcoeff-A", (t,n) -> [-A*n[1] - B(n[1]/4.03e21) * n[1]*n[1] ; n[1]],[],[]), # SRH 'A' term from above; Pooya B(n)
    Model("Bconst-A", (t,n) -> [-A*n[1] - Bconst * n[1]*n[1]; n[1]],[],[]), # SRH 'A' term and bimolecular fit from above
    Model("SRH-A",    (t,n) -> [-A*n[1]; n[1] ], [], []) # SRH 'A' term only
]

using Plots

function plotsoln(t,xv)
	
    plot(t,xv[:,1]) # Time on x-axis, versus n[1] (density) on Y axis
	yaxis!("Density n (cm^-3)")
	xaxis!("Time (s)")
	png("density_linear.png")
	yaxis!(:log10)
	png("density_log.png")

    plot(t,xv[:,2]) # Time on x-axis, versus n[2] (integrating dn/dt, emission) on Y axis
	yaxis!("Integrated Emission (???)")
	xaxis!("Time (s)")
	png("integrated_emission.png")

    I(n) = -B(n/4.03e21) * n*n
	intensity=[-I(x) for x in xv[:,1]] #extract intensity as a function of time, by feeding the solved densities (from the ODE) into the solver
	# Nb: probably not the most elegant way to do this (!)
	intensity=hcat(intensity...).'

	# calculates intensity reusing the same functional toolkit - NB: assumes all recombination is emissive
	plot(t,intensity[:,1])
	yaxis!("Emission Intensity")
	xaxis!("Time (s)")
	png("emission.png")
	yaxis!(:log10)
	png("emission_log.png")
end

# Initial vector; density is first part
initial=[0.08*4.03e21, 0.] # Start at n=0.08 Pooyas

for model in Models
	println("Simulating model: ",model.label)
    model.t,model.xv=ode23(model.ODE,initial,[0.;2e-8]) # Numerically Integrate from 0 to ... seconds
	model.xv=hcat(model.xv...).'
    #plotsoln(model.t,model.xv)
end

println("Plotting charge density vs. time")
plot()
yaxis!("Density n (cm^-3)")
xaxis!("Time (s)")
for model in Models
    plot!(model.t,model.xv[:,1],label=model.label)
end
png("density_linear.png")

yaxis!(:log10)
png("density_log.png")

println("Plotting integrated emission.")
plot()
for model in Models
    plot(model.t,model.xv[:,2],label=model.label) # Time on x-axis, versus n[2] (integrating dn/dt, emission) on Y axis
end
yaxis!("Integrated Emission (???)")
xaxis!("Time (s)")
png("integrated_emission.png")

println("Simulating radiative emission (TRPL) and plotting...")
plot()

I(n) = -B(n/4.03e21) * n*n # Emission model
for model in Models
	intensity=[-I(x) for x in model.xv[:,1]] #extract intensity as a function of time, by feeding the solved densities (from the ODE) into the solver
	# Nb: probably not the most elegant way to do this (!)
	intensity=hcat(intensity...).'

	# calculates intensity reusing the same functional toolkit - NB: assumes all recombination is emissive
	plot!(model.t,intensity[:,1],label=model.label)
end

yaxis!("Emission Intensity")
xaxis!("Time (s)")
png("emission.png")
yaxis!(:log10)
png("emission_log.png")

