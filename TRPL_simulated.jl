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
    #df = readdlm("B.cgs.kT-.0259.nh3ch3")
    #df=readdlm("jagged.dat") # Only 8 points, very sawtoothy
    pts=df[:,1] # Points
    vals=df[:,2] #Values at these points

    # For ...(this)... case, make sure `length(pts) >> n`.
    n=13
    V=vandermonde(c,n,pts)
    # Are you ready for the magic?
    af=Fun(V\vals,c) # Approximate Function (af)
    # me is now an ApproxFun representation of the tabulated data. 
    # As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.
    return af,df
end

B,df=loadB("B.cgs.kT-.0259.nh3ch3")

function graphB(af,df)
    # Logscale version...
    plot(af,label="Chebyshev B")
    plot!(df[:,1],df[:,2],label="Raw B")
    yaxis!("B coeff")
    xaxis!("Density",:log10)

    # Linear version...
    plot(af,label="Chebyshev B")
    plot!(df[:,1],df[:,2],label="Raw B")
    yaxis!("B coeff")
    xaxis!("Density")
end

#using Plots
#graphB(B,df)

### OK, solving the ODE
#i=Fun([0.,20.])
#i0=0.i
#B0=B(0) # Constant value approx
#I=n->[n(0.)-B0,n(20.), n' + B0*n^2 ]
#n=newton(I,i0)

BScale=1e10 # scale from Pooya's internal units to SI

using ODE
function I(t, n)
  [n[2]; - BScale*B(n[2]) * n[2]*n[2] ]
end

initial=[0., 0.8]
T,xv=ode23(I,initial,[0.;1e-6]) # Integrate from 0 to 1 microsecond
xv=hcat(xv...).'

using Plots
plot(xv[:,1],xv[:,2]) # Probably a more elegant way of getting an XY plot!
