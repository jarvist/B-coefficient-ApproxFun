using ApproxFun

#x = Fun(identity,[0.,10.])
#f = sin(x^2)
#g = cos(x)

c=Chebyshev([0,0.08]) #Define Chebyshev domain in this range (to match data imported)

#pts=points(c,100) #This is a built in function; Approxfun would then do eg. vals=cos(pts) to get the values
#show(pts)

# Standard two column data form
df = readdlm("B.cgs.kT-.0259.nh3ch3")
#df=readdlm("jagged.dat") # Only 8 points, very sawtoothy
pts=df[:,1] # Points
vals=df[:,2] #Values at these points

# From: https://github.com/ApproxFun/ApproxFun.jl/issues/275 , courtesy of private communication with Sheehan Olver
# Least squares approximation of data on an evenly spaced grid with Chebyshev series

function vandermonde(S,n,x::AbstractVector)
    V=Array(Float64,length(x),n)
    for k=1:n
        V[:,k]=Fun([zeros(k-1);1],S)(x)
    end
    V
end

#x=collect(0:0.001:1.)
#S=Chebyshev([0.,1.])
#V=vandermonde(S,20,x)
#v=cos(cos(x))
#f=Fun(V\v,S)

# For ...(this)... case, make sure `length(pts) >> n`.
n=13
V=vandermonde(c,n,pts)
me=Fun(V\vals,c)


#cfs=ApproxFun.transform(c,vals)
# Are you ready for the magic?
#me=Fun(cfs,c) # me is now an ApproxFun representation of the tabulated data. 
# As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.

using Plots

# Logscale version...
plot(me,label="Chebyshev B")
plot!(df[:,1],df[:,2],label="Raw B")
yaxis!("B coeff")
xaxis!("Density",:log10)

# Linear version...
#plot(me,label="Chebyshev B")
#plot!(df[:,1],df[:,2],label="Raw B")
#yaxis!("B coeff")
#xaxis!("Density)
