using ApproxFun

#x = Fun(identity,[0.,10.])
#f = sin(x^2)
#g = cos(x)

c=Chebyshev([0,0.08]) #Define Chebyshev domain in this range (to match data imported)

#pts=points(c,100) #This is a built in function; Approxfun would then do eg. vals=cos(pts) to get the values
#show(pts)

# Standard two column data form
#df = readdlm("B.cgs.kT-.0259.nh3ch3")
df=readdlm("jagged.dat")
pts=df[:,1] # Points
vals=df[:,2] #Values at these points

cfs=ApproxFun.transform(c,vals)
# Are you ready for the magic?
me=Fun(cfs,c) # me is now an ApproxFun representation of the tabulated data. 
# As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.

using Plots
plot(me)
plot(df[:,1],df[:,2])
