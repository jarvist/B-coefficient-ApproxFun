using ApproxFun

#x = Fun(identity,[0.,10.])
#f = sin(x^2)
#g = cos(x)

c=Chebyshev([0,0.09])

pts=points(c,100)
show(pts)

# Standard two column data form
df = readdlm("B.cgs.kT-.0259.nh3ch3")
pts=df[:,1] # Points
vals=df[:,2] #Values at these points

cfs=ApproxFun.transform(c,vals)
# Are you ready for the magic?
me=Fun(cfs,c) # me is now an ApproxFun representation of the tabulated data. 
# As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.


