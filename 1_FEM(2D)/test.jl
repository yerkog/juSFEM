include("LinearTriangleAssemble.jl")
include("LinearTriangleElementStiffness.jl")
include("LinearTriangleElementStresses.jl")

using LinearAlgebra
using Gadfly

#E = 2.06e11
E = 21e10
NU = 0.3
t = 0.01

stifflike1 = LinearTriangleElementStiffness(E,NU,t,0.16,0.04,0.12,0.04,0.12,0.02)
stifflike2 = LinearTriangleElementStiffness(E,NU,t,0.16,0.04,0.12,0.02,0.16,0.02)

KE = zeros(30,30)
for j = 1:3:10
    global KE=LinearTriangleAssemble(KE,stifflike1,j,j+3,j+4)
    global KE=LinearTriangleAssemble(KE,stifflike2,j,j+4,j+1)
end
for j = 2:3:11
    global KE=LinearTriangleAssemble(KE,stifflike1,j,j+3,j+4)
    global KE=LinearTriangleAssemble(KE,stifflike2,j,j+4,j+1)
end

Q = zeros(30,1)
Q[2,1] = -1000
K = zeros(30,30)
K[1:24,1:24] = KE[1:24,1:24]
K[25:30,25:30] = Matrix{Float16}(I,6,6)

d = K\Q

xx = 0:0.04:0.16;
yy = d[28:-6:4]
plot(x=xx,y=yy,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"))
