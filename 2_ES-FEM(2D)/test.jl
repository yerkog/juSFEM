include("AssembleQuad.jl")
include("AssembleTriangle.jl")
include("BInner.jl")
include("BOuter.jl")
include("Stiffness.jl")

using LinearAlgebra
using Gadfly
#using Cairo
#using Fontconfig

E = 21e10
NU = 0.3
t = 0.01

KE = zeros(30,30)

edge1 = BT(0.16,0.04,0.12,0.04,0.12,0.02)
s1 = Stiffness(edge1,0.0008/6)
for j = 1:3:10
    global KE = AssembleTriangle(KE,s1,j,j+3,j+4)
end

edge2 = BT(0.12,0,0.16,0,0.16,0.02)
s2 = Stiffness(edge2,0.0008/6)
for j = 6:3:15
    global KE = AssembleTriangle(KE,s2,j,j-3,j-4)
end

edge3 = BT(0.16,0.02,0.16,0.04,0.12,0.02)
s3 = Stiffness(edge3,0.0008/6)
for j = 2:1:3
    global KE = AssembleTriangle(KE,s3,j,j-1,j+3)
end

edge4 = BT(0,0.04,0,0.02,0.04,0.04)
s4 = Stiffness(edge4,0.0008/6)
for j = 13:1:14
    global KE = AssembleTriangle(KE,s4,j,j+1,j-3)
end

edge5 = BQ(0.16,0.04,0.12,0.04,0.12,0.02,0.16,0.02)
s5 = Stiffness(edge5,0.0008/3)
for j = 1:3:10
    global KE = AssembleQuad(KE,s5,j,j+3,j+4,j+1)
end
for j = 2:3:11
    global KE = AssembleQuad(KE,s5,j,j+3,j+4,j+1)
end

edge6 = BQ(0.12,0.02,0.12,0,0.16,0.02,0.16,0.04)
s6 = Stiffness(edge6,0.0008/3)
for j = 5:3:14
    global KE = AssembleQuad(KE,s6,j,j+1,j-3,j-4)
end

edge7 = BQ(0.12,0.04,0.08,0.02,0.12,0.02,0.16,0.04)
s7 = Stiffness(edge7,0.0008/3)
for j = 4:3:10
    global KE = AssembleQuad(KE,s7,j,j+4,j+1,j-3)
end
for j = 5:3:11
    global KE = AssembleQuad(KE,s7,j,j+4,j+1,j-3)
end

Q = zeros(30,1)
Q[2,1] = -1000
K = zeros(30,30)
K[1:24,1:24] = KE[1:24,1:24]
K[25:30,25:30] = Matrix{Float16}(I,6,6)

d = K\Q

#set_default_plot_size(20cm,8cm)
xx = 0:0.04:0.16;
yy = d[28:-6:4]
yyy = [0.0 -4.762389143492468e-6 -1.595175961665792e-5 -3.132983984866875e-5 -4.8661741894599886e-5]
c = map(i -> (((-1000*i^2)/(6*E*(t*0.04^3/12)))*(3*0.16-i)),0:0.04:0.16;)
line1 = layer(x=xx,y=yy,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="green"))
line2 = layer(x=xx,y=yyy,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="blue"))
line3 = layer(x=xx,y=c,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="red"))

set_default_plot_size(20cm,8cm)
Gadfly.with_theme(:dark) do
plot(line1,line2,line3,Guide.manual_color_key("Methods", ["ES-FEM", "FEM", "analytical solution"],["green", "blue", "red"]),Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"),Guide.title("RESULT"))
    #,Coord.cartesian(xmin=0, xmax=0.18)
    #draw(PDF("foo.pdf",20cm,8cm,dpi=1000),f)
end
