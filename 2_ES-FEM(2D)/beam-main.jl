include("AssembleQuad.jl")
include("AssembleTriangle.jl")
include("BInner.jl")
include("BOuter.jl")
include("Stiffness.jl")

using LinearAlgebra
using Gadfly
using CSV
using DataFrames
#using Cairo
#using Fontconfig

E = 21e10
NU = 0.3
t = 0.01

KE = zeros(110,110)

edge1 = BT(0.4,0.08,0.36,0.08,0.36,0.06)
s1 = Stiffness(edge1,0.0008/6)
for j = 1:5:46
    global KE = AssembleTriangle(KE,s1,j,j+5,j+6)
end
edge2 = BT(0.36,0,0.4,0,0.4,0.02)
s2 = Stiffness(edge2,0.0008/6)
for j = 10:5:55
    global KE = AssembleTriangle(KE,s2,j,j-5,j-6)
end
edge3 = BT(0.4,0.06,0.4,0.08,0.36,0.06)
s3 = Stiffness(edge3,0.0008/6)
for j = 2:1:5
    global KE = AssembleTriangle(KE,s3,j,j-1,j+5)
end
edge4 = BT(0,0.08,0,0.06,0.04,0.08)
s4 = Stiffness(edge4,0.0008/6)
for j = 51:1:54
    global KE = AssembleTriangle(KE,s4,j,j+1,j-5)
end

edge5 = BQ(0.4,0.08,0.36,0.08,0.36,0.06,0.4,0.06)
s5 = Stiffness(edge5,0.0008/3)
for j = 1:5:46
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
for j = 2:5:47
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
for j = 3:5:48
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
for j = 4:5:49
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
edge6 = BQ(0.36,0.06,0.36,0.04,0.4,0.06,0.4,0.08)
s6 = Stiffness(edge6,0.0008/3)
for j = 7:5:52
    global KE = AssembleQuad(KE,s6,j,j+1,j-5,j-6)
end
for j = 8:5:53
    global KE = AssembleQuad(KE,s6,j,j+1,j-5,j-6)
end
for j = 9:5:54
    global KE = AssembleQuad(KE,s6,j,j+1,j-5,j-6)
end
edge7 = BQ(0.36,0.08,0.32,0.06,0.36,0.06,0.4,0.08)
s7 = Stiffness(edge7,0.0008/3)
for j = 6:5:46
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end
for j = 7:5:47
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end
for j = 8:5:48
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end
for j = 9:5:49
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end

#=
TEMPQ = zeros(110,1)
TEMPQ[2,1] = -4000

TEMPMatrix = zeros(110,110)
TEMPMatrix[1:100,1:100] = KE[1:100,1:100]
TEMPMatrix[101:110,101:110] = Matrix{Float16}(I,10,10)

d=TEMPMatrix\TEMPQ

XX = collect(0:0.04:0.4)
Y1 = d[106:-10:6]
Y2 = [0.0 -1.0822437722258944e-5 -4.035503570515363e-5 -8.588051757808588e-5 -0.00014535720223828857 -0.00021680303958405118 -0.0002982248351603411 -0.00038762920996542776 -0.0004830328435555055 -0.0005824293893356926 -0.0006833792103817733]
Y3 = map(i -> (((-4000*i^2)/(6*E*(t*0.08^3/12)))*(3*0.4-i)),0:0.04:0.4;)
line1 = layer(x=XX,y=Y1,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="green"))
line2 = layer(x=XX,y=Y2,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="blue"))
line3 = layer(x=XX,y=Y3,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="red"))

plot(x=XX,y=YY,Geom.point,Geom.smooth(method=:loess,smoothing=0.9),
Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"))
set_default_plot_size(20cm,8cm)
Gadfly.with_theme(:default) do
f = plot(line1,line2,line3,Guide.manual_color_key("Methods", ["ES-FEM", "FEM", "analytical solution"],["green", "blue", "red"]),Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"),Guide.title("RESULT"))
    #,Coord.cartesian(xmin=0, xmax=0.18)
    draw(PDF("foo.pdf",20cm,8cm,dpi=1000),f)
end
=#
TEMPQ = zeros(110,1)
TEMPQ[2,1] = -200
for j = 12:10:92
    TEMPQ[j,1] = -400
end

TEMPMatrix = zeros(110,110)
TEMPMatrix[1:100,1:100] = KE[1:100,1:100]
TEMPMatrix[101:110,101:110] = Matrix{Float16}(I,10,10)

df = DataFrame(TEMPMatrix)
CSV.write("F://MasterTask//2_ES-FEM(2D)//KE.csv",df)

d=TEMPMatrix\TEMPQ

XX = collect(0:0.04:0.4)
Y1 = d[106:-10:6]

plot(x=XX,y=Y1,Geom.point,Geom.smooth(method=:loess,smoothing=0.5),
Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"))
#set_default_plot_size(20cm,8cm)
