include("LinearTriangleAssemble.jl")
include("LinearTriangleElementStiffness.jl")
include("LinearTriangleElementStresses.jl")

using LinearAlgebra
using DataFrames
using CSV

E = 1
NU = 0.25
t = 1

s1 = LinearTriangleElementStiffness(E,NU,t,0,0,1,0,0,1)
s2 = LinearTriangleElementStiffness(E,NU,t,1,1,0,1,1,0)

KE = zeros(8,8)

KE=LinearTriangleAssemble(KE,s1,1,2,4)
KE=LinearTriangleAssemble(KE,s2,3,4,2)

TEMPQ = zeros(8,1)
TEMPQ[6,1] = -10

TEMPMatrix = zeros(8,8)
TEMPMatrix[3:6,3:6] = KE[3:6,3:6]
TEMPMatrix[1:2,1:2] = Matrix{Float16}(I,2,2)
TEMPMatrix[7:8,7:8] = Matrix{Float16}(I,2,2)

d=TEMPMatrix\TEMPQ
df = DataFrame(d)
s = (d[4]+d[6])/2