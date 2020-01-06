include("AssembleQuad.jl")
include("AssembleTriangle.jl")
include("BInner.jl")
include("BOuter.jl")
include("Stiffness.jl")

using LinearAlgebra
using CSV
using DataFrames

E = 1
NU = 0.25
t = 1

edge1 = BT(0,0,1,0,0,1)
ke1 = Stiffness(edge1,1/6)
edge2 = BT(1,0,1,1,0,1)
ke2 = Stiffness(edge2,1/6)
edge3 = BT(1,1,0,1,1,0)
ke3 = Stiffness(edge3,1/6)
edge4 = BT(0,1,0,0,1,0)
ke4 = Stiffness(edge4,1/6)
edge5 = BQ(0,1,0,0,1,0,1,1)
ke5 = Stiffness(edge5,1/3)

KE = zeros(8,8)
AssembleTriangle(KE,ke1,1,2,4)
AssembleTriangle(KE,ke2,2,3,4)
AssembleTriangle(KE,ke3,3,4,2)
AssembleTriangle(KE,ke4,4,1,2)
AssembleQuad(KE,ke5,4,1,2,3)

F = zeros(8,1)
F[6] = -10

K = zeros(8,8)
K[3:6,3:6] = KE[3:6,3:6]
K[1:2,1:2] = Matrix{Float16}(I,2,2)
K[7:8,7:8] = Matrix{Float16}(I,2,2)
df = DataFrame(K)

d = K \ F
s = (d[4]+d[6])/2
print("位移值为：")
print(s)