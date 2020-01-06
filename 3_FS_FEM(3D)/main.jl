using .SFEM
using DelimitedFiles
#using BenchmarkTools

A = "E://SMESH//1576.cen"
B = "E://SMESH//1576.ele"
C = "E://SMESH//1576.facein"
D = "E://SMESH//1576.faceout"
E = "E://SMESH//1576.node"
#F = "E://MESH//1576.face"

EE, NU, FF = 2E8, 0.3, -2500
Node, Ele, Cen, Face_out, Face_in = IOsfem(A, B, C, D, E)
KE, Q = Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)
#@time begin
result1 = MKL(KE, Q)
#end
#Draw_mesh(Node,Disp,B)
#Draw_disp(Node,Disp,F)
#Compare(Node,result,EE,FF)
#AA = WriteDISP(Node,Ele,result)
function dd(Node,result)
    value = Vector{Float64}()
    for i = 1:size(Node, 1)
        if Node[i,2] == 0 && Node[i,3] == 0
            append!(value,result[i,3])
        end
    end
    return value
end
