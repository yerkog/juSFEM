using .SFEM
using BenchmarkTools

A = "G://SMESH//9.cen"
B = "G://SMESH//9.ele"
C = "G://SMESH//9.facein"
D = "G://SMESH//9.faceout"
E = "G://SMESH//9.node"

const EE, NU, FF = 3E7, 0.3, -1000.0
Node, Ele, Cen, Face_out, Face_in = IOsfem(A, B, C, D, E)
@time begin
    KE, Q = Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)
end
@time begin
    Disp = MKL(KE, Q)
end
