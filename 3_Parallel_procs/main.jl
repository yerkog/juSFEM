using .SFEM
using BenchmarkTools

A = "E://SMESH//1576.cen"
B = "E://SMESH//1576.ele"
C = "E://SMESH//1576.facein"
D = "E://SMESH//1576.faceout"
E = "E://SMESH//1576.node"

const EE, NU, FF = 2E8, 0.3, -2500.0
Node, Ele, Cen, Face_out, Face_in = IOsfem(A, B, C, D, E)
@time begin
    KE, Q = Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)
end

@time begin
    Disp = MKL(KE, Q)
end
