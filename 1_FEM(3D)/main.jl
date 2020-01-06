using .FEM3D

A = "G://MESH//test.node"
B = "G://MESH//test.face"
C = "G://MESH//test.ele"
EE, NU, FF = 2E8, 0.3, -2500

Node, Ele = IOFEM3D(A,B,C)
result = KFEM3D(EE,NU,FF,Node,Ele)
#Disp = WriteDISP(Node,Ele,result)
#Draw_mesh(Node,Disp,B)
#Draw_disp(Node,Disp,B)
#WriteSTRESS(EE,NU,Node,Ele,result)
function dd(Node,result)
    value = Vector{Float64}()
    for i = 1:size(Node, 1)
        if Node[i,2] == 0 && Node[i,3] == 0
            append!(value,result[i,3])
        end
    end
    return value
end
