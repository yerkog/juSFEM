using Distributed
addprocs(4)

@everywhere module SMESH

    using LinearAlgebra
    using SharedArrays
    using DelimitedFiles
    using Distributed

    include("E://Project-SFEM//pre-treated//faceout.jl")
    include("E://Project-SFEM//pre-treated//facein.jl")
    include("E://Project-SFEM//pre-treated//center.jl")

    export IOsfem

        function IOsfem(A,B,C)
            #==============================================================
                                     @点坐标文件读取
            ==============================================================#
            node = readdlm(A)
            node_num = Int(node[1,1])
            Node = zeros(Float64,node_num,3)
            @inbounds for i = 1:node_num
                Node[i,1] = node[i+1,2]
                Node[i,2] = node[i+1,3]
                Node[i,3] = node[i+1,4]
            end
            printstyled("[STATUS] --> 点坐标读取完成\n", color=:green)
            #==============================================================
                                     @四面体序列文件读取
            ==============================================================#
            ele = readdlm(C)
            ele_num = Int(ele[1,1])
            Ele = zeros(Int64,ele_num,4)
            @inbounds for i = 1:ele_num
                Ele[i,1] = ele[i+1,2]+1
                Ele[i,2] = ele[i+1,3]+1
                Ele[i,3] = ele[i+1,4]+1
                Ele[i,4] = ele[i+1,5]+1
            end
            printstyled("[STATUS] --> 四面体序列读取完成\n", color=:green)

            Cen = zeros(Float64,ele_num,3)
            @inbounds for i = 1:ele_num
                point = center(Node[Ele[i,1],1],Node[Ele[i,1],2],Node[Ele[i,1],3],
                               Node[Ele[i,2],1],Node[Ele[i,2],2],Node[Ele[i,2],3],
                               Node[Ele[i,3],1],Node[Ele[i,3],2],Node[Ele[i,3],3],
                               Node[Ele[i,4],1],Node[Ele[i,4],2],Node[Ele[i,4],3])
                Cen[i,1] = point[1]
                Cen[i,2] = point[2]
                Cen[i,3] = point[3]
            end
            printstyled("[STATUS] --> 体心序列计算完成\n", color=:green)
            #==============================================================
                                     @三角形面序列文件读取
            ==============================================================#
            face = readdlm(B)
            face_num = Int(face[1,1])
            printstyled("[STATUS] --> 面序列读取完成，正在划分边界\n", color=:green)

            Fx = Vector{Int64}(); Fy = Vector{Int64}(); Fz = Vector{Int64}()
            Fa = Vector{Int64}(); Fb = Vector{Int64}(); Fc = Vector{Int64}()

            @inbounds for i = 1:face_num
                if face[i+1,5] == 1
                    Fx = append!(Fx,Int(face[i+1,2])+1)
                    Fy = append!(Fy,Int(face[i+1,3])+1)
                    Fz = append!(Fz,Int(face[i+1,4])+1)
                elseif face[i+1,5] == 0
                    Fa = append!(Fa,Int(face[i+1,2])+1)
                    Fb = append!(Fb,Int(face[i+1,3])+1)
                    Fc = append!(Fc,Int(face[i+1,4])+1)
                else
                    printstyled("[ERROR] --> 面边界读取错误\n", color=:red)
                end
            end

            face_out = hcat(Fx,Fy,Fz)
            face_out_num = size(face_out,1)
            face_in = hcat(Fa,Fb,Fc)
            face_in_num = size(face_in,1)

            a1 = Matrix(transpose(face_out))
            b1 = Matrix(transpose(Ele))
            Face_out = hcat(face_out, Faceout(a1, b1))
            printstyled("[STATUS] --> 外部面索引完毕\n", color=:green)

            a2 = Matrix(transpose(face_in))
            b2 = Matrix(transpose(Ele))
            Face_in = hcat(face_in, Facein(a2, b2))
            printstyled("[STATUS] --> 内部面索引完毕\n", color=:green)
            println("==================================================")
            println("[INFO] --> 四面体共有 $ele_num 个")
            println("[INFO] --> 三角形面共有 $face_num 个")
            println("[INFO] --> 内部面共有 $face_in_num 个")
            println("[INFO] --> 外部面共有 $face_out_num 个")
            println("==================================================")
            open("E://SMESH//1576.node","w") do io
            writedlm(io, Node)
            end
            open("E://SMESH//1576.ele","w") do io
            writedlm(io, Ele)
            end
            open("E://SMESH//1576.cen","w") do io
            writedlm(io, Cen)
            end
            open("E://SMESH//1576.faceout","w") do io
            writedlm(io, Face_out)
            end
            open("E://SMESH//1576.facein","w") do io
            writedlm(io, Face_in)
            end
            printstyled("[STATUS] --> 索引完毕\n", color=:green)
        end

end
