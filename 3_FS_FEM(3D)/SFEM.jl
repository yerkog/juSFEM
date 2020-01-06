module SFEM

    using LinearAlgebra
    using DelimitedFiles
    using Pardiso
    using SparseArrays
    using WriteVTK
    #using Gadfly
    #using Cairo
    #using Fontconfig
    #using DataFrames
    #using CSV
    #using PlotlyJS
    #using ORCA

    include("E://Project-SFEM//3_FS_FEM(3D)//area.jl")
    include("E://Project-SFEM//3_FS_FEM(3D)//vectorout.jl")
    include("E://Project-SFEM//3_FS_FEM(3D)//vectorin.jl")
    include("E://Project-SFEM//3_FS_FEM(3D)//volume.jl")
    include("E://Project-SFEM//3_FS_FEM(3D)//assemblehex.jl")
    include("E://Project-SFEM//3_FS_FEM(3D)//assembletet.jl")


    export IOsfem
    export Ksfem
    export MKL
    #export Draw_mesh
    #export Draw_disp
    export WriteDISP
    #export Compare

#======================================================================
函数一：文件读取函数
======================================================================#
    function IOsfem(A, B, C, D, E)

        Cen = readdlm(A)
        Cen = Array{Float64,2}(Cen)

        Ele = readdlm(B)
        Ele = Array{Int64,2}(Ele)

        Face_in = readdlm(C)
        Face_in = Array{Int64,2}(Face_in)

        Face_out = readdlm(D)
        Face_out = Array{Int64,2}(Face_out)

        Node = readdlm(E)
        Node = Array{Float64,2}(Node)

        return Node, Ele, Cen, Face_out, Face_in
    end
#======================================================================
函数二：单元循环计算
======================================================================#
    function Ksfem(Node, Ele, Cen, Face_out, Face_in, EE, NU, FF)
        #==============================================================
                              @外力荷载条件（添加外荷载）
        ==============================================================#
        printstyled("[STATUS] --> 开始添加荷载条件\n", color=:green)
        Q = zeros(Float64,size(Node,1)*3,1)
        max = maximum(Node,dims = 1)
        min = minimum(Node,dims = 1)
        force_point = Vector{Int64}()

        for i = 1:size(Node,1)
            if Node[i,3] == max[3] && Node[i,1] != min[1]
                force_point = append!(force_point,i)
            end
        end

        for i = force_point
            Q[i*3] = FF/size(force_point,1)
        end
        printstyled("[STATUS] --> 荷载条件添加完毕\n", color=:green)
        #==============================================================
                                 @外部面
        ==============================================================#
        IK = Vector{Int64}(); JK = Vector{Int64}(); VK = Vector{Float64}()

        for i = 1:size(Face_out,1)
            #println("===========================================================================================")
            #printstyled("[STATUS] --> 开始计算第 $i 个外部面单元\n", color=:green)
            v = Volume(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                       Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                       Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                       Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
            #println("[INFO] --> 第 $i 个四面体的体积为 $v ")
            #@show v
            A1 = area(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                      Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                      Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3])
            #@show A1
                      #println("[INFO] --> 第1个面的面积为 $A1 ")
            A2 = area(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                      Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                      Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
            #@show A2
                      #println("[INFO] --> 第2个面的面积为 $A2 ")
            A3 = area(Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                      Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                      Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
            #@show A3
                      #println("[INFO] --> 第3个面的面积为 $A3 ")
            A4 = area(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                      Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                      Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
            #@show A4
                      #println("[INFO] --> 第4个面的面积为 $A4 ")
            n = Vectorout(Node[Face_out[i,1],1],Node[Face_out[i,1],2],Node[Face_out[i,1],3],
                          Node[Face_out[i,2],1],Node[Face_out[i,2],2],Node[Face_out[i,2],3],
                          Node[Face_out[i,3],1],Node[Face_out[i,3],2],Node[Face_out[i,3],3],
                          Cen[Face_out[i,4],1],Cen[Face_out[i,4],2],Cen[Face_out[i,4],3])
            #@show n
                          #printstyled("[STATUS] --> 开始计算应变矩阵\n", color=:green)
            bAx = (1/v)*(n[1,1]*(1/3)*A1+n[2,1]*(5/12)*A2+n[3,1]*(1/12)*A3+n[4,1]*(5/12)*A4)
            #@show bAx
            bAy = (1/v)*(n[1,2]*(1/3)*A1+n[2,2]*(5/12)*A2+n[3,2]*(1/12)*A3+n[4,2]*(5/12)*A4)
            #@show bAy
            bAz = (1/v)*(n[1,3]*(1/3)*A1+n[2,3]*(5/12)*A2+n[3,3]*(1/12)*A3+n[4,3]*(5/12)*A4)
            #@show bAz
            bA = [bAx 0 0; 0 bAy 0; 0 0 bAz; bAy bAx 0; 0 bAz bAy; bAz 0 bAx]
            bBx = (1/v)*(n[1,1]*(1/3)*A1+n[2,1]*(5/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(1/12)*A4)
            bBy = (1/v)*(n[1,2]*(1/3)*A1+n[2,2]*(5/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(1/12)*A4)
            bBz = (1/v)*(n[1,3]*(1/3)*A1+n[2,3]*(5/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(1/12)*A4)
            bB = [bBx 0 0; 0 bBy 0; 0 0 bBz; bBy bBx 0; 0 bBz bBy; bBz 0 bBx]
            bCx = (1/v)*(n[1,1]*(1/3)*A1+n[2,1]*(1/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(5/12)*A4)
            bCy = (1/v)*(n[1,2]*(1/3)*A1+n[2,2]*(1/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(5/12)*A4)
            bCz = (1/v)*(n[1,3]*(1/3)*A1+n[2,3]*(1/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(5/12)*A4)
            bC = [bCx 0 0; 0 bCy 0; 0 0 bCz; bCy bCx 0; 0 bCz bCy; bCz 0 bCx]
            bDx = (1/v)*(n[2,1]*(1/12)*A2+n[3,1]*(1/12)*A3+n[4,1]*(1/12)*A4)
            bDy = (1/v)*(n[2,2]*(1/12)*A2+n[3,2]*(1/12)*A3+n[4,2]*(1/12)*A4)
            bDz = (1/v)*(n[2,3]*(1/12)*A2+n[3,3]*(1/12)*A3+n[4,3]*(1/12)*A4)
            bD = [bDx 0 0; 0 bDy 0; 0 0 bDz; bDy bDx 0; 0 bDz bDy; bDz 0 bDx]
            B = hcat(bA,bB,bC,bD)
            #@show B
            #break
            #println("[INFO] --> 第 $i 个四面体的应变矩阵为 $B ")
            D = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0; NU 1-NU NU 0 0 0; NU NU 1-NU 0 0 0; 0 0 0 (1-2*NU)/2 0 0; 0 0 0 0 (1-2*NU)/2 0; 0 0 0 0 0 (1-2*NU)/2]
            ke = B'*D*B*v
            #@show ke
            #break
            #println("[INFO] --> 第 $i 个四面体的刚度矩阵为 $ke ")
            #printstyled("[STATUS] --> 开始组装刚度矩阵\n", color=:green)
            IK, JK, VK = Assembletet(IK,JK,VK,ke,Face_out[i,1],Face_out[i,2],Face_out[i,3],Face_out[i,5])
            #printstyled("[STATUS] --> 第 $i 个刚度矩阵组装完毕\n", color=:green)
        end
        printstyled("[STATUS] --> 外部面刚度矩阵组装完毕\n", color=:green)
        #==============================================================
                                 @内部面
        ==============================================================#
        for i = 1:size(Face_in,1)
            #println("===========================================================================================")
            #printstyled("[STATUS] --> 开始计算第 $i 个内部面单元\n", color=:green)
            v1 = Volume(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                        Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                        Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                        Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])
            v2 = Volume(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                        Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                        Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                        Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])
            v = v1 + v2
            #@show v
            #println("[INFO] --> 第 $i 个六面体的体积为 $v ")
            A1 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                      Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                      Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])
            #println("[INFO] --> 第1个面的面积为 $A1 ")
            #@show A1
            A2 = area(Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                      Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                      Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])
            #println("[INFO] --> 第2个面的面积为 $A2 ")
            #@show A2
            A3 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                      Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                      Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3])
            #println("[INFO] --> 第3个面的面积为 $A3 ")
            #@show A3
            A4 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                      Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                      Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])
            #println("[INFO] --> 第4个面的面积为 $A3 ")
            #@show A4
            A5 = area(Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                      Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                      Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])
            #println("[INFO] --> 第5个面的面积为 $A3 ")
            #@show A5
            A6 = area(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                      Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                      Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])
            #println("[INFO] --> 第6个面的面积为 $A4 ")
            #@show A6
            n = Vectorin(Node[Face_in[i,1],1],Node[Face_in[i,1],2],Node[Face_in[i,1],3],
                         Node[Face_in[i,2],1],Node[Face_in[i,2],2],Node[Face_in[i,2],3],
                         Node[Face_in[i,3],1],Node[Face_in[i,3],2],Node[Face_in[i,3],3],
                         Cen[Face_in[i,4],1],Cen[Face_in[i,4],2],Cen[Face_in[i,4],3],
                         Cen[Face_in[i,5],1],Cen[Face_in[i,5],2],Cen[Face_in[i,5],3])
            #printstyled("[STATUS] --> 开始计算应变矩阵\n", color=:green)
            #@show n
            bCx = (1/v)*(n[1,1]*(5/12)*A1+n[2,1]*(1/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(5/12)*A4+n[5,1]*(1/12)*A5+n[6,1]*(5/12)*A6)
            #@show bCx
            bCy = (1/v)*(n[1,2]*(5/12)*A1+n[2,2]*(1/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(5/12)*A4+n[5,2]*(1/12)*A5+n[6,2]*(5/12)*A6)
            #@show bCy
            bCz = (1/v)*(n[1,3]*(5/12)*A1+n[2,3]*(1/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(5/12)*A4+n[5,3]*(1/12)*A5+n[6,3]*(5/12)*A6)
            #@show bCz
            bC = [bCx 0 0; 0 bCy 0; 0 0 bCz; bCy bCx 0; 0 bCz bCy; bCz 0 bCx]
            bDx = (1/v)*(n[1,1]*(5/12)*A1+n[2,1]*(5/12)*A2+n[3,1]*(1/12)*A3+n[4,1]*(5/12)*A4+n[5,1]*(5/12)*A5+n[6,1]*(1/12)*A6)
            bDy = (1/v)*(n[1,2]*(5/12)*A1+n[2,2]*(5/12)*A2+n[3,2]*(1/12)*A3+n[4,2]*(5/12)*A4+n[5,2]*(5/12)*A5+n[6,2]*(1/12)*A6)
            bDz = (1/v)*(n[1,3]*(5/12)*A1+n[2,3]*(5/12)*A2+n[3,3]*(1/12)*A3+n[4,3]*(5/12)*A4+n[5,3]*(5/12)*A5+n[6,3]*(1/12)*A6)
            bD = [bDx 0 0; 0 bDy 0; 0 0 bDz; bDy bDx 0; 0 bDz bDy; bDz 0 bDx]
            bBx = (1/v)*(n[1,1]*(1/12)*A1+n[2,1]*(5/12)*A2+n[3,1]*(5/12)*A3+n[4,1]*(1/12)*A4+n[5,1]*(5/12)*A5+n[6,1]*(5/12)*A6)
            bBy = (1/v)*(n[1,2]*(1/12)*A1+n[2,2]*(5/12)*A2+n[3,2]*(5/12)*A3+n[4,2]*(1/12)*A4+n[5,2]*(5/12)*A5+n[6,2]*(5/12)*A6)
            bBz = (1/v)*(n[1,3]*(1/12)*A1+n[2,3]*(5/12)*A2+n[3,3]*(5/12)*A3+n[4,3]*(1/12)*A4+n[5,3]*(5/12)*A5+n[6,3]*(5/12)*A6)
            bB = [bBx 0 0; 0 bBy 0; 0 0 bBz; bBy bBx 0; 0 bBz bBy; bBz 0 bBx]
            bAx = (1/v)*(n[1,1]*(1/12)*A1+n[2,1]*(1/12)*A2+n[3,1]*(1/12)*A3)
            bAy = (1/v)*(n[1,2]*(1/12)*A1+n[2,2]*(1/12)*A2+n[3,2]*(1/12)*A3)
            bAz = (1/v)*(n[1,3]*(1/12)*A1+n[2,3]*(1/12)*A2+n[3,3]*(1/12)*A3)
            bA = [bAx 0 0; 0 bAy 0; 0 0 bAz; bAy bAx 0; 0 bAz bAy; bAz 0 bAx]
            bEx = (1/v)*(n[4,1]*(1/12)*A4+n[5,1]*(1/12)*A5+n[6,1]*(1/12)*A6)
            bEy = (1/v)*(n[4,2]*(1/12)*A4+n[5,2]*(1/12)*A5+n[6,2]*(1/12)*A6)
            bEz = (1/v)*(n[4,3]*(1/12)*A4+n[5,3]*(1/12)*A5+n[6,3]*(1/12)*A6)
            bE = [bEx 0 0; 0 bEy 0; 0 0 bEz; bEy bEx 0; 0 bEz bEy; bEz 0 bEx]
            B = hcat(bC,bD,bB,bA,bE)
            #@show B
            #break
            #println("[INFO] --> 第 $i 个六面体的应变矩阵为 $B ")
            D = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0; NU 1-NU NU 0 0 0; NU NU 1-NU 0 0 0; 0 0 0 (1-2*NU)/2 0 0; 0 0 0 0 (1-2*NU)/2 0; 0 0 0 0 0 (1-2*NU)/2]
            ke = B'*D*B*v
            #@show ke
            #break
            #println("[INFO] --> 第 $i 个六面体的刚度矩阵为 $ke ")
            #printstyled("[STATUS] --> 开始组装刚度矩阵\n", color=:green)
            IK, JK, VK = Assemblehex(IK,JK,VK,ke,Face_in[i,1],Face_in[i,2],Face_in[i,3],Face_in[i,6],Face_in[i,7])
            #printstyled("[STATUS] --> 第 $i 个刚度矩阵组装完毕\n", color=:green)
        end
        printstyled("[STATUS] --> 内部面刚度矩阵组装完毕\n", color=:green)
        #@show IK
        #@show JK
        #@show VK
        KE = sparse(IK, JK, VK)
        #df = DataFrame(KE)
        #CSV.write("F://MasterTask//3_FS-FEM(3D)//KE.csv",df)
        #==============================================================
                                @边界条件修改（添加固定端）
        ==============================================================#
        printstyled("[STATUS] --> 确定边界索引\n", color=:green)
        min = minimum(Node, dims = 1)
        node_dex = Vector{Int64}()
        @inbounds for i = 1:size(Node, 1)
            if Node[i,1] == min[1]
                append!(node_dex, i)
            end
        end
        printstyled("[STATUS] --> 生成固定点序列\n", color=:green)

        printstyled("[STATUS] --> 正在写入边界修改后的刚度矩阵\n", color=:green)
        @inbounds for i = 1:size(node_dex, 1)
            KE[node_dex[i]*3, node_dex[i]*3] = KE[node_dex[i]*3, node_dex[i]*3]*1E20
            KE[node_dex[i]*3-1, node_dex[i]*3-1] = KE[node_dex[i]*3-1, node_dex[i]*3-1]*1E20
            KE[node_dex[i]*3-2, node_dex[i]*3-2] = KE[node_dex[i]*3-2, node_dex[i]*3-2]*1E20
        end

        return KE, Q
    end
        #==============================================================
                                 @解方程
        ==============================================================#
    function MKL(KE,Q)
        printstyled("[STATUS] --> 开始解方程\n", color=:green)
        ps = MKLPardisoSolver()
        set_nprocs!(ps,24)
        result = similar(Q)
        result = solve(ps, KE, Q)
        println(get_nprocs(ps))

        Disp = zeros(Float64,Int(size(Q,1)/3),3)
        for i = 1:Int(size(Q,1)/3)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        return Disp
    end
#======================================================================
函数三：绘制悬臂梁长度方向位移示意图
======================================================================#
    function Draw_mesh(Node,Disp,B)
        show_point = Vector{Int64}()
        for i = 1:size(Node,1)
            if Node[i,2] == 0 && Node[i,3] == 0
                show_point = append!(show_point,i)
            end
        end
        println("[INFO] --> 绘制的点序号为: $show_point ")
        #==============================================================
                    @根据需要绘制的点的x坐标，由小到大重新排列点序号
        ==============================================================#
        show_point_xvalue = Vector{Float64}()
        for i = show_point
            show_point_xvalue = append!(show_point_xvalue,Node[i,1])
        end
        sort_point = sortperm(show_point_xvalue)
        ax = Vector{Float64}()
        for i = sort_point
            ax = append!(ax,Node[show_point[i],1])
        end
        ay = Vector{Float64}()
        for i = sort_point
            ay = append!(ay,Disp[show_point[i],3])
        end

        face = readdlm(B)
        face_num = Int(face[1,1])
        F1 = Vector{Int64}();F2 = Vector{Int64}();F3 = Vector{Int64}()
        for i = 1:face_num
                F1 = append!(F1,Int(face[i+1,2])+1)
                F2 = append!(F2,Int(face[i+1,3])+1)
                F3 = append!(F3,Int(face[i+1,4])+1)
        end
        Face = hcat(F1,F2,F3)

        data = GenericTrace[]
        for i = 1:size(Face,1)
            X = [Node[Face[i,1],1],Node[Face[i,2],1],Node[Face[i,3],1],Node[Face[i,1],1]]
            Y = [Node[Face[i,1],2],Node[Face[i,2],2],Node[Face[i,3],2],Node[Face[i,1],2]]
            Z = [Node[Face[i,1],3],Node[Face[i,2],3],Node[Face[i,3],3],Node[Face[i,1],3]]
            trace = scatter3d(;x=X,y=Y,z=Z,mode="lines",showlegend=false,line=attr(color="#1f77b4", width=1))
            push!(data,trace)
        end

        layout = Layout(title="Beam",
    					autosize="false",
    					width=1280,
    					height=720,
    					margin=attr(l=0,r=0,b=50,t=50),
                        scene=attr(aspectmode="manual",
                                   aspectratio=attr(x=1,y=0.17,z=0.17),
                                   xaxis=attr(title="X",range=[-1,6],showgrid=false,zeroline=false,showline=false,showticklabels=false),
                                   yaxis=attr(title="Y",range=[-1,1.5],showgrid=false,zeroline=false,showline=false,showticklabels=false),
                                   zaxis=attr(title="Z",range=[-1,1.5],showgrid=false,zeroline=false,showline=false,showticklabels=false))
                        )

        plot(data,layout)
        #savehtml(p,"1.html")
    end


    function Draw_disp(Node,Disp,F)
        face = readdlm(F)
        face_num = Int(face[1,1])
        F1 = Vector{Int64}();F2 = Vector{Int64}();F3 = Vector{Int64}()
        for i = 1:face_num
                F1 = append!(F1,Int(face[i+1,2])+1)
                F2 = append!(F2,Int(face[i+1,3])+1)
                F3 = append!(F3,Int(face[i+1,4])+1)
        end
        Face = hcat(F1,F2,F3)

        for i = 1:size(Node,1)
            Node[i,1] = Node[i,1] + Disp[i,1]
            Node[i,2] = Node[i,2] + Disp[i,2]
            Node[i,3] = Node[i,3] + Disp[i,3]
        end

        data = GenericTrace[]
        for i = 1:size(Face,1)
            X = [Node[Face[i,1],1],Node[Face[i,2],1],Node[Face[i,3],1],Node[Face[i,1],1]]
            Y = [Node[Face[i,1],2],Node[Face[i,2],2],Node[Face[i,3],2],Node[Face[i,1],2]]
            Z = [Node[Face[i,1],3],Node[Face[i,2],3],Node[Face[i,3],3],Node[Face[i,1],3]]
            trace = scatter3d(;x=X,y=Y,z=Z,mode="lines",showlegend=false,line=attr(color="#1f77b4", width=1.5))
            push!(data, trace)
        end

        layout = Layout(title="Beam",
                        autosize="false",
                        width=2000,
                        height=1080,
                        margin=attr(l=0,r=0,b=0,t=0),
                        scene=attr(aspectmode="manual",
                                   aspectratio=attr(x=1,y=1/6,z=1/6),
                                   xaxis=attr(title="X",range=[0,6.2],showgrid=true,zeroline=false,showline=true,showticklabels=true,
                                   yaxis=attr(title="Y",range=[0,1.0],showgrid=true,zeroline=false,showline=true,showticklabels=true),
                                   zaxis=attr(title="Z",range=[-0.03,1.0],showgrid=true,zeroline=false,showline=true,showticklabels=true))
                                   )
                        )
        plot(data,layout)
        #savefig(p,"C://Users//huoze//Desktop//1.png")
    end
#======================================================================
函数四：编写VTK格式文件
======================================================================#
    function WriteDISP(Node,Ele,result)

        Disp = zeros(Float64,size(Node,1),3)
        for i = 1:size(Node,1)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        pts = Node
        for i = 1:size(pts, 1)
            pts[i,1] = Node[i,1] + Disp[i,1]
            pts[i,2] = Node[i,2] + Disp[i,2]
            pts[i,3] = Node[i,3] + Disp[i,3]
        end

        pts = transpose(pts)
        pdata = Array{Float64}(undef,size(pts, 2))
        for i = 1:size(pts, 2)
            pdata[i] = Disp[i,3]
        end

        celltype = VTKCellTypes.VTK_TETRA
        cells = MeshCell[]

        for i = 1:size(Ele, 1)
            inds = Array{Int64}(undef, 4)
            inds[1] = Ele[i,1]
            inds[2] = Ele[i,2]
            inds[3] = Ele[i,3]
            inds[4] = Ele[i,4]

            c = MeshCell(celltype, inds)
            push!(cells, c)
        end

        vtk = vtk_grid("3D_disp", pts, cells,compress=3)
        vtk_point_data(vtk, pdata, "Z Disp (m)")
        outfiles =  vtk_save(vtk)
        return outfiles::Vector{String}, Disp
    end
#======================================================================
函数五：与解析解的对比图
======================================================================#
    function Compare(Node,result,EE,FF)
        Disp = zeros(Float64,size(Node,1),3)
        for i = 1:size(Node,1)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        l1 = Vector{Int64}(); l2 = Vector{Int64}(); l3 = Vector{Int64}(); l4 = Vector{Int64}()
        max = maximum(Node, dims = 1)
        min = minimum(Node, dims = 1)
        println(max)
        for i = 1:size(Node, 1)
            if Node[i,2] == min[2] && Node[i,3] == min[3]
                append!(l1, i)
            end
            if Node[i,2] == max[2] && Node[i,3] == min[3]
                append!(l2,i)
            end
            if Node[i,2] == max[2] && Node[i,3] == max[3]
                append!(l3,i)
            end
            if Node[i,2] == min[2] && Node[i,3] == max[3]
                append!(l4,i)
            end
        end
        println(l1)

        XX = Array{Float64}(undef,size(l1, 1))
        for i = 1:size(l1, 1)
            XX[i] = Node[i,1]
        end
        XX = sort(XX)

        v1 = Vector{Float64}(); v2 = Vector{Float64}(); v3 = Vector{Float64}(); v4 = Vector{Float64}()
        for i = size(l1,1)
            append!(v1,Disp[l1[i],3])
        end
        v1 = sort(v1,rev=true)

        for i = size(l2,1)
            append!(v2,Disp[l2[i],3])
        end
        v2 = sort(v2,rev=true)

        for i = size(l3,1)
            append!(v3,Disp[l3[i],3])
        end
        v3 = sort(v3,rev=true)

        for i = size(l4,1)
            append!(v4,Disp[l4[i],3])
        end
        v4 = sort(v4,rev=true)

        Y1 = (v1+v2+v3+v4)/4
        Y2 = map(i -> (FF*i^2)/(24*EE*((max[2]-min[2])*(max[3]-min[3])^3/12))*(i^2-4*(max[1]-min[1])*i+6*(max[1]-min[1])^2),XX)
        line1 = layer(x=XX,y=Y1,Geom.point,Geom.smooth,Theme(default_color="green"))
        line2 = layer(x=XX,y=Y2,Geom.point,Geom.smooth,Theme(default_color="blue"))
        plot(line1,line2)
        #set_default_plot_size(20cm,8cm)
        #plot(line1,line2,Guide.manual_color_key("Methods", ["FS-FEM", "analytical solution"],["blue", "red"]),Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"),Guide.title("RESULT"))
        #,Coord.cartesian(xmin=0, xmax=0.18)


    end

end
