module FEM3D

    using LinearAlgebra
    using DelimitedFiles
    using WriteVTK
    #using DataFrames
    #using CSV
    #using PlotlyJS

    include("volume.jl")
    include("assembletet.jl")


    export IOFEM3D
    export KFEM3D
    export Draw_disp
    export WriteDISP
    export WriteSTRESS

    function IOFEM3D(A,B,C)
        #==============================================================
                                 @点坐标文件读取
        ==============================================================#
        node = readdlm(A)
        node_num = Int(node[1,1])
        Node = zeros(Float64,node_num,3)
        for i = 1:node_num
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
        for i = 1:ele_num
            Ele[i,1] = ele[i+1,2]+1
            Ele[i,2] = ele[i+1,3]+1
            Ele[i,3] = ele[i+1,4]+1
            Ele[i,4] = ele[i+1,5]+1
        end
        printstyled("[STATUS] --> 四面体序列读取完成\n", color=:green)

        return Node, Ele
    end

    function KFEM3D(EE,NU,FF,Node,Ele)
        KE = zeros(Float64,size(Node,1)*3,size(Node,1)*3)
        for i = 1:size(Ele,1)
            x1 = Node[Ele[i,1],1]
            x2 = Node[Ele[i,2],1]
            x3 = Node[Ele[i,3],1]
            x4 = Node[Ele[i,4],1]
            y1 = Node[Ele[i,1],2]
            y2 = Node[Ele[i,2],2]
            y3 = Node[Ele[i,3],2]
            y4 = Node[Ele[i,4],2]
            z1 = Node[Ele[i,1],3]
            z2 = Node[Ele[i,2],3]
            z3 = Node[Ele[i,3],3]
            z4 = Node[Ele[i,4],3]
            mbeta1 = [1 y2 z2; 1 y3 z3; 1 y4 z4]
            mbeta2 = [1 y1 z1; 1 y3 z3; 1 y4 z4]
            mbeta3 = [1 y1 z1; 1 y2 z2; 1 y4 z4]
            mbeta4 = [1 y1 z1; 1 y2 z2; 1 y3 z3]
            mgamma1 = [1 x2 z2; 1 x3 z3; 1 x4 z4]
            mgamma2 = [1 x1 z1; 1 x3 z3; 1 x4 z4]
            mgamma3 = [1 x1 z1; 1 x2 z2; 1 x4 z4]
            mgamma4 = [1 x1 z1; 1 x2 z2; 1 x3 z3]
            mdelta1 = [1 x2 y2; 1 x3 y3; 1 x4 y4]
            mdelta2 = [1 x1 y1; 1 x3 y3; 1 x4 y4]
            mdelta3 = [1 x1 y1; 1 x2 y2; 1 x4 y4]
            mdelta4 = [1 x1 y1; 1 x2 y2; 1 x3 y3]
            beta1 = -1*det(mbeta1)
            beta2 = det(mbeta2)
            beta3 = -1*det(mbeta3)
            beta4 = det(mbeta4)
            gamma1 = det(mgamma1)
            gamma2 = -1*det(mgamma2)
            gamma3 = det(mgamma3)
            gamma4 = -1*det(mgamma4)
            delta1 = -1*det(mdelta1)
            delta2 = det(mdelta2)
            delta3 = -1*det(mdelta3)
            delta4 = det(mdelta4)
            B1 = [beta1 0 0; 0 gamma1 0; 0 0 delta1; gamma1 beta1 0; 0 delta1 gamma1; delta1 0 beta1]
            B2 = [beta2 0 0; 0 gamma2 0; 0 0 delta2; gamma2 beta2 0; 0 delta2 gamma2; delta2 0 beta2]
            B3 = [beta3 0 0; 0 gamma3 0; 0 0 delta3; gamma3 beta3 0; 0 delta3 gamma3; delta3 0 beta3]
            B4 = [beta4 0 0; 0 gamma4 0; 0 0 delta4; gamma4 beta4 0; 0 delta4 gamma4; delta4 0 beta4]
            v = Volume(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
            B = [B1 B2 B3 B4]/(6*v)
            D = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0; NU 1-NU NU 0 0 0; NU NU 1-NU 0 0 0;
                 0 0 0 (1-2*NU)/2 0 0; 0 0 0 0 (1-2*NU)/2 0; 0 0 0 0 0 (1-2*NU)/2]
            ke = v*B'*D*B
            KE = Assembletet(KE,ke,Ele[i,1],Ele[i,2],Ele[i,3],Ele[i,4])
        end
        #df = DataFrame(KE)
        #CSV.write("F://MasterTask//1_FEM(3D)//KE.csv",df)
        #==============================================================
                            @边界条件修改（添加固定端）
        ==============================================================#
        printstyled("[STATUS] --> 开始添加边界条件\n", color=:green)
        min = minimum(Node,dims = 1)
        for i = 1:size(Node,1)
            if Node[i,1] == min[1]
                for j = 1:size(Node,1)*3
                    KE[i*3,j] = 0
                    KE[i*3-1,j] = 0
                    KE[i*3-2,j] = 0
                    KE[j,i*3] = 0
                    KE[j,i*3-1] = 0
                    KE[j,i*3-2] = 0
                    KE[i*3,i*3] = 1
                    KE[i*3-1,i*3-1] = 1
                    KE[i*3-2,i*3-2] = 1
                end
            end
        end
        #==============================================================
                             @外力荷载条件（添加固定端）
        ==============================================================#
        printstyled("[STATUS]--> 开始添加荷载条件\n", color=:green)
        Q = zeros(Float64,size(Node,1)*3,1)
        max = maximum(Node,dims = 1)
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
                                 @解方程
        ==============================================================#
        result = KE \ Q

        Disp = zeros(Float64,size(Node,1),3)
        for i = 1:size(Node,1)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        return Disp
    end
    #==============================================================
                             @后处理
    ==============================================================#
    function Draw_disp(Node,Disp,B)
        face = readdlm(B)
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
                        width=1280,
                        height=720,
                        margin=attr(l=0,r=0,b=0,t=0),
                        scene=attr(aspectmode="manual",
                                   aspectratio=attr(x=1,y=0.3,z=0.3),
                                   xaxis=attr(title="X",range=[-1,6],showgrid=false,zeroline=false,showline=false,showticklabels=false),
                                   yaxis=attr(title="Y",range=[-1,1.5],showgrid=false,zeroline=false,showline=false,showticklabels=false),
                                   zaxis=attr(title="Z",range=[-2.5,1.5],showgrid=false,zeroline=false,showline=false,showticklabels=false))
                        )

        plot(data,layout)
        #savefig(p,"F://MasterTask//3_FS-FEM(3D)//1.png")
    end
    #==============================================================
                             @VTK文件
    ==============================================================#
    function WriteDISP(Node,Ele,result)

        Disp = zeros(Float64,size(Node,1),3)
        for i = 1:size(Node,1)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        for i = 1:size(Node, 1)
            Node[i,1] = Node[i,1] + Disp[i,1]
            Node[i,2] = Node[i,2] + Disp[i,2]
            Node[i,3] = Node[i,3] + Disp[i,3]
        end

        pts = transpose(Node)
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
        vtk_point_data(vtk, pdata, "my_point_data")
        outfiles =  vtk_save(vtk)
        return outfiles::Vector{String}, Disp
    end

    function WriteSTRESS(EE,NU,Node,Ele,result)

        Disp = zeros(Float64,size(Node,1),3)
        for i = 1:size(Node,1)
            Disp[i,1] = result[i*3-2,1]
            Disp[i,2] = result[i*3-1,1]
            Disp[i,3] = result[i*3,1]
        end

        Misis = zeros(Float64, size(Ele, 1))
        #s = zeros(Float64, 6, size(Ele, 1))
        for i = 1:size(Ele,1)
            x1 = Node[Ele[i,1],1]
            x2 = Node[Ele[i,2],1]
            x3 = Node[Ele[i,3],1]
            x4 = Node[Ele[i,4],1]
            y1 = Node[Ele[i,1],2]
            y2 = Node[Ele[i,2],2]
            y3 = Node[Ele[i,3],2]
            y4 = Node[Ele[i,4],2]
            z1 = Node[Ele[i,1],3]
            z2 = Node[Ele[i,2],3]
            z3 = Node[Ele[i,3],3]
            z4 = Node[Ele[i,4],3]
            mbeta1 = [1 y2 z2; 1 y3 z3; 1 y4 z4]
            mbeta2 = [1 y1 z1; 1 y3 z3; 1 y4 z4]
            mbeta3 = [1 y1 z1; 1 y2 z2; 1 y4 z4]
            mbeta4 = [1 y1 z1; 1 y2 z2; 1 y3 z3]
            mgamma1 = [1 x2 z2; 1 x3 z3; 1 x4 z4]
            mgamma2 = [1 x1 z1; 1 x3 z3; 1 x4 z4]
            mgamma3 = [1 x1 z1; 1 x2 z2; 1 x4 z4]
            mgamma4 = [1 x1 z1; 1 x2 z2; 1 x3 z3]
            mdelta1 = [1 x2 y2; 1 x3 y3; 1 x4 y4]
            mdelta2 = [1 x1 y1; 1 x3 y3; 1 x4 y4]
            mdelta3 = [1 x1 y1; 1 x2 y2; 1 x4 y4]
            mdelta4 = [1 x1 y1; 1 x2 y2; 1 x3 y3]
            beta1 = -1*det(mbeta1)
            beta2 = det(mbeta2)
            beta3 = -1*det(mbeta3)
            beta4 = det(mbeta4)
            gamma1 = det(mgamma1)
            gamma2 = -1*det(mgamma2)
            gamma3 = det(mgamma3)
            gamma4 = -1*det(mgamma4)
            delta1 = -1*det(mdelta1)
            delta2 = det(mdelta2)
            delta3 = -1*det(mdelta3)
            delta4 = det(mdelta4)
            B1 = [beta1 0 0; 0 gamma1 0; 0 0 delta1; gamma1 beta1 0; 0 delta1 gamma1; delta1 0 beta1]
            B2 = [beta2 0 0; 0 gamma2 0; 0 0 delta2; gamma2 beta2 0; 0 delta2 gamma2; delta2 0 beta2]
            B3 = [beta3 0 0; 0 gamma3 0; 0 0 delta3; gamma3 beta3 0; 0 delta3 gamma3; delta3 0 beta3]
            B4 = [beta4 0 0; 0 gamma4 0; 0 0 delta4; gamma4 beta4 0; 0 delta4 gamma4; delta4 0 beta4]
            v = Volume(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
            B = [B1 B2 B3 B4]/(6*v)
            D = (EE/((1+NU)*(1-2*NU)))*[1-NU NU NU 0 0 0; NU 1-NU NU 0 0 0; NU NU 1-NU 0 0 0;
                 0 0 0 (1-2*NU)/2 0 0; 0 0 0 0 (1-2*NU)/2 0; 0 0 0 0 0 (1-2*NU)/2]
            q = [Disp[Ele[i,1],1],Disp[Ele[i,1],2],Disp[Ele[i,1],3],
                 Disp[Ele[i,2],1],Disp[Ele[i,2],2],Disp[Ele[i,2],3],
                 Disp[Ele[i,3],1],Disp[Ele[i,3],2],Disp[Ele[i,3],3],
                 Disp[Ele[i,4],1],Disp[Ele[i,4],2],Disp[Ele[i,4],3]]
            stress = D*B*q
            #=
            s[1,i] = stress[1,1]
            s[2,i] = stress[2,1]
            s[3,i] = stress[3,1]
            s[4,i] = stress[4,1]
            s[5,i] = stress[5,1]
            s[6,i] = stress[6,1]
            =#
            m = (1/sqrt(2))*sqrt((stress[1,1]-stress[2,1])^2+(stress[2,1]-stress[3,1])^2+(stress[3,1]-stress[1,1])^2+6*((stress[4,1])^2+(stress[5,1])^2+(stress[5,1])^2))
            #m = stress[1,1]
            #Misis[Ele[i,1]] = Misis[Ele[i,1]]+m/4
            #Misis[Ele[i,2]] = Misis[Ele[i,2]]+m/4
            #Misis[Ele[i,3]] = Misis[Ele[i,3]]+m/4
            #Misis[Ele[i,4]] = Misis[Ele[i,4]]+m/4
            Misis[i] = m
        end
        #println(Misis)


        for i = 1:size(Node, 1)
            Node[i,1] = Node[i,1] + Disp[i,1]
            Node[i,2] = Node[i,2] + Disp[i,2]
            Node[i,3] = Node[i,3] + Disp[i,3]
        end

        pts = transpose(Node)
        cdata = Misis
#=
        cdata = Array{Float64}(undef,size(Ele, 1))
        for i = 1:size(pts, 2)
            cdata[i] = s[1,i]
        end
=#
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

        vtk = vtk_grid("3D_stress", pts, cells,compress=3)
        vtk_cell_data(vtk, cdata, "Mises")
        outfiles =  vtk_save(vtk)
        return outfiles::Vector{String}
    end

end
