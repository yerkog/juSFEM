using LinearAlgebra
function Vectorin(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5)
    #分别为C、D、B、H、I五个点
    v = zeros(Float64,6,3)
    #第一个面(x1,y1,z1),(x2,y2,z2),(x4,y4,z4)
    CD = [x2-x1;y2-y1;z2-z1]
    CH = [x4-x1;y4-y1;z4-z1]
    na = LinearAlgebra.cross(CD,CH)
    BC = [x1-x3;y1-y3;z1-z3]
    j1 = dot(na,BC)
    if j1 > 0
        n1 = na
    #elseif j1 == 0
        #printstyled("[ERROR] --> 法向量计算错误：1st\n", color=:red)
    else
        n1 = -na
        #printstyled("[WARN] --> 法向量方向转换：1st\n", color=:yellow)
    end
    union1 = n1/sqrt(n1[1]^2+n1[2]^2+n1[3]^2)
    #println("[INFO] --> 第1个面单位法向量: $union1 ")
    v[1,1] = union1[1]
    v[1,2] = union1[2]
    v[1,3] = union1[3]

    #第二个面(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)
    DH = [x4-x2;y4-y2;z4-z2]
    DB = [x3-x2;y3-y2;z3-z2]
    nb = LinearAlgebra.cross(DH,DB)
    CD = [x2-x1;y2-y1;z2-z1]
    j2 = dot(nb,CD)
    if j2 > 0
        n2 = nb
    #elseif j2 == 0
        #printstyled("[ERROR] --> 法向量计算错误：2nd\n", color=:red)
    else
        n2 = -nb
        #printstyled("[WARN] --> 法向量方向转换：2nd\n", color=:yellow)
    end
    union2 = n2/sqrt(n2[1]^2+n2[2]^2+n2[3]^2)
    #println("[INFO] --> 第2个面单位法向量: $union2 ")
    v[2,1] = union2[1]
    v[2,2] = union2[2]
    v[2,3] = union2[3]

    #第三个面(x1,y1,z1),(x3,y3,z3),(x4,y4,z4)
    HB = [x3-x4;y3-y4;z3-z4]
    HC = [x1-x4;y1-y4;z1-z4]
    nc = LinearAlgebra.cross(HB,HC)
    DH = [x4-x2;y4-y2;z4-z2]
    j3 = dot(nc,DH)
    if j3 > 0
        n3 = nc
    #elseif j3 == 0
        #printstyled("[ERROR] --> 法向量计算错误：3rd\n", color=:red)
    else
        n3 = -nc
        printstyled("[WARN] --> 法向量方向转换：3rd\n", color=:yellow)
    end
    union3 = n3/sqrt(n3[1]^2+n3[2]^2+n3[3]^2)
    #println("[INFO] --> 第3个面单位法向量: $union3 ")
    v[3,1] = union3[1]
    v[3,2] = union3[2]
    v[3,3] = union3[3]

    #第四个面(x1,y1,z1),(x2,y2,z2),(x5,y5,z5)
    CD = [x2-x1;y2-y1;z2-z1]
    CI = [x5-x1;y5-y1;z5-z1]
    nd = LinearAlgebra.cross(CD,CI)
    BC = [x1-x3;y1-y3;z1-z3]
    j4 = dot(nd,BC)
    if j4 > 0
        n4 = nd
    #elseif j4 == 0
        #printstyled("[ERROR] --> 法向量计算错误：4th\n", color=:red)
    else
        n4 = -nd
        #printstyled("[WARN] --> 法向量方向转换：4th\n", color=:yellow)
    end
    union4 = n4/sqrt(n4[1]^2+n4[2]^2+n4[3]^2)
    #println("[INFO] --> 第4个面单位法向量: $union4 ")
    v[4,1] = union4[1]
    v[4,2] = union4[2]
    v[4,3] = union4[3]

    #第五个面(x2,y2,z2),(x3,y3,z3),(x5,y5,z5)
    BD = [x2-x3;y2-y3;z2-z3]
    BI = [x5-x3;y5-y3;z5-z3]
    ne = LinearAlgebra.cross(BD,BI)
    CB = [x3-x1;y3-y1;z3-z1]
    j5 = dot(ne,CB)
    if j5 > 0
        n5 = ne
    #elseif j5 == 0
        #printstyled("[ERROR] --> 法向量计算错误：5th\n", color=:red)
    else
        n5 = -ne
        #printstyled("[WARN] --> 法向量方向转换：5th\n", color=:yellow)
    end
    union5 = n5/sqrt(n5[1]^2+n5[2]^2+n5[3]^2)
    #println("[INFO] --> 第5个面单位法向量: $union5 ")
    v[5,1] = union5[1]
    v[5,2] = union5[2]
    v[5,3] = union5[3]

    #第六个面(x1,y1,z1),(x3,y3,z3),(x5,y5,z5)
    IB = [x3-x5;y3-y5;z3-z5]
    IC = [x1-x5;y1-y5;z1-z5]
    nf = LinearAlgebra.cross(IB,IC)
    DI = [x5-x2;y5-y2;z5-z2]
    j6 = dot(nf,DI)
    if j6 > 0
        n6 = nf
    #elseif j6 == 0
        #printstyled("[ERROR] --> 法向量计算错误：6th\n", color=:red)
    else
        n6 = -nf
        #printstyled("[WARN] --> 法向量方向转换：6th\n", color=:yellow)
    end
    union6 = n6/sqrt(n6[1]^2+n6[2]^2+n6[3]^2)
    #println("[INFO] --> 第6个面单位法向量: $union6 ")
    v[6,1] = union6[1]
    v[6,2] = union6[2]
    v[6,3] = union6[3]

    return v
end
