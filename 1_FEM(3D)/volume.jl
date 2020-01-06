using LinearAlgebra
function Volume(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)

    a = [1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]
    v = 1/6*det(a)

    if v < 0
        printstyled("[STATUS] --> 体积正负转换\n", color=:green)
    elseif v == 0
        printstyled("[ERROR] --> 体积计算错误\n", color=:red)
    end

    u = abs(v)
    return u
end
