include("Vector3.jl")
function BT(x1,y1,x2,y2,x3,y3)
    xi = (x1+x2+x3)/3
    yi = (y1+y2+y3)/3
    lp = [sqrt((x1-x2)^2+(y1-y2)^2) sqrt((x2-xi)^2+(y2-yi)^2) sqrt((xi-x1)^2+(yi-y1)^2)]
     n = vector3(x1,y1,x2,y2,xi,yi)
    Gp1 = [1/2 1/6 4/6]
    b11 = n[1]*Gp1[1]*lp[1]
    b12 = n[3]*Gp1[2]*lp[2]
    b13 = n[5]*Gp1[3]*lp[3]
    bx1 = (3/0.0004)*(b11+b12+b13)
    b14 = n[2]*Gp1[1]*lp[1]
    b15 = n[4]*Gp1[2]*lp[2]
    b16 = n[6]*Gp1[3]*lp[3]
    by1 = (3/0.0004)*(b14+b15+b16)
    B1 = [bx1 0; 0 by1; by1 bx1]
    Gp2 = [1/2 4/6 1/6]
    b21 = n[1]*Gp2[1]*lp[1]
    b22 = n[3]*Gp2[2]*lp[2]
    b23 = n[5]*Gp2[3]*lp[3]
    bx2 = (3/0.0004)*(b21+b22+b23)
    b24 = n[2]*Gp2[1]*lp[1]
    b25 = n[4]*Gp2[2]*lp[2]
    b26 = n[6]*Gp2[3]*lp[3]
    by2 = (3/0.0004)*(b24+b25+b26)
    B2 = [bx2 0; 0 by2; by2 bx2]
    Gp3 = [0 1/6 1/6]
    b31 = n[1]*Gp3[1]*lp[1]
    b32 = n[3]*Gp3[2]*lp[2]
    b33 = n[5]*Gp3[3]*lp[3]
    bx3 = (3/0.0004)*(b31+b32+b33)
    b34 = n[2]*Gp3[1]*lp[1]
    b35 = n[4]*Gp3[2]*lp[2]
    b36 = n[6]*Gp3[3]*lp[3]
    by3 = (3/0.0004)*(b34+b35+b36)
    B3 = [bx3 0; 0 by3; by3 bx3]
    B = [B1 B2 B3]
    return B
end
