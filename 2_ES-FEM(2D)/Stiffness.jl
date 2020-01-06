function Stiffness(B,A)
    D = (E/(1-NU*NU))*[1 NU 0; NU 1 0; 0 0 (1-NU)/2]
    ke = t*A*B'*D*B
    return ke
end
#test pass
