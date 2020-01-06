function LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym)
    A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2
    βi = yj-ym
    βj = ym-yi
    βm = yi-yj
    γi = xm-xj
    γj = xi-xm
    γm = xj-xi
    B = [βi 0 βj 0 βm 0; 0 γi 0 γj 0 γm; γi βi γj βj γm βm]/(2*A)
    D = (E/(1-NU*NU))*[1 NU 0; NU 1 0; 0 0 (1-NU)/2]
    ke = t*A*B'*D*B
    return ke
end
