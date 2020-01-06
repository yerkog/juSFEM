function vector3(x1,y1,x2,y2,x3,y3)
    A = y2-y1
    B = x1-x2
    v1 = [A/sqrt(A^2+B^2) B/sqrt(A^2+B^2)]
    D = y3-y2
    E = x2-x3
    v2 = [D/sqrt(D^2+E^2) E/sqrt(D^2+E^2)]
    F = y1-y3
    G = x3-x1
    v3 = [F/sqrt(F^2+G^2) G/sqrt(F^2+G^2)]
    return [v1 v2 v3]
end
