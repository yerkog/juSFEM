function LinearTriangleAssemble(KE,ke,i,j,m)
    DOF = [2*i-1 2*i 2*j-1 2*j 2*m-1 2*m]
    for n1 = 1:6
        for n2  = 1:6
            KE[DOF[n1],DOF[n2]] = KE[DOF[n1],DOF[n2]] + ke[n1,n2]
        end
    end
#=    KE[2*i-1,2*i-1] = KE[2*i-1,2*i-1] + ke[1,1]
    KE[2*i-1,2*i] = KE[2*i-1,2*i] + ke[1,2]
    KE[2*i-1,2*j-1] = KE[2*i-1,2*j-1] + ke[1,3]
    KE[2*i-1,2*j] = KE[2*i-1,2*j] + ke[1,4]
    KE[2*i-1,2*m-1] = KE[2*i-1,2*m-1] + ke[1,5]
    KE[2*i-1,2*m] = KE[2*i-1,2*m] + ke[1,6]
    KE[2*i,2*i-1] = KE[2*i,2*i-1] + ke[2,1]
    KE[2*i,2*i] = KE[2*i,2*i] + ke[2,2]
    KE[2*i,2*j-1] = KE[2*i,2*j-1] + ke[2,3]
    KE[2*i,2*j] = KE[2*i,2*j] + ke[2,4]
    KE[2*i,2*m-1] = KE[2*i,2*m-1] + ke[2,5]
    KE[2*i,2*m] = KE[2*i,2*m] + ke[2,6]
    KE[2*j-1,2*i-1] = KE[2*j-1,2*i-1] + ke[3,1]
    KE[2*j-1,2*i] = KE[2*j-1,2*i] + ke[3,2]
    KE[2*j-1,2*j-1] = KE[2*j-1,2*j-1] + ke[3,3]
    KE[2*j-1,2*j] = KE[2*j-1,2*j] + ke[3,4]
    KE[2*j-1,2*m-1] = KE[2*j-1,2*m-1] + ke[3,5]
    KE[2*j-1,2*m] = KE[2*j-1,2*m] + ke[3,6]
    KE[2*j,2*i-1] = KE[2*j,2*i-1] + ke[4,1]
    KE[2*j,2*i] = KE[2*j,2*i] + ke[4,2]
    KE[2*j,2*j-1] = KE[2*j,2*j-1] + ke[4,3]
    KE[2*j,2*j] = KE[2*j,2*j] + ke[4,4]
    KE[2*j,2*m-1] = KE[2*j,2*m-1] + ke[4,5]
    KE[2*j,2*m] = KE[2*j,2*m] + ke[4,6]
    KE[2*m-1,2*i-1] = KE[2*m-1,2*i-1] + ke[5,1]
    KE[2*m-1,2*i] = KE[2*m-1,2*i] + ke[5,2]
    KE[2*m-1,2*j-1] = KE[2*m-1,2*j-1] + ke[5,3]
    KE[2*m-1,2*j] = KE[2*m-1,2*j] + ke[5,4]
    KE[2*m-1,2*m-1] = KE[2*m-1,2*m-1] + ke[5,5]
    KE[2*m-1,2*m] = KE[2*m-1,2*m] + ke[5,6]
    KE[2*m,2*i-1] = KE[2*m,2*i-1] + ke[6,1]
    KE[2*m,2*i] = KE[2*m,2*i] + ke[6,2]
    KE[2*m,2*j-1] = KE[2*m,2*j-1] + ke[6,3]
    KE[2*m,2*j] = KE[2*m,2*j] + ke[6,4]
    KE[2*m,2*m-1] = KE[2*m,2*m-1] + ke[6,5]
    KE[2*m,2*m] = KE[2*m,2*m] + ke[6,6]=#
    return KE
end
