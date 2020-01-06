function Faceout(a,b)
    c = zeros(Int64,2,size(a,2))
    c = SharedArray(c)
    @inbounds @sync @distributed for i = 1:size(a,2)
        for j = 1:size(b,2)
            is_1 = a[1,i] == b[1,j] || a[2,i] == b[1,j] || a[3,i] == b[1,j]
            is_2 = a[1,i] == b[2,j] || a[2,i] == b[2,j] || a[3,i] == b[2,j]
            is_3 = a[1,i] == b[3,j] || a[2,i] == b[3,j] || a[3,i] == b[3,j]
            is_4 = a[1,i] == b[4,j] || a[2,i] == b[4,j] || a[3,i] == b[4,j]

            is_123 = is_1 & is_2 & is_3
            is_124 = is_1 & is_2 & is_4
            is_134 = is_1 & is_3 & is_4
            is_234 = is_2 & is_3 & is_4

            if is_123
                c[1,i] = j
                c[2,i] = b[4,j]
            elseif is_124
                c[1,i] = j
                c[2,i] = b[3,j]
            elseif is_134
                c[1,i] = j
                c[2,i] = b[2,j]
            elseif is_234
                c[1,i] = j
                c[2,i] = b[1,j]
            end
        end
    end
    return Matrix(transpose(c))
end
