function Facein(a, b)
    c = zeros(Int64, 4, size(a, 2))
    c = SharedArray(c)
    @inbounds @sync @distributed for i=1:size(a, 2)
        tmp = 1
        for j=1:size(b, 2)
            tmp > 2 && continue
            is_1 = a[1, i] == b[1, j] || a[2, i] == b[1, j] || a[3, i] == b[1, j]
            is_2 = a[1, i] == b[2, j] || a[2, i] == b[2, j] || a[3, i] == b[2, j]
            is_3 = a[1, i] == b[3, j] || a[2, i] == b[3, j] || a[3, i] == b[3, j]
            is_4 = a[1, i] == b[4, j] || a[2, i] == b[4, j] || a[3, i] == b[4, j]
            # This line could potentially save time with larger matrix sizes.
            #!is_1 & !is_2 & !is_3 & !is_4 && continue
            is_123 = is_1 & is_2 & is_3
            is_124 = is_1 & is_2 & is_4
            is_134 = is_1 & is_3 & is_4
            is_234 = is_2 & is_3 & is_4
            if is_123
                c[tmp, i] = j
                c[tmp+2, i] = b[4, j]
                tmp += 1
            elseif is_124
                c[tmp, i] = j
                c[tmp+2, i] = b[3, j]
                tmp += 1
            elseif is_134
                c[tmp, i] = j
                c[tmp+2, i] = b[2, j]
                tmp += 1
            elseif is_234
                c[tmp, i] = j
                c[tmp+2, i] = b[1, j]
                tmp += 1
            end
        end
    end
    return Matrix(transpose(c))
end
