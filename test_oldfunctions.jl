LAPLACIAN_KERNEL = [1   4  1
                    4 -20  4
                    1   4  1] / 6.

LAPLACIAN_KERNEL_RADIUS = size(LAPLACIAN_KERNEL, 1)÷2


function h2laplacian(array, z, x)
    window = view(array, i-LAPLACIAN_KERNEL_RADIUS:i+LAPLACIAN_KERNEL_RADIUS,
                         j-LAPLACIAN_KERNEL_RADIUS:j+LAPLACIAN_KERNEL_RADIUS)
    return(sum(LAPLACIAN_KERNEL .* window))
end


function new_p(old_P, cur_P, v, dtoh2, z, x)
    return(2*cur_P[z, x] - old_P[z, x] +
           dtoh2*v[z, x]^2 * h2laplacian(cur_P, z, x))
end


function new_p!(old_P, cur_P, new_P!, v, dtoh2, z, x)
    new_P![z,x] = 2*cur_P[z, x] - old_P[z, x] +
                  dtoh2*v[z, x]^2 * h2laplacian(cur_P, z, x)
end


function pad_extremes(v, padding)
    _v = Array{Real}(undef, (size(v, 1)+2*padding,
                             size(v, 2)+2*padding))
    # center
    _v[1+padding:size(v,1)+padding, 1+padding:size(v,2)+padding] .= v
    # extremes
    _v[1:padding, 1+padding:end-padding] .= repeat(v[1:1,1:end], padding,1)
    _v[1+padding:end-padding, 1:padding] .= repeat(v[1:end,1:1], 1, padding)
    _v[1+end-padding:end, 1+padding:end-padding] .= repeat(v[end:end,1:end], padding,1)
    _v[1+padding:end-padding, 1+end-padding:end] .= repeat(v[1:end,end:end], 1, padding)
    # edges
    _v[1:padding, 1:padding] .= v[1,1]
    _v[1+end-padding:end, 1:padding] .= v[end,1]
    _v[1+end-padding:end, 1+end-padding:end] .= v[end,end]
    _v[1:padding, 1+end-padding:end] .= v[1,end]
    return(_v)
end

# This function was made for odd sided kernels
function old_valid_convolve(array, kernel, convolved)
    padding = (size(kernel,1) ÷ 2, size(kernel,2) ÷ 2)
    _array = PaddedView(0, array, (1-padding[1]:size(array,1)+padding[1],
                                   1-padding[2]:size(array,2)+padding[2]))
    for j = 1:size(array, 2)
        @async for i = 1:size(array, 1)
            @async begin
                window = view(_array, i-padding[1]:i + padding[1],
                                      i-padding[2]:i + padding[2])
                convolved[i,j] = sum(kernel .* window)
            end
        end
    end
    display(convolved)
end

# This function was made for odd sided kernels
function valid_convolve!(array, kernel, convolved!)
    #@assert (size(kernel, 1)%2≢0) && (size(kernel, 2)%2≢0)
    padding = (size(kernel,1) ÷ 2, size(kernel,2) ÷ 2)
    _array = PaddedView(0, array, (1-padding[1]:size(array,1)+padding[1],
                                   1-padding[2]:size(array,2)+padding[2]))
    @sync for j = 1:size(array, 2)
        @async for i = 1:size(array, 1)
            @async begin
                window = view(_array, i-padding[1]:i + padding[1],
                                      j-padding[2]:j + padding[2])
                convolved![i,j] = sum(kernel .* window)
            end
        end
    end
end
