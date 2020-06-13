h²∇² = [1.0   4.0  1.0
        4.0 -20.0  4.0
        1.0   4.0  1.0] / 6.

∇²rad = size(h²∇², 1)÷2
ABS = 80
#ABS = 2


function pad_extremes(v, padding)
    _v = Array{Real}(undef, (size(v, 1)+2*padding,
                             size(v, 2)+2*padding))
    # center
    _v[1+padding:size(v,1)+padding, 1+padding:size(v,2)+padding] .= v
    # borders
    _v[1:padding, 1+padding:end-padding] .= repeat(v[1:1,1:end], padding,1)
    _v[1+padding:end-padding, 1:padding] .= repeat(v[1:end,1:1], 1, padding)
    _v[1+end-padding:end, 1+padding:end-padding] .= repeat(v[end:end,1:end], padding,1)
    _v[1+padding:end-padding, 1+end-padding:end] .= repeat(v[1:end,end:end], 1, padding)
    # diagonal edges
    _v[1:padding, 1:padding] .= v[1,1]
    _v[1+end-padding:end, 1:padding] .= v[end,1]
    _v[1+end-padding:end, 1+end-padding:end] .= v[end,end]
    _v[1:padding, 1+end-padding:end] .= v[1,end]
    return(_v)
end

function reduce_pointwise_product(A::AbstractArray, B::AbstractArray)
    prod_sum = 0
    for i in eachindex(A)
        prod_sum += A[i] .* B[i]
    end
    return(prod_sum)
end


function new_p(cur_P, old_P, v, dtoh², h²∇², zP, xP, zv, xv)
    win = @view cur_P[zP-1:zP+1, xP-1:xP+1]
    h²∇²cur_P = reduce_pointwise_product(h²∇², win)
    return(2*cur_P[zP, xP] - old_P[zP, xP] + v[zv, xv]^2*dtoh² * h²∇²cur_P)
end


"""
    ricker(t::Real, σ::Real=1)
Compute the normalized ricker function of t with standard deviation σ.
"""
function norm_ricker(t::Real, σ::Real=1)
    return(2/(√(3*σ)π^(1/4)) * (1-(t/σ)^2) * exp(-t^2/(2σ^2)))
end


"""
    rickerwave(ν::Real, Δt::Real)
Compute a ricker wave with main frequency ν, sampled in time intervals Δt.

For a stable signal assert ``ν << \\frac{1}{2 dt}``
"""
function rickerwave(ν::Real, Δt::Real)
    @assert ν < 0.2*(1.0/(2.0*Δt))
    function ricker(t::Real)
        return((1-2*t^2) * exp(-t^2))
    end
    nw = 2*Integer(floor(2.2/ν/Δt/2))
    nc = Integer(floor(nw/2))
    t = (π*ν*Δt) .* (nc-nw:nc)
    rickerwave = ricker.(t)
end

#plt.plot(rickerwave(60, .0001))
#plt.show()
