using Base.Threads
using PyPlot
using InteractiveUtils
include("./reff.jl")


# parameters
NX = 9
NZ = 9
NT = 1

dt = .1
h = 1

# three layered model (temporary)
H1 = (2*NZ)÷6
H2 = (3*NZ)÷6
H3 = (4*NZ)÷6
V1 = 3
V2 = 3
V3 = 3

# reflectors position
z1 = H1
z2 = z1+H2

# actually defining velocity field
v0 = Array{Real}(undef, (NZ,NX))
v0[1:z1, 1:end] .= V1
v0[z1+1:z2, 1:end] .= V2
v0[z2+1:end, 1:end] .= V3

# starting from pressure field P0
P0 = zero(v0)
#P0[1,1] = 1     # endon test
#P0[1, NX÷2] = 1  # split test


v = pad_extremes(v0, ABS)

P = zeros((size(P0,1)+2*(ABS+1),
           size(P0,2)+2*(ABS+1),
           3))
P[1+ABS+1:end-ABS-1, 1+ABS+1:end-ABS-1, 1] .= P0

σ = .5
ricker = map(t->2/(√(3*σ)π^(1/4))*(1-(t/σ)^2)*exp(-t^2/(2σ^2)), collect(-20:20)/8)

function propagate(P, v, h²∇², nz, nx, nt, h, dt, ABS)
    dtoh2 = dt/h^2
    for timeiter in eachindex(1:nt)
        new_t = 1 + (timeiter)%3    # 2
        cur_t = 1 + (timeiter+2)%3  # 1
        old_t = 1 + (timeiter+1)%3  # 3

        #println(timeiter, " ", size(ricker, 1))
        #timeiter < size(ricker,1) ? P[100, 100, cur_t] = ricker[timeiter] : nothing
        #println(timeiter, " ", size(ricker, 1), " ", ricker[timeiter])
        #timeiter < size(ricker,1) ? P[2+ABS, 2+ABS+nx/2, cur_t] = ricker[timeiter] : nothing

        #println(new_t, " ", cur_t, " ", old_t)
        @threads for spciter in CartesianIndices((2:nz+2*ABS+1, 2:nx+2*ABS+1))
            begin
                zP, xP = spciter[1], spciter[2]
                zv, xv = zP-1, xP-1
                @views P[zP, xP, new_t] = new_p(P[:,:,cur_t], P[:,:,old_t],
                                                 v, dtoh2, h²∇², zP, xP, zv, xv)
            end
        end
    end
end

@time propagate(P, v, h²∇², NZ, NX, NT, h, dt, ABS)

#imshow(P[:,:,1])
#plt.show()
