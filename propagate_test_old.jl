using Base.Threads
using PaddedViews
using OffsetArrays
using InteractiveUtils
include("./reff.jl")


# parameters
NX = 5
NZ = 5
NT = 100

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
P0[1,1] = 1     # endon test
#P0[1, NX÷2] = 1  # split test


v = OffsetArray(pad_extremes(v0, ABS),
                1-ABS:size(v0, 1)+ABS,
                1-ABS:size(v0, 2)+ABS)

P = zeros((size(P0,1)+2*ABS, size(P0,2)+2*ABS, 3))
P[1+ABS:end-ABS, 1+ABS:end-ABS, 1] .= P0




function propagate(P, v, h²∇², nz, nx, nt, h, dt, ABS)
    dtoh2 = dt/h^2
    padded_P = PaddedView(0, P, (size(P,1)+2, size(P,2)+2, 3))
    _P = OffsetArray(padded_P, 1-ABS-1:size(padded_P,1)-ABS-1, 1-ABS-1:size(padded_P,2)-ABS-1, 1:3)
    _v = OffsetArray(v, 1-ABS:size(v,1)-ABS, 1-ABS:size(v,2)-ABS)

    # display(_P)# -2:8
    # display(_v)# -1:7

    for timeiter in eachindex(1:nt)
        new_t = 1 + (timeiter)%3    # 2
        cur_t = 1 + (timeiter+2)%3  # 1
        old_t = 1 + (timeiter+1)%3  # 3
        #println(new_t, " ", cur_t, " ", old_t)
        @threads for spciter in CartesianIndices((2-ABS:nz+ABS, 1-ABS:nx+ABS))
#        #@threads for spciter in CartesianIndices(P[2:end-1, 2:end-1, 1])
#        #@threads for spciter in CartesianIndices((2:nz+2*ABS-2, 2:nx+2*ABS-2))
            begin
                z, x = spciter[1], spciter[2]
                @views P[z, x, new_t] = new_p(_P[:,:,cur_t], _P[:,:,old_t],
                                         _v, dtoh2, h²∇², z, x)
            end
        end
    end
end

propagate(P, v, h²∇², NZ, NX, NT, h, dt, ABS)


#function prop(
#              P,
#              v,
#              h²∇²,
#              nt::Integer=NT,
#              nx::Integer=NX,
#              nz::Integer=NZ,
#              ABS::Integer=ABS,
#             )
#
#    dtoh2 = dt/h^2
#    for timeiter in eachindex(1:nt)
#        new_t = 1 + (timeiter)%3    # 2
#        cur_t = 1 + (timeiter+2)%3  # 1
#        old_t = 1 + (timeiter+1)%3  # 3
#        #println(new_t, " ", cur_t, " ", old_t)
#        @threads for spciter in CartesianIndices((1-ABS:nz+ABS, 1-ABS:nx+ABS))
#        #@threads for spciter in CartesianIndices(P[2:end-1, 2:end-1, 1])
#        #@threads for spciter in CartesianIndices((2:nz+2*ABS-2, 2:nx+2*ABS-2))
#            begin
#                z, x = spciter[1], spciter[2]
#                @views P[z, x, new_t] = new_p(P[:,:,cur_t], P[:,:,old_t],
#                                         v, dtoh2, h²∇², z, x)
#                println(z," ",x)
#            end
#        end
#    end
#end
#
##@time @code_warntype prop(_P, h²∇²)
#@time prop(_P,_v, h²∇²)
##::AbstractArray{Real, 2}=h²∇²
#display(_P)
##imshow(_P[:,:,1])
##plt.show()
