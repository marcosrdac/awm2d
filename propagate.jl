using Base.Threads
using OffsetArrays
using InteractiveUtils
#using PyPlot
#using Plots
include("./reff.jl")
#using OffsetArrays
#using LinearAlgebra
#using ImageFiltering


# parameters
NX = 5
NZ = 5
NT = 100

dt = .1
h = 1
dtoh2 = dt/h^2

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
v = Array{Real}(undef, (NZ,NX))
v[1:z1, 1:end] .= V1
v[z1+1:z2, 1:end] .= V2
v[z2+1:end, 1:end] .= V3

# starting from pressure field P0
P0 = zero(v)
#P0[1,1] = 1     # endon test
P0[1, NX÷2] = 1  # split test


# defining padded pressure fields for 3 times
P0_PADDING = ABS+∇²rad
P = Array{Real}(undef, (size(P0,1)+2*P0_PADDING, size(P0,2)+2*P0_PADDING, 3))
for i = 1:3
    P[:,:,i] .= zeros(Real, size(P,1), size(P,2))
end
P[1+P0_PADDING:end-P0_PADDING, 1+P0_PADDING:end-P0_PADDING, 1] .= P0


# defining padded velocity field _v (with indexes that make send)
_v = OffsetArray(pad_extremes(v, ABS),
                 1-ABS:size(v, 1)+ABS,
                 1-ABS:size(v, 2)+ABS)

# creating a padded version of pressure field (with indexes that make send)
_P = OffsetArray(P, 1-P0_PADDING:size(v,1)+P0_PADDING,
                    1-P0_PADDING:size(v,2)+P0_PADDING,
                    1:3)

function prop(
              #P::OffsetArray{Real,3,Array{Real,3}},
              _P,
              _v,
              h²∇²,
              nt::Integer=NT,
              nx::Integer=NX,
              nz::Integer=NZ,
              ABS::Integer=ABS,
             )
    for timeiter in eachindex(1:nt)
        new_t = 1 + (timeiter)%3    # 2
        cur_t = 1 + (timeiter+2)%3  # 1
        old_t = 1 + (timeiter+1)%3  # 3
        #println(new_t, " ", cur_t, " ", old_t)
        @threads for spciter in CartesianIndices((1-ABS:nz+ABS, 1-ABS:nx+ABS))
        #@threads for spciter in CartesianIndices(P[2:end-1, 2:end-1, 1])
        #@threads for spciter in CartesianIndices((2:nz+2*ABS-2, 2:nx+2*ABS-2))
            begin
                z, x = spciter[1], spciter[2]
                @views P[z, x, new_t] = new_p(P[:,:,cur_t], P[:,:,old_t],
                                         v, dtoh2, h²∇², z, x)
                println(z," ",x)
            end
        end
    end
end

#@time @code_warntype prop(_P, h²∇²)
@time prop(_P,_v, h²∇²)
#::AbstractArray{Real, 2}=h²∇²
display(_P)
#imshow(_P[:,:,1])
#plt.show()
