include("./reff.jl")

using .Propagate

P_file = "/mnt/hdd/home/tmp/awp_data/P.bin"
direct_seis_file = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
P_sub_direct_file = "/mnt/hdd/home/tmp/awp_data/P_sub_direct.bin"
reversed_P_file="/mnt/hdd/home/tmp/awp_data/reversed_P.bin"
migrated_file = "/mnt/hdd/home/tmp/awp_data/migrated.bin"

# mesh and solving parameters
begin
    h  = 1.0 # km
    Δt = .001 # s
    #NX = 321
    #NZ = 321
    NX = 100
    NZ = 100
    #NT = 1900
    NT = 100
    grid = FDM_Grid(h, Δt, NZ, NX, NT)
end


# signal parameters
begin
    ν = 6 # Hz
    array = "split"
    signature = rickerwave(ν, Δt)

    if array === "split"      position = CartesianIndex(1, NX÷2+1)
    elseif array === "endon"  position = CartesianIndex(1, 1)
    elseif array === "center" position = CartesianIndex(NZ÷2+1, NX÷2+1)
    end
    signal = Signal1D(signature, position)
end


# three layered model parameters
begin
    H1 = (1*NZ)÷3
    H2 = (1*NZ)÷3
    V1 = 3. # km/s
    V2 = 5. # km/s
    V3 = 9. # km/s

    # reflectors position
    z1 = H1
    z2 = z1+H2

    # actually defining velocity field
    v = Array{Float64}(undef, (NZ,NX))
    v[   1:z1,  1:end] .= V1
    v[z1+1:z2,  1:end] .= V2
    v[z2+1:end, 1:end] .= V3
end


# starting pressure field
begin
    P0 = zero(v)
end
