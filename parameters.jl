# mesh and solving parameters
begin
    h  = 1.0 # km
    Δt = .001 # s
    NX = 321
    NZ = 321
    NT = 3000
    grid = FDM_Grid(h, Δt, NZ, NX, NT)
end


# three layered model parameters
begin
    H1 = (1*NZ)÷8
    H2 = (1*NZ)÷8
    V1 = 3. # km/s
    V2 = 4. # km/s
    V3 = 5. # km/s

    # reflectors position
    z1 = H1
    z2 = z1+H2

    # actually defining velocity field
    v = Array{Float64}(undef, (NZ,NX))
    v[   1:z1,  1:end] .= V1
    v[z1+1:z2,  1:end] .= V2
    v[z2+1:end, 1:end] .= V3
end
