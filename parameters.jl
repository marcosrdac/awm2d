P_file = "/mnt/hdd/home/tmp/awp_data/P.bin"
direct_seis_file = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
reversed_P_file="/mnt/hdd/home/tmp/awp_data/reversed_P.bin"
migrated_file = "/mnt/hdd/home/tmp/awp_data/migrated.bin"

# mesh definition
h  = 1.0 # km
Δt = .001 # s
NX = 321
NZ = 321
NT = 4500
#NX = 100
#NZ = 100
#NT = 100

# signal parameters
ν = 6 # Hz
array = "split"

# three layered model parameters
H1 = (1*NZ)÷3
H2 = (1*NZ)÷3
V1 = 3. # km/s
V2 = 5. # km/s
V3 = 9. # km/s
