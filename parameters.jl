P_file = "/mnt/hdd/home/tmp/awp_data/P.bin"
shots_file = "/mnt/hdd/home/tmp/awp_data/shots.bin"
v_file = "/mnt/hdd/home/tmp/awp_data/v.bin"
source_signature_file = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
direct_seis_file = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
reversed_P_file="/mnt/hdd/home/tmp/awp_data/reversed_P.bin"
migrated_file = "/mnt/hdd/home/tmp/awp_data/migrated.bin"
multi_seis_file = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"

# mesh definition
Δz = 1.0 # km
Δx = 1.0 # km
Δt = 0.001 # s
NX = 321
NZ = 321
# NT = 4500
NT = 3000
# NX = 50
# NZ = 50
# NT = 500

# stencil order
stencil_order = 8

# signal parameters
array = "split"
nrec = 10
Δxrec = 10
