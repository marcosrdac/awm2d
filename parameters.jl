Pfile = "/mnt/hdd/home/tmp/awp_data/P.bin"
# Pfile = "/home/marcosrdac/P.bin" # 5 s faster w/ @time
shotsfile = "/mnt/hdd/home/tmp/awp_data/shots.bin"
vfile = "/mnt/hdd/home/tmp/awp_data/v.bin"
sourcesignaturefile = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
seisfile = "/mnt/hdd/home/tmp/awp_data/seis.bin"
seisgainfile = "/mnt/hdd/home/tmp/awp_data/seis_gain.bin"
directseisfile = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
reversedPfile="/mnt/hdd/home/tmp/awp_data/reversed_P.bin"
migratedfile = "/mnt/hdd/home/tmp/awp_data/migrated.bin"
migratedgainfile = "/mnt/hdd/home/tmp/awp_data/migrated_gain.bin"
multiseisfile = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"
multiseisfile = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"
multiseisgainfile = "/mnt/hdd/home/tmp/awp_data/multi_seis_gain.bin"

# mesh definition
Δz = 1.0 # km
Δx = 1.0 # km
# Δt = 0.001 # s
Δt = 0.001 # s
# nt = 4500
# nt = 3000
# nx = 50
# nz = 50
# nt = 800
# nt = 2600
nt = 300
# nt = 300
# nt = 100

# stencil order
stencil_order = 8

# signal parameters
array = "split"
nrec = 10
Δxrec = 10
