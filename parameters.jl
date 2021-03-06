# Modeling parameters

# === Mesh settings === #
Δz = 8.0  # m
Δx = 25.0  # m 
Δt = 0.002  # s
nt = 2400

# === Modeling kind === #
kind = "single shot modeling"
# kind = "multiple shot modeling"
# kind = "reverse time migration"

# === Simple direct modeling source position === #
(sx, sz) = (1, 1)

# === PDE Solver settings === #
# --- time --- #
# use REM or fall back to FDM
rem = true
# --- space --- #
# use PSM or fall back to FDM
pseudo = false


##### PATH SETTINGS #####

if kind == "single shot modeling"

    # === Direct modeling path settings === #
    # --- INPUT --- #
    sourcesignaturefile = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
    vfile = "/mnt/hdd/home/tmp/awp_data/v.bin"
    # --- OUTPUT --- #
    # 3D out
    Pfile = "/mnt/hdd/home/tmp/awp_data/P.bin"
    # seismogram out (alternative)
    seisfile = "/mnt/hdd/home/tmp/awp_data/seis.bin"

elseif kind == "multiple shot modeling"

    # Multiple-shot modeling parameters
    # --- INPUT --- #
    sourcesignaturefile = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
    # --- AUXILIAR --- #
    seisfile = "/mnt/hdd/home/tmp/awp_data/seis.bin"
    Pfile = "/mnt/hdd/home/tmp/awp_data/P.bin"
    # --- OUTPUT --- #
    multiseisfile = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"

elseif kind == "reverse time migration"

    # === RTM path settings === #
    # --- INPUT --- #
    # Multiple-shot 
    sourcesignaturefile = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
    multiseisfile = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"
    # --- AUXILIAR --- #
    Pfile = "/mnt/hdd/home/tmp/awp_data/P.bin"
    reversedPfile="/mnt/hdd/home/tmp/awp_data/reversed_P.bin"
    # --- OUTPUT --- #
    migratedfile = "/mnt/hdd/home/tmp/awp_data/migrated.bin"

else
    println("Invalid modeling kind: '", kind ,"'.")

end
