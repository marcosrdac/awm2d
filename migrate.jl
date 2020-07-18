using InteractiveUtils
include("./reff.jl")


# signal parameters
begin
    P_file = "/mnt/hdd/home/tmp/awp_data/P.bin"
    reversed_P_file = "/mnt/hdd/home/tmp/awp_data/reversed_P.bin"
    migrated_file = "/mnt/hdd/home/tmp/awp_data/migrated.bin"
end

@time image_condition(P_file, reversed_P_file, migrated_file)
# 102 s
# 149 s
