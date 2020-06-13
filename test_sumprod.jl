include("./functions.jl")


function sumprod1(A, B)
    window = @view B[1:m,1:m]
    return(sum(A .* window))
end

function sumprod2(A, B)
    window = B[1:m,1:m]
    return(sum(A .* window))
end

function sumprod3(A, B)
    return(sum(A .* @view B[1:m,1:m]))
end

function sumprod4(A, B)
    return(sum(A .* view(B,1:m,1:m)))
end

function sumprod5(A, B)
    return(sum(A .* B))
end

function sumprod6(A, B)
    sum = 0
    for i in 1:m*m
        sum += A[i] .* B[i]
    end
    return(sum)
end

function sumprod7(A, B)
    sum = 0
    for i in eachindex(@view A[1:m,1:m])
        @inbounds sum += A[i] .* B[i]
    end
    return(sum)
end

function sumprod8(A, B)
    sum = 0
    for i in eachindex(A[1:m,1:m])
        @inbounds sum += A[i] .* B[i]
    end
    return(sum)
end

function sumprod9(A, B)
    sum = 0
    for i in eachindex(A)
        @inbounds sum += A[i] .* B[i]
    end
    return(sum)
end

function sumprod10(A, B)
    sum = 0
    for i in eachindex(A)
        sum += A[i] .* B[i]
    end
    return(sum)
end

function sumprod11(A, B)
    sum = 0
    for i in eachindex(copy(B))
        sum += A[i] .* B[i]
    end
    return(sum)
end

function sumprod12(B)
    return(sum(arr .* @view B[1:m,1:m]))
end

function sumprod13(A,B)
    mapreduce((A,B)->A.*B,+, A, copy(B))
end
    
n=1000
m=3
a = reshape(collect(1:n^2), (n,n))
arr = a[2:m+1,2:m+1]
arr2 = view(a,2:m+1, 2:m+1)
arr3 = view(a,3:m+2, 3:m+2)

# while compiling
println("before compilation")
println(@time sumprod1(arr,arr2))
println(@time sumprod2(arr,arr2))
println(@time sumprod3(arr,arr2))
println(@time sumprod4(arr,arr2))
println(@time sumprod5(arr,arr2))
println(@time sumprod6(arr,arr2))
println(@time sumprod7(arr,arr2))
println(@time sumprod8(arr,arr2))
println(@time sumprod9(arr,arr2))
println(@time sumprod10(arr,arr2))
println(@time sumprod11(arr,arr2))
println(@time sumprod12(arr2))
println(@time sumprod13(arr,arr2))

println()
println("after compilation")
println(@time sumprod1(arr,arr3))
println(@time sumprod2(arr,arr3))
println(@time sumprod3(arr,arr3))
println(@time sumprod4(arr,arr3))
println(@time sumprod5(arr,arr3))
println(@time sumprod6(arr,arr3))
println(@time sumprod7(arr,arr3))
println(@time sumprod8(arr,arr3))
println(@time sumprod9(arr,arr3))
println(@time sumprod10(arr,arr3))
println(@time sumprod11(arr,arr3))
println(@time sumprod12(arr3))
println(@time sumprod13(arr,arr3))


