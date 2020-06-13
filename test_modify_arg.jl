include("./functions.jl")

arr = ones(3,3,2)
arr[:,:,1] = zeros(3,3)

function modify!(array, result!)
  result![1,1] = array[1,1]
end

modify!(arr[1:end,1:end,1], view(arr,:,:,2))
display(arr)
