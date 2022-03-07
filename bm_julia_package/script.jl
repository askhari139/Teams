include("/mnt/d/Github/Projects/Ongoing/Relative_boolean/bm_julia_package/bmodel.jl")
using Base.Threads

fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end

println(Threads.nthreads())


# Threads.@threads for topoFile in topoFiles
# 	y1 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = -1)
# 	y2 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = 0)
# 	println(topoFile, " - ", y1, "and", y2, " seconds.")
# end
for topoFile in topoFiles
	y3 = @elapsed x = edgeWeightPert(topoFile; nPerts = 100, nInit = 10000, types = [0])
	println(topoFile, " - ", y3, " seconds.")
end