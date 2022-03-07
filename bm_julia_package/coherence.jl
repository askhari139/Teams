include("dependencies.jl")
include("utils.jl")
include("bmodel.jl")

function asyncInit(update_matrix2,
    nIter::Int, state::Array{Int,1})
    # state = rand(stateVec, n_nodes) #pick random state
    init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
    flag = 0
    n_nodes = length(state)
    for j in 1:nIter
        s1 = sign.(update_matrix2*state)
        u = rand(1:n_nodes, 1)
        if iszero(j%2) # check after every two steps,hopefully reduce the time
            if s1 == state
                flag = 1
                break
            end
        end
        state[u] = s1[u]
    end
    return state, flag
end

function coherenceIter(update_matrix,
    nIter::Int, state::Array{Int,1}, pertNodes::Union{Int, Array{Int,1}}, 
    nSim::Int)
    init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
    state[pertNodes] = -1*state[pertNodes]
    pert = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])

    sInit = copy(state)
    collection = reshape([], 0, 3)
    for i in 1:nSim
        sFin, flag = asyncInit(update_matrix, nIter, sInit)
        fin = join(["'", join(replace(x -> x == -1 ? 0 : x, sFin)), "'"])
        collection = vcat(collection, [init fin flag])
    end
    cols = [:init, :fin, :flag]
    collection = DataFrame(collection, cols)
    df = dfFreqGen(collection, cols)
    return df
end

function coherence(topoFile::String; nIter::Int=1000, simulate::Bool=false,
    bmodelArgs::Dict=Dict(), nPert::Int=1, randPert::Bool=false, nInit::Int=100, nSim::Int=10)
    if simulate
        x = bmodel_reps(topoFile, bmodelArgs...)
    end
    update_matrix,Nodes = topo2interaction(topoFile)
    n_nodes = length(Nodes)
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_matrix2 = sparse(update_matrix2')
    fileName = replace(topoFile, ".topo" => "_finFlagFreq.csv")
    states = CSV.read(fileName, DataFrame) 
    states = filter(:flag => x -> x== 1, states)
    if typeof(states.Avg0) == Vector{String}
        states = filter(:Avg0 => x -> x!="NA", states)
    else
        states = filter(:Avg0 => x -> !ismissing(x), states)
    end
    states = states[:, :states]
    nNodes = length(states[1]) - 2
    if nPert > 1
        randPert = true
    end
    choice = reshape(rand(1:nNodes, nPert*nInit), nInit, nPert)
    if !randPert
        choice = collect(1:nNodes)
        nInit = nNodes
    end
    collectionLarge = DataFrame(init=[], fin=[], flag=[], frequency=Float64[])
    for st in states
        st = st[2:(nNodes+1)]
        st = [x=='0' ? -1 : 1 for x in st]
        for i in 1:nInit
            pertNodes = choice[i,:]
            collectionLarge = vcat(collectionLarge, coherenceIter(update_matrix2, nIter, st, pertNodes, nSim))
        end
    end
    nm = "singleNode"
    if !randPert
        nm = "allNode"
    end
    if nPert > 1
        nm = "multiNode"
    end
    nm = join([nm, "_coherence.csv"])
    rootName = replace(topoFile, ".topo" => nm)
    CSV.write(rootName, collectionLarge)
end

