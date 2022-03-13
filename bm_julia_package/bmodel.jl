include("dependencies.jl")
include("utils.jl")
include("async_update.jl")

#=
Author : Kishore Hari
function name : bmodel
Inputs :
    topoFile : string; path to the topo file for the network
    nInit : integer; number of initial conditions to use for simulation
    nIter : integer; maximum number of time steps for simulation
    mode : string; ['Async', 'Sync']; mode of update
    stateRep : integer; [-1, 0]; whether to use -1 or 0 to represent low-expression
    rep : integer; replicate number. 
    csv : boolean; whether to write the output table into csv
    type : integer; [0, 1, 2]; the form of update rule to be used. 0 is equal weightage, 1 is activation dominant and 2 is inhibition dominant
Active outputs :
    state_df : DataFrame; 3 columns - [init, fin, flag]:
        init : string; Initial state
        fin : string; Final state
        flag : integer; 1 - fin is a steady state, 0 - not
    Nodes : string array; List of nodes with the order in which state is listed
Passive outputs :
    The function writes states_df to a csv file if csv = true in input
=#
function bmodel(topoFile::String; nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, type::Int=0, randSim::Bool = false,
    randVec::Array{Float64, 1}=[0.0])
    update_matrix,Nodes = topo2interaction(topoFile, type)
    if mode == "Async"
        if stateRep == -1
            if randSim
                state_df, frust_df = asyncRandUpdate(update_matrix, nInit, nIter, randVec)
            else 
                state_df, frust_df = asyncUpdate(update_matrix, nInit, nIter)
            end
        else
            state_df, frust_df = asyncUpdate2(update_matrix, nInit, nIter)
        end
    else
        print("Method under construction.")
    end
    # file_name = join([replace(topoFile, ".topo" => "_"), repl])
    # CSV.write(join(name,"_bmRes.csv"]), state_df)
    return state_df,Nodes, frust_df
end

#=
Author : Kishore Hari
function name : bmodel_reps
Inputs :
    topoFile : string; path to the topo file for the network
    nInit : integer; number of initial conditions to use for simulation
    nIter : integer; maximum number of time steps for simulation
    mode : string; ['Async', 'Sync']; mode of update
    stateRep : integer; whether to use -1 or 0 to represent low-expression
    reps : integer; number of replicates 
    csv : boolean; whether to write the output table of the function bmodel into csv
    types : integer array; subset of [0, 1, 2]; the forms of update rule to be used. 0 is equal weightage, 1 is activation dominant and 2 is inhibition dominant
    init : bool; Whether or not to include the initial conditions in the output. 
        If checked true, output will contain the frequency of all unique initial condition-steady state pair. For asynchronous boolean, this increases the table size by a lot
Active outputs :
    finFreqFinal_df : DataFrame; [states, Avg0, SD0, frust0,...]:
        states : string; Steady states
        Avg0 : float; Mean frequency (from 'reps' replicates) of the steady states obtained using rule 0
        SD0 : float; Standard deviation of the frequency
        frust0 : float; Frustration of the steady state
        Similarly we have Avg1, SD1, frust1, Avg2, SD2 and frust 2 depending upon the types argument for the function
    finFlagFreqFinal_df : DataFrame; [states, flag, Avg0,...]
        flag : integer; 1 - states is a steady state, 0 - not
    initFinFlagFreqFinal_df : DataFrame; [init, states, flag, Avg0,...]
        init : string; initial state
Passive outputs :
    The function writes finFreqFinal_df, finFlagFreqFinal_df, initFinFlagFreqFinal_df to files.
=#

function bmodel_reps(topoFile::String; nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3, csv::Bool=false, 
    types::Array{Int, 1} = [0],init::Bool=false, randSim::Bool=false, root::String="", 
    randVec::Array{Float64,1}=[0.0])
    update_matrix,Nodes = topo2interaction(topoFile)
    if length(Nodes)>60
        print("Network is too big")
        return 0
    end
    finFreqFinal_df_list_list = []
    finFlagFreqFinal_df_list_list = []
    initFinFlagFreqFinal_df_list_list = []
    frust_df_list = []


    for type in types

        finFreqFinal_df_list = []
        finFlagFreqFinal_df_list = []
        initFinFlagFreqFinal_df_list = []
        frust_df_list = []

        for rep in 1:reps
            states_df, Nodes, frust_df = bmodel(topoFile, nInit = nInit, 
                nIter = nIter, mode = mode, stateRep = stateRep, type = type, 
                randSim = randSim, randVec = randVec)
            # state_df = dropmissing(state_df, disallowmissing = true)
            push!(frust_df_list, frust_df)
            # Frequnecy table 
            finFreq_df = dfFreq(states_df, [:fin])
            #final states with flag
            finFlagFreq_df = dfFreq(states_df, [:fin, :flag])

            # all counts
            if init
                initFinFlagFreq_df = dfFreq(states_df, [:fin, :flag, :init])
                push!(initFinFlagFreqFinal_df_list, initFinFlagFreq_df)
            end
            push!(finFreqFinal_df_list, finFreq_df)
            push!(finFlagFreqFinal_df_list, finFlagFreq_df)
        end

        finFreqFinal_df = finFreqFinal_df_list[1]
        # println(typeof(finFreqFinal_df))
        finFlagFreqFinal_df = finFlagFreqFinal_df_list[1]
        if init
            initFinFlagFreqFinal_df = initFinFlagFreqFinal_df_list[1]
        end
        for i in 2:reps
            finFreqFinal_df = outerjoin(finFreqFinal_df, finFreqFinal_df_list[i], 
                on = [:states], makeunique = true)
            finFlagFreqFinal_df = outerjoin(finFlagFreqFinal_df, 
                finFlagFreqFinal_df_list[i], 
                on = [:states, :flag], makeunique=true)
            if init
                initFinFlagFreqFinal_df = outerjoin(initFinFlagFreqFinal_df, 
                    initFinFlagFreqFinal_df_list[i],
                    on = [:init, :states, :flag], makeunique = true)
            end
        end

        frust_df = DataFrame()
        for i in frust_df_list
            frust_df = vcat(frust_df, i)
        end
        frust_df = unique(frust_df, :fin)
        finFreqFinal_df = meanSD(finFreqFinal_df, "frequency")
        finFreqFinal_df = outerjoin(finFreqFinal_df, frust_df, 
            on = :states => :fin, makeunique =true)
        finFreqFinal_df = rename(finFreqFinal_df, 
            Dict(:Avg => Symbol(join(["Avg", type])), 
                :SD => Symbol(join(["SD", type])),
                :frust => Symbol(join(["frust", type]))))
        push!(finFreqFinal_df_list_list, finFreqFinal_df)


        finFlagFreqFinal_df = meanSD(finFlagFreqFinal_df, "frequency")
        finFlagFreqFinal_df = outerjoin(finFlagFreqFinal_df, frust_df, 
            on = :states => :fin, makeunique =true)
        finFlagFreqFinal_df = rename(finFlagFreqFinal_df, 
            Dict(:Avg => Symbol(join(["Avg", type])), 
                :SD => Symbol(join(["SD", type])),
                :frust => Symbol(join(["frust", type]))))
        push!(finFlagFreqFinal_df_list_list, finFlagFreqFinal_df)


        if init
            initFinFlagFreqFinal_df = meanSD(initFinFlagFreqFinal_df, "frequency")
            initFinFlagFreqFinal_df = outerjoin(initFinFlagFreqFinal_df, frust_df, 
                on = :states => :fin, makeunique =true)
            initFinFlagFreqFinal_df = rename(initFinFlagFreqFinal_df, 
                Dict(:Avg => Symbol(join(["Avg", type])), 
                :SD => Symbol(join(["SD", type])),
                :frust => Symbol(join(["frust", type]))))
            push!(initFinFlagFreqFinal_df_list_list, initFinFlagFreq_df)

        end
    end

    finFreqFinal_df = finFreqFinal_df_list_list[1]
        # println(typeof(finFreqFinal_df))
    finFlagFreqFinal_df = finFlagFreqFinal_df_list_list[1]
    if init
        initFinFlagFreqFinal_df = initFinFlagFreqFinal_df_list_list[1]
    end
    n = length(types)
    if n > 1
        for i in 2:n
            finFreqFinal_df = outerjoin(finFreqFinal_df, finFreqFinal_df_list_list[i], 
                on = [:states], makeunique = true)
            finFlagFreqFinal_df = outerjoin(finFlagFreqFinal_df, 
                finFlagFreqFinal_df_list_list[i], 
                on = [:states, :flag], makeunique=true)
            if init
                initFinFlagFreqFinal_df = outerjoin(initFinFlagFreqFinal_df, 
                    initFinFlagFreqFinal_df_list_list[i],
                    on = [:init, :states, :flag], makeunique = true)
            end
        end
    end
    

    rootName = replace(topoFile, ".topo" => "")
    if root !=""
        rootName = join([rootName, "_",root])
    end
    # println(rootName)
    if stateRep == 0
        rootName = join([rootName, "0"])
    end
    finFreqName = join([rootName, "_finFreq.csv"])
    finFlagFreqName = join([rootName, "_finFlagFreq.csv"])

    CSV.write(finFreqName, finFreqFinal_df)
    CSV.write(finFlagFreqName, finFlagFreqFinal_df)


    if init
        initFinFlagFreqName = join([rootName, "_initFinFlagFreq.csv"])
        CSV.write(initFinFlagFreqName, initFinFlagFreqFinal_df)
    end

    # write csv files
    if !randSim
        nodesName = join([replace(topoFile, ".topo" => ""), "_nodes.txt"])
        update_matrix,Nodes = topo2interaction(topoFile)
        io = open(nodesName, "w")
        for i in Nodes
            println(io, i)
        end
        close(io);
    end
    # return the dataframes
    if init
        return(finFreqFinal_df, finFlagFreqFinal_df, 
            initFinFlagFreqFinal_df)
    else
        return(finFreqFinal_df, finFlagFreqFinal_df)
    end
end

function edgeWeightPert(topoFile::String; nPerts::Int=10000, nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3, csv::Bool=false, 
    types::Array{Int, 1} = [0,1,2],init::Bool=false, randSim::Bool=true)
    updMat, nodes = topo2interaction(topoFile)
    nZ = length(findall(updMat.!=0))
    nRand = nZ*nPerts
    rands = reshape(rand(nRand), nPerts, nZ)
    randFold = replace(topoFile, ".topo" => "_rand")
    d1 = pwd()
    mkpath(randFold)
    cpPath = join([randFold, "/", topoFile])
    cp(topoFile, cpPath, force = true)
    cd(randFold)
    Threads.@threads for i in 1:nPerts
        # println(string(i))
        bmodel_reps(topoFile; nInit = nInit, nIter = nIter, mode = mode, stateRep = stateRep, randSim=true, root = string(i), 
        randVec = rands[i,:], types = types)
    end
    nodesName = join([replace(topoFile, ".topo" => ""), "_nodes.txt"])
        update_matrix,Nodes = topo2interaction(topoFile)
        io = open(nodesName, "w")
        for i in Nodes
            println(io, i)
        end
    close(io);
    cd(d1)

end



function stg(topoFile::String, mode::String="Async")
    print(topoFile)
    update_matrix,Nodes = topo2interaction(topoFile)
    if mode == "Async"
        stg = async_stg(update_matrix)
    else
        stg = sync_stg(update_matrix)
    end
    CSV.write(join([replace(topoFile, ".topo" => ""), "_stg.csv"]), stg)
    return stg,Nodes
end





