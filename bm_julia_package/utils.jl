## Converts topo file to interaction matrix
function topo2interaction(topoFile::String, type::Int=0)
    df = DataFrame(CSV.File(topoFile))
    dropmissing!(df)
    Nodes = sort(unique(vcat(df.Source, df.Target)))
    n_nodes = length(Nodes)
    update_matrix = zeros(Int64, n_nodes, n_nodes)
    for i in 1:size(df, 1)
        if df[i, 3] == 2
            df[i,3] = -1
        end
        j = findfirst(x->x==df[i,1], Nodes)
        k = findfirst(x->x==df[i,2], Nodes)
        update_matrix[j,k] = Int64(df[i,3])
    end
    if type == 1
        replace!(x -> x == 1 ? 100 : x, update_matrix)
    end

    if type == 2
        replace!(x -> x == -1 ? -100 : x, update_matrix)
    end
    return update_matrix,Nodes
end

## create an exhaustive list of binary states of a given length
function listStates(nStates::Int)
    x = collect(Iterators.product([[-1,1] for i = 1:nStates]...))
    y = []
    for i in x
        j = collect(i)
        push!(y,j)
    end
    return y
end

## table function in R
function freqCalc(x::Array{String,1})
    y = unique(x)
    d = Dict([(i,count(k->k==i,x)) for i in y])
    df = DataFrame(d)
end

## converts binary state (-1s) to decimal
function stateToDec(state::Array{Int,1}, binVec::Array{Int,1})
    y = sum([(i != 1 && i != -1) for i in state])
    if y != 0
        print("Invalid state")
        return
    end
    state = replace(x -> x == -1 ? 0 : x, state)
    if size(state) != size(binVec)
        state = state'
    end
    s = sum(state.*binVec)
    return s
end

## calculates frustration of a state
function frustration(state::Array{Int,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        if s<0
            frustration = frustration + 1
        end
    end
    frustration = frustration/nEdges
    return frustration
end
function frustration(state::Array{Int,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        if s<0
            frustration = frustration + 1
        end
    end
    frustration = frustration/nEdges
    return frustration
end

## 
function frustration0(state::Array{Int,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        p = state[x] == state[y]
        if ((p == true) & (v == 1)) | ((p == false) & (v == -1))
            frustration = frustration + 1
        end
    end
    frustration = frustration/nEdges
    return frustration
end

## get frequency for select columns in a dataframe
function dfFreq(state_df::DataFrame, cols::Array{Symbol, 1})
    df = @> begin
        state_df
        DataFrames.groupby(cols)
        DataFrames.combine(nrow => :Count)
        @transform(frequency = :Count./sum(:Count))
        select(push!(cols, :frequency))
        rename(:fin => :states)
    end
    return df
end

function dfFreqGen(state_df::DataFrame, cols::Array{Symbol, 1})
    df = DataFrames.groupby(state_df, cols)
    df = DataFrames.combine(df, nrow => :Count)
    df.frequency = df.Count./sum(df.Count)
    df = select(df, push!(cols, :frequency))
    return df
end

## calculate average
function avg(x)
    m = sum(x)/length(x)
    return m
end

## calculate standard deviation
function SD(x)
    m = avg(x)
    v = sum((x.-m).^2)/length(x)
    s = sqrt(v)
    return s
end

## apply a function to each row of a dataframe
function rowWise(df::DataFrame, func::Function)
    v = []
    for i in eachrow(df)
        push!(v, func(Array(i)))
    end
    return v
end

## get mean and SD rowwise of columns containing a keyword in their name
function meanSD(df::DataFrame, keyword::String)
    cols = names(df)
    cols = cols[[occursin(keyword, i) for i in cols]]
    df_new = df[:, cols]
    m = rowWise(df_new, avg)
    s = rowWise(df_new, SD)
    d = @> begin
        df
        @transform(Avg = m, SD = s)
        select(Not(cols))
    end
    return d
end

## convert a state from -1 to 0
function zeroConv(state::Array{Int64,1})
    replace(x -> x == -1 ? 0 : x, sign.(state))
end

function zeroConv(state::Array{Float64,1})
    replace(x -> x == -1 ? 0 : x, Int.(sign.(state)))
end

  