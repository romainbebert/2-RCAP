function load2RCAP(fname)
    f = open(fname)
    
    m, lim = parse.(Int, split(readline(f)) )

    C1 = zeros(Int, m, m)
    C2 = zeros(Int, m, m)
    weights = zeros(Int, m, m)
#=
    for i = 1:m
        j = 1 
        for valeur in split(readline(f))
            C1[i,j] = parse(Int, valeur)
            j+=1
        end
    end

    for i = 1:m
        j = 1 
        for valeur in split(readline(f))
            C2[i,j] = parse(Int, valeur)
            j+=1
        end
    end
    
    for i = 1:m
        j = 1 
        for valeur in split(readline(f))
            weights[i,j] = parse(Int, valeur)
            j+=1
        end
    end
    =#

    C1 = parse(Array{Int,m}, readline(f))
    C2 = parse(Array{Int,m}, readline(f))
    weights = parse(Array{Int,m}, readline(f))
    

    return m,lim,C1,C2,weights
end