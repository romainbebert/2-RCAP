function load2RCAP(fname)
    f = open(fname)

    m, lim = parse.(Int, split(fname[15:end],"_") )

    C1 = zeros(Int, m, m)
    C2 = zeros(Int, m, m)
    weights = zeros(Int, m, m)

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

    #=
    println(m, " ", lim)
    println(C1)
    println(C2)
    println(weights)
    =#
    return m,lim,C1,C2,weights
end