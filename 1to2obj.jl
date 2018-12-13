using Random
using DelimitedFiles

function createInstance(fname)
    f = open(fname)
    m = parse.(Int,readline(f))

    C1 = zeros(Int, m, m)
    C2 = zeros(Int, m, m)
    weights = zeros(Int, m, m)
    limit = 0

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
            weights[i,j] = parse(Int, valeur)
            j+=1
        end
    end

    C2 = shuffle(C1)
    limit::Int = floor(sum(C1)/2)

    #=
    println(m, " ", limit)
    println(C1)
    println(C2)
    println(weights)
    =#
    writedlm(string("Instances2obj/",m,"_",limit),[C1;C2;weights], " ")
end