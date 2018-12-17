#= 	Author : Romain Bernard
	Problem : 2-RCAP

	Genetic parameters :
		(crossover chance : xchance) -> pas applicable à cette implémentation
		mutation chance : mchance
		generation size : gen_size
		stop condition : stop (time ms)

	Steps :
		create first gen (at random with a portion from GRASP or a solver)
		apply crossover using elite parents
		apply mutation to the childrens randomly
		rince and repeat until stop condition and return best solution found

	Strategy :
		elite
		First gen created partially at random and partially with GRASP (or vOPT)

=#

#----------------------	  Packages	------------------------------	
using Random 

#------------------- Type(s) definitions ---------------------------

mutable struct Generation
    nbInd::Int
	people::Array{BitArray} #Population of solutions
	ranking::Array{Int} #Rank for each solution
	fronts::Array{Array{Int}} #List of pareto fronts (fronts[i] is rank i)
	mchance::AbstractFloat #Mutation chance
end

struct Problem_Variables
	n::Int #nb of items
	c1::Array{Int} #cost matrix for z1
	c2::Array{Int} #cost matrix for z2
	weights::Array{Int} 
	limit::Int
end


#Number of individuals, problem specific variables
function firstGen(nbInd::Int, vars::Problem_Variables)	
	startingPop::Array{BitArray} = []

    for i = 1:nbInd 
    	push!(startingPop, bitrand(vars.n,vars.n)) 
    end

end

#Sorting using non-dominance as an order relation (fast sort algorithm) 
function NDsort(gen::Generation, vars::Problem_Variables)
	S = [Int[] for i=1:gen.nbInd] #Initializes to nbInd empty arrays
	n = zeros(Int,gen.nbInd)
	for i = 1:gen.nbInd
		for j = 1:gen.nbInd
			obji = getObjectiveValues(gen.people[i])
			objj = getObjectiveValues(gen.people[j]) 
			if obji < objj
				push!(S[i],j)
			elseif obji > objj
				n[i] += 1
			end
		end

		if n[i] == 0
			ranking[i] = 1
			push!(gen.fronts[1],i)
		end
	end

	i = 1 #initiating front counter

	while size(gen.fronts[i],1) != 0
		Q = [] #Temporary array for the next front
		
		for p in gen.fronts[i]
			for q in S[p]
				n[q] -= 1
				if n[q] == 0
					ranking[q] = i+1
					push!(Q,q)
				end
			end
		end
		i+=1
		gen.fronts[i] = Q
	end

end


function getObjectiveValues(x::Array{BitArray}, vars::Problem_Variables)
	z1 = 0
	z2 = 0

	for i = 1:vars.n
		for j = 1:vars.n 
			z1 += vars.c1[i][j]*x[i][j]
			z2 += vars.c2[i][j]*x[i][j]
		end
	end

	return z1,z2
end


function isValid(x::Array{BitArray}, vars::Problem_Variables)
	i = 1
	j = 1

	while i <= vars.n && total <= vars.limit 
		total = 0
		while j <= vars.n && total <= vars.limit 
			total += w[i][j]*x[i][j]
			j += 1
		end
		i += 1
	end

	return total <= vars.limit
end


function nsga2(nbInd::Int, vars::Problem_Variables, stopTime = 10000)

end