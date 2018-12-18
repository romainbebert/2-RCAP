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
	mchance::AbstractFloat #Mutation chance (from 0 to 1)
end

struct Problem_Variables
	n::Int #nb of items
	limit::Int
	c1::Array{Int} #cost matrix for z1
	c2::Array{Int} #cost matrix for z2
	weights::Array{Int} 
end


#Number of individuals, problem specific variables
function firstGen(nbInd::Int, vars::Problem_Variables)	
	startingPop::Array{BitArray} = []

    while size(startingPop,1) < nbInd
    	tmp = bitrand(vars.n,vars.n)
    	if isValid(tmp,vars)
    		push!(startingPop,tmp)
    	end
    end

    return startingPop
end

#Generate offsprings from a generation (and the obligatory problem variables)
function offspring(gen::Generation, vars::Problem_Variables)

end
#nbInd = 5; gen = Generation(nbInd,firstGen(nbInd,vars),zeros(Int,nbInd),[Int[] for i=1:nbInd],0.1)
#Sorting using non-dominance as an order relation (fast sort algorithm) 
function NDsort(gen::Generation, vars::Problem_Variables)
	S = [Int[] for i=1:gen.nbInd] #Initializes to nbInd empty arrays
	n = zeros(Int,gen.nbInd)
	gen.fronts = [Int[] for i=1:gen.nbInd]
	for i = 1:gen.nbInd
		for j = 1:gen.nbInd
			obji = getObjectiveValues(gen.people[i],vars)
			objj = getObjectiveValues(gen.people[j],vars) 
			if obji[1] < objj[1] && obji[2] < objj[2]
				#println(i, " strictly dominates ", j)
				push!(S[i],j)
			elseif obji[1] > objj[1] && obji[2] > objj[2]
				n[i] += 1
			end
		end

		if n[i] == 0
			gen.ranking[i] = 1
			push!(gen.fronts[1],i)
			#println("Added ", i, " to the first front")
		end
	end

	i = 1 #initiating front counter

	while size(gen.fronts[i],1) != 0 && i < gen.nbInd
		Q = [] #Temporary array for the next front
		for p in gen.fronts[i]
			for q in S[p]
				n[q] -= 1
				if n[q] == 0
					gen.ranking[q] = i+1
					push!(Q,q)
					#println("Added ", q, " to the ", i+1 ," front")
				end
			end
		end
		i+=1
		gen.fronts[i] = Q
	end
	println("Fronts : ")
	println(gen.fronts)

end


function getObjectiveValues(x::BitArray{2}, vars::Problem_Variables)
	z1 = 0
	z2 = 0
	for i = 1:vars.n
		for j = 1:vars.n 
			z1 += vars.c1[i,j]*x[i,j]
			z2 += vars.c2[i,j]*x[i,j]
		end
	end

	#println("z1 = ",z1,"; z2 = ", z2)
	return z1,z2
end


function isValid(x::BitArray{2}, vars::Problem_Variables)
	i = 1
	j = 1
	total = 0

	while i <= vars.n && total <= vars.limit 
		total = 0
		while j <= vars.n && total <= vars.limit 
			total += vars.weights[i,j]*x[i,j]
			j += 1
		end
		i += 1
	end

	return total <= vars.limit
end

# Core function, args : number of individuals, problem variables 
# and optionnal stop condition (time in ms)
function nsga2(nbInd::Int, vars::Problem_Variables, stopTime = 10000)

end