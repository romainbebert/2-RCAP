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
    nbInd::Int64
	people::Array{BitArray} #Population of solutions
	ranking::Array{Int} #Rank for each solution
	mchance::AbstractFloat #Mutation chance
	gen_mean::AbstractFloat #Fitness mean for this generation
	gen_best::Int #Best fitnesses of the solutions (to replace with the pareto front ?)
end

struct Problem_Variables
	n::Int
	c1::Array{Int}
	c2::Array{Int}
	weights::Array{Int}
	limit::Int
end


#Number of individuals, problem specific variables
function firstGen(nbInd::Int, vars::Problem_Variables)
	
	startingPop::Array{BitArray} = []

    for i = 1:nbInd 
    	push!(startingPop, bitrand(vars.n)) 
    end

    println(startingPop)
end


function updateRanks(gen::Generation, vars::Problem_Variables)

end