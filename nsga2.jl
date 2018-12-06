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

#------------------- Type(s) definitions ---------------------------

type Generation {
    nbInd::Int64
	people::Array{Array{Int64,1},1} #Population of solutions
	ranking::Array{Int64} #Rank for each solution
	mchance::Float64 #Mutation chance
	gen_mean::Float64 #Fitness mean for this generation
	gen_best::Int64 #Best fitnesses of the solutions
	fittest::Array{Int64,1} #Best solutions of the generation (probable duplicate with gen_best but it makes some things easier to code)
}


#Number of individuals, optional starting pop, percent of randomised individuals, problem specific variables
function firstGen(nbInd, startingPop = [], randpercent, costs, duedates) {
    if startingPop == []
        for i in 1:nbInd
            push!(startingPop, rand(bounds.inf, bounds.sup))
        end
    else 
        for i in 1:(nbInd-startingPop.size)
            push!(startingPop, rand(bounds.inf, bounds.sup))
        end
    end
}