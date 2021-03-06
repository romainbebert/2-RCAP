#= 	Author : Romain Bernard
	Problem : 2-RCAP

	Genetic parameters :
		(crossover chance : xchance) -> pas applicable à cette implémentation
		mutation chance : mchance
		generation size : nbInd
		stop condition : nbGen (number of generation, Int)

	Steps :
		create first gen (at random for now)
		apply crossover using elite parents
		apply mutation to the childrens randomly
		rince and repeat until stop condition and return best solution found

	Strategy :
		elite
		First gen created partially at random and partially with GRASP (or vOPT)

=#

#----------------------	  Packages	------------------------------	
using Random 
using PyPlot
using vOptSpecific
include("load2RCAP.jl")

#------------------- Type(s) definitions ---------------------------

mutable struct Generation
    nbInd::Int
	people::Array{Array{Int}} #Population of solutions
	ranking::Array{UInt16} #Rank for each solution (might be useless)
	fronts::Array{Array{Int}} #List of pareto fronts (fronts[i] is rank i)
	crowddist::Array{Float64} #List of crowding distances (crowddist[i] is rank i)
	mchance::Float64 #Mutation chance (from 0 to 1)
	#Ideal and Anti-Ideal coordinates
	best1::Int
	best2::Int
	worst1::Int
	worst2::Int
end

struct Problem_Variables
	n::Int #nb of items
	limit::Int
	c1::Array{Int} #cost matrix for z1
	c2::Array{Int} #cost matrix for z2
	weights::Array{Int} 
end

#-------------------------------------------------------------------

#--------------- Generation fabrication functions ------------------


#Number of individuals, problem specific variables
function firstGen(nbInd::Int, vars::Problem_Variables)	
	startingPop = []
	println("Creating first Gen...")
	while size(startingPop,1) < nbInd
    	push!(startingPop,shuffle(1:vars.n))
    end
    println("Done !")
    return startingPop
end


#-------------------------------------------------------------------


#------------------------ NSGA-II Operators ------------------------


#Sorting using non-dominance as an order relation (fast sort algorithm) 
function NDsort(gen::Generation, vars::Problem_Variables)
	S = [Int[] for i=1:size(gen.people,1)] #Initializes to nbInd empty arrays
	n = zeros(Int,size(gen.people,1))
	gen.fronts = [Int[] for i=1:size(gen.people,1)]
	gen.ranking = zeros(Int, size(gen.people,1))
	for i = 1:size(gen.people,1)
		obji = getObjectiveValues(gen.people[i],vars)
		istate = isValid(gen.people[i], vars)
		for j = 1:size(gen.people,1)
			objj = getObjectiveValues(gen.people[j],vars)
			jstate =  isValid(gen.people[j], vars)
			
			#Dominance test (validity is included in the test)
			if istate != jstate
				if istate > jstate
					push!(S[i],j)
				else 
					n[i] += 1
				end
			elseif obji[1] <= objj[1] && obji[2] <= objj[2]
				push!(S[i],j)
			elseif obji[1] >= objj[1] && obji[2] >= objj[2]
				n[i] += 1
			end
		end

		#Update best/worst values for objective 1
		if obji[1] < gen.best1
			gen.best1 = obji[1]
		elseif obji[1] > gen.worst1
			gen.worst1 = obji[1]
		end
		#Update best/worst values for objective 2
		if obji[2] < gen.best2
			gen.best2 = obji[2]
		elseif obji[2] > gen.worst2
			gen.worst2 = obji[2]
		end

		if n[i] == 0
			gen.ranking[i] = 1
			push!(gen.fronts[1],i)
			#println("Added ", i, " to the first front")
		end
	end

	i = 1 #initiating front counter

	while size(gen.fronts[i],1) != 0 && i < size(gen.people,1)
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

end

function crowdingDistance(gen::Generation, vars::Problem_Variables) 

	i = 1 #initiating front counter
	gen.crowddist = zeros(Float64,size(gen.people,1))

	while size(gen.fronts[i],1) != 0 && i < size(gen.people,1)

		f1 = sort(gen.fronts[i], by = x -> getObjectiveValues(gen.people[x],vars)[1], alg=Base.Sort.QuickSort) 
		f2 = sort(gen.fronts[i], by = x -> getObjectiveValues(gen.people[x],vars)[2], alg=Base.Sort.QuickSort) 

		gen.crowddist[f1[1]] = Base.Inf
		gen.crowddist[f1[size(f1,1)]] = Base.Inf
		gen.crowddist[f2[1]] = Base.Inf
		gen.crowddist[f2[size(f2,1)]] = Base.Inf

		for j = 2:(size(f1,1)-1)
			gen.crowddist[f1[j]] += (getObjectiveValues(gen.people[f1[j+1]], vars)[1] - 
									 getObjectiveValues(gen.people[f1[j-1]], vars)[1])/(gen.worst1-gen.best1)
			#println(gen.crowddist[f1[j]])
			gen.crowddist[f2[j]] += (getObjectiveValues(gen.people[f2[j+1]], vars)[2] - 
									 getObjectiveValues(gen.people[f2[j-1]], vars)[2])/(gen.worst2-gen.best2)
			#println(gen.crowddist[f2[j]])
		end

		i += 1
	end
end

#-------------------------------------------------------------------


#------------------------ Genetic Operators ------------------------

#Selection by binary tournament (returns indexes of the solutions)
function tournament(gen::Generation)
	p1,p2 = rand(1:size(gen.people,1)), rand(1:size(gen.people,1))
	while p1 == p2 
		p2 = rand(1:size(gen.people,1))
	end
	if gen.ranking[p1] < gen.ranking[p2]
		return p1 
	elseif gen.ranking[p1] > gen.ranking[p2]
		return p2 
	else 
		return gen.crowddist[p1] > gen.crowddist[p2] ? p1 : p2
	end
end


function two_point_crossover(p1::Array{Int}, p2::Array{Int}, vars::Problem_Variables)
	
	cut_a = cut_b = rand(2:vars.n-1)
	c1 = c2 = 1:vars.n

	while cut_a == cut_b 
		cut_b = rand(2:vars.n-1)
	end

	cut_a,cut_b = minmax(cut_a,cut_b)

	copyto!(c1, 1, p1, 1, cut_a-1)
	copyto!(c1, cut_a, p2, cut_a, cut_b-cut_a+1)
	copyto!(c1, cut_b+1, p1, cut_b+1, vars.n-cut_b)

	copyto!(c2, 1, p2, 1, cut_a-1)
	copyto!(c2, cut_a, p1, cut_a, cut_b-cut_a+1)
	copyto!(c2, cut_b+1, p2, cut_b+1, vars.n-cut_b)

	return c1,c2
end

#= Etapes Ordered Crossover (OX) :
	- Choisir deux points de coupe
	- Enfant 1 prend la partie coupée de parent1
	- Enfant 1 récupère les éléments de parent2 ne se trouvant pas dans subseq1
	- Idem pour Enfant 2 en inversant les parents
=#
function OX_crossover(p1::Array{Int},p2::Array{Int},vars::Problem_Variables)
	nbVar = size(p1,1)
	num1 = rand(1:nbVar)
	num2 = rand(1:nbVar)
	deb = min(num1, num2)
	fin = max(num1, num2)
	c1 = zeros(Int64,nbVar); c2 = zeros(Int64,nbVar)

	subseq1 = p1[deb:fin]
	subseq2 = p2[deb:fin]
	c1[deb:fin] = subseq1
	c2[deb:fin] = subseq2
	remaining1 = filter(x -> !in(x,subseq1), p2)
	remaining2 = filter(x -> !in(x,subseq2), p1)

	for i in 1:nbVar
		if c1[i] == 0
			c1[i] = popfirst!(remaining1)
		end
		if c2[i] == 0
			c2[i] = popfirst!(remaining2)
		end
	end

	return c1,c2
end

#Mutation by random swapping
function mutation(x::Array{Int}, vars::Problem_Variables)
	nbVar = size(x,1)
	oldX = x

	num1 = rand(1:nbVar)
	num2 = rand(1:nbVar)
	#Pas indispensable mais ça serait dommage de gâcher une occurence de mutation
	while num1 == num2
		num2 = rand(1:nbVar)
	end

	swap1 = min(num1, num2)
	swap2 = max(num1, num2)
	tmp = x[swap1]
	x[swap1] = x[swap2]
	x[swap2] = tmp

	return x
end


#-------------------------------------------------------------------


#---------------------- PB specific Operators ----------------------


function getObjectiveValues(x::Array{Int}, vars::Problem_Variables)
	z1 = 0
	z2 = 0
	for i = 1:vars.n
		z1 += vars.c1[i,x[i]]
		z2 += vars.c2[i,x[i]]
	end

	#println("z1 = ",z1,"; z2 = ", z2)
	return z1,z2
end


function isValid(x::Array{Int}, vars::Problem_Variables)
	i = 1
	total = 0 #weight of the tasks

	while i <= vars.n && total <= vars.limit
		total += vars.weights[i,x[i]]
		i += 1
	end

	return total <= vars.limit
end

#local search in the 1-swap neighbourghood (to improve, takes way too long)
function localsearch(x::Array{Int}, vars::Problem_Variables)
	i = 1
	acc1, acc2 = getObjectiveValues(x, vars)
	tmp1, tmp2 = acc1, acc2

	while i < size(x,1) && (tmp1 >= acc1 || tmp2 >= acc2)
		j = 1
		while j < size(x,1) && (tmp1 >= acc1 || tmp2 >= acc2)
			tmpSol = deepcopy(x)
			tmpSol[i] = x[j]
			tmpSol[j] = x[i]
			tmp1, tmp2 = getObjectiveValues(tmpSol,vars)
			if isValid(tmpSol,vars) && tmp1 < acc1 && tmp2 < acc2
				x = tmpSol
			end
			j += 1
		end
		i += 1
	end

	return x
end


#-------------------------------------------------------------------


#-------------------------- Core Functions -------------------------

# args : number of individuals, problem variables and optionnal stop condition (number of generations)
function nsga2(fname, nbInd::Int = 50, nbGen::Int = 100, mchance::Float64 = 0.1, progress = false)
	
	#Getting Problem variables
	n, lim, c1, c2, weights = load2RCAP(fname)
	vars = Problem_Variables(n, lim, c1, c2, weights)

	#initializing plotting
	figure("Graph",figsize=(6,6)) # Create a new figure
	title("Espace des objectifs")
	xlabel("z1")
	ylabel("z2")
	show()
	xlim(0, sum(c1))
	ylim(0, sum(c2))

	#initiating first generation
	gen = Generation(nbInd,firstGen(nbInd,vars),zeros(Int,nbInd),[Int[] for i=1:nbInd], zeros(Float64,nbInd),mchance,0,0,0,0)
	gen.best1, gen.best2 = getObjectiveValues(gen.people[1], vars)
	gen.worst1, gen.worst2 = gen.best1, gen.best2
	
	NDsort(gen, vars)
	crowdingDistance(gen, vars)
	for i = 1:gen.nbInd
		gen.people[i] = localsearch(gen.people[i], vars)
	end
	printgraph(gen, vars,0)
	println(gen.people)
	
	for i = 1:nbGen
		println("----------", i,"th Generation Starting ----------")
		offsprings = []
		#Generating offsprings through crossover and mutation operators
		for j = 1:2:gen.nbInd
			p1 = tournament(gen)
			p2 = tournament(gen)

			c1,c2 = OX_crossover(gen.people[p1],gen.people[p2], vars)

			if rand() < gen.mchance 
				c1 = localsearch(c1, vars)
			end
			if rand() < gen.mchance
				c2 = localsearch(c2, vars)
			end

			push!(offsprings, c1)
			push!(offsprings, c2)
		end	

		#Here we use "unique" to remove multiple instances of a solution that might erase some pareto values
		gen.people = unique(cat(gen.people,offsprings,dims=1))
		NDsort(gen, vars)
		crowdingDistance(gen, vars)
		println("Ranked everything")
		k=1
		newPop = []

		while size(newPop,1)+size(gen.fronts[k],1) <= gen.nbInd && k < gen.nbInd
			l = 1
			while l <= size(gen.fronts[k],1) && size(newPop,1) <= gen.nbInd
				push!(newPop,gen.people[gen.fronts[k][l]])
				l+=1
			end
			k+=1
		end

		#If a front is too big to fit, we order it by CD and insert what we can
		if size(newPop,1) < gen.nbInd && size(gen.fronts[k],1) != 0
			sort!(gen.fronts[k],by = x -> gen.crowddist[x])
			l = 1
			while size(newPop,1) < gen.nbInd
				push!(newPop, gen.people[gen.fronts[k][l]])
				l += 1
			end
		end

		gen.people = newPop
		

		NDsort(gen, vars)
		crowdingDistance(gen, vars)
		println("Ranked the new set !")

		println("Starting localsearch")
		for i = 1:2
			map(x -> localsearch(gen.people[x], vars),gen.fronts[i])
		end
		println("Localsearch done")

		if i%10 == 0 && progress
			printgraph(gen, vars,i)
		end
		#readline()
		println("Pareto front : ", gen.fronts[1])
		println("Objective values : ", map(x -> getObjectiveValues(gen.people[x],vars), gen.fronts[1]),"\n")
	end

	printgraph(gen, vars,nbGen)
	#vOptSolution(vars)
end

function printgraph(gen::Generation, vars::Problem_Variables, genNb::Int)
	clf()
	plot(map(x -> getObjectiveValues(gen.people[x], vars)[1], gen.fronts[1]), 
		 map(x -> getObjectiveValues(gen.people[x], vars)[2], gen.fronts[1]), "go")
	for i = 2:gen.nbInd
		plot(map(x -> getObjectiveValues(gen.people[x], vars)[1], gen.fronts[i]), 
			 map(x -> getObjectiveValues(gen.people[x], vars)[2], gen.fronts[i]), "r.")
	end
	show()
	savefig(string("plots/",vars.n,"_",vars.limit,"/gen",genNb,".svg"))
end


#Résoud un 2LAP sans limite de poids
function vOptSolution(vars::Problem_Variables)
	figure("Graphe vOpt",figsize=(6,6)) # Create a new figure
	title("Espace des objectifs")
	xlabel("z1")
	ylabel("z2")
	clf()
	#Solver is already Przybylski2008 by default
	z1,z2, S = vSolve(set2LAP(vars.n,vars.c1,vars.c2))
	plot(z1,z2, "go")
	savefig(string("plots/",vars.n,"_",vars.limit,"/vOptSpecific.svg"))
	show()	
end

using vOptGeneric
using GLPK,GLPKMathProgInterface
#Résoud le problème du sujet et affiche Y_N
function vOptGenericSolution(vars::Problem_Variables)

	model = vModel(solver = GLPKSolverMIP())

	@variable(model, x[1:vars.n, 1:vars.n], Bin)
	@addobjective(model, Min, sum(vars.c1[i,j]*x[i,j] for i=1:vars.n, j=1:vars.n))
	@addobjective(model, Min, sum(vars.c2[i,j]*x[i,j] for i=1:vars.n, j=1:vars.n))
	@constraint(model, sum(vars.weights[i,j]*x[i,j] for i=1:vars.n, j=1:vars.n) <= vars.limit)
	@constraint(model, cols[i=1:vars.n], sum(x[i,j] for j=1:vars.n) == 1 )
	@constraint(model, rows[j=1:vars.n], sum(x[i,j] for i=1:vars.n) == 1 )

	solve(model, method = :epsilon, step = 0.5)
	
	Y_N = getY_N(model)
	f1 = map(y -> y[1], Y_N)
	f2 = map(y -> y[2], Y_N)
	figure("Graphe vOpt",figsize=(6,6)) # Create a new figure
	title("Espace des objectifs")
	xlabel("z1")
	ylabel("z2")
	plot(f1,f2, "go")
	show()
end

#=
include("load2RCAP.jl")
include("nsga2.jl")
n, lim, c1, c2, weights = load2RCAP("Instances2obj/5_11")
vars = Problem_Variables(n, lim, c1, c2, weights)
nbInd = 10; gen = Generation(nbInd,firstGen(nbInd,vars),zeros(Int,nbInd),[Int[] for i=1:nbInd], zeros(Float64,nbInd),0.01,0,0,0,0)
NDsort(gen,vars)
crowdingDistance(gen, vars)
=#
