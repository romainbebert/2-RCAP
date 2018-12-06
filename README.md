# 2-RCAP

Bi-objective RCAP implementation using NSGA-II

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
