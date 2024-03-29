## Tree network structure generation for heat conduction by Cellular Automaton

This repo contains the codes used in the following paper entitled:
[Tree-network structure generation for heat conduction by cellular automaton, by R. Boichot, L. Luo, Y. Fan - Energy Conversion and Management, 2009](https://doi.org/10.1016/j.enconman.2008.09.003)

More exactly, this is a Matlab port of the original code made in Visual Basic for application ([very first version here](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/tree/main/Initial%20code)) which took forever to converge (I was far from being a numerical method expert at that time). The Matlab code itself is hopefully easy to use, much faster and very optimized. Enter the filling ratio and the ratio of conductivity of two materials on a heating surface linked to a localized heat sink and it makes the conductive matter (in dark) evolve following temperature gradients. The shape obtained is not stricto sensu optimal as there is no objective function in the code but presents a very efficient design to cool a distributed heated surfaces like a computer chips, battery stacks, some parts of fuel cells, etc. Cooling effectiveness could be even more increased by using a mixed criterium between thermal gradients and temperature at the interface to guide the cellular automaton. The theoretical optimal solution to the problem must have equalized temperatures along the adiabatic borders which is not the case here, so it falls in some (beautiful) local optimum. Considering laws of Nature, it is probably an "enough" solution regarding its law cost and local aspect for explaining shape of plant and tree roots.

**Code free to use, please cite the author according to the license !**

The code is based on an finite difference approximation of the temperature equation solved on arbitrary domains. It uses a sparse direct solver and can converges on moderately powerfull computer within an hour. The cellular automaton algorithm part by itself is a pure fabrication of this study and has no anteriority. 

## Step 1: arbitrary starting shape
We start with an arbitrary cooling shape with gradients at the interface between cooling material (dark) and heating surface (white). The heat can flow only through a localized heat sink (blue).
![Step1](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP1.png)

## Step 2: Cellular automaton driven by gradients
An evolutionnary algorithm (Cellular Automaton) is used to evolve the shape by attracting conductive matter along temperature gradients. It mimicks the chemotaxis of roots in a nutrient medium that follows concentration gradients. Conductive matter evolves like a growing organism trying to maximize food extraction (heat) in a closed space (heating pad). 
![Step2](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP2.png)

## Step 3: convergence
The shape ceases to evolve once thermal gradients are equalized around conductive matter. The solution is not optimal as no objective function is minimized but the topology evolved is an efficient local minima and most of all: it converges very fast ! It can outperform most of the other methods based on pure mathematics to attack the problem. It is however an average performer compared to other evolutionnaty algorithms (both gradient based or global).
![Step3](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP3.png)

## Code output during calculation
Science must be fancy so the code outputs cool figures while converging. As the problem is mainly symmetrical, only half geometries are considered for calculations.
![output](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/Code_Output.png)

## Exemple of convergence path for kp/k0= 10 and filling ratio = 0.3 (10 steps per frame)
![convergence](Pictures/CA_Output.gif)

## Epilogue and funfacts

This code and some of my following papers have more or less killed the [constructal theory](https://en.wikipedia.org/wiki/Adrian_Bejan#Constructal_law) which was very famous in the early 2000s as a "theory of everything" in engineering. As killing [bullshit theories](https://cambridgescholars.com/product/978-1-5275-3839-9) is my favorite activity, I consider the job as done as this algorithm was more or less the first nail in the coffin. The author of the now forgotten constructal theory first tried to prevent my papers to be published and then tried to claim my results to be own part of his theory, which I always reject. This particular algorithm was also attacked for plagiarism without any success and the strikers went crying as soon as I took out my lab notebook with proof of anteriority. Finally, the code principle was copied (without citing me) and patented (citing me, yes this is particularly stupid as it renders the patent immediately obsolete !) by a noticeable number of assholes publishing in countries from all around the world were respecting deontology seems to be elective. 

I had clearly imagined a more peaceful begining of scientific career but hey, plagiarism is just a form of flattery after all...
That said, enjoy the code !
