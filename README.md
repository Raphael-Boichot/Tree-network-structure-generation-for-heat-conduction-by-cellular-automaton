# Tree network structure generation for heat conduction by Cellular Automaton

This repo contains the codes used in the following paper entitled:
[Tree-network structure generation for heat conduction by cellular automaton, by R. Boichot, L. Luo, Y. Fan - Energy Conversion and Management, 2009](https://doi.org/10.1016/j.enconman.2008.09.003)

More exactly, this is a Matlab update of the original code made in Visul Basic 6 (and lost since). The code is easy to use. Enter the [filling ratio and the ratio of conductivity of two materials](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Codes/Main_CA.m) on a heating surface linked to a localized heat sink and it makes the conductive matter (in dark) evolve following temperature gradients. The shape obtained is not stricto sensu optimal as there is no objective function in the code but presents a very efficient design to cool a distributed heated surfaces like a computer chips, battery stacks, some parts of fuel cells, etc. Cooling effectiveness could be even more increased by using a mixed criterium between thermal gradients and temperature at the interface to guide the cellular automaton ([modify the code here](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/421a51fb88f8051d2978c1b49b94973a7481aa89/Codes/automate_cell_direct.m#L105) to activate this feature).

**Code free to use, please cite the author according to the license !**

The code is based on an finite difference approximation of the temperature equation solved on arbitrary domains. It uses a sparse direct solver and can converges on moderately powerfull computer within an hour. The cellular automaton algorithm part by itself is a pure fabrication of this study and has no anteriority. Funfacts: this study was attacked for plagiarism without any success and the strikers went crying as soon as I took out my lab notebook. The code principle was copied (without citing me) and patented (citing me !) by a noticeable number of assholes publishing in countries were deontology seems to be a big joke.

# Step 1: arbitrary starting shape
We start with an arbitrary cooling shape with gradients at the interface between cooling material (dark) and heating surface (white). The heat can flow only through a localized heat sink (blue).
![Step1](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP1.png)

# Step 2: Cellular automaton driven by gradients
An evolutionnary algorithm (Cellular Automaton) is used to evolve the shape by attracting conductive matter along temperature gradients. It mimicks the chemotaxis of roots in a nutrient medium that follows concentration gradients. Conductive matter evolves like a growing organism trying to maximize food extraction (heat) in a closed space (heating pad). 
![Step2](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP2.png)

# Step 3: convergence
The shape ceases to evolve once thermal gradients are equalized around conductive matter. The solution is not optimal as no objective function is minimized but the topology evolved is an efficient local minima and most of all: it converges very fast ! It can outperform most of the other methods based on pure mathematics to attack the problem. It is however an average performer compared to other evolutionnaty algorithms (both gradient based or global).
![Step3](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP3.png)

# Code output during calculation
Science must be fancy so the code outputs cool figures while converging. As the problem is mainly symmetrical, only half geometries are considered for calculations. Everything is in French.
![output](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/Code_Output.png)

# Exemple of shape at convergence with a 800x800 square elements topology, filling ratio=30% and kp/k0=20
![convergence](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/Converged_shape.png)
