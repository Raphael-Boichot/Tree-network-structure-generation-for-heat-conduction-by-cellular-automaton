# Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton

This repo contains the codes used in the following paper entitled:
[Tree-network structure generation for heat conduction by cellular automaton, by R. Boichot, L. Luo, Y. Fan - Energy Conversion and Management, 2009](https://doi.org/10.1016/j.enconman.2008.09.003)

More exactly, this is a Matlab update of the original code made in Visul Basic 6 (and lost since). The code is easy to use. Enter the [filling ratio and the ratio of conductivity of two materials](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Codes/Main_CA.m) on a heating surface linked to a localized heat sink and it makes the conductive matter (in dark) evolve following temperature gradients. The shape obtained is not stricto sensu optimal as there is no objective function in the code but presents a very efficient design to cool a distributed heated surface. Cooling effectiveness could be even more increased by using a mixed criterium between thermal gradients and temperature at the interface to guide the cellular automaton ([modify the code here](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/421a51fb88f8051d2978c1b49b94973a7481aa89/Codes/automate_cell_direct.m#L105) to activaye this feature).

Code free to use, please cite the author according to the license !

# Step 1: An arbitrary cooling shape with gradients at the interface between cooling material (dark) and heating surface (white)
![Step1](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP1.png)

# Step 2: Evolutionnary algorithm (Cellular Automaton) to evolve the shape by attracting conductive matter along temperature gradients
![Step2](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP2.png)

# Step 3: The shape ceases to evolve once thermal gradients are equalized around itself
![Step3](https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton/blob/main/Pictures/STEP3.png)

# Example of shape at convergence with a 800x800 square elements topology
![convergence]()
