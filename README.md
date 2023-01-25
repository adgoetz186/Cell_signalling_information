# Cell_signalling_information
Code used to generate data for paper "Information Transmission by Individual Cells"

# Navigation

## Project Organization

The project was structured around generating plots for a paper. As a result programs which are directly used for making a specific figure are stored in a directory which indicates the figure and subfigure. For example, the program contained in Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/C/ generates figure 2C, which can be viewed in Mutual_Information_Final_Version/Figures/Figure_2/C/.

Custom function libraries are stored in Mutual_Information_Final_Version/functions/, this encompasses the bulk of the machinery of the project with Mutual_Information_Final_Version/functions/mutual_information_functions/MI_Calculation_Functions.py in particular being a central function library as it contains functions for obtaining the mutual information as well as finding the optimal input distribution and corresponding channel capacity.

# Additional Notes

## Context and Motivation of Work

At the current time the paper manuscript is under construction. It will be linked here to provide context for the programs when completed. Until that point, a poster on the work titled "Mutual Information Poster.png" will be provided to give context and scientific motivation for many of the algorithms and figures contained in this project.

## Arguement Files

Programs are stored along with corresponding argument files which tells the program exactly what to do. This convention was chosen to allow for multiple argument configurations to be ran and saved (enabling full reproducibility) while also cutting down on bulky path strings in program files.
