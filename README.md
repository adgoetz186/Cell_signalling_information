# Cell Signalling Information

This code was used to generate part of the data for the paper "[The ability to sense the environment is heterogeneously distributed in cell populations](https://www.biorxiv.org/content/10.1101/2023.03.07.531554v2.abstract)" by Andrew Goetz*, Hoda Akl*, Purushottam D Dixit *contributed equally

# Project Organization

Programs which are directly used for making a specific figure in the paper are stored in a directory which indicates the figure and subfigure. For example, the [program](Mutual_Information_Final_Version/Figure_Generating_Programs/Figure_2/C/single_cell_responses.py) which generates [figure 2C](Mutual_Information_Final_Version/Figures/Figure_2/C/response_distributions.png) is contained in the [figure generating section](Mutual_Information_Final_Version/Figure_Generating_Programs) while figure 2C is stored in the [figure](Mutual_Information_Final_Version/Figures/) section.

The [directory](Mutual_Information_Final_Version/functions/) containing custom function libraries encompasses the bulk of the machinery of the project with the [mutual information functions](Mutual_Information_Final_Version/functions/mutual_information_functions/MI_Calculation_Functions.py) in particular being a central function library as it contains functions for obtaining the mutual information as well as finding the optimal input distribution and corresponding channel capacity.

# Additional Notes

## Context and Motivation of Work

At the current time the paper manuscript is under construction. It will be linked here to provide context for the programs when completed. Until that point, a [poster](Mutual_Information_Final_Version/Mutual_Information_Poster.png) on the work will be provided to give context and scientific motivation for many of the algorithms and figures contained in this project.
