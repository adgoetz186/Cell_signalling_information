# Cell Signalling Information

This code was used to generate part of the data for the paper "[The ability to sense the environment is heterogeneously distributed in cell populations](https://www.biorxiv.org/content/10.1101/2023.03.07.531554v2.abstract)" by Andrew Goetz*, Hoda Akl*, Purushottam D Dixit *contributed equally

doi: [https://doi.org/10.1101/2023.03.07.531554](https://doi.org/10.1101/2023.03.07.531554)

# Project Organization

This project is structured around the generation of figures for the paper. As a result, programs which are directly used for making a specific figure in the paper are stored in a directory reflects the figure and subfigure they are used to generate. "[Figure_Generating_Programs](Mutual_Information_Main/Figure_Generating_Programs)" is the main directory which contains all figure generating programs. Each figure generating program is stored in a directory along with a readme file which provides instructions for generating the data used to obtain the figure. All generated figure images are stored in "[Figures](Mutual_Information_Main/Figures)".

The "[functions](Mutual_Information_Main/functions/)" directory contains custom function libraries and encompasses the bulk of the machinery of the project.

The "[General_Data_Generation](Mutual_Information_Main/General_Data_Generation/)" directory contains programs used to generate essential data for plots but do not directly generate plots themselves.

The "[Initialization_File](Mutual_Information_Main/functions/)" directory contains the ".pth" file which may be used to add the project to sys.path.

The "[Mathematica_Moment_Code](Mutual_Information_Main/Mathematica_Moment_Code/)" directory contains mathematica code used to generate and process the moment closure equations.

# Additional Notes

## Context and Motivation of Work

At the current time the paper manuscript is under construction. It will be linked here to provide context for the programs when completed. Until that point, a [poster](Mutual_Information_Final_Version/Mutual_Information_Poster.png) on the work will be provided to give context and scientific motivation for many of the algorithms and figures contained in this project.
