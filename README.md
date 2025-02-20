# General readme

Land-use intensification is considered a major driver of biodiversity loss. However, most evidence for land use effects come from space-for-time substitutions, comparing biodiversity in areas differing in current land use intensity. Despite their frequent use, the robustness of space-for-time substitutions in capturing temporal changes in biodiversity remains unclear. We compared spatial and temporal responses of plant and arthropod communities to changes in land-use intensity, using 150, ten-year time series from grasslands in three German regions (Schwäbische Alb, Hainich-Dün, Schorfheide-Chorin) from sampling years 2008-2018. 
The data have been collected within the long-term biodiversity monitoring project 1374 “The Biodiversity Exploratories” funded by the German Research Foundation. The files in this project provide the R code for dataset assembly and analysis of spatial and temporal pair-wise comparisons of plant- and arthropod-community and land-use data, as well as for the graphical visualisation of results. This R code was used to obtain the results and their graphical visualisation of the research publication. 
The source data used in the provided scripts are archived in the project-owned information system BEXIS (https://www.bexis.uni-jena.de/). See the dataset description document (“dataset_description_GitHub. docx”) for details on the used datasets and their accessibility in BEXIS. First, this script provides the source code that assembles the data used in the analyses and is provided in its assembled form already in BEXIS (dataset id: 31921, 31924, 31925, 31926). Second, the files in this project provide the R code for the analytical steps for estimating spatial and temporal diversity responses of plant- and arthropod-communities to land-use changes. Third, the files contain the R code for creating and assembling the figures used in the abovementioned publication.
The scripts in this project shall be run in the following order (see the readme of each section for details):


1)	Dataset creation 
2)	Data analysis: general dissimilarity models (GDM)
3)	Data analysis: linear mixed effects models (LM)
4)	Plotting of results
