# GoM_LC_disease

This manuscript is currently in review at Global Change Biology Communications

Below is a breakdown of the files and how they were used in the manuscript

# /code
FGB_processing_final.m - This is the primary analysis script. It takes all output from the CMS (traj files) each year and creates euclidian distance matrices and density plots w/ non-convex hulls

LC_HYCOM_final.m - This script was used to download daily HYCOM SSH and determine the loop current

connectivity.m - This script used a shapefile of coral reef areas, buffer it to match model resolution, and determine the proportion of maticles that reach reef areas.

Chapter_3_GoM.Rmd - R script to organize and plot data

Fig1.m - matlab script to create figure 1 

# /data
Yearly folder with all trajectory output files from the CMS. Used in FGB_processing_final.m

Coral polygons

# /figures
PDF files of all manuscript figures
