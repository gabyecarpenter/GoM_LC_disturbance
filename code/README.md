This folder contains all of the scipts used in the manuscript

Matlab scripts (.m)
-Used to format and analyze output from the Connectivity Modeling System (CMS)

FGB_processing_final.m
-Ran this script for every year of CMS output
1. Formats data into matlab structure
2. Creates figures for particle density
3. Calculates distance traveled
4. Looks at connectivity

Fig1.m
- Create basemap for figure 1 - mesoscale features and regions added in Adobe Illustrator

LC_HYCOM_final.m
1. Download daily HYCOM SSH files and extract LC contour based on 0.17m SSH
2. Refine and check countours based on lat/lon
3. Calculate contour length
4. Plot daily contours based on length on a map
5. Plot histogram and calculate disturbance period medians
6. Plot LC lenth through time as line

connectivity.m
- Creates shapefiles for buffered coral reef areas
- Assigns connectivty to regions in the GoM


R markdown files (.Rmd)
- Used to create figures in the manuscript
- Imports text files of data created in matlab and re-formats for figures
- Some final aesthetic figure edits may have been made in Illustrator

