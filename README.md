This repository contains code to accompany the manuscript titled

### **Coral color and depth drive symbiosis ecology of Montipora capitata in Kāne'ohe Bay, O'ahu, Hawai'i**

by **Teegan Innis, Ross Cunning, Raphael Ritson-Williams, Christopher Wall and Ruth Gates**

In this manuscript, we describe the spatial distribution of Symbiodinium among colonies of M. capitata, a dominant coral across Kāne‘ohe Bay, Hawai‘i. We demonstrate a striking depth-distribution, suggesting multiple symbiont types may compete within individual corals with outcomes determined interactively by attributes of the host (color) and the abiotic environment (light).

**Repository contents:**

* **Data/:** Contains field collection data and qPCR data used in this study. Also contains NOAA daily tide data used in analysis. Sub-folders are as follows:
  * *Collection Data*: Field collection data and color scores from each of 5 independent observers
  * *Tide_Tables*: NOAA daily tide tables for Coconut Island, Kāne'ohe Bay, O'ahu, Hawai'i used for depth correction in analysis
  * *coast_n83.shp*: Shape file outlining O'ahu for use in Fig. 4
  * *qPCR_data*: qPCR data for 707 sample colonies 
  * *KBMap.Rdata*: R code to call satellite image of Kāne'ohe Bay when conducting spatial analyses and formulating Fig. 4
* **Figures/:** Contains png files for each figure included in the manuscript and 2 supplementary figures. Code to generate these figures can be found in Analysis/Markdown.Rmd.
* **Analysis/:** Contains R Markdown documents with commented code to reproduce all data analysis and figures presented in the manuscript. Knitting these documents produces the HTML output. 
  * *Markdown.Rmd*: R Markdown with commented code to reproduce all data analysis and figures presented in manuscript
  * *Markdown.html*: HTML output of Markdown.Rmd
  * *Color_Score.Rmd*: R Markdown with commented code to seaparetely analyze color morph of samples colonies given 5 independent observers' scores
  * *Color_Score.html*: HTML output of Color_Score.Rmd
