# Precipitation isotopes data analysis
Data analysis and visualization of $\delta^{18}O$ and $\delta^{2}H$ values for the 40-year record of precipitation isotopes in Ljubljana. The results obtained using this analysis are presented in Gačnik et al. (2026), "Climate change reflected in 40-year isotopic composition trends of precipitation in Slovenia", _Environmental Research_

## Instructions:
1) Open the _"Data analysis and visualization.R"_ in the _"scripts"_ directory, preferably with RStudio as the text editor.
2) The first line of code in the opened script should be ran only once and then uncommented, as it installs all the needed R packages for running the script.
3) Set path parameters (_path_scripts_, _path_data_, and _path_figures_) to match the paths of your downloaded GitHub repository. The path parameters are located in the lines 11-13 of the _"Data analysis and visualization.R"_ script.
Pay attention to structure paths with "/" and not with "\\"  (the default for Windows paths).
4) Run the code line by line with CTRL + ENTER from line 11 onward.
    Notes: 
    - The part of the script where Hilbert-Huang transform (HHT) and Ensemble Empirical Mode Decomposition (EEMD) are performed (starting at line 441) will require up to 20 min to complete the HHT & EEMD for all variables and all stations.
    (3 stations with 6 variables each). Consider running the script on high-performance computing clusters to decrease the computing time.
    - The plot viewer in Rstudio will not show the plots in the intended aspect ratio and size. The actual figures in the intended format are saved into the directory specified by the "path_figures" parameter at the start of the script
    - The description of how the HYSPLIT timetrends were obtained for all stations is available in the main body of the Gačnik et al. (2026) manuscript.

Prepared by Jan Gačnik, Klara Žagar, István Gábor Hatvani, Zoltán Kern, and Polona Vreča. 
