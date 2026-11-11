# Precipitation isotopes data analysis (SZGG figures)
Simplified visualization of $\delta^{18}O$ and $\delta^{2}H$ values for the 40-year record of precipitation isotopes in Ljubljana, prepared for the contribution in the SZGG proceedings (in Slovenian). The analysis is based on Gačnik et al. (2025), "Climate change reflected in 40-year isotopic composition trends of precipitation in Slovenia", Environmental Research. The full code for that article is available on the _main_ branch of this repository.

## Instructions:
1) Open the _"SZGG_figures.R"_ in the _"scripts"_ directory, preferably with RStudio as the text editor.
2) The first line of code in the opened script should be ran only once and then uncommented, as it installs all the needed R packages for running the script.
3) Set path parameters (_path_scripts_, _path_data_, and _path_figures_) to match the paths of your downloaded GitHub repository. The path parameters are located in the lines 11-13 of the _"SZGG_figures.R"_ script.
Pay attention to structure paths with "/" and not with "\\"  (the default for Windows paths).
4) Run the code line by line with CTRL + ENTER from line 11 onward.
    Notes: 
    - The Ensemble Empirical Mode Decomposition (EEMD) uses the same settings as in Gačnik et al. (2025), but without the Monte Carlo significance test, so the whole script runs in less than a minute.
    - The plot viewer in Rstudio will not show the plots in the intended aspect ratio and size. The actual figures in the intended format are saved into the directory specified by the "path_figures" parameter at the start of the script
    - The description of how the HYSPLIT timetrends were obtained is available in the main body of the Gačnik et al. (2025) manuscript.

Prepared by Jan Gačnik, Klara Žagar, István Gábor Hatvani, Zoltán Kern, and Polona Vreča. 
