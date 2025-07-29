# Precipitation isotopes data analysis
Data analysis and visualization of $\delta$^{18}O and $\delta$^{2}H values for the 40-year record of precipitation isotopes in Ljubljana. The results obtained using this analysis are presented in Gačnik et al. (_in preparation_), "Climate change reflected in 40-year isotopic composition trends of precipitation in Slovenia"

## Preparation prior use:
Prior to running R scripts, downloading shapefile data is required for generating map plots. The shapefiles used in in Gačnik et al. (_in preparation_) can be found on the data portal of Republic of Slovenia, Ministry of Natural Resources and Spatial Planning, via the following link: https://ipi.eprostor.gov.si/jgp/data. After entering the portal, follow _"Državni topografski sistem" -> "Državna pregledna karta" -> "Državna pregledna karta merila 250.000 - rastrski podatki" -> "Prenos podatkov"_ and download the zipfiles _"naselja in objekti.zip"_, _"komunikacije.zip"_, and _"hidrografija.zip"_. The zipfiles should be extracted to the directory _"shapes"_ which exists the directory where this GitHub repository is downloaded. For detailed desription of the contents and the name codes of shapefile data, consult the "shapefile_codes_description_GURS.docx" file located in the  _"shapes"_ directory.

## Instructions:
1) Open the _"Data analysis and visualization.R"_ script.
2) The first line of code should be ran only once and then deleted, as it installs all the needed R packages for running the script.
3) Set path parameters to match the paths of your downloaded GitHub repositury. The path parameters are located in the lines 11-14 of the _"Data analysis and visualization.R"_ script.
Pay attention to structure paths with "/" and not with "\"  (the default for Windows paths).
3) Run the code line by line with CTRL + ENTER.
    Notes: 
    - The part of the script where Hilbert-Huang transform (HHT) and Ensemble Empirical Mode Decomposition (EEMD) are performed will require up to 20min to complete the HHT & EEMD for all variables and all stations
    (3 stations with 6 variables each). Consider running the script on high-performance computing clusters to decrease the computing time
    - If correct shapefiles are NOT downloaded to the _"shapes"_ directory, the script will report errors when it comes to reading shapefiles and generating map plots.
    However, the rest of the script is not dependent on these shapefiles, so the code can still be used even without needed shapefiles.