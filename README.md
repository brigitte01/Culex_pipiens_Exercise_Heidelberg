\# Culex pipiens Habitat Suitability Modelling — Heidelberg 2026



\## Overview



This repository contains R scripts, datasets, and outputs for modelling the habitat suitability of \*Culex pipiens\* mosquitoes in Italy and projecting the predictions across Europe.



A MaxEnt model was trained using Italian mosquito surveillance data and environmental predictors, then applied to European environmental data to produce a habitat suitability map.



\## Data



Five environmental predictors were used:



\- Built-up land proportion

\- Cropland proportion

\- Air temperature

\- Relative humidity

\- Total precipitation



\*\*Italy:\*\* Land cover derived from MODIS (2008–2020) and climate variables from ERA5-Land (2014–2022).



\*\*Europe:\*\* Land cover derived from MODIS (2024) and climate variables from CERRA (2024–2025).



The Italian mosquito surveillance data were processed to obtain 626 unique presence raster cells.



\## Repository Structure



\- `Data/Italy/` — Italian environmental predictors, mosquito presence data, and background points.

\- `Data/Europe/` — European environmental predictors.

\- `Scripts/` — R scripts for data preparation, model training, and prediction.

\- `Outputs/Italy/` — MaxEnt model, evaluation results, and plots.

\- `Outputs/Europ/` — European habitat suitability prediction.



\## Analysis Workflow



The analysis follows five main steps:



1\. Prepare and align environmental predictors.

2\. Process Italian mosquito surveillance data.

3\. Generate 2,000 KDE-weighted background points.

4\. Tune and train the MaxEnt model using spatial block cross-validation.

5\. Project the selected Italian model across Europe.



The corresponding R scripts are available in the `Scripts/` folder.



\## Main Output



`Outputs/Europ/Culex\_pipiens\_MaxEnt\_Europe\_CERRA\_2024\_2025.tif`



This GeoTIFF contains predicted habitat suitability across Europe, with values ranging from 0 to 1.



