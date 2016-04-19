# HSI-CBI-validation
Validating spatial species distribution models with presence-only data

© 2016, Jens G. Froese (jens.froese@uq.net.au)

## Description
This R script describes a methodology for validating a spatially-explicit numerical model prediction (habitat suitability index, probability of occurrence etc.) against presence-only species records using the *Continuous Boyce Index (CBI)* method (Hirzel et al. 2006 after Boyce et al. 2002)

Our implementation allows validation of multiple model predictions (here: habitat suitability index (HSI)) against multiple validation data sets.

The R script requires two `.TXT` files for each model/validation data combination:

1. Expected HSI across validation background:
   + define justifiable validation background (representative of data)
   + mask raster layer describing the spatial distribution of the model-predicted value by the validation background
   + export masked raster attribute table to `.TXT` with 3 columns: [ID], [predicted value, e.g. HSI], [pixel count]
2. Predicted HSI at presence records:
   + convert presence records into raster layer
   + combine masked raster layers of the model-predicted value and presence records
   + export combined raster attribute table to `.TXT` with 4 columns: [ID], [pixel count], [predicted value, e.g. HSI], [number of presence records per pixel (esp. if coarse grain, multiple records may fall within one pixel]

For each model prediction, the R script returns:
   + the *CBI* validation metric
   + the proportion of the validation background with HSI above a user-defined threshold
   + a plot of the predicted-to-expected ratio, i.e. the proportion of presence records with a given HSI (or range of HSI) divided by the proportion of the validation background covered by that HSI (or range of HSI)

## References
* Boyce, M.S. et al. 2002. Evaluating resource selection functions. — Ecological Modelling 157: 281-300.
* Broennimann, O. 2015. Package 'ecospat': spatial ecology miscellaneous methods. URL http://cran.r-project.org/web/packages/ecospat/.
* Dowle, M. et al. 2015. Package 'data.table': extension of data.frame. URL https://github.com/Rdatatable/data.table/wiki/.
* Hijmans, R.J. 2015. Package 'raster': geographic data analysis and modeling. URL http://cran.r-project.org/web/packages/raster/.
* Hirzel, A.H. et al. 2006. Evaluating the ability of habitat suitability models to predict species presences. — Ecological Modelling 199: 142-152.
* RCoreTeam 2015. R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
* Zeileis, A. et al. 2015. Package 'zoo': S3 infrastructure for regular and irregular time series. URL http://zoo.R-Forge.R-project.org/.
