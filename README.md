# Climate drivers of enterovirus transmission

**System requirements and installation guide:**

All code for this analysis was developed in R version 4.2.1 (2022-06-23) "Funny-Looking Kid", R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

**Main scripts:**

Fig1.R develops calculations for the center of gravity (mean timing of cases) and intensity of the epidemic across locations in China using the data published in Takahashi 2016 (reproduced here) https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001958. 

Fig2_make_empbetas_**.R constructs the dataset of empirical transmission rates for EVA71 and CVA16. 

Fig2_Regressions.R fits the regression model and plots the coefficients as displayed in Fig 2A. 

Fig3_Projections.R used the coefficients of the regression models and province-level climate projections to simulate model outbreaks under current and future climates to estimate impact on outbreak size for CVA16 and EVA71

**Data files:**

Temperature averages for Fig. 1: china_adm1_levelwt.RData

Province level weekly climate data for merging with empirical beta time series to test the effect of these variables: era5.t2m.daily.China.subregion.2009-2014.K.csv, era5.q2m.daily.China.subregion.2009-2014.kg_per_kg.csv, chirps.precip.daily.China.subregion.2009-2014.mm_per_day.csv

Climate historic and SSP585 projections at the province-level: historic_daily_proj_China.csv and ssp585_daily_proj_China.csv

Long-run temperature data for calculating impact of variability: era5.t2m.daily.China.subregion.1990-2020.K.csv

Calculated empirical betas for EVA71, CVA16 and polio: China (eva71.RData and cva16.RData), Japan (cvajp.RData and evajp.RData), USA (polioeb.RData)

Main regression model for input into projections:jp_ch_eva_tempschooling.RData, jp_ch_cva_tempschooling.RData










