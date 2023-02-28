# lmmpred
Linear Mixed Model Prediction
This project is to create biomarker predictions from linear mixed models.

Programs here will be reading in data from the ACCORD clinical trial.  

The population of subjects is those in the Lipid half of the ACCORD trial.

Evaluable subjects will be those in the population with >= 4 observations for each of three biomarkers, HBA1C, LDL, and SBP.

Evaluable subjects also must have a measured biomarker value at time S+TAU.

For each evaluable subject, the program should create predictions for each biomarker at S+TAU and compare those predictions to the observed values.

Predictions will be based on fitting certain linear mixed models to the population of subjects biomarker data with the current evaluable subjects data censored at time S.  
