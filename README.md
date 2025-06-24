# Plant_Health_sens_spec_workflow
Code and data repository for paper 'Unlocking plant health survey data: an approach to quantify the sensitivity and specificity of visual inspections'


# (1) WorkflowSimulationStudyScript_feb25.R
Script for simulating observation datasets of plant health inspectors surveying a given number of trees across two sites (with different true disease prevalence).
Bayesian modelling is applied to this observation dataset using 'Plummer M. Just Another Gibbs Sampler (JAGS) Version 4.3.0. 2017.' to estimate the sensitivity and specificity
of surevyors in the absence of a gold-standard dataset. Output files of model estimates and observation datasets are saved. Uses input file 'AODprocesseddata.csv' to simulate surveyor sensitivity and specificity based on AOD stem bleed dataset.

# (2) Processing simulation output_Feb25.R
Script to process the output from the modelling in (1) to provide model estimates of sensitivity and specificity values of individual simulated surveyors and compare with the actual sensitivity and specificity values
of simulated surveyors.
Script also contains the code to combine the individual model estimates of simulated surveyor sensitivity and specificity values to provide estimates of the distribution of all surveyors' sensitivity and specicity. This model estimate is compared against the distributions that surveyor sensitivity and specificity were simulated from.
Script for plotting both these outputs is provided.

# (3) S3 File.xlsx 
Output statistics from model runs from (1) prsented in the paper (individuals and distributions) with mean squared error values, and details of the prior probability distributions. This now contains data on the runs presented in the final manuscript with misspecified
ture disease prevalence prior estimates.

# (4) Visualising output errors_Feb25.R
Script to visualise the model error estimates by model priors - reads file output_errors_feb25.csv

# (5) Feb25_AOD analyses code.R
Script to analyse the sensitivity and specificity of real-world trained citizen science observation data for the sensitivity and specificity of detecting three acute oak decline symptoms (stem bleeds, agrilus beetle exit holes, bark cracking)
using expert data as a 'gold-standard' dataset. Uses files 'Citizen scientist observations.csv' and 'Expert Observations.csv'

# (6) S2 File.xlsx
Summary of citizen scientist surveyor sensitivity and specificity values from (5)
