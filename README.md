# Plant_Health_sens_spec_workflow
Code and data repository for paper 'Unlocking plant health survey data: an approach to quantify the sensitivity and specificity of visual inspections'


# (1) WorkflowSimulationStudyScript_feb25.R
Script for simulating observation datasets of plant health inspectors surveying a given number of trees across two sites (with different true disease prevalence) with simulated sensitivity and specificity values.
Bayesian modelling is applied to this observation dataset using 'Plummer M. Just Another Gibbs Sampler (JAGS) Version 4.3.0. 2017.' Output files of model estimates and observation datasets are saved.

# (2) Processing simulation output_Feb25.R
Script to process the output from the modelling in (1) to provide model estimates of sensitivity and specificity values of individual simulated surveyors and compare with the actual sensitivity and specificity values
of simulated surveyors for calculation of mean squared error values.
Script also contains the code to combine the individual estimates of simulated surveyor sensitivity and specificity to provide estimates of the distribution of all surveyors' sensitivity and specicity, and compare against the distributions surveyor sensitivity and specificity were simulated from.
Script for plotting both these outputs is provided.

# (3) Visualising output errors_Feb25.R
Script to visualise the model error estimates by model priors - reads file output_errors_feb25.csv

# (4) Feb25_AOD analyses code.R
Script to analyse the sensitivity and specificity of real-world citizen science observation data for the sensitivity and specificity of detecting three acute oak decline symptoms (stem bleeds, agrilus beetle exit holes, bark cracking)
using expert data as a 'gold-standard' dataset. Uses files 'Citizen scientist observations.csv' and 'Expert Observations.csv'
