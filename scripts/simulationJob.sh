#!bin/bash
# This is a sample PBS script. 
# Specify which shell to use
#$ -S /bin/bash
#   
#
#   Request 20 hours of walltime
#
#$ -l h_vmem=4G
#$ -l h_rt=20:00:00
#
# Maintain environment variables
##$ -V
#
#   The following is the body of the script. By default,
#   PBS scripts execute in your home directory, not the
#   directory from which they were submitted. The following
#   line places you in the directory from which the job
#   was submitted.
#
#
#   Send mail when job begins
#$ -m b
#   Send mail when job ends
#$ -m e
#   Send mail when job aborts
#$ -m a

DATA_FOLDER=$1
RESULTS_FOLDER=$2
BASE_FOLDER=$3
SIMULATION_OUTPUT=$4
SIMULATION_YEARS=$5
echo "Reading the data from "${DATA_FOLDER}" and writing the results into "${RESULTS_FOLDER}". Base folder is "${BASE_FOLDER}". Output file is "${SIMULATION_OUTPUT}". Simulation years is "${SIMULATION_YEARS}
Rscript ${BASE_FOLDER}/MetapopulationTernModel.R ${RESULTS_FOLDER} ${DATA_FOLDER} ${SIMULATION_YEARS} &> ${SIMULATION_OUTPUT}