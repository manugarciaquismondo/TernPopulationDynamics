BASE_DIR=${PWD}/../
SIMULATION_YEARS=50
RESULTS_DIRECTORY=new_results
mkdir ${BASE_DIR}/errors ${BASE_DIR}output ${BASE_DIR}data ${BASE_DIR}${RESULTS_DIRECTORY} ${BASE_DIR}sim_outputs
for simulation in `seq -f "%03g" 1 100`
do
  mkdir ${BASE_DIR}${RESULTS_DIRECTORY}/simulation_${simulation}
	qsub -e ${BASE_DIR}errors/error_${simulation}.txt -o ${BASE_DIR}output/output_${simulation}.txt ./simulationJob.sh ${BASE_DIR}data ${BASE_DIR}${RESULTS_DIRECTORY}/simulation_${simulation} ${BASE_DIR}r-scripts ${BASE_DIR}sim_outputs/output_${simulation}.txt ${SIMULATION_YEARS}
#	echo ${DATASET} ${CLUSTERS}
done