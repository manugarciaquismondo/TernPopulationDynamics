BASE_DIR=${PWD}/../
mkdir ${BASE_DIR}/errors ${BASE_DIR}output ${BASE_DIR}data ${BASE_DIR}results ${BASE_DIR}sim_outputs
for simulation in `seq -f "%03g" 1 100`
do
  mkdir ${BASE_DIR}results/simulation_${simulation}
	qsub -e ${BASE_DIR}errors/error_${simulation}.txt -o ${BASE_DIR}output/output_${simulation}.txt ./simulationJob.sh ${BASE_DIR}data ${BASE_DIR}results/simulation_${simulation} ${BASE_DIR}src ${BASE_DIR}sim_outputs/output_${simulation}.txt
#	echo ${DATASET} ${CLUSTERS}
done