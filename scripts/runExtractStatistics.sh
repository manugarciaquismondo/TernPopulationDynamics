BASE_DIR=${PWD}/../
RESULTS_DIRECTORY=new_results
cd ${BASE_DIR}${RESULTS_DIRECTORY}
for directory in `ls -d */`
do
  Rscript ${BASE_DIR}r-scripts/extractMetrics.R ${BASE_DIR}${RESULTS_DIRECTORY}/${directory} ${BASE_DIR}statistics/${directory} ${BASE_DIR}data
done