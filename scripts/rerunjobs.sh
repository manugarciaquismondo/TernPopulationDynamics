qdel -u $(whoami)
rm -f ../output/*.txt ../errors/*.txt ../sim_outputs/*.txt
./runSimulations.sh
