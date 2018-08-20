v=v001

for i in `seq 0 4`
do
python parameter_estimation_mem.py mock_mem_noMCMC_${v}_${i}b.csv params_chain_noMCMC_${v}_${i}b.csv mclog_noMCMC_${v}_${i}b.csv --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > noMCMC_${v}_${i}b.log.txt &
done
