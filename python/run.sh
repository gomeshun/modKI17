v=v001

for i in `seq 0 49`
do
python parameter_estimation.py mock_tot_noMCMC_${v}_${i}.csv mock_tot_noMCMC_${v}_${i}.csv params_chain_noMCMC_tot_${v}_${i}.csv mclog_noMCMC_tot_${v}_${i}.csv --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > noMCMC_tot_${v}_${i}.log.txt &
done
