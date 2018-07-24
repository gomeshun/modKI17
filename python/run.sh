v=v001

for i in `seq 0 9`
do
fnmem=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
fnfg=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
pchain=params_chain_noMCMC_plummerlike_tot_${v}_${i}.csv
mclog=mclog_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
logtxt=modKI17_noMCMC_plummerlike_tot_${v}_${i}.log.txt
python parameter_estimation.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
done
wait
for i in `seq 10 19`
do
fnmem=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
fnfg=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
pchain=params_chain_noMCMC_plummerlike_tot_${v}_${i}.csv
mclog=mclog_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
logtxt=modKI17_noMCMC_plummerlike_tot_${v}_${i}.log.txt
python parameter_estimation.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
done
wait
for i in `seq 20 29`
do
fnmem=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
fnfg=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
pchain=params_chain_noMCMC_plummerlike_tot_${v}_${i}.csv
mclog=mclog_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
logtxt=modKI17_noMCMC_plummerlike_tot_${v}_${i}.log.txt
python parameter_estimation.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
done
wait
for i in `seq 30 39`
do
fnmem=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
fnfg=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
pchain=params_chain_noMCMC_plummerlike_tot_${v}_${i}.csv
mclog=mclog_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
logtxt=modKI17_noMCMC_plummerlike_tot_${v}_${i}.log.txt
python parameter_estimation.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
done
wait
for i in `seq 40 49`
do
fnmem=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
fnfg=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
pchain=params_chain_noMCMC_plummerlike_tot_${v}_${i}.csv
mclog=mclog_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
logtxt=modKI17_noMCMC_plummerlike_tot_${v}_${i}.log.txt
python parameter_estimation.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
done
wait
