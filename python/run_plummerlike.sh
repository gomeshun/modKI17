v=v002

for n in `seq 0 4`
do
  for j in `seq 0 4`
  do
    i=$((10*n+j))
    echo "calc start for ${i}th mock."
    fnmem=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
    fnfg=mock_tot_plummerlike_noMCMC_${v}_${i}.csv
    pchain=params_chain_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
    mclog=mclog_modKI17_noMCMC_plummerlike_tot_${v}_${i}.csv
    logtxt=modKI17_noMCMC_plummerlike_tot_${v}_${i}.log.txt
    python parameter_estimation_mod.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
  done
  wait
done
