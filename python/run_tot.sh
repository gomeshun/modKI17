v=v003
n=4
kmax=3

for j in `seq 0 $((n-1))`
do
  for k in `seq 0 $((kmax-1))`
  do
    i=$((n*j+k))
    echo "calc start for ${i}th mock."
    fnmem=mock_tot_noMCMC_${v}_${i}.csv
    fnfg=mock_tot_noMCMC_${v}_${i}.csv
    pchain=params_chain_KI17_noMCMC_tot_${v}_${i}.csv
    mclog=mclog_KI17_noMCMC_tot_${v}_${i}.csv
    logtxt=KI17_noMCMC_tot_${v}_${i}.log.txt
    python parameter_estimation.py ${fnmem} ${fnfg} ${pchain} ${mclog} --no-confirmation --burnin 5000 --burnin-repeat 3 --sampling 20000 > ${logtxt} &
  done
  wait
done
