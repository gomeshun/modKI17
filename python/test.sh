v=v001
for n in `seq 0 4`
do
  for i in `seq 0 9`
  do
  j=$((10*n+i))
    echo hoge${j}
  done
done
