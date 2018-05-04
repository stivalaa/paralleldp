echo threads iter ms dcmiss  > mundara.dcmiss.rtab ;
( cd mundara_output >/dev/null ;
n=1; while [ $n -le 32 ]; do i=1; while  [ $i -le 10 ]; do dcmiss=`grep ' exit' dcmisstest_${n}_${i}.out | awk '{print $4}'`; etime=`grep 'elapsed time' dcmisstest_${n}_${i}.out | awk '{print $3}'` ; echo $n $i $etime $dcmiss >> ../mundara.dcmiss.rtab ; i=`expr $i + 1`; done ; n=`expr $n + 1`; done
)
