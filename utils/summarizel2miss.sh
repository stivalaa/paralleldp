echo threads iter ms l2dmiss > mundara.l2miss.rtab ;
( cd mundara_output ;
n=1; while [ $n -le 32 ]; do i=1; while  [ $i -le 10 ]; do l2miss=`grep ' exit' test_${n}_${i}.out | awk '{print $4}'`; etime=`grep 'elapsed time' test_${n}_${i}.out | awk '{print $3}'` ; echo $n $i $etime $l2miss >> ../mundara.l2miss.rtab ; i=`expr $i + 1`; done ; n=`expr $n + 1`; done
)
