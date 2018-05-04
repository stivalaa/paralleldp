# make actions file for oahttslftest -r from mosp_oahttslf debug output
awk '/insert/{print "insert " $2,$3} /lookup/{print "lookup "$2}' ../mosp/out | sed 's/)//g'|sed 's/,//g' > actions
