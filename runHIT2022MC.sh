
#!/bin/bash
# how to run:
# running with command (for example):
#./runHIT2022MC.sh

echo "$@"

runNumber='100'
#runNumber='100 140 200 220'

echo Running Macro \in Batch

for x in $runNumber 
    do
echo submit MC energy $x
./runMacroBatchT1MC.sh -i /mnt/d/DecodedMC_HIT2022_MC_"$x".root -o /mnt/c/Users/Lorenzo/Desktop/shoe/OutputMacro/MC/TW -p AnaFOOT_TW -m 1 -n 1000 -r AnalyzeTWFragMC
done
echo All the jobs for HIT2022_MC submitted
