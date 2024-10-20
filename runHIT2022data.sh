#!/bin/bash
# how to run:
# running with command (for example):
#./runHIT2022data.sh

echo "$@"

runNumber='4742'
#runNumber='4742 4743 4744 4745 4828'
#Runs without target:
#runNumber='4742 4743 4744 4745 4766 4801 4828 4830 4837 '
#Runs in local:
#runNumber='4742 4743 4744 4745 4828 4895 4896 4897 4898 4899 4900'
#Runs total:
#runNumber='4742 4743 4744 4745 4766 4801 4828 4830 4837 4895 4896 4897 4898 4899 4900 4901 4902 4903 4904 4905 4906'

echo Running Macro \in Batch
    
for x in $runNumber
    do
echo submit run number $x

#./runMacroBatchT1.sh -i /mnt/c/Users/Lorenzo/Desktop/shoe/DataDecoded/Decoded_HIT2022_"$x".root -o /mnt/c/Users/Lorenzo/Desktop/shoe/OutputMacro/ -p AnaFOOT -m 1 -n 1000 -r AnalyzeFOOT

#./runMacroBatchT1.sh -i /mnt/c/Users/Lorenzo/Desktop/shoe/DataDecoded/Decoded_HIT2022_"$x".root -o /mnt/c/Users/Lorenzo/Desktop/shoe/OutputMacro/ -p AnaFOOT_TW -m 1 -n 1000 -r AnalyzeTWChargeTime

./runMacroBatchT1.sh -i /mnt/c/Users/Lorenzo/Desktop/shoe/DataDecoded/Decoded_HIT2022_"$x".root -o /mnt/c/Users/Lorenzo/Desktop/shoe/OutputMacro/ -p AnaFOOT_Calo -m 1 -n 1000 -r AnalyzeCalo

done
echo  All the files needed for HIT2022 analysis are submitted!
