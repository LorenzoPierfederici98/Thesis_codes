#!/bin/bash
# how to run (example)
#./runMacroBatchT1.sh -i /storage/gpfs_data/foot/mtoppi/DataDecoded/fullGSI/Merge_GSI2021_4306.root -o /storage/gpfs_data/foot/mtoppi/OutputMacro/ -p AnaFOOT -m 0 -n 100

INPUT_BASE_PATH="."
OUTPUT_BASE_PATH="."
SHOE_PATH="shoe"

while getopts i:o:p:m:n:r: flag
do
    case "${flag}" in
        i) inFile=${OPTARG};;
        o) outFolder=${OPTARG};;
        p) prefix=${OPTARG};;
        m) max=${OPTARG};;
        n) nev=${OPTARG};;
        r) runMacro=${OPTARG};;
    esac
done

echo "$@"

if [ -z ${runMacro+x} ]; then
    echo "runMacro is unset"
    runMacro="AnalyzeFOOT"
    echo "set to '$runMacro'"
else
    echo "Name of Macro to run is set to '$runMacro'"
fi

if [ -z ${prefix+x} ]; then
    echo "prefix is unset"
    prefix="Ana"
    echo "set to '$prefix'"
else
    echo "prefix is set to '$prefix'"
fi

if [ -z ${max+x} ]; then
    echo "max is unset"
    max=1
    echo "set to $max"
else
    echo "max is set to '$max'"
fi

if [ -z ${nev+x} ]; then
    echo "nev is unset"
    nev=100
    echo "set to $nev, while max is $max"
else
    echo "nev is set to '$nev'"
fi

if [ ! -f "$inFile" ]; then
    echo "Input file ${inFile} not found!"
    exit 1
fi

# Remove slash from input file if present
if [[ ${inFile: -1} == "/" ]]; then
    inFile=${inFile::-1}
fi

if [ ! -d "$outFolder" ]; then
    mkdir -p $outFolder
    if [ $? -ne 0 ]; then
        echo "Failed to create output directory. Exiting"
        exit 1
    fi
    echo "Directory ${outFolder} did not exist, created now!"
fi

# Remove slash from output folder if present
if [[ ${outFolder: -1} == "/" ]]; then
    outFolder=${outFolder::-1}
fi

echo
echo "-----------------------------------------------------"
echo "Running on file = $inFile"
echo "Output folder = $outFolder"
echo "-----------------------------------------------------"
echo

# Create folder for HTC auxiliary files if not present
HTCfolder="${outFolder}/HTCfiles"

if [ ! -d $HTCfolder ]; then
    mkdir -p $HTCfolder
    if [ $? -ne 0 ]; then
        echo "Failed to create HTC files directory. Exiting"
        exit 1
    fi
    echo "Directory ${HTCfolder} did not exist, created now!"
fi

outFile_base="${outFolder}/${prefix}_$(basename ${inFile})"
outFile="${prefix}_$(basename ${inFile})"
echo $outFile
echo $outFile_base

# Create executable
outFile=${outFile::-5}
jobExec="${HTCfolder}/runMacroInBatch_${outFile}.sh"
jobExec_base=${jobExec::-3}
echo $jobExec
echo $jobExec_base

# Create executable file for jobs
cat <<EOF > $jobExec
#!/bin/bash

export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/lib64/:/lib/

export HOME=\$(pwd)

cd ${SHOE_PATH}/build/Reconstruction
source ${SHOE_PATH}/build/setupFOOT.sh

root -l -b -q "$runMacro".cc++g\(\"$inFile\",$max,$nev,\"$prefix\",\"$outFolder\"\)

if [ \$? -eq 0 ]; then
    echo "Job done successfully..."
else
    echo "Unexpected error in processing of file"
fi
EOF

# Make the job executable
chmod 754 ${jobExec}

# Create a list of jobs to run in parallel
jobList="${HTCfolder}/jobList.txt"

for i in $(seq 1 $max); do
    echo "${jobExec}" >> ${jobList}
done

# Run jobs in parallel
cat ${jobList} | parallel -j 4
rm -f $jobList
