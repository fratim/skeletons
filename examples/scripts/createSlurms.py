import os
import numpy as np
import sys
from skeletons.utilities import dataIO

template ='''#!/bin/bash
#
# add all other SBATCH directives here
#
#SBATCH -p {PARTITION}                                       # use the COX partition
#SBATCH -n 1                                                 # Number of cores
#SBATCH -N 1                                                 # Ensure that all cores are on one matching
#SBATCH --mem={MEMORY}                                       # CPU memory in MBs
#SBATCH -t 0-{HOURS}:00                                      # time in dd-hh:mm to run the code for
#SBATCH --mail-type=NONE                                     # send all email types (start, end, error, etc.)
#SBATCH --mail-user=tfranzmeyer@g.harvard.edu                # email address to send to
#SBATCH -o {OUTPUT_PATH}/{JOBNAME}.out                       # where to write the log files
#SBATCH -e {ERROR_PATH}/{JOBNAME}.err                        # where to write the error files
#SBATCH -J fillholes_{JOBNAME}                               # jobname given to job

module load Anaconda3/5.0.1-fasrc02
module load cuda/9.0-fasrc02 cudnn/7.1_cuda9.0-fasrc01

source activate fillholes

export PYTHONPATH=$PYTHONPATH:/n/home12/tfranzmeyer/Code/

cd {RUNCODEDIRECTORY}

python scripts/{COMMAND}

echo "DONE"

'''


def makeFolder(folder_path):
    if os.path.exists(folder_path):
        raise ValueError("Folderpath " + folder_path + " already exists!")
    else:
        os.mkdir(folder_path)

def writeFile(filename, data):
    if os.path.exists(filename):
        raise ValueError("File " + filename + " already exists!")
    else:
        with open(filename, 'w') as fd:
            fd.write(data)

if(len(sys.argv))!=5:
    raise ValueError(" Scripts needs 4 cluster partitions as input, put 0 if not less desired")
else:
    n_part = 0
    partitions = ["0","0","0","0"]

    if sys.argv[1]!="0":
        partitions[0] = sys.argv[1]
        n_part +=1
    if sys.argv[2]!="0":
        partitions[1] = sys.argv[2]
        n_part +=1
    if sys.argv[3]!="0":
        partitions[2] = sys.argv[3]
        n_part +=1
    if sys.argv[4]!="0":
        partitions[3] = sys.argv[4]
        n_part +=1

files_written = 0
code_run_path = "/n/home12/tfranzmeyer/Code/skeletons/examples/"
run_hours = "3"
slurm_path = "/n/home12/tfranzmeyer/slurms_skeletons/"

prefix = "Zebrafinch"

ID_max = 410
refinement_chunksize = 50

error_path = dataIO.OutputDirectory(prefix) + "error_files/"
output_path = dataIO.OutputDirectory(prefix) + "output_files/"
template = template.replace('{RUNCODEDIRECTORY}', code_run_path)
template = template.replace('{HOURS}', run_hours)
memory = str(39000)
# sizex * sizey * sizez * 3 * 8 bytes

start_blocks = dataIO.StartBlocks(prefix)
n_blocks = dataIO.NBlocks(prefix)


SLURM_OUTPUT_FOLDER = slurm_path

step01folderpath = SLURM_OUTPUT_FOLDER+"step01/"
step02folderpath = SLURM_OUTPUT_FOLDER+"step02/"
step03folderpath = SLURM_OUTPUT_FOLDER+"step03/"
step04folderpath = SLURM_OUTPUT_FOLDER+"step04/"
step05folderpath = SLURM_OUTPUT_FOLDER+"step05/"

makeFolder(step01folderpath)
makeFolder(step02folderpath)
makeFolder(step03folderpath)
makeFolder(step04folderpath)
makeFolder(step05folderpath)

# write slurm for step one
command = "step1.py"
jobname = "S1"

t = template
t = t.replace('{JOBNAME}', jobname)
t = t.replace('{COMMAND}', command)
t = t.replace('{ERROR_PATH}', error_path)
t = t.replace('{OUTPUT_PATH}', output_path)
t = t.replace('{MEMORY}', memory)
t = t.replace('{PARTITION}', partitions[np.random.randint(0,n_part)])
filename = step01folderpath + jobname + ".slurm"
writeFile(filename, t)
files_written += 1

# write slurm for step two
for bz in range(start_blocks[0], start_blocks[0] + n_blocks[0]):

    command = "step2.py" + " " + str(bz)
    jobname = "S2"+"_" +"z"+str(bz).zfill(2)

    t = template
    t = t.replace('{JOBNAME}', jobname)
    t = t.replace('{COMMAND}', command)
    t = t.replace('{ERROR_PATH}', error_path)
    t = t.replace('{OUTPUT_PATH}', output_path)
    t = t.replace('{MEMORY}', memory)
    t = t.replace('{PARTITION}', partitions[np.random.randint(0,n_part)])

    filename = step02folderpath + jobname + ".slurm"
    writeFile(filename, t)
    files_written += 1

# write slurm for step three
for bz in range(start_blocks[0], start_blocks[0] + n_blocks[0]):

    command = "step3.py" + " " + str(bz)
    jobname = "S3"+"_" +"z"+str(bz).zfill(2)

    t = template
    t = t.replace('{JOBNAME}', jobname)
    t = t.replace('{COMMAND}', command)
    t = t.replace('{ERROR_PATH}', error_path)
    t = t.replace('{OUTPUT_PATH}', output_path)
    t = t.replace('{MEMORY}', memory)
    t = t.replace('{PARTITION}', partitions[np.random.randint(0,n_part)])

    filename = step03folderpath + jobname + ".slurm"
    writeFile(filename, t)
    files_written += 1

# write slurm for step four
for bz in range(start_blocks[0], start_blocks[0] + n_blocks[0]):
    for by in range(start_blocks[1], start_blocks[1] + n_blocks[1]):
        for bx in range(start_blocks[2], start_blocks[2] + n_blocks[2]):

            command = "step4.py" + " " + str(bz) + " " + str(by) + " " + str(bx)
            jobname = "S4"+"_"+ "z"+str(bz).zfill(2)+"y"+str(by).zfill(2)+"x"+str(bx).zfill(2)

            t = template
            t = t.replace('{JOBNAME}', jobname)
            t = t.replace('{COMMAND}', command)
            t = t.replace('{ERROR_PATH}', error_path)
            t = t.replace('{OUTPUT_PATH}', output_path)
            t = t.replace('{MEMORY}', memory)
            t = t.replace('{PARTITION}', partitions[np.random.randint(0,n_part)])

            filename = step04folderpath + jobname + ".slurm"
            writeFile(filename, t)
            files_written += 1


ID_range = np.arange(1,ID_max)

for ID_start in ID_range[::refinement_chunksize]:

    ID_end = ID_start+refinement_chunksize-1

    command = "step5.py" + " " + str(ID_start) + " " + str(ID_end)
    jobname = "S5" + "_" + str(ID_start) + "_" + str(ID_end)

    t = template
    t = t.replace('{JOBNAME}', jobname)
    t = t.replace('{COMMAND}', command)
    t = t.replace('{ERROR_PATH}', error_path)
    t = t.replace('{OUTPUT_PATH}', output_path)
    t = t.replace('{MEMORY}', str(memory))
    t = t.replace('{PARTITION}', partitions[np.random.randint(0,n_part)])

    filename = step05folderpath + jobname + ".slurm"
    writeFile(filename, t)
    files_written += 1


print ("Files written: " + str(files_written))
