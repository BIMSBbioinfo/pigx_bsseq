numjobs=12
RUNTIME="32:00:00" 
MEM="8G"

#===== PATHS TO BE LINKED TO SYMBOLICALLY FOR THE SNAKEMAKE SCRIPT ======
PATH_IN="/scratch/AG_Akalin/bosberg/SRA_files/"
PATH_OUT="/scratch/AG_Akalin/bosberg/SRA_files/SRA_deconv_pipe_run_out/"
PATH_GENOME="/scratch/AG_Akalin/Base/Genomes/hg19/"
PATH_DEPS="/home/bosberg/.guix-profile/bin/"

#===== PATHS SPECIFIC TO THE R-SCRIPT FOR DECONVOLUTION:
R_WORKDIR="/home/bosberg/bs/pigx_bsseq/"
R_SIGMAT_PATH="/home/bosberg/bs/pigx_out/HK_Sun_data/"
R_PATH_DATA="/scratch/AG_Akalin/bosberg/SRA_files/SRA_deconv_pipe_run_out/06_deduped/"
R_PATHOUT="/scratch/AG_Akalin/bosberg/SRA_files/SRA_deconv_pipe_run_out/07_deconv/"

# ================================================================


if [ -d "path_links" ]; then
rm -r path_links
fi
#--- start clean ---

mkdir path_links
#--- set symbolic links to the various paths.
ln -s ${PATH_IN}     path_links/in
ln -s ${PATH_OUT}    path_links/out
ln -s ${PATH_GENOME} path_links/genome_ref
ln -s ${PATH_DEPS}   path_links/dependencies

# ================================================================
# create a unique log file for this run:
i=1
LOG="./nohup_SM_submit_"${i}".log"
while [ -f ${LOG} ]
do
  i=$((i+1))
  LOG="./nohup_SM_submit_"${i}".log"
done
# ================================================================

echo "starting Snakemake session on " $(date) >  ${LOG}
echo "" >>${LOG}

echo ""                                       >> ${LOG}
echo "------ using samples : ------"          >> ${LOG}
grep -i "files" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ from folder : --------"          >> ${LOG}
grep   "PATHIN" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ to folder : ----------"          >> ${LOG}
grep PATHOUT config.json                      >> ${LOG}

echo ""                                       >> ${LOG}

#========================================================================================
#----------  NOW START DEFINING THE VARIABLES NECESSARY FOR THE DCONVOLUTION SCRIPT -----

echo "R_WORKDIR="        ${R_WORKDIR}          > Define_BSvars.R
echo "R_SIGMAT_PATH="    ${R_SIGMAT_PATH}    >> Define_BSvars.R
echo "R_PATH_DATA="      ${R_PATH_DATA}      >> Define_BSvars.R
echo "R_PATHOUT="        ${R_PATHOUT}        >> Define_BSvars.R

#========================================================================================
#----------  NOW START RUNNING SNAKEMAKE:  ----------------------------------------------

echo "------ Commencing Snakemake : ------"   >> ${LOG}

if [ $# -eq 0 ]
  then
    echo "No arguments supplied. Exiting."
    exit
fi

if [ $1 == "dry" ] 
then
snakemake -s BSseq_pipeline.py -n
elif [ $1 == "cluster" ] 
then
snakemake -s BSseq_pipeline.py --jobs ${numjobs}   --cluster "qsub -V -l h_vmem=${MEM} -pe smp 1  -l h_rt=${RUNTIME} -l h_stack=128k"       >> ${LOG} &
elif [ $1 == "single_core" ] 
then
nohup  snakemake -s BSseq_pipeline.py    >> ${LOG} &
else
echo "command line arg not understood. Exiting without submission"
fi

