numjobs=6

WORKDIR="/home/bosberg/bs/pigx_bsseq/"
SIGMAT_PATH="/home/bosberg/bs/pigx_out/HK_Sun_data/"
PATH_DATA="/data/local/bosberg/SRA_cell_files/deconv_pipe_run_out/06_sorted/"
PATHOUT="/data/local/bosberg/SRA_cell_files/deconv_pipe_run_out/07_deconv/"

i=1
LOG="./nohup_SM_submit_"${i}".log"
# ================================================================

while [ -f ${LOG} ]
do
  i=$((i+1))
  LOG="./nohup_SM_submit_"${i}".log"
done

echo "starting Snakemake session on " $(date) >  ${LOG}
echo "" >>${LOG}

echo ""                                       >> ${LOG}
echo "------ using samples : ------"          >>${LOG}
grep -i "files" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ from folder : --------"          >>${LOG}
grep   "PATHIN" config.json                   >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ to folder : ----------"          >>${LOG}
grep PATHOUT config.json                      >> ${LOG}

echo ""                                       >> ${LOG}
echo "------ Commencing Snakemake : ------"   >>${LOG}

#========================================================================================
#----------  NOW START DEFINING THE VARIABLES NECESSARY FOR THE DCONVOLUTION SCRIPT -----

echo "WORKDIR="        ${WORKDIR}          > Define_BSvars.R 
echo "SIGMAT_PATH="    ${SIGMAT_PATH}    >> Define_BSvars.R   
echo "PATH_DATA="      ${PATH_DATA}      >> Define_BSvars.R   
echo "PATHOUT="        ${PATHOUT}        >> Define_BSvars.R   

#========================================================================================
#----------  NOW START RUNNING SNAKEMAKE:  ----------------------------------------------

nohup  snakemake -s BSseq_pipeline.py --jobs ${numjobs}        >> ${LOG} &

