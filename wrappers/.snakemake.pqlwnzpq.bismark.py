
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/kwreczy/miniconda3/lib/python3.5/site-packages/snakemake-3.12.0-py3.5.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05X=\x00\x00\x00/home/kwreczy//tmp/my_output/04_mapped_n_deduped/sampleB1.bamq\x06a}q\x07(X\x01\x00\x00\x00oq\x08h\x06X\x06\x00\x00\x00_namesq\t}q\nh\x08K\x00N\x86q\x0bsubX\x06\x00\x00\x00configq\x0c}q\r(X\x07\x00\x00\x00PATHOUTq\x0eX\x1d\x00\x00\x00/home/kwreczy//tmp/my_output/q\x0fX\n\x00\x00\x00CHROM_INFOq\x10X@\x00\x00\x00/home/kwreczy/repositories/pigx_bsseq/test_dataset/chromInfo.txtq\x11X\x10\x00\x00\x00trim_galore_argsq\x12X\x00\x00\x00\x00q\x13X\x08\x00\x00\x00GTOOLBOXq\x14X\x1c\x00\x00\x00/home/kwreczy/tmp/ggtoolbox/q\x15X\x03\x00\x00\x00LOGq\x16X\x13\x00\x00\x00/home/kwreczy/logs/q\x17X\x05\x00\x00\x00PROGSq\x18}q\x19(X\x08\x00\x00\x00CUTADAPTq\x1aX\x08\x00\x00\x00cutadaptq\x1bX\x07\x00\x00\x00BISMARKq\x1cX\x07\x00\x00\x00bismarkq\x1dX\x0e\x00\x00\x00BISMARK2REPORTq\x1eX\x0e\x00\x00\x00bismark2reportq\x1fX\x1d\x00\x00\x00BISMARK_METHYLATION_EXTRACTORq X\x1d\x00\x00\x00bismark_methylation_extractorq!X\x07\x00\x00\x00BOWTIE2q"X\x07\x00\x00\x00bowtie2q#X\x13\x00\x00\x00DEDUPLICATE_BISMARKq$X\x13\x00\x00\x00deduplicate_bismarkq%X\x08\x00\x00\x00SAMTOOLSq&X\x08\x00\x00\x00samtoolsq\'X\x1a\x00\x00\x00BISMARK_GENOME_PREPARATIONq(X\x1a\x00\x00\x00bismark_genome_preparationq)X\x06\x00\x00\x00FASTQCq*X\x06\x00\x00\x00fastqcq+X\n\x00\x00\x00TRIMGALOREq,X\x0b\x00\x00\x00trim_galoreq-X\x0b\x00\x00\x00DEDUPLICATEq.X\x13\x00\x00\x00deduplicate_bismarkq/uX\x06\x00\x00\x00PATHINq0X6\x00\x00\x00/home/kwreczy/repositories/pigx_bsseq/test_dataset/in/q1X\n\x00\x00\x00GENOMEPATHq2X\x18\x00\x00\x00/home/kwreczy//tmp/ce10/q3X\x0e\x00\x00\x00GENOME_VERSIONq4X\x04\x00\x00\x00ce10q5X\x0b\x00\x00\x00fastqc_argsq6h\x13X\x07\x00\x00\x00SAMPLESq7}q8X\x08\x00\x00\x00sampleB1q9}q:(X\x08\x00\x00\x00SampleIDq;h9X\t\x00\x00\x00fastq_extq<]q=X\x05\x00\x00\x00fq.gzq>aX\x0b\x00\x00\x00Fastqc_argsq?h\x13X\x08\x00\x00\x00ReadTypeq@X\x04\x00\x00\x00WGBSqAX\t\x00\x00\x00TreatmentqBX\x01\x00\x00\x00BqCX\x05\x00\x00\x00fastqqD]qEX\x10\x00\x00\x00pe1.single.fq.gzqFaX\n\x00\x00\x00fastq_nameqG]qHX\n\x00\x00\x00pe1.singleqIaX\x10\x00\x00\x00Trim_galore_argsqJh\x13X\x0b\x00\x00\x00Bismark.argqKh\x13usX\x0c\x00\x00\x00bismark_argsqLX\x0b\x00\x00\x00 -N 1 -L 2 qMX\n\x00\x00\x00NUMTHREADSqNX\x01\x00\x00\x002qOuX\x05\x00\x00\x00inputqPcsnakemake.io\nInputFiles\nqQ)\x81qRX@\x00\x00\x00/home/kwreczy//tmp/my_output/02_trimmed/pe1.single_trimmed.fq.gzqSa}qT(h\t}qUX\x06\x00\x00\x00infileqVK\x00K\x01\x86qWshVcsnakemake.io\nNamedlist\nqX)\x81qYhSa}qZh\t}q[sbubX\t\x00\x00\x00wildcardsq\\csnakemake.io\nWildcards\nq])\x81q^X\x08\x00\x00\x00sampleB1q_a}q`(h\t}qaX\x06\x00\x00\x00sampleqbK\x00N\x86qcsX\x06\x00\x00\x00sampleqdh_ubX\x06\x00\x00\x00paramsqecsnakemake.io\nParams\nqf)\x81qg(hMX>\x00\x00\x00--output_dir /home/kwreczy//tmp/my_output/04_mapped_n_deduped/qhX(\x00\x00\x00--temp_dir /home/kwreczy//tmp/my_output/qiX\n\x00\x00\x00--bowtie2 qjX,\x00\x00\x00--path_to_bowtie /home/kwreczy/tmp/ggtoolboxqkX#\x00\x00\x00/home/kwreczy/tmp/ggtoolbox/bismarkqlh3X+\x00\x00\x00--samtools_path /home/kwreczy/tmp/ggtoolboxqme}qn(X\r\x00\x00\x00genome_folderqoh3X\x05\x00\x00\x00extraqphMX\x07\x00\x00\x00sampathqqhmX\n\x00\x00\x00useBowtie2qrhjh\t}qs(hoK\x06N\x86qthpK\x00N\x86quhqK\x07N\x86qvX\x0b\x00\x00\x00bowtie2pathqwK\x04N\x86qxhrK\x03N\x86qyX\x07\x00\x00\x00tempdirqzK\x02N\x86q{X\x06\x00\x00\x00outdirq|K\x01N\x86q}X\x07\x00\x00\x00programq~K\x05N\x86q\x7fuhzhihwhkh|hhh~hlubX\x07\x00\x00\x00threadsq\x80K\x01X\t\x00\x00\x00resourcesq\x81csnakemake.io\nResources\nq\x82)\x81q\x83(K\x01K\x01e}q\x84(X\x06\x00\x00\x00_nodesq\x85K\x01h\t}q\x86(X\x06\x00\x00\x00_coresq\x87K\x00N\x86q\x88h\x85K\x01N\x86q\x89uh\x87K\x01ubX\x03\x00\x00\x00logq\x8acsnakemake.io\nLog\nq\x8b)\x81q\x8cXM\x00\x00\x00/home/kwreczy//tmp/my_output/04_mapped_n_deduped/sampleB1_bismark_mapping.logq\x8da}q\x8e(h\t}q\x8fh\x8aK\x00N\x86q\x90sh\x8ah\x8dubX\x04\x00\x00\x00ruleq\x91X\n\x00\x00\x00bismark_seq\x92ub.')
######## Original script #########
__author__ = "Katarzyna Wreczycka"
__copyright__ = "Copyright 2017, KW"


from snakemake.shell import shell



bismark_args = " ".join(list(snakemake.params))
print("bismark_args")
print(bismark_args)


log = snakemake.log_fmt_shell(stdout=True, stderr=True)

n = len(snakemake.input.infile)
assert n == 1 or n == 2, "input->sample must have 1 (single-end) or 2 (paired-end) elements."

if n == 1:
    reads = "{}".format(*snakemake.input.infile)
else:
    reads = "-1 {} -2 {}".format(*snakemake.infile)


command = bismark_args + reads + " <2 " + log

print('command')
print(command)

shell(command)


# shell(
#     "(bowtie2 --threads {snakemake.threads} {snakemake.params.extra} "
#     "-x {snakemake.params.index} {reads} "
#     "| samtools view -Sbh -o {snakemake.output[0]} -) {log}")
#     
# 
# print('--------bismark-------infile---------')
# print(input.infile)
# 
#         cmd1 = " ".join(
#         [BISMARK,
#          params.extra,
#          params.outdir,
#          params.tempdir,
#          params.sampath,
#          params.bowtie2path,
#          params.useBowtie2,
#          params.genome_folder,
#          '--multicore', str(threads)
#          ])
# 
#         if len(input.infile)==1:
#           cmd2 = " ".join(
#             [" ", input.infile[0],
#              '2>',log.log
#              ])
#         elif len(input.infile)==2:
#           cmd2 = " ".join(
#             [" ", " -1 ", input.infile[0],
#              " -2 ", input.infile[1],
#              '2>',log.log
#              ])
# 
#         #real_list_files_bismark(input.infile)
#         command = cmd1 + cmd2 + "; touch " +  output.o
#         print("----------commad bismark--------")
#         print(command)
#         shell(command)
# def get_fq_after_trimming(wc):
# 
#    print('---------------wc------get_fq_after_trimming-----2222')
#    print(list(wc))
# 
#    samps = config['SAMPLES'][wc.sample]['fastq_name']
# 
#    if type(samps) is str:
#         samps = [samps]
# 
#    if(len(samps)==2):
#      d=[os.path.join(DIR_trimmed, samps[0]+'_val_1.fq.gz'), os.path.join(DIR_trimmed, samps[1]+'_val_2.fq.gz')]
#    else:
#      d=[os.path.join(DIR_trimmed, samps[0]+'_trimmed.fq.gz')]
