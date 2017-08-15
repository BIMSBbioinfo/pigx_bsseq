
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/kwreczy/miniconda3/lib/python3.5/site-packages/snakemake-3.12.0-py3.5.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05(X+\x00\x00\x00--samtools_path /home/kwreczy/tmp/ggtoolboxq\x06X,\x00\x00\x00--path_to_bowtie /home/kwreczy/tmp/ggtoolboxq\x07X\x0b\x00\x00\x00 -N 1 -L 2 q\x08X>\x00\x00\x00--output_dir /home/kwreczy//tmp/my_output/04_mapped_n_deduped/q\tX\x18\x00\x00\x00/home/kwreczy//tmp/ce10/q\nX#\x00\x00\x00/home/kwreczy/tmp/ggtoolbox/bismarkq\x0bX(\x00\x00\x00--temp_dir /home/kwreczy//tmp/my_output/q\x0cX\n\x00\x00\x00--bowtie2 q\re}q\x0e(X\x07\x00\x00\x00sampathq\x0fh\x06X\x06\x00\x00\x00_namesq\x10}q\x11(h\x0fK\x00N\x86q\x12X\x0b\x00\x00\x00bowtie2pathq\x13K\x01N\x86q\x14X\x05\x00\x00\x00extraq\x15K\x02N\x86q\x16X\r\x00\x00\x00genome_folderq\x17K\x04N\x86q\x18X\x06\x00\x00\x00outdirq\x19K\x03N\x86q\x1aX\x07\x00\x00\x00programq\x1bK\x05N\x86q\x1cX\x07\x00\x00\x00tempdirq\x1dK\x06N\x86q\x1eX\n\x00\x00\x00useBowtie2q\x1fK\x07N\x86q uh\x13h\x07h\x15h\x08h\x17h\nh\x19h\th\x1bh\x0bh\x1dh\x0ch\x1fh\rubX\x05\x00\x00\x00inputq!csnakemake.io\nInputFiles\nq")\x81q#X@\x00\x00\x00/home/kwreczy//tmp/my_output/02_trimmed/pe1.single_trimmed.fq.gzq$a}q%(X\x06\x00\x00\x00infileq&csnakemake.io\nNamedlist\nq\')\x81q(h$a}q)h\x10}q*sbh\x10}q+h&K\x00K\x01\x86q,subX\x07\x00\x00\x00threadsq-K\x01X\x06\x00\x00\x00outputq.csnakemake.io\nOutputFiles\nq/)\x81q0X=\x00\x00\x00/home/kwreczy//tmp/my_output/04_mapped_n_deduped/sampleB1.bamq1a}q2(X\x01\x00\x00\x00oq3h1h\x10}q4h3K\x00N\x86q5subX\x03\x00\x00\x00logq6csnakemake.io\nLog\nq7)\x81q8XM\x00\x00\x00/home/kwreczy//tmp/my_output/04_mapped_n_deduped/sampleB1_bismark_mapping.logq9a}q:(h\x10}q;h6K\x00N\x86q<sh6h9ubX\x06\x00\x00\x00configq=}q>(X\x03\x00\x00\x00LOGq?X\x13\x00\x00\x00/home/kwreczy/logs/q@X\x08\x00\x00\x00GTOOLBOXqAX\x1c\x00\x00\x00/home/kwreczy/tmp/ggtoolbox/qBX\n\x00\x00\x00NUMTHREADSqCX\x01\x00\x00\x002qDX\x10\x00\x00\x00trim_galore_argsqEX\x00\x00\x00\x00qFX\x05\x00\x00\x00PROGSqG}qH(X\x07\x00\x00\x00BOWTIE2qIX\x07\x00\x00\x00bowtie2qJX\x0b\x00\x00\x00DEDUPLICATEqKX\x13\x00\x00\x00deduplicate_bismarkqLX\x1a\x00\x00\x00BISMARK_GENOME_PREPARATIONqMX\x1a\x00\x00\x00bismark_genome_preparationqNX\x08\x00\x00\x00SAMTOOLSqOX\x08\x00\x00\x00samtoolsqPX\x0e\x00\x00\x00BISMARK2REPORTqQX\x0e\x00\x00\x00bismark2reportqRX\n\x00\x00\x00TRIMGALOREqSX\x0b\x00\x00\x00trim_galoreqTX\x1d\x00\x00\x00BISMARK_METHYLATION_EXTRACTORqUX\x1d\x00\x00\x00bismark_methylation_extractorqVX\x06\x00\x00\x00FASTQCqWX\x06\x00\x00\x00fastqcqXX\x13\x00\x00\x00DEDUPLICATE_BISMARKqYX\x13\x00\x00\x00deduplicate_bismarkqZX\x07\x00\x00\x00BISMARKq[X\x07\x00\x00\x00bismarkq\\X\x08\x00\x00\x00CUTADAPTq]X\x08\x00\x00\x00cutadaptq^uX\x07\x00\x00\x00PATHOUTq_X\x1d\x00\x00\x00/home/kwreczy//tmp/my_output/q`X\x07\x00\x00\x00SAMPLESqa}qbX\x08\x00\x00\x00sampleB1qc}qd(X\x05\x00\x00\x00fastqqe]qfX\x10\x00\x00\x00pe1.single.fq.gzqgaX\t\x00\x00\x00fastq_extqh]qiX\x05\x00\x00\x00fq.gzqjaX\t\x00\x00\x00TreatmentqkX\x01\x00\x00\x00BqlX\n\x00\x00\x00fastq_nameqm]qnX\n\x00\x00\x00pe1.singleqoaX\x0b\x00\x00\x00Fastqc_argsqphFX\x08\x00\x00\x00ReadTypeqqX\x04\x00\x00\x00WGBSqrX\x0b\x00\x00\x00Bismark.argqshFX\x08\x00\x00\x00SampleIDqthcX\x10\x00\x00\x00Trim_galore_argsquhFusX\n\x00\x00\x00GENOMEPATHqvh\nX\x0e\x00\x00\x00GENOME_VERSIONqwX\x04\x00\x00\x00ce10qxX\x0b\x00\x00\x00fastqc_argsqyhFX\x0c\x00\x00\x00bismark_argsqzh\x08X\x06\x00\x00\x00PATHINq{X6\x00\x00\x00/home/kwreczy/repositories/pigx_bsseq/test_dataset/in/q|X\n\x00\x00\x00CHROM_INFOq}X@\x00\x00\x00/home/kwreczy/repositories/pigx_bsseq/test_dataset/chromInfo.txtq~uX\t\x00\x00\x00resourcesq\x7fcsnakemake.io\nResources\nq\x80)\x81q\x81(K\x01K\x01e}q\x82(X\x06\x00\x00\x00_coresq\x83K\x01h\x10}q\x84(h\x83K\x00N\x86q\x85X\x06\x00\x00\x00_nodesq\x86K\x01N\x86q\x87uh\x86K\x01ubX\t\x00\x00\x00wildcardsq\x88csnakemake.io\nWildcards\nq\x89)\x81q\x8aX\x08\x00\x00\x00sampleB1q\x8ba}q\x8c(X\x06\x00\x00\x00sampleq\x8dh\x8bh\x10}q\x8eX\x06\x00\x00\x00sampleq\x8fK\x00N\x86q\x90subX\x04\x00\x00\x00ruleq\x91X\n\x00\x00\x00bismark_seq\x92ub.')
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


command = bismark_args + reads + " <2 {log}"

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
