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


command = bismark_args + reads + " <2 " + log_fmt_shell

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
