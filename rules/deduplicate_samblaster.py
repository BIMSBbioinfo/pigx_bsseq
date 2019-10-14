# ==========================================================================================
# Deduplication by samblaster 
#
# https://github.com/GregoryFaust/samblaster
#
rule samblaster_markdup_sort:
    input:
        DIR_mapped+"{sample}.bwameth.bam"
    output:
        bam = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam",
        index = DIR_sorted+"{sample}.bwameth.sorted.markdup.bam.bai"
    params:
        threads=config['execution']['rules']['samtools_sort_bam']['threads'],
        memory=config['execution']['rules']['samtools_sort_bam']['memory'],
        tmpdir=DIR_sorted+"{sample}/"
    log:
        DIR_sorted+"{sample}_markdups.log"
    shell:
        nice("samtools", 
        ["view -h {input}"," | ", 
        tool("samblaster"),toolArgs("samblaster"),"2> {log}","|",
        tool("samtools"),"sort","-T={params.tmpdir}",
         "-o {output.bam}", "-@ {params.threads}", 
         "-m {params.memory}", "-l 9","2> {log}",";",
         tool("samtools"),"index {output.bam}"],("{log}"))
