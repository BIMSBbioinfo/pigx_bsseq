# ---------------------------------------------------------------------------- #


def multiqc_files(branch):
    files = files_for_sample(list_files_rawQC)
    files += files_for_sample(list_files_TG)
    files += files_for_sample(list_files_posttrim_QC)
    if branch == "bismark":
        files += files_for_sample(list_files_bismark)
        files += files_for_sample(list_files_methyldackel_mbias_bismark)
    elif branch == "bwameth":
        files += files_for_sample(list_files_bwameth)
        files += files_for_sample(list_files_bwamethMappingStats)
        files += files_for_sample(list_files_markdup)
        files += files_for_sample(list_files_methyldackel_mbias_bwameth)
    files = list(chain.from_iterable(files))
    return files


rule multiqc:
    input:
        files=lambda wc: multiqc_files(wc.branch),
    output:
        os.path.join(
            DIR_final,
            "multiqc",
            "{branch}_multiqc_report.html",
        ),
    params:
        trim_galore=DIR_trimmed,
        raw_qc=DIR_rawqc,
        posttrim_qc=DIR_posttrim_QC,
    log:
        os.path.join(
            DIR_final,
            "multiqc",
            "{branch}_multiqc.log",
        ),
    message:
        fmt("Generating multi-sample QC report for {wildcards.branch} branch.")
    shell:
        nice(
            "multiqc",
            [
                "-f",
                "-n {output}",
                "{input.files}",
                "{params}",
            ],
            "{log}",
        )
