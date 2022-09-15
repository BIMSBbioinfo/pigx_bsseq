# PiGx BSseq Pipeline.
#
# Copyright 2022 Alexander Blume <alexander.blume@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ==========================================================================================
# Deduplication with picard:
#
# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates


rule picard_MarkDuplicates:
    input:
        DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
    output:
        outfile=DIR_sorted + "{sample}.bwameth.sorted.picard_MarkDuplicates.bam",
        metrics=DIR_sorted + "{sample}.bwameth.sorted.picard_MarkDuplicates.metrics.txt",
    log:
        DIR_sorted + "{sample}.bwameth.sorted.picard_MarkDuplicates.log",
    message:
        fmt("Marking duplicate reads using Picard for sample {wildcards.sample}")
    shell:
        nice(
            "java",
            [
                "-jar",
                tool("picard"),
                "MarkDuplicates",
                "I={input}",
                "O={output.outfile}",
                "M={output.metrics}",
                "2>> {log}",
            ],
            "{log}",
        )
