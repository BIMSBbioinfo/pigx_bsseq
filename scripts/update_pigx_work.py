# PiGx BSseq Pipeline.
#
# Copyright © 2022 Alexander blume <alexander.blume@mdc-berlin.de>
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

import os
from os import path
from func_defs import bail
import argparse
import shutil

# --------------------------------------
# copy and link files to pigx-work directory
# --------------------------------------
# Create symbolic links to the inputs and reference genome

# Create links within the output folder that point directly to the
# reference genome, as well as to each sample input file so that it's
# clear where the source data came from.

# N.B. Any previously existing links will be kept in place, and no
# warning will be issued if this is the case.


def makelink(src, target):
    if not path.isfile(src):
        bail("Refusing to link non-existent file %s" % src)
    elif not path.isdir(path.dirname(target)):
        bail(
            "%s or subdirectory does not exist for linking %s"
            % path.dirname(target),
            target,
        )
    else:
        try:
            os.symlink(src, target)
        except FileExistsError:
            pass


def update_pigx_work(config):
    os.makedirs(
        path.join(config["locations"]["output-dir"], "pigx_work/input"), exist_ok=True
    )

    # copy documentation file to the output work directory.
    if os.getenv("PIGX_UNINSTALLED"):
        contents = path.join(config["locations"]["pkgdatadir"], "etc/CONTENTS.txt")
    else:
        contents = path.join(config["locations"]["pkgdatadir"], "CONTENTS.txt")
    shutil.copyfile(
        contents, path.join(config["locations"]["output-dir"], "pigx_work/CONTENTS.txt")
    )

    # Link the reference genome
    try:
        os.symlink(
            os.path.dirname(config["locations"]["genome-fasta"]),
            path.join(config["locations"]["output-dir"], "pigx_work/refGenome"),
        )
    except FileExistsError:
        pass

    # Create file links
    for sample in config["SAMPLES"]:
        flist = config["SAMPLES"][sample]["files"]
        single_end = len(flist) == 1

        for idx, f in enumerate(flist):
            if not f.endswith(".gz"):
                # FIXME: Future versions should handle unzipped .fq or .bz2.
                bail("Input files must be gzipped: %s." % f)

            tag = "" if single_end else "_" + str(idx + 1)
            linkname = config["SAMPLES"][sample]["SampleID"] + tag + ".fq.gz"
            makelink(
                path.join(config["locations"]["input-dir"], f),
                path.join(
                    config["locations"]["output-dir"], "pigx_work/input/", linkname
                ),
            )
