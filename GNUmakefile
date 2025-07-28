# -*- GNUMakefile -*-
# PiGx Developer GNUMakefile
#
# Copyright © 2021 Aexander Blume <alexander.blume@mdc-berlin.de> 
#
# This file is provided to help the developers of the PiGx Pipelines.  
# Change History
# 08/06/2021 Alexander Blume    Update help message.
#                               Add release and sign.
#
# 11/05/2024 Alexander Blume    Convert bash script into GNUMakefile.
# 28/07/2025 Alexander Blume    Add more targets and improve help.
#
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

.PHONY: all clean test release sign

# Delegate any unspecified target to the original Makefile
%:
	@$(MAKE) -f Makefile $@

all: 
	make init
	make build


##? help: Show usage and available commands
# this target finds helpstrings beginning with '##' in the Makefile 
help:
	@echo "Usage: $(MAKE) [subcommand] [-v]"
	@echo
	@echo "Available subcommands are:"
	@grep -h "##" $(MAKEFILE_LIST) | grep -v grep  | sed -e 's/##//' | sed -e 's/^ /    /' | sed -e 's/: /\t\t/'
	@exit 1


##? init: Initialize the project and build the executable
init: pigx-bsseq
	@echo "Initializing the PiGx project..."
	@echo "Fetching submodules..."
	git submodule update --init --recursive
	
## build: Build the executable
build: pigx-bsseq
	./bootstrap.sh
	./configure

## build-guix: Build the executable in a pure environment using guix shell
#  --pure:        unset existing environment variables
#  -D:            include the development inputs of the next package
#  -f guix.scm:   use the given file as the build manifest.
build-guix: pigx-bsseq
	guix shell --pure -D -f guix.scm -- ./bootstrap.sh
	guix shell --pure -D -f guix.scm -- ./configure PYTHONPATH='${GUIX_PYTHONPATH}'

## test: Run tests with sample configuration
test:
	PIGX_UNINSTALLED=1 ./pigx-bsseq -s tests/settings.yaml tests/sample_sheet.csv

## dry: Run a dry-run of the pipeline
dry:
	PIGX_UNINSTALLED=1 ./pigx-bsseq -s tests/settings.yaml tests/sample_sheet.csv -n --force --printshellcmds

## tarball: Create a distribution tarball
tarball:
	make distcheck

VERSION = $(shell cat VERSION)
## sign: Sign tag and release (requires gpg)
sign:
	git tag --sign v$(VERSION)
	gpg --detach-sign pigx_*-$(VERSION).tar.gz

## upload-release: Upload the release to GitHub (requires gh)
upload-release:
	git push --tags
	gh release create v$(VERSION) $(shell ls pigx_*-$(VERSION).tar.gz{,.sig}) --draft 

## release: Create a release (requires gpg and gh)
release: 
	make tarball
	make sign
	make upload-release

format:
	@echo "Formatting code with snakefmt..."
	snakefmt snakefile.py rules/*.py
	@echo "Formatting code with air..."	
	air format scripts/
	@echo "Formatting complete."

check-format:
	@echo "Formatting code with snakefmt..."
	snakefmt --check snakefile.py rules/*.py
	@echo "Formatting code with air..."	
	air format --check scripts/
	@echo "Formatting complete."

# Execute the specified command
$(eval $(ARGS) : ; @true)

