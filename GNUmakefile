# -*- GNUMakefile -*-
# PiGx Developer GNUMakefile
#
# Copyright © 2021-2026 Aexander Blume <alexander.blume@mdc-berlin.de> 
#
# This file is provided to help the developers of the PiGx Pipelines.  
# Change History
# 08/06/2021 Alexander Blume    Update help message.
#                               Add release and sign.
#
# 11/05/2024 Alexander Blume    Convert bash script into GNUMakefile.
# 28/07/2025 Alexander Blume    Add more targets and improve help.
# 29/07/2025 Alexander Blume    Add more targets and update default target.
#								Add clean target to run maintainer-clean.
# 								Add target to create Makefile.
# 02/06/2026 Alexander Blume    Rename init to fetch-submodules.
# 								Add target for pigx-common/common/pigx-runner.in.
# 								Make build targets depend on pigx-runner
# 								
#
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

PIPELINE := bsseq
PIPELINE_RUNNER := pigx-$(PIPELINE)
PIGX_RUNNER := pigx-common/common/pigx-runner.in
PIPELINE_TEMPLATES := etc/settings.yaml.in report_templates/index.Rmd.in report_templates/diffmeth.Rmd.in

BUILD_TARGET := $(if $(GUIX_PYTHONPATH),build-guix,build)
BUILD_DEPS := $(PIPELINE_RUNNER) $(PIGX_RUNNER) $(subst .in,,$(PIPELINE_TEMPLATES))

.PHONY: all
all: $(BUILD_TARGET)


##? help: Show usage and available commands
# this target finds helpstrings beginning with '##' in the Makefile 
help:
	@echo "Usage: $(MAKE) [subcommand] [-v]"
	@echo
	@echo "Available subcommands are:"
	@grep -E '^[ \t]*## [a-zA-Z_-]+:' $(MAKEFILE_LIST) | sed -E 's/^[ \t]*## //' | sed -E 's/: */\t/' | column -s $$'\t' -t
	@exit 0


define init-submodules
@if git submodule status | grep -q -E '^[-+]' ; then \
	echo "INFO: Need to reinitialize git submodules"; \
	git submodule update --init; \
	fi
	endef

.PHONY: init-submodules
## init-submodules: Initialize the submodules
init-submodules:
	$(init-submodules)

$(PIGX_RUNNER):
	$(init-submodules)

$(PIPELINE_RUNNER): $(PIGX_RUNNER)
	./configure

# -- Autoconf bootstrap dependencies --

# Tell make how to turn any .in file into its configured counterpart
%: %.in config.status
	./config.status --file=$@

config.status: $(CONFIGURE_DEPS)
	./configure


.PHONY: build
## build: Build the executable
build: $(BUILD_DEPS)
	@:

.PHONY: build-guix
## build-guix: Build the executable in a pure environment using guix shell
#  --pure:        unset existing environment variables
#  -D:            include the development inputs of the next package
#  -f guix.scm:   use the given file as the build manifest.
build-guix: $(BUILD_DEPS) guix.scm
	guix shell --pure -D -f guix.scm -- ./bootstrap.sh
	guix shell --pure -D -f guix.scm -- sh -c './configure PYTHONPATH="$$GUIX_PYTHONPATH"'

.PHONY: clean
# https://www.gnu.org/prep/standards/html_node/Standard-Targets.html
## clean: Delete almost everything that can be reconstructed with the Makefile. 
clean:
	$(MAKE) maintainer-clean

.PHONY: test
## test: Run tests with sample configuration
test: $(BUILD_DEPS)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings.yaml tests/sample_sheet.csv

.PHONY: dry
## dry: Run a dry-run of the pipeline
dry: $(BUILD_DEPS)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings.yaml tests/sample_sheet.csv -n --force --printshellcmds

.PHONY: tarball
## tarball: Create a distribution tarball
tarball:
	$(MAKE) distcheck

VERSION := $(shell cat VERSION 2>/dev/null || echo "0")

TARBALL :=  pigx_$(PIPELINE)-$(VERSION).tar.gz
SIGNED_TAR :=  pigx_$(PIPELINE)-$(VERSION).tar.gz.sig

$(TARBALL):
	$(MAKE) tarball

.PHONY: sign
## sign: Sign tag and release (requires gpg)
sign: $(TARBALL)
	git tag --sign v$(VERSION)
	gpg --detach-sign $(TARBALL)

.PHONY: upload-release
## upload-release: Upload the release to GitHub (requires gh)
upload-release: $(TARBALL) $(SIGNED_TAR)
	git push --tags
	gh release create v$(VERSION) $(shell ls $(TARBALL) $(SIGNED_TAR)) --draft 

.PHONY: release
## release: Create a release (requires gpg and gh)
release: 
	$(MAKE) tarball
	$(MAKE) sign
	$(MAKE) upload-release

.PHONY: format
format:
	@echo "Formatting code with snakefmt..."
	snakefmt snakefile.py rules/*.py
	@echo "Formatting code with air..."	
	air format scripts/
	@echo "Formatting complete."

.PHONY: check-format
check-format:
	@echo "Formatting code with snakefmt..."
	snakefmt --check snakefile.py rules/*.py
	@echo "Formatting code with air..."	
	air format --check scripts/
	@echo "Formatting complete."


# The local Makefile is generated by autoconf.
Makefile: $(PIGX_RUNNER)
	./bootstrap.sh
	./configure

# Delegate any unspecified target to the original Makefile (only when it exists)
ifneq ($(wildcard Makefile),)
	%: Makefile
	@$(MAKE) -f Makefile $@
endif


# Execute the specified command
$(eval $(ARGS) : ; @true)

