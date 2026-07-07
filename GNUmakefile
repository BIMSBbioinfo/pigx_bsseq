# -*- GNUMakefile -*-
# PiGx Developer GNUMakefile
#
# Copyright © 2021-2026 Alexander Blume <alexander.blume@mdc-berlin.de>
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
PIPELINE_TEMPLATES := \
	etc/settings.yaml.in \
	report_templates/index.Rmd.in \
	report_templates/diffmeth.Rmd.in

BUILD_TARGET := $(if $(GUIX_PYTHONPATH),build-guix,$(if $(CONDA_PREFIX),build-conda,build))
BUILD_DEPS := $(PIPELINE_RUNNER) $(subst .in,,$(PIPELINE_TEMPLATES))
DEV_ENV_TARGET := $(if $(GUIX_PYTHONPATH),dev-guix,dev-conda)

$(if $(and $(GUIX_PYTHONPATH),$(CONDA_PREFIX)),$(warning Both GUIX_PYTHONPATH and CONDA_PREFIX are set; using Guix environment))

# Keep generated config outputs synchronized with their sources.
CONFIGURE_DEPS := \
	configure.ac \
	aclocal.m4 \
	Makefile.am

ENV_FILES := guix.scm requirements.yaml

CONDA_ENV := $(CURDIR)/.conda
CONDA_LOCK := $(CONDA_ENV)/requirements.lock

VERSION_FILE := VERSION
VERSION_TAG := $(shell cat $(VERSION_FILE) 2>/dev/null || echo "0")
TARBALL := pigx_$(PIPELINE)-$(VERSION_TAG).tar.gz
SIGNED_TARBALL := $(TARBALL).sig

# keep list of version controlled source files to add them at the bottom as
# explicit rules so that the catch-all pattern below cannot match them.
SRC_FILES := \
	$(PIGX_RUNNER) \
	$(PIPELINE_TEMPLATES) \
	$(CONFIGURE_DEPS) \
	$(ENV_FILES) \
	$(VERSION_FILE)

all: $(BUILD_TARGET)

.PHONY: help
## help: Show usage and available commands
# this target finds helpstrings beginning with '##' in the Makefile
help:
	@echo "Usage: $(MAKE) [subcommand] [-v]"
	@echo
	@echo "Available subcommands are:\n"
	@grep -E '^[ \t]*## [a-zA-Z_-]+:' $(MAKEFILE_LIST) | sed -E 's/^[ \t]*## //' | awk -F': ' '{printf "%-24s %s\n", $$1, $$2}'
	@echo
	@echo "More available via tab-completion. Try 'make <TAB><TAB>' to see all available targets."
	@exit 0

# ---------------------------------------------------------------------------
# == Bootstrap and configure helpers ==
# ---------------------------------------------------------------------------


pigx-common/common/m4: 
	@if git submodule status | grep -q -E '^[-+]' ; then \
		echo "INFO: Need to reinitialize git submodules"; \
		git submodule update --init; \
	fi

## init-submodules: Initialize the submodules
init-submodules: pigx-common/common/m4
	@echo "INFO: Submodules are initialized."

configure: $(CONFIGURE_DEPS) pigx-common/common/m4
	@[ -f build-aux/install-sh ] || ./bootstrap.sh

$(PIGX_RUNNER): | pigx-common/common/m4

$(PIPELINE_RUNNER): $(PIGX_RUNNER) | configure
	./configure

# ---------------------------------------------------------------------------
# == Build targets ==
# ---------------------------------------------------------------------------

.PHONY: build build-guix build-conda

## build: Build the executable using the current environment
build: $(PIPELINE_RUNNER)

## build-guix: Build in a pure Guix environment
build-guix: require-guix guix.scm configure $(PIGX_RUNNER)
	@echo "Building in guix environment: $(GUIX_PYTHONPATH)"
	@guix shell --pure -D -f guix.scm -- sh -c '\
		export PYTHONPATH="$$GUIX_PYTHONPATH"; \
		./configure; \
	'

## build-conda: Build in a conda/micromamba environment
build-conda: $(CONDA_LOCK) configure $(PIGX_RUNNER)
	@echo "Building in conda environment: $(CONDA_ENV)"
	@$(MAMBA_EXE) run -p $(CONDA_ENV) --clean-env bash -c '\
		export R_LIBS_SITE="$${CONDA_PREFIX}/lib/R/library"; \
		export PYTHONPATH="$$(python -c '\''import sysconfig; print(sysconfig.get_paths()["purelib"])'\'')"; \
		./configure; \
	'

# ---------------------------------------------------------------------------
# == Install targets ==
# ---------------------------------------------------------------------------

.PHONY: install-local

## build: Build the executable using the current environment
install-local:
	./bootstrap.sh
	./configure --prefix=$$HOME/.local
	$(MAKE) install

# ---------------------------------------------------------------------------
# == Test and pipeline run targets ==
# ---------------------------------------------------------------------------

.PHONY: test-local test-unlock test-cluster test-slurm-dry test-qsub-dry test-dry

## test-local: Run the sample pipeline locally
test-local: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings.yaml tests/sample_sheet.csv

## test-unlock: Unlock the sample pipeline run
test-unlock: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings.yaml tests/sample_sheet.csv --unlock

## test-cluster: Dry-run the sample pipeline with Slurm cluster settings
test-cluster: test-slurm-dry

## test-slurm-dry: Dry-run the sample pipeline with Slurm cluster settings
test-slurm-dry: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings_slurm.yaml tests/sample_sheet.csv -n --force --printshellcmds

## test-qsub-dry: Dry-run the sample pipeline with qsub cluster settings
test-qsub-dry: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings_qsub.yaml tests/sample_sheet.csv -n --force --printshellcmds

## test-gpu: Run the sample pipeline with GPU settings
test-gpu: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings_gpu.yaml tests/sample_sheet.csv --force --printshellcmds

## test-gpu-dry: Dry-run the sample pipeline with GPU settings
test-gpu-dry: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings_gpu.yaml tests/sample_sheet.csv -n --force --printshellcmds

## test-slurm-gpu-dry: Dry-run the sample pipeline with Slurm and GPU settings
test-slurm-gpu-dry: $(PIPELINE_RUNNER)
	PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings_slurm_gpu.yaml tests/sample_sheet.csv -n --force --printshellcmds

## test-dry: Run a dry-run of the pipeline
test-dry: $(PIPELINE_RUNNER)
	@PIGX_UNINSTALLED=1 ./$(PIPELINE_RUNNER) -s tests/settings.yaml tests/sample_sheet.csv -n --force --printshellcmds > dry-run.log && echo "Dry-run successful. No errors detected." || (echo "Dry-run failed. Check dry-run.log for details." && exit 1)

TEST_CONFIG_FILE := config.json
$(TEST_CONFIG_FILE): | test-dry

# ---------------------------------------------------------------------------
# == Release targets ==
# ---------------------------------------------------------------------------

.PHONY: release-dist sign tag upload-release release check-version

## check-version: Ensure VERSION file is committed and matches the working tree
check-version:
	@git ls-files --error-unmatch $(VERSION_FILE) >/dev/null 2>&1 || \
		(echo "Error: $(VERSION_FILE) is not tracked by git." >&2; exit 1)
	@git diff --quiet -- $(VERSION_FILE) || \
		(echo "Error: $(VERSION_FILE) has uncommitted changes. Commit it before tagging." >&2; exit 1)
	@git diff --quiet --cached -- $(VERSION_FILE) || \
		(echo "Error: $(VERSION_FILE) has staged but uncommitted changes. Commit it before tagging." >&2; exit 1)
	@committed_version=$$(git show HEAD:$(VERSION_FILE) 2>/dev/null); \
	if [ "$$committed_version" != "$(VERSION_TAG)" ]; then \
		echo "Error: committed $(VERSION_FILE) ($$committed_version) does not match working copy ($(VERSION_TAG))." >&2; \
		exit 1; \
	fi

$(TARBALL): $(VERSION_FILE) Makefile
	$(MAKE) -f Makefile distcheck

## release-dist: Create a distribution tarball for release
release-dist: $(TARBALL)

$(SIGNED_TARBALL): $(TARBALL)
	gpg --detach-sign $(TARBALL)

## sign: Sign the release tarball
sign: $(SIGNED_TARBALL)

## tag: Tag the release in git
tag: sign
	git tag --sign v$(VERSION_TAG)

## upload-release: Upload release artifacts to GitHub
upload-release: tag
	git push --tags
	gh release create v$(VERSION_TAG) $(TARBALL) $(SIGNED_TARBALL) --draft

## release: Create a complete release
release: release-dist sign upload-release

# ---------------------------------------------------------------------------
# == Lint and formatting helpers ==
# ---------------------------------------------------------------------------

.PHONY: lint format format-check

## lint-check: Lint the snakefile using Snakemake (non-fatal report)
lint-check: require-snakemake $(TEST_CONFIG_FILE)
	snakemake -s snakefile.py --configfile $(TEST_CONFIG_FILE) --lint || { echo "Linting failed. Please list the issues with 'make lint'." && exit 1; }

## lint: Lint the snakefile using Snakemake
lint: require-snakemake $(TEST_CONFIG_FILE)
	snakemake -s snakefile.py --configfile $(TEST_CONFIG_FILE) --lint

## format: Format rules and scripts
format: require-snakefmt require-air
	snakefmt snakefile.py rules/*.py;
	air format scripts/

## format-check: Check formatting of rules and scripts
format-check: require-snakefmt require-air
	snakefmt --check snakefile.py rules/*.py;\
	air format --check scripts/ || { echo "Formatting check failed. Please run 'make format' to fix formatting issues." && exit 1; }

# ---------------------------------------------------------------------------
# == Development environment helpers ==
# ---------------------------------------------------------------------------

.PHONY: dev dev-guix dev-conda

## dev: Enter the preferred development environment
dev:
	$(MAKE) $(DEV_ENV_TARGET)

## dev-guix: Enter the development environment with all tools available
dev-guix: require-guix
	@echo "Entering development environment with Guix..."
	guix shell -D -f guix.scm

$(CONDA_LOCK): requirements.yaml
	@echo "Creating/updating conda environment..."
	@if test -d "$(CONDA_ENV)/conda-meta"; then \
		$(MAMBA_EXE) install -y -p $(CONDA_ENV) -f requirements.yaml; \
	else \
		$(MAMBA_EXE) create -y -p $(CONDA_ENV) -f requirements.yaml; \
	fi
	@echo "Generating locked requirements from requirements.yaml..."
	@$(MAMBA_EXE) env export -p $(CONDA_ENV) --explicit > $(CONDA_LOCK)

## dev-conda: Enter the development environment using micromamba
dev-conda: $(CONDA_LOCK)
	@echo "Entering conda environment..."
	@$(MAMBA_EXE) run -p $(CONDA_ENV) bash -c '\
		export R_LIBS_SITE="$${CONDA_PREFIX}/lib/R/library"; \
		export PYTHONPATH="$$(python -c '\''import sysconfig; print(sysconfig.get_paths()["purelib"])'\'')"; \
		export HOME="$${HOME}"; \
		bash \
	'

# ---------------------------------------------------------------------------
# == Cleanup and fallback rules ==
# ---------------------------------------------------------------------------

.PHONY: clean
## clean: Delete almost everything that can be reconstructed with the Makefile.
clean: Makefile
	@$(MAKE) -f Makefile maintainer-clean

require-%:
	@command -v $* >/dev/null 2>&1 || \
		{ echo "Error: '$*' not found. Run 'make dev-guix' or 'make dev-conda' to enter the development environment."; exit 1; }

# The generated Makefile is created by autoconf.
Makefile: | configure
	@echo "Generating Makefile..."
	@./configure -q


# Source files tracked by version control.
$(SRC_FILES): ;

# Delegate any unspecified target to the generated Makefile.
%:
	@if [ -f Makefile ]; then \
		$(MAKE) --no-print-directory -f Makefile $@ || { echo "Error: Target '$@' not found in Makefile. Run 'make help' to see available targets."; exit 1; }; \
	else \
		echo "Error: Makefile not found. Run 'make' in a development environment with all dependencies installed to generate the Makefile."; exit 1; \
	fi

# Execute the specified command
$(eval $(ARGS) : ; @true)
