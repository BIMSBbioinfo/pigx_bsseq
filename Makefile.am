ACLOCAL_AMFLAGS = -I pigx-common/common/m4

nodist_bin_SCRIPTS = \
  pigx-bsseq

dist_pkglibexec_SCRIPTS =			\
  snakefile.py						\
  qsub-template.sh

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =	\
  scripts/generate_report.R			\
  scripts/fetch_procedures.R		\
  scripts/func_defs.py				\
  scripts/update_pigx_work.py				\
  scripts/validate_config.py				\
  scripts/ideoDMC.R					\
  scripts/export_tbx2bw.R         	\
  scripts/methCall.R          		\
  scripts/makeTabix.R          		\
  scripts/methSeg.R           		\
  scripts/methUnite.R          		\
  scripts/methDiff.R				\
  scripts/find_dmrs.R				

pkglibexec_rulesdir = $(pkglibexecdir)/rules
dist_pkglibexec_rules_SCRIPTS =			\
  rules/Align_bwameth_rules.py 			\
  rules/deduplicate_samblaster.py 		\
  rules/preprocessing_methyldackel.py 	\
  rules/mapping_stats.py 				\
  rules/multiqc.py 						\
  rules/export_tabix_bigwig.py 			\
  rules/perform_diffmeth.py

pkgdata_report_templatesdir = $(pkgdatadir)/report_templates
dist_pkgdata_report_templates_SCRIPTS = 			\
  report_templates/index.Rmd						\
  report_templates/diffmeth.Rmd						\
  report_templates/reports.bib

dist_pkgdata_DATA =							\
  etc/sample_sheet.csv.example				\
  etc/settings.yaml							\
  etc/pretty.txt 							\
  etc/CONTENTS.txt            \
  images/Logo_PiGx.png

dist_doc_DATA = 					\
  README.md

EXTRA_DIST =										\
  tests/genome/sample.fasta	\
  tests/genome/refGene.hg19.bed.gz                       \
  tests/genome/cpgIslandExt.hg19.bed.gz                  \
  tests/in/PE_1.fq.gz								\
  tests/in/PE_2.fq.gz								\
  tests/in/SE_techrep1.fq.gz							\
  tests/in/SE_techrep2.fq.gz							\
  tests/sample_sheet.csv								\
  tests/settings.yaml

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)"

TESTS = \
  test.sh											\
  tests/test-check-sample_sheet-header.sh			\
  tests/test-write-access-to-out.sh

EXTRA_DIST += $(TESTS)

EXTRA_DIST += META.in VERSION

clean-local:
	-rm -rf .snakemake
	-rm -rf tests/out
	-rm -rf tests/out/.snakemake

CLEANFILES = $(nodist_bin_SCRIPTS) config.json


# Build docker image with Guix
docker: dist
	guix pack -e '(load "guix.scm")' --with-source=pigx_bsseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils

# For developers: fetch and convert bed files for hg19 assembly.  Requires bedops.
fetch-cpg:
	wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz \
	   | gunzip -c \
	   | awk 'BEGIN{ OFS="\t"; }{ print $$2, $$3, $$4, $$5$$6, substr($$0, index($$0, $$7)); }' \
	   | sort-bed - | gzip > cpgIslandExt.hg19.bed.gz

fetch-ref:
	wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz \
	   | gunzip -c \
	   | awk 'BEGIN{ OFS="\t"; }{ print $$3, $$5, $$6, $$2, $$12, $$4, $$7, $$8; }' \
	   | sort-bed - | gzip > refGene.hg19.bed.gz
