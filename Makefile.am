ACLOCAL_AMFLAGS = -I m4

nodist_bin_SCRIPTS = \
  pigx-rnaseq

dist_pkglibexec_SCRIPTS = \
  snakefile.py			\
  qsub-template.sh

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS =					\
  scripts/runDeseqReport.R							\
  scripts/translate_sample_sheet_for_report.R		\
  scripts/deseqReport.Rmd		\
  scripts/counts_matrix_from_SALMON.R		\
  scripts/validate_input.py	\
  scripts/norm_counts_deseq.R	\
  scripts/count_reads.R	\
  scripts/collate_read_counts.R \
  scripts/collate_deseq_results.R \
  scripts/validate_input_annotation.R

dist_pkgdata_DATA =									\
  etc/sample_sheet.csv.example								\
  etc/settings.yaml									\
  etc/pretty.txt									\
  images/Logo_PiGx.png

dist_doc_DATA = 									\
  README.md

EXTRA_DIST =										\
  VERSION											\
  META.in											\
  tests/sample_sheet.csv							\
  tests/sample_data/sample.gtf						\
  tests/sample_data/reads/HBR_Rep1.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep1.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep1.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep3.read1.fastq.gz	\
  tests/sample_data/reads/HBR_Rep2.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep2.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep2.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep3.read2.fastq.gz	\
  tests/sample_data/reads/UHR_Rep1.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep3.read1.fastq.gz	\
  tests/sample_data/reads/UHR_Rep3.read2.fastq.gz	\
  tests/sample_data/reads/HBR_Rep2.read2.fastq.gz	\
  tests/sample_data/sample.fasta					\
  tests/sample_data/sample.cdna.fasta				\
	tests/sample_data/preprocessed/colData.tsv     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/cmd_info.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/quant.sf     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/quant.genes.sf     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/lib_format_counts.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/logs/salmon_quant.log     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/libParams/flenDist.txt     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/fld.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/expected_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/observed_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/observed_bias_3p.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/obs5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/obs3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/exp5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/exp3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/obs_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/exp_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/meta_info.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep1/aux_info/ambig_info.tsv     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/cmd_info.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/quant.sf     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/quant.genes.sf     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/lib_format_counts.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/logs/salmon_quant.log     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/libParams/flenDist.txt     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/fld.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/expected_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/observed_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/observed_bias_3p.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/obs5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/obs3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/exp5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/exp3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/obs_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/exp_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/meta_info.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep3/aux_info/ambig_info.tsv     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/cmd_info.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/quant.sf     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/quant.genes.sf     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/lib_format_counts.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/logs/salmon_quant.log     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/libParams/flenDist.txt     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/fld.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/expected_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/observed_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/observed_bias_3p.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/obs5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/obs3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/exp5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/exp3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/obs_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/exp_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/meta_info.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep2/aux_info/ambig_info.tsv     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/cmd_info.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/quant.sf     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/quant.genes.sf     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/lib_format_counts.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/logs/salmon_quant.log     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/libParams/flenDist.txt     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/fld.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/expected_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/observed_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/observed_bias_3p.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/obs5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/obs3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/exp5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/exp3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/obs_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/exp_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/meta_info.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep2/aux_info/ambig_info.tsv     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/cmd_info.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/quant.sf     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/quant.genes.sf     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/lib_format_counts.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/logs/salmon_quant.log     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/libParams/flenDist.txt     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/fld.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/expected_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/observed_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/observed_bias_3p.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/obs5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/obs3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/exp5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/exp3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/obs_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/exp_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/meta_info.json     \
	tests/sample_data/preprocessed/salmon_output/UHR_Rep3/aux_info/ambig_info.tsv     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/cmd_info.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/quant.sf     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/quant.genes.sf     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/lib_format_counts.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/logs/salmon_quant.log     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/libParams/flenDist.txt     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/fld.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/expected_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/observed_bias.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/observed_bias_3p.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/obs5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/obs3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/exp5_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/exp3_seq.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/obs_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/exp_gc.gz     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/meta_info.json     \
	tests/sample_data/preprocessed/salmon_output/HBR_Rep1/aux_info/ambig_info.tsv     \
	tests/sample_data/preprocessed/feature_counts/raw_counts/salmon/counts_from_SALMON.transcripts.tsv     \
	tests/sample_data/preprocessed/feature_counts/raw_counts/salmon/counts_from_SALMON.genes.tsv     \
	tests/sample_data/preprocessed/feature_counts/raw_counts/hisat2/counts.tsv     \
	tests/sample_data/preprocessed/feature_counts/normalized/salmon/TPM_counts_from_SALMON.transcripts.tsv     \
	tests/sample_data/preprocessed/feature_counts/normalized/salmon/TPM_counts_from_SALMON.genes.tsv     \
	tests/sample_data/preprocessed/feature_counts/normalized/hisat2/deseq_size_factors.txt     \
	tests/sample_data/preprocessed/feature_counts/normalized/hisat2/deseq_normalized_counts.tsv     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep2_Aligned.sortedByCoord.out.bam.bai     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep2_Aligned.sortedByCoord.out.bam     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep3_Aligned.sortedByCoord.out.bam     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep1_Aligned.sortedByCoord.out.bam.bai     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep1_Aligned.sortedByCoord.out.bam     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep2_Aligned.sortedByCoord.out.bam     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep2_Aligned.sortedByCoord.out.bam.bai     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep1_Aligned.sortedByCoord.out.bam     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep3_Aligned.sortedByCoord.out.bam     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep3_Aligned.sortedByCoord.out.bam.bai     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep3_Aligned.sortedByCoord.out.bam.bai     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep1.read_counts.csv     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep1_Aligned.sortedByCoord.out.bam.bai     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep2.read_counts.csv     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep2.read_counts.csv     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/UHR_Rep3.read_counts.csv     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep3.read_counts.csv     \
	tests/sample_data/preprocessed/mapped_reads/hisat2/HBR_Rep1.read_counts.csv     \
	tests/sample_data/preprocessed/trimmed_reads/HBR_Rep1.trimmed.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/HBR_Rep2.trimmed.R1.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/HBR_Rep2.trimmed.R2.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/HBR_Rep3.trimmed.R1.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/HBR_Rep3.trimmed.R2.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/UHR_Rep1.trimmed.R1.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/UHR_Rep1.trimmed.R2.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/UHR_Rep2.trimmed.R1.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/UHR_Rep2.trimmed.R2.fq.gz \
	tests/sample_data/preprocessed/trimmed_reads/UHR_Rep3.trimmed.fq.gz \
  tests/settings.yaml \
  tests/settings_no_de.yaml

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/test_genome_coverage/test.sh					\
  tests/test_deseq_reports/test.sh					\
  tests/test_multiqc/test.sh						\
  tests/test_hisat2_counts/test.sh					\
  tests/test_salmon_counts/test.sh					\
  tests/test_salmon/test_salmon_index.sh				\
  tests/test_salmon/test_salmon_quant.sh				\
  tests/test_hisat2/test.sh						\
  test.sh

clean-local:
	-rm -rf .snakemake
	-rm -rf tests/output
	-rm -rf tests/output/.snakemake
	-rm -rf tests/test_deseq_reports/config.json
	-rm -rf tests/test_deseq_reports/logs
	-rm -rf tests/test_deseq_reports/pigx_work
	-rm -rf tests/test_deseq_reports/report
	-rm -rf tests/test_deseq_reports/.snakemake
	-rm -rf tests/test_deseq_reports/test.sh
	-rm -rf tests/test_multiqc/config.json
	-rm -rf tests/test_multiqc/logs
	-rm -rf tests/test_multiqc/multiqc
	-rm -rf tests/test_multiqc/pigx_work
	-rm -rf tests/test_multiqc/.snakemake
	-rm -rf tests/test_multiqc/test.sh
	-rm -rf tests/test_salmon_counts/colData.tsv
	-rm -rf tests/test_salmon_counts/logs/
	-rm -rf tests/test_salmon_counts/pigx_work/
	-rm -rf tests/test_salmon_counts/salmon_output/
	-rm -rf tests/test_salmon_counts/.snakemake
	-rm -rf tests/test_salmon_counts/test.sh
	-rm -rf tests/test_salmon_counts/test.sh.log
	-rm -rf tests/test_salmon_counts/test.sh.trs
	-rm -rf tests/test_salmon/logs/
	-rm -rf tests/test_salmon/pigx_work/
	-rm -rf tests/test_salmon/salmon_index/
	-rm -rf tests/test_salmon/salmon_output/
	-rm -rf tests/test_salmon/.snakemake
	-rm -rf tests/test_salmon/test_salmon_index.sh
	-rm -rf tests/test_salmon/test_salmon_index.sh
	-rm -rf tests/test_salmon/test_salmon_index.sh.log
	-rm -rf tests/test_salmon/test_salmon_index.sh.trs
	-rm -rf tests/test_salmon/test_salmon_quant.sh
	-rm -rf tests/test_salmon/test_salmon_quant.sh
	-rm -rf tests/test_salmon/test_salmon_quant.sh.log
	-rm -rf tests/test_salmon/test_salmon_quant.sh.trs
	-rm -rf tests/test_hisat2_counts/config.json
	-rm -rf tests/test_hisat2_counts/pigx_work
	-rm -rf tests/test_hisat2_counts/preprocessed_data
	-rm -rf tests/test_hisat2_counts/.snakemake
	-rm -rf tests/test_hisat2_counts/test.sh
	-rm -rf tests/test_hisat2/Log.out
	-rm -rf tests/test_hisat2/logs/
	-rm -rf tests/test_hisat2/mapped_reads/
	-rm -rf tests/test_hisat2/pigx_work/
	-rm -rf tests/test_hisat2/.snakemake
	-rm -rf tests/test_hisat2/star_index/
	-rm -rf tests/test_hisat2/test.sh
	-rm -rf tests/test_hisat2/test.sh.log
	-rm -rf tests/test_hisat2/test.sh.trs
CLEANFILES = $(nodist_bin_SCRIPTS) config.json
# Build docker image with Guix
docker: dist
	guix pack -C none -e '(load "guix.scm")' --with-source=pigx_rnaseq-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
