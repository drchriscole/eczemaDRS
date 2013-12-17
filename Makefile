## Makefile for autogenerating analysis and figures for paper...

# Source data
ALL_EXPR = data/all_gene_expression.tsv
CTRL_GENO = data/control_FLG_genotypes.dat
CASE_GENO = data/eczema_FLG_genotypes.dat
CASE_EXPR = data/eczema_gene_expression.tsv

# Output to be generated
FIGS = sample_correlation_heatmap.pdf sample_correlation_dendrogram.pdf FLG_boxplot.pdf
CASE_CNTRL_DATA = edgeR_ctrl_wt_vs_cmpdhet.csv edgeR_ctrl_wt_vs_wt.csv edgeR_ctrl_wt_vs_het.csv
CASE_ONLY_DATA = eczema_wt_vs_cmpdhet.csv eczema_wt_vs_het.csv
SIMPLE_DATA = simple_ctrl_vs_eczema_edgeR.csv
FLG_COR = FLG_correlation_eczema.csv

.PHONY: clean figures analysis all

all: analysis figures

analysis: $(CASE_CNTRL_DATA) $(CASE_ONLY_DATA) $(SIMPLE_DATA) $(FLG_COR)
	@echo "Analysis - Done!"

figures: $(FIGS)
	@echo "Creating figures - Done!"

## Plot sample correlation data
sample_correlation_heatmap.pdf sample_correlation_dendrogram.pdf: bin/sample_correlation.R $(ALL_EXPR) $(CTRL_GENO) $(CASE_GENO)
	Rscript bin/sample_correlation.R $(ALL_EXPR) $(CTRL_GENO) $(CASE_GENO)

## Plot filaggrin gene expression across cases and controls stratified by genotype
FLG_boxplot.pdf: bin/plot_expression_boxplot.R $(ALL_EXPR) $(CTRL_GENO) $(CASE_GENO)
	Rscript bin/plot_expression_boxplot.R $(ALL_EXPR) $(CTRL_GENO) $(CASE_GENO) FLG

## Perform edgeR case-control analysis with genotype stratification
$(CASE_CNTRL_DATA): bin/edgeR_case_control_analysis.R $(ALL_EXPR) $(CTRL_GENO) $(CASE_GENO)
	Rscript bin/edgeR_case_control_analysis.R $(ALL_EXPR) $(CTRL_GENO) $(CASE_GENO)

## Perform edgeR case-case analysis
$(CASE_ONLY_DATA): bin/edgeR_eczema_cases_only.R $(CASE_EXPR) $(CASE_GENO)
	Rscript bin/edgeR_eczema_cases_only.R $(CASE_EXPR) $(CASE_GENO)

## Perform simple case-control analysis with no genotype
$(SIMPLE_DATA): bin/edgeR_case_control_analysis_simple.R $(ALL_EXPR)
	Rscript bin/edgeR_case_control_analysis_simple.R $(ALL_EXPR)

## Perform filaggrin gene expression correlation analysis 
$(FLG_COR): bin/FLG_correlation_analysis.R $(CASE_EXPR) $(CASE_GENO)
	Rscript bin/FLG_correlation_analysis.R $(CASE_EXPR) $(CASE_GENO)

## remove output files
clean:
	rm -f *.pdf *.csv