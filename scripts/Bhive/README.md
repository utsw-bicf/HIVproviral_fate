# An integrated genomics approach deciphering human genome codes shaping HIV proviral transcription and fate

## BHIVE Data processing
BHIVE data was processed with a sinuglarity container https://github.com/gui11aume/BHIVE_for_single_provirus_transcriptomics
Update the expr.nf file to updated_expr.nf before running the data.





## Scripts
updated_expr.nf

### Figures 2B
bhive_annotate_expression.sh (which uses exp_txt2bed.pl) (creates pie chart)
hiv_expression_annotated_bedfiles.sh (makes files for rest of analysis)

### Figures 2C-2H
make_circos.sh (which uses config files in circos_config_files)

### Figures 3B - 3F
BHIVE_expression_by_group_scatterplots_closestTSS.R

### Figures 4B - 4F
BHIVE_expression_by_group_scatterplots_geneexpression.R

BHIVE_expression_by_group_closest_enhancers.sh
BHIVE_expression_by_group_scatterplots_closestEnhancers.R
BHIVE_expression_by_group_scatterplots_closestSuperEnhancers.R

