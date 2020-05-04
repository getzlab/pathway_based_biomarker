This folder contains the following files:

1) Drug_Pw_RankSum_Tissue_Specific.m - The main function that recieves as input pathway activity levels or gene expression and drug response AUCs
2) Data.xlsx - This is the input data that was used in this paper (CTD2 pathway activity levels, gene expression and AUCs Z-score)
3) Print_significant_results.m - a sub-function used in the main code(1) for printing results with FDR < 0.25.
4) QQ_PLOT.m - calculation of qq plot for evaluation of the results per tissue type.
5) Down_sampling.m - down samplinh analysis that was performed as disscused in the MS.

The PathOlogist is an open source tool that was developed in the NIH and can be access via the following paper: - PMID: 21542931.