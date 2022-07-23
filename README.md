# COSMOS-TPP_paper

## Repository to summarise Master Thesis for paper

For more information please check doc/Notes.md

## Quick summary

Hormone-dependent tumours like ovarian and breast cancers belong to the leading causes of cancer-related deaths in women and are associated with very high relapse and resistance rates. The comprehensive adaption and mutation capacities of these cancer cells challenge biomedical research to improve existing and develop new treatment options. Cancer specific characteristics are exploited in therapy by for example PARP inhibitors, which affect deregulated DNA repair mechanisms. Application of PARP inhibitors requires a deep knowledge of the underlying biology of DNA repair processes, the molecular mechanism of PARP inhibitors and the biological effect of the mutational status on the patient's phenotype.\\

The rapid technological progress and decreasing costs of high throughput experiments result in an increasing number of studies exploiting multi-omics data to understand complex diseases such as cancer. COSMOS is a recently published multi-omics pipeline, which extracts enriched sub-networks for a given functional context via causal reasoning and was used to draw valuable mechanistic hypotheses from multi-omics data analysed in first studies. In this project, we combined COSMOS with Thermal Proteome Profiling (TPP), a mass-spectrometry based technique assessing changes in thermal protein stability, which can reveal a series of functional changes of the proteome upon perturbations. The unique insights TPP provides into the state of the proteome might contribute to COSMOS and complement its current modeling approach. \\

Here, we show how transcription factors, kinases and proteins with altered thermal stability which characterize the response of ovarian cancer cells to PARP inhibition can be identified and integrated in a network-based approach. We successfully derived protein activities for TPP data via the integration with phosphoproteomics measurements and adapted the COSMOS workflow to model these proteins together with kinases and transcription factors. This allowed us to elucidate known consequences on cell cycle and DNA damage response in detail and to suggest additional signaling via the hippo pathway in response to the PARP inhibitor Olaparib.\\

We anticipate that the introduced workflow can be used and adapted for further studies to integrate TPP data with other data types. It holds great potential to be extended with existing knowledge about TPP in the context of metabolism or protein complexes. Also, we expect that applying this workflow to more PARP inhibition data sets will further improve the analysis, interpretation and exploration of personalized PARP inhibition therapies in the future.

![](doc/Overview.jpg)