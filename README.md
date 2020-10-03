# OBmap
PhD thesis project for Kevin Zhu
Advisor: Hiro Matsunami, PhD
Institution: Duke University
Collaborators: Justin Silverman MD PhD, Matt Wachowiak PhD, Shawn Burton PhD, Luis Saraiva PhD, Antonio Scialdone PhD, Mayra Ruiz

Goal: Determine the position of olfactory receptor (OR) glomeruli

Problem: Since 1993, the field has mapped glomeruli positions for only 3% of the 1100 ORs with most of these having incomparable positional information.  By mapping more OR glomeruli, we can gain a better understanding of how olfactory information is organized in the brain.

Approach: 100 micron serial sections of the mouse olfactory bulb (OB) are taken along the anterior-posterior, medial-lateral, and ventral-dorsal axis (1 mouse per dimension). Sections from each dimension are pooled and enriched for low abundance OR sequences using target enrichment probes. Reads are aligned using STAR and quantified with RSEM. All analysis is performed in R. 

Git repo:
inputs/make_TPMmtx/make_obmtx.Rmd performs the merging of TPM columns from multiple RSEM output files.
mri_to_R/mri_OB_cubed.Rmd constructs a voxel representation of a mouse OB from 3D modeling files that originated from MRI scans of the mouse brain.
heatmaps/ob_heatmaps.Rmd creates custom sorted heatmap representations of OR expression across dimension position.  Additional analysis looks into OR covariance and glomeruli symmetry line position.
sc_style/obmap_scstyle.Rmd performs single-cell cluster analysis by treating sections from the same dimensional position as a cluster. UMAP dimensional reduction indicates spatially close samples are transcriptionally similar.
3dimOB/3d_obmapping.Rmd generates a statistical 3D model from single dimension targeted transcriptomics data. OR glomeruli positions are assigned with the help of cluster and dorsal-ventral constraint functions.    
other_analysis/ contains .Rmd files that investigate the processes of read alignment and capture enrichment.

