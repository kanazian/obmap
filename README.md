# OBmap
PhD thesis project for Kevin Zhu
Advisor: Hiro Matsunami, PhD
Institution: Duke University
Collaborators: Justin Silverman MD PhD, Matt Wachowiak PhD, Shawn Burton PhD, Luis Saraiva PhD, Antonio Scialdone PhD, Mayra Ruiz

### Goal: 
Determine the position of olfactory receptor (OR) glomeruli

### Problem: 
Since 1993, the field has mapped glomeruli positions for ~3% of the 1100 ORs with most of these having incomparable positional information.  By mapping more OR glomeruli, we can gain a better understanding of how olfactory information is organized in the brain.

### Approach: 
100 micron serial sections of the mouse olfactory bulb (OB) are taken along the anterior-posterior, medial-lateral, and ventral-dorsal axis (1 mouse per dimension). Sections from each dimension are pooled and enriched for low abundance OR sequences using target enrichment probes. Reads are aligned using STAR and quantified with RSEM. All analysis is performed in R. 
