## GOALS
The idea of this is to make a matrix of features for each cell line out of a binned genome
This way, when we call dots on Hi-C/Micro-C/HiChIP (eventually), we can filter them based on features present on the linear genome
Sort of like a compendium for by-celltype epigenomic features
This gives meaning to dots called on a matrix  - the nature of these methods obscures the meaning of the loops, as they are commonly presented as dots on the interaction matrix

### Todo:
* Suite for analyzing called dots - script that will do actual analysis of dots
  * Anchors treated as features
  * Percent occupancy of anchor spaces by arch. proteins (both anchors, 1 anchor, 0 anchors)
  * Types of loops with respect to enhancers and promoters (E-E / E-P / P-P)
  * Visualization - % occupied, % type of interaction
  * Retain identity of loop when it is in feature matrix (do not mark with 1, instead make a "loop ID")
  * Where to get promoter list from?
  * Gene expression levels?
* How can I implement gene regulation aspects ... annotate promoter bins w/ gene name?
