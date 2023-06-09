## GOALS
The idea of this is to make a matrix of features for each cell line out of a binned genome
This way, when we call dots on Hi-C/Micro-C/HiChIP (eventually), we can filter them based on features present on the linear genome
Sort of like a compendium for by-celltype epigenomic features
This gives meaning to dots called on a matrix  - the nature of these methods obscures the meaning of the loops, as they are commonly presented as dots on the interaction matrix

### 3/26/23
Clustering of genomic bins marked for presence of multiple eipgenomic features
Sample clustering algorithms that work well with binary (1 or 0) data

### K-Modes
https://cran.r-project.org/web/packages/klaR/index.html

### Bernoulli Mixture model 
https://bayespy.org/examples/bmm.html#results

### Proximus
https://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Clustering/Proximus

### Hidden Markov Model (HMM)
https://medium.com/@postsanjay/hidden-markov-models-simplified-c3f58728caab

### Apriori 
https://github.com/search?q=apriori%20algorithm

### Frequent Pattern (FP) Growth
https://www.javatpoint.com/fp-growth-algorithm-in-data-mining
