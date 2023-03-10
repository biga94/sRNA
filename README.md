# MSc project of thesis: *In silico* target prediction for small noncoding RNAs in *Mycobacterium tuberculosis*

**Summary:**


Understanding how bacteria resist antimicrobial treatments, by regulating gene expression in response to antibiotics is a crucial step towards fighting resistance. Also, one of the infections at higher risk to become untreatable is multidrug resistance tuberculosis, due to the lack of new drugs available against Mycobacterium tuberculosis (MTB).  

In this thesis project, we focused on a promising class of agents that intervene at the post- transcriptional level, known as smallRNAs, which alter gene expression in the face of a stressful event for the microorganism, such as the effect of an antibiotic.
We used three target prediction tools - [IntaRNA](https://github.com/BackofenLab/IntaRNA), [IntaRNAsTar](https://github.com/BackofenLab/IntaRNA#IntaRNAsTar) and [sRNARFTarget](https://github.com/BioinformaticsLabAtMUN/sRNARFTarget) - to obtain separate predictions set of interactions between 647 sRNA sequences, coming from previous work, and the entire MTB transcriptome (3976 genes).  
To handle the amount of prediction hits obtained (more than 1.5 million hits per each output) we decided to arbitrarily reduce the datasets to the top 40 interactions per each sRNA candidate, normalize hybridization energy by sRNA length and test the predictive ability by cross-validating of algorithm’s output. IntaRNA from cross-validation results and literature confirmation, resulted being the best prediction algorithm, so we decided to go ahead exploiting only its output in downstream analysis. In parallel, we fitted a generalized extreme value to estimate p-value for each energy interaction. Taking only interactions with p < 0.05, we obtained 21 interactions, 20 of which were not previously known in literature and from previous work. This list of interactions clearly needs biological validation with the most accurate laboratory assays.
Furthermore, we combined sRNA length information with hybridization energy from corresponding mRNA target interaction to develop a k-means clustering and produced
principal component analysis (PCA) graphs. However, results from PCA didn’t suggest distinction of clusters of sRNA with different energy values.
Briefly, candidate 1405 was predicted to target genes involved in fitness recovery upon antibiotic exposure (Rv0098, rpoC, fbpC, hemL, rpsA). Despite candidate 2023 was found to be more likely implicated in virulence pathways, we found a predicted role in regulating bedaquiline response (espI, atpE). Target prediction for candidate 2223 suggested a redundant regulatory role of a network linked to linezolid resistance (pknG, rplC) and cell wall integrity (mgtA, glmU). Candidate 63 was predicted to target mmpL4, a gene reported to be relevant for virulence and pathogenic features rather than drug resistance. Candidate 899 was found to target gyrB encoding a key protein (DNA gyrase) for the mechanism of action of fluoroquinolones. Notably, our mutation analysis based on WGS data and associated phenotypic drug susceptibility testing for nearly 2000 clinical isolates supports the associations outlined above.
To further characterize our sRNA panel, we first run multiple local alignments, also accounting for the secondary structure, then supported by the production of coverage and region plots. Results didn’t suggest similarities between sRNA sequences, apart from a small region in 5’ between candidates 2023 and 2223 that requires further investigation.
An in silico compensatory point mutations analysis has been taken out on our sRNA panel in order to evaluate minimal free energy variations in seed region with complementary mRNA target that could lead to a compromise of the interaction. However, any relevant mutation emerged from the simulations. A future endeavor would be the evaluation of multiple-pairs mutations along with the only evaluation of single-point mutations. Apart from in silico assays, a wet validation is required to confirm biologically our sRNA panel and testing effectively the predictive ability of the tools used.


In conclusion, in this thesis we have performed and determined the following results: 
1. We tested and compared 3 different tools to study RNA-RNA target prediction;
2. We determined and characterized 5 new sRNA, involved in drug resistance mechanism;
3. We established a new in silico protocol to help in targeting RNA-RNA binding.

For the first time, our analysis suggests the existence of more complex post-transcriptional regulatory networks around key drug resistant genes commonly considered to identify mutations marker of resistance to be included in molecular diagnostic assays. Our study shed lights on putative involvement of post-transcriptional regulation for genes relevant for drug resistance in M. tuberculosis, thus highlighting the importance to investigate sRNAs and their post-transcriptional regulatory role in antibiotic resistance.  

## Workflow
![Workflow](https://github.com/biga94/sRNA/blob/main/workflow.png)
