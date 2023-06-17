## COIMMR
These files are the key script for the COIMMR
Code used for the analysis of the paper:

COIMMR ：a computational framework to reveal the contribution of herbal ingredients against human cancer via immune microenvironment and metabolic reprogramming.

Saisai Tian, Yanan Li, Jia Xu, Lijun Zhang, Jinbo Zhang, Jinyuan Lu，Xike Xu, Xin Luan, Jing Zhao, Weidong Zhang

## Syetem requirement
R studio 4.0.4

## STEP 1 GSEA analysis
Description: CalculatING the NES values of corresponding immune signetures and metabolism signatures to establish the immunological and metabolic landscapes.
input 
     Cancer Gene Expression Profiles
     Gene expression profiles of herbal ingredients
     Immune signatures
     Metabolic signatures
     
## STEP 2 Correlation analysis
Description: Calculation of the immunological similarity score (Iscore) or metabolic similarity score (Mscore).
input
     Enrichment analysis results for cancer
     Enrichment analysis results for herbal ingredient
     
## STEP 3 permutation test
input
     Randomly generated immune/metabolic signature
     Cancer Gene Expression Profiles
     Gene expression profiles of herbal ingredients
     
## STEP 4 
Description: Calculating the direct modulating effect of herbal ingredients on cancer(OCES).
input
     cancer signature (Cancer up Gene ; Cancer down Gene)
     Gene expression profiles of herbal ingredients
     
## STEP 5
Description:The correlation between OCES and Iscore or Mscore was determined to quantify the contribution of the regulatory effects of herbal ingredients on the immune microenvironment and reprogrammed metabolic processes to their anti-human cancer activities.
input
     OCES values
     immunological similarity score/metabolic similarity score
     
## STEP 6
Description:Screening and ranking of drug candidates




     
    
