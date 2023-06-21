# COIMMR
These files are the key script for the COIMMR
Code used for the analysis of the paper:

COIMMR ：a computational framework to reveal the contribution of herbal ingredients against human cancer via immune microenvironment and metabolic reprogramming.

Saisai Tian, Yanan Li, Jia Xu, Lijun Zhang, Jinbo Zhang, Jinyuan Lu，Xike Xu, Xin Luan, Jing Zhao, Weidong Zhang

# STEP 1 GSEA analysis
Description: CalculatING the NES values of corresponding immune signetures and metabolism signatures to establish the immunological and metabolic landscapes.
## input 
     Cancer Gene Expression Profiles
     Gene expression profiles of herbal ingredients
     Immune signatures
     Metabolic signatures
     
# STEP 2 Correlation analysis
Description: Calculation of the immunological similarity score (Iscore) or metabolic similarity score (Mscore).
## input
     Enrichment analysis results for cancer
     Enrichment analysis results for herbal ingredient
     
# STEP 3 permutation test
## input
     Randomly generated immune/metabolic signature (step3.1)
     Cancer Gene Expression Profiles 
     Gene expression profiles of herbal ingredients
     
# STEP 4 
Description: Calculating the direct modulating effect of herbal ingredients on cancer(OCES).
## input
     cancer signature (Cancer up Gene ; Cancer down Gene)
     Gene expression profiles of herbal ingredients
     
# STEP 5
Description:The correlation between OCES and Iscore or Mscore was determined to quantify the contribution of the regulatory effects of herbal ingredients on the immune microenvironment and reprogrammed metabolic processes to their anti-human cancer activities.
## input
     OCES values
     immunological similarity score/metabolic similarity score
     
# STEP 6
Description:Screening and ranking of drug candidates

# System requirement
R studio 4.0.4 R version 4.0.4 (2021-02-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrepel_0.9.1         cmapR_1.6.0           htmltools_0.5.4      
 [4] shinyjs_2.1.0         shinyWidgets_0.5.7    plotly_4.9.3         
 [7] leaflet_2.1.1         DT_0.17               shinythemes_1.2.0    
[10] shinyBS_0.61          forcats_0.5.1         stringr_1.5.0        
[13] purrr_1.0.1           readr_1.4.0           tibble_3.1.8         
[16] tidyverse_1.3.0       tidyr_1.2.0           dplyr_1.0.8          
[19] ggplot2_3.3.5         pROC_1.17.0.1         memoise_2.0.1        
[22] shiny_1.7.4          

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.4.1            
  [3] plyr_1.8.6                  lazyeval_0.2.2             
  [5] flowCore_2.6.0              crosstalk_1.1.1            
  [7] listenv_0.8.0               GenomeInfoDb_1.26.4        
  [9] digest_0.6.31               fansi_1.0.4                
 [11] magrittr_2.0.3              openxlsx_4.2.3             
 [13] globals_0.16.1              modelr_0.1.8               
 [15] RcppParallel_5.0.3          matrixStats_0.58.0         
 [17] cytolib_2.6.2               anytime_0.3.9              
 [19] colorspa



     
    
