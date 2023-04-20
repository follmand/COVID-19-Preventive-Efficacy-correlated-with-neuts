# COVID-19 preventive efficacy orrelated with neutralizing antibodies

In this repository include two R files, example_code_for_COV-2069.r and example_code_for_COVE.r. The COV-2069 file provides the code used to fit the 3PL model + the antibody prediction from concentration to neutralization titer. The COV file provides the code used to fit the log-linear model and also antibody decay and generates simulated data.

The code has been run on R 4.2.2 and requires recent versions of the libraries listed at the top of the files. All required libraries can be downloaded from CRAN via the standard install.packages() function in R. Users should be aware that installation of the brms package in the COVE example may require additional setup to enable compiled code to run. Additional instructions may be found in the brms documentation. The session information for the machine on which these scripts were developed and run is given at the bottom of this Readme file. The total run time for the COVE example should be under one minutes on most modern computers while the COV-2069 file should take roughly 15 minutes. 

You can run example_code_for_COVE.r by opening up R and executing the source code for each script in the console. Both examples are self contained and can be run without special instructions once the required libraries are installed. 


## Session information
sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.3.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets 
[6] methods   base     

other attached packages:
[1] tidyr_1.3.0    survival_3.4-0 mgcv_1.8-41   
[4] nlme_3.1-160   dplyr_1.1.0    brms_2.18.0   
[7] Rcpp_1.0.10   

loaded via a namespace (and not attached):
 [1] matrixStats_0.63.0   xts_0.12.2          
 [3] threejs_0.3.3        rstan_2.21.7        
 [5] tensorA_0.36.2       tools_4.2.2         
 [7] backports_1.4.1      utf8_1.2.3          
 [9] R6_2.5.1             DT_0.26             
[11] DBI_1.1.3            colorspace_2.1-0    
[13] tidyselect_1.2.0     gridExtra_2.3       
[15] prettyunits_1.1.1    processx_3.8.0      
[17] Brobdingnag_1.2-9    emmeans_1.8.3       
[19] compiler_4.2.2       cli_3.6.0           
[21] shinyjs_2.1.0        sandwich_3.0-2      
[23] colourpicker_1.2.0   posterior_1.3.1     
[25] scales_1.2.1         dygraphs_1.1.1.6    
[27] checkmate_2.1.0      mvtnorm_1.1-3       
[29] callr_3.7.3          stringr_1.5.0       
[31] digest_0.6.31        StanHeaders_2.21.0-7
[33] base64enc_0.1-3      pkgconfig_2.0.3     
[35] htmltools_0.5.4      fastmap_1.1.0       
[37] htmlwidgets_1.5.4    rlang_1.0.6         
[39] rstudioapi_0.14      shiny_1.7.3         
[41] farver_2.1.1         generics_0.1.3      
[43] zoo_1.8-11           crosstalk_1.2.0     
[45] gtools_3.9.4         distributional_0.3.1
[47] inline_0.3.19        magrittr_2.0.3      
[49] loo_2.5.1            bayesplot_1.10.0    
[51] Matrix_1.5-3         munsell_0.5.0       
[53] fansi_1.0.4          abind_1.4-5         
[55] lifecycle_1.0.3      stringi_1.7.8       
[57] multcomp_1.4-20      MASS_7.3-58.1       
[59] pkgbuild_1.4.0       plyr_1.8.8          
[61] grid_4.2.2           parallel_4.2.2      
[63] promises_1.2.0.1     crayon_1.5.2        
[65] miniUI_0.1.1.1       lattice_0.20-45     
[67] splines_4.2.2        ps_1.7.2            
[69] pillar_1.8.1         igraph_1.3.5        
[71] markdown_1.4         estimability_1.4.1  
[73] shinystan_2.6.0      reshape2_1.4.4      
[75] codetools_0.2-18     stats4_4.2.2        
[77] rstantools_2.2.0     glue_1.6.2          
[79] RcppParallel_5.1.5   vctrs_0.5.2         
[81] httpuv_1.6.6         purrr_1.0.1         
[83] gtable_0.3.1         ggplot2_3.4.0       
[85] mime_0.12            xtable_1.8-4        
[87] coda_0.19-4          later_1.3.0         
[89] tibble_3.1.8         shinythemes_1.2.0   
[91] TH.data_1.1-1        ellipsis_0.3.2      
[93] bridgesampling_1.1-2
