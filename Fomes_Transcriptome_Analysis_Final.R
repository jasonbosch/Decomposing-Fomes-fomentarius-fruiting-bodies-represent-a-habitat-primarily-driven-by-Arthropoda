#Analysis for: Decomposing Fomes fomentarius fruiting bodies, unlike fresh ones, represent a bacteria-rich habitat primarily driven by Arthropoda
#Jason Bosch

#In order to run the script and repeat the analysis, it is necessary to download the data files.
#All relevant, processed data files can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.10081831

#####SET UP#####

###Reproducibility###

library(groundhog)

groundhog_day <- "2023-09-01"

meta.groundhog(groundhog_day)

###Libraries###

#CRAN libraries
groundhog.library(c("ggplot2", #Plotting
                    "RColorBrewer", #For colours
                    "pheatmap", #Generate heatmap
                    "ggrepel", # Keep PCA labels from overlapping
                    "cowplot", #Joining plots together
                    "cluster", #For clustering
                    "NbClust", #For clustering
                    "shipunov", #For bootstrapping clustering
                    "ggdendro", #For plotting bootstrapped clusters
                    "stringr", #For manipulating taxonomy strings
                    "ggplotify", #Convert pheatmap into GGplot
                    "ggnewscale"), #Allows multiple scales to be used in a single ggplot
                  groundhog_day) 

#BioConductor libraries
library(Rsubread) #Doing counting reads

###Directory###

setwd("~/PostDoc/02_Projects/05_Fomes_transcription/04_Analysis_Results/R_Final/")

###Functions###

#Combine all the low abundance taxa (below 1%)
Combine_Low_Abundance <- function(INPUT_DF,COMBINE_COL,COUNT_COL,THRESHOLD) {
  #Function will take a data frame and join everything below a certain abundance threshold.
  #Find the abundances of the targets and whether they are less than the threshold
  Combine_Abundances <- aggregate(as.formula(paste(".~",as.name(COMBINE_COL),sep = "")),INPUT_DF[,c(COMBINE_COL,COUNT_COL)],sum)
  Combine_Below_Threshold <- Combine_Abundances[Combine_Abundances[,COUNT_COL]/sum(Combine_Abundances[,COUNT_COL])<THRESHOLD,COMBINE_COL]
  #If nothing is below then nothing needs to be done
  if (length(Combine_Below_Threshold)==0) {
    return(INPUT_DF)
  }
  #Otherwise combine
  else {
    Result_DF <- INPUT_DF
    Result_DF[which(Result_DF[,COMBINE_COL]%in%Combine_Below_Threshold),COMBINE_COL] <- paste("Below ",THRESHOLD*100,"%",sep = "")
    return(Result_DF) 
  }
}

#Takes a list of KOs and returns a dataframe with different KEGG modules and whether the complete pathway is possible from the
#given KOs. Due to the way data was stored, it evaluates KOs within square brackets.
EvaluatePathways <- function(KOs) {
  Res_DF <- as.data.frame(matrix(nrow = 0, ncol = 4,dimnames = list(NULL,c("Energy Metabolism","Module","Name","Complete"))))
  
  #The definition of the module is a logical expression of K numbers. Comma separated K numbers indicate alternatives. 
  #Plus signs are used to represent a complex and a minus sign denotes a non-essential component in the complex.
  #Non-essential components are removed.
  
  #M00175
  #Nitrogen fixation, nitrogen => ammonia
  #K02588+K02586+K02591-K00531,K22896+K22897+K22898+K22899
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Nitrogen metabolism","M00175","Nitrogen fixation, nitrogen => ammonia",
      ("[K02588]"%in%KOs&"[K02586]"%in%KOs&"[K02591]"%in%KOs)|("[K22896]"%in%KOs&"[K22897]"%in%KOs&"[K22898]"%in%KOs&"[K22899]"%in%KOs)
    )
  
  #Assimilatory nitrate reduction, nitrate => ammonia
  #M00531
  #(K00367,K10534,K00372-K00360) (K00366,K17877,K26139+K26138,K00361)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Nitrogen metabolism","M00531","Assimilatory nitrate reduction, nitrate => ammonia",
      ("[K00367]"%in%KOs|"[K10534]"%in%KOs|"[K00372]"%in%KOs)&("[K00366]"%in%KOs|"[K17877]"%in%KOs|("[K26139]"%in%KOs&"[K26138]"%in%KOs)|"[K00361]"%in%KOs)
    )
  
  #Dissimilatory nitrate reduction, nitrate => ammonia
  #M00530
  #(K00370+K00371+K00374,K02567+K02568) (K00362+K00363,K03385+K15876)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Nitrogen metabolism","M00530","Dissimilatory nitrate reduction, nitrate => ammonia",
      (("[K00370]"%in%KOs&"[K00371]"%in%KOs&"[K00374]"%in%KOs)|("[K02567]"%in%KOs&"[K02568]"%in%KOs))&(("[K00362]"%in%KOs&"[K00363]"%in%KOs)|("[K03385]"%in%KOs&"[K15876]"%in%KOs))
    )
  
  #Denitrification, nitrate => nitrogen
  #M00529
  #(K00370+K00371+K00374,K02567+K02568) (K00368,K15864) (K04561+K02305) K00376
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Nitrogen metabolism","M00529","Denitrification, nitrate => nitrogen",
      (("[K00370]"%in%KOs&"[K00371]"%in%KOs&"[K00374]"%in%KOs)|("[K02567]"%in%KOs&"[K02568]"%in%KOs))&("[K00368]"%in%KOs|"[K15864]"%in%KOs)&("[K04561]"%in%KOs&"[K02305]"%in%KOs)&"[K00376]"%in%KOs
    )
  
  #Nitrification, ammonia => nitrite
  #M00528
  #K10944+K10945+K10946 K10535
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Nitrogen metabolism","M00528","Nitrification, ammonia => nitrite",
      ("[K10944]"%in%KOs&"[K10945]"%in%KOs&"[K10946]"%in%KOs)&"[K10535]"%in%KOs
    )
  
  #Complete nitrification, comammox, ammonia => nitrite => nitrate
  #M00804
  #K10944+K10945+K10946 K10535 K00370+K00371
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Nitrogen metabolism","M00804","Complete nitrification, comammox, ammonia => nitrite => nitrate",
      ("[K10944]"%in%KOs&"[K10945]"%in%KOs&"[K10946]"%in%KOs)&"[K10535]"%in%KOs&("[K00370]"%in%KOs&"[K00371]"%in%KOs)
    )
  
  #M00176
  #Assimilatory sulfate reduction, sulfate => H2S
  #(K13811,K00958+K00860,K00955+K00957,K00956+K00957+K00860) K00390 (K00380+K00381,K00392)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Sulfur metabolism","M00176","Assimilatory sulfate reduction, sulfate => H2S",
      ("[K13811]"%in%KOs|("[K00958]"%in%KOs&"[K00860]"%in%KOs)|("[K00955]"%in%KOs&"[K00957]"%in%KOs)|("[K00956]"%in%KOs&"[K00957]"%in%KOs&"[K00860]"%in%KOs))&"[K00390]"%in%KOs&(("[K00380]"%in%KOs&"[K00381]"%in%KOs)|"[K00392]"%in%KOs)
    )
  
  #M00596
  #Dissimilatory sulfate reduction, sulfate => H2S
  #K00958 (K00394+K00395) (K11180+K11181)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Sulfur metabolism","M00596","Dissimilatory sulfate reduction, sulfate => H2S",
      "[K00958]"%in%KOs&("[K00394]"%in%KOs&"[K00395]"%in%KOs)&("[K11180]"%in%KOs&"[K11181]"%in%KOs)
    )
  
  #M00595
  #Thiosulfate oxidation by SOX complex, thiosulfate => sulfate
  #K17222+K17223+K17224-K17225-K22622+K17226+K17227
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Sulfur metabolism","M00595","Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
      "[K17222]"%in%KOs&"[K17223]"%in%KOs&"[K17224]"%in%KOs&"[K17226]"%in%KOs&"[K17227]"%in%KOs
    )
  
  #M00161
  #Photosystem II
  #K02703+K02706+K02705+K02704+K02707+K02708
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Photosynthesis","M00161","Photosystem II",
      "[K02703]"%in%KOs&"[K02706]"%in%KOs&"[K02705]"%in%KOs&"[K02704]"%in%KOs&"[K02707]"%in%KOs&"[K02708]"%in%KOs
    )
  
  #M00163
  #Photosystem I
  #K02689+K02690+K02691+K02692+K02693+K02694
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Photosynthesis","M00163","Photosystem I",
      "[K02689]"%in%KOs&"[K02690]"%in%KOs&"[K02691]"%in%KOs&"[K02692]"%in%KOs&"[K02693]"%in%KOs&"[K02694]"%in%KOs
    )
  
  #M00597
  #Anoxygenic photosystem II
  #K08928+K08929
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Photosynthesis","M00597","Anoxygenic photosystem II",
      "[K08928]"%in%KOs&"[K08929]"%in%KOs
    )
  
  #M00598
  #Anoxygenic photosystem I
  #K08940+K08941+K08942+K08943
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Photosynthesis","M00598","Anoxygenic photosystem I",
      "[K08940]"%in%KOs&"[K08941]"%in%KOs&"[K08942]"%in%KOs&"[K08943]"%in%KOs
    )
  
  #M00567
  #Methanogenesis, CO2 => methane
  #(K00200+K00201+K00202+K00203-K11261+(K00205,K11260,K00204)) K00672 K01499 (K00319,K13942) K00320 (K00577+K00578+K00579+K00580+K00581-K00582-K00583+K00584) (K00399+K00401+K00402) (K22480+K22481+K22482,K03388+K03389+K03390,K08264+K08265,K03388+K03389+K03390+K14127+(K14126+K14128,K22516+K00125))
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00567","Methanogenesis, CO2 => methane",
      ("[K00200]"%in%KOs&"[K00201]"%in%KOs&"[K00202]"%in%KOs&"[K00203]"%in%KOs&("[K00205]"%in%KOs|"[K11260]"%in%KOs|"[K00204]"%in%KOs))&"[K00672]"%in%KOs&"[K01499]"%in%KOs&("[K00319]"%in%KOs|"[K13942]"%in%KOs)&"[K00320]"%in%KOs&("[K00577]"%in%KOs&"[K00578]"%in%KOs&"[K00579]"%in%KOs&"[K00580]"%in%KOs&"[K00581]"%in%KOs&"[K00584]"%in%KOs)&("[K00399]"%in%KOs&"[K00401]"%in%KOs&"[K00402]"%in%KOs)&(("[K22480]"%in%KOs&"[K22481]"%in%KOs&"[K22482]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs)|("[K08264]"%in%KOs&"[K08265]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs&"[K14127]"%in%KOs&(("[K14126]"%in%KOs&"[K14128]"%in%KOs)|("[K22516]"%in%KOs&"[K00125]"%in%KOs))))
    )
  
  #M00357 
  #Methanogenesis, acetate => methane
  #((K00925 K00625),K01895) (K00193+K00197+K00194) (K00577+K00578+K00579+K00580+K00581-K00582-K00583+K00584) (K00399+K00401+K00402) (K22480+K22481+K22482,K03388+K03389+K03390,K08264+K08265,K03388+K03389+K03390+K14127+(K14126+K14128,K22516+K00125))
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00357","Methanogenesis, acetate => methane",
      (("[K00925]"%in%KOs&"[K00625]"%in%KOs)|"[K01895]"%in%KOs)&("[K00193]"%in%KOs&"[K00197]"%in%KOs&"[K00194]"%in%KOs)&("[K00577]"%in%KOs&"[K00578]"%in%KOs&"[K00579]"%in%KOs&"[K00580]"%in%KOs&"[K00581]"%in%KOs&"[K00584]"%in%KOs)&("[K00399]"%in%KOs&"[K00401]"%in%KOs&"[K00402]"%in%KOs)&(("[K22480]"%in%KOs&"[K22481]"%in%KOs&"[K22482]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs)|("[K08264]"%in%KOs&"[K08265]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs&"[K14127]"%in%KOs&(("[K14126]"%in%KOs&"[K14128]"%in%KOs)|("[K22516]"%in%KOs&"[K00125]"%in%KOs))))
    )
  
  #M00356
  #Methanogenesis, methanol => methane
  #K14080+K04480+K14081 K00399+K00401+K00402 (K22480+K22481+K22482,K03388+K03389+K03390,K08264+K08265,K03388+K03389+K03390+K14127+(K14126+K14128,K22516+K00125))
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00356","Methanogenesis, methanol => methane",
      ("[K14080]"%in%KOs&"[K04480]"%in%KOs&"[K14081]"%in%KOs)&("[K00399]"%in%KOs&"[K00401]"%in%KOs&"[K00402]"%in%KOs)&(("[K22480]"%in%KOs&"[K22481]"%in%KOs&"[K22482]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs)|("[K08264]"%in%KOs&"[K08265]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs&"[K14127]"%in%KOs&(("[K14126]"%in%KOs&"[K14128]"%in%KOs)|("[K22516]"%in%KOs&"[K00125]"%in%KOs))))
    )
  
  #M00563
  #Methanogenesis, methylamine/dimethylamine/trimethylamine => methane
  #K14082 ((K16177-K16176),(K16179-K16178),(K14084-K14083)) K00399+K00401+K00402 (K22480+K22481+K22482,K03388+K03389+K03390,K08264+K08265,K03388+K03389+K03390+K14127+(K14126+K14128,K22516+K00125))
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00563","Methanogenesis, methylamine/dimethylamine/trimethylamine => methane",
      "[K14082]"%in%KOs&(("[K16177]"%in%KOs)|("[K16179]"%in%KOs)|("[K14084]"%in%KOs))&("[K00399]"%in%KOs&"[K00401]"%in%KOs&"[K00402]"%in%KOs)&(("[K22480]"%in%KOs&"[K22481]"%in%KOs&"[K22482]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs)|("[K08264]"%in%KOs&"[K08265]"%in%KOs)|("[K03388]"%in%KOs&"[K03389]"%in%KOs&"[K03390]"%in%KOs&"[K14127]"%in%KOs&(("[K14126]"%in%KOs&"[K14128]"%in%KOs)|("[K22516]"%in%KOs&"[K00125]"%in%KOs))))
    )
  
  #M00358
  #Coenzyme M biosynthesis
  #K08097 K05979 K05884 K13039+K06034
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00358","Coenzyme M biosynthesis",
      "[K08097]"%in%KOs&"[K05979]"%in%KOs&"[K05884]"%in%KOs&("[K13039]"%in%KOs&"[K06034]"%in%KOs)
    )
  
  #M00608
  #2-Oxocarboxylic acid chain extension, 2-oxoglutarate => 2-oxoadipate => 2-oxopimelate => 2-oxosuberate
  #K10977 K16792+K16793 K10978
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00608","2-Oxocarboxylic acid chain extension, 2-oxoglutarate => 2-oxoadipate => 2-oxopimelate => 2-oxosuberate",
      "[K10977]"%in%KOs&("[K16792]"%in%KOs&"[K16793]"%in%KOs)&"[K10978]"%in%KOs
    )
  
  #M00174
  #Methane oxidation, methanotroph, methane => formaldehyde
  #((K10944+K10945+K10946),(K16157+K16158+K16159+K16160+K16161+K16162)) ((K14028-K14029),K23995)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00174","Methane oxidation, methanotroph, methane => formaldehyde",
      (("[K10944]"%in%KOs&"[K10945]"%in%KOs&"[K10946]"%in%KOs)|("[K16157]"%in%KOs&"[K16158]"%in%KOs&"[K16159]"%in%KOs&"[K16160]"%in%KOs&"[K16161]"%in%KOs&"[K16162]"%in%KOs))&(("[K14028]"%in%KOs)|"[K23995]"%in%KOs)
    )
  
  #M00346
  #Formaldehyde assimilation, serine pathway
  #K00600 K00830 K00018 K11529 K01689 K01595 K00024 K08692+K14067 K08691
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00346","Formaldehyde assimilation, serine pathway",
      "[K00600]"%in%KOs&"[K00830]"%in%KOs&"[K00018]"%in%KOs&"[K11529]"%in%KOs&"[K01689]"%in%KOs&"[K01595]"%in%KOs&"[K00024]"%in%KOs&("[K08692]"%in%KOs&"[K14067]"%in%KOs)&"[K08691]"%in%KOs
    )
  
  #M00345
  #Formaldehyde assimilation, ribulose monophosphate pathway
  #(((K08093,K13812) K08094),K13831) (K00850,K16370) K01624
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00345","Formaldehyde assimilation, ribulose monophosphate pathway",
      ((("[K08093]"%in%KOs|"[K13812]"%in%KOs)&"[K08094]"%in%KOs)|"[K13831]"%in%KOs)&("[K00850]"%in%KOs|"[K16370]"%in%KOs)&"[K01624]"%in%KOs
    )
  
  #M00344
  #Formaldehyde assimilation, xylulose monophosphate pathway
  #K17100 K00863 K01624 K03841
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00344","Formaldehyde assimilation, xylulose monophosphate pathway",
      "[K17100]"%in%KOs&"[K00863]"%in%KOs&"[K01624]"%in%KOs&"[K03841]"%in%KOs
    )
  
  #M00378
  #F420 biosynthesis, archaea
  #K11780+K11781
  #K14941
  #K11212 K12234
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00378","F420 biosynthesis, archaea",
      ("[K11780]"%in%KOs&"[K11781]"%in%KOs)|"[K14941]"%in%KOs|("[K11212]"%in%KOs&"[K12234]"%in%KOs)
    )
  
  #M00935
  #Methanofuran biosynthesis
  #K09733 K19793 K07144
  #K18933 K06914
  #K07072
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00935","Methanofuran biosynthesis",
      ("[K09733]"%in%KOs&"[K19793]"%in%KOs&"[K07144]"%in%KOs)|("[K18933]"%in%KOs&"[K06914]"%in%KOs)|"[K07072]"%in%KOs
    )
  
  #M00422
  #Acetyl-CoA pathway, CO2 => acetyl-CoA
  #K00192+K00195 K00193+K00197+K00194
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Methane metabolism","M00422","Acetyl-CoA pathway, CO2 => acetyl-CoA",
      ("[K00192]"%in%KOs&"[K00195]"%in%KOs)&("[K00193]"%in%KOs&"[K00197]"%in%KOs&"[K00194]"%in%KOs)
    )
  
  #M00165
  #Reductive pentose phosphate cycle (Calvin cycle)
  #K00855 (K01601-K01602) K00927 (K05298,K00150,K00134) (K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01623,K01624) (K01100,K11532,K01086) K00615 (K01807,K01808)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00165","Reductive pentose phosphate cycle (Calvin cycle)",
      "[K00855]"%in%KOs&("[K01601]"%in%KOs)&"[K00927]"%in%KOs&("[K05298]"%in%KOs|"[K00150]"%in%KOs|"[K00134]"%in%KOs)&("[K01623]"%in%KOs|"[K01624]"%in%KOs)&("[K03841]"%in%KOs|"[K02446]"%in%KOs|"[K11532]"%in%KOs|"[K01086]"%in%KOs)&"[K00615]"%in%KOs&("[K01623]"%in%KOs|"[K01624]"%in%KOs)&("[K01100]"%in%KOs|"[K11532]"%in%KOs|"[K01086]"%in%KOs)&"[K00615]"%in%KOs&("[K01807]"%in%KOs|"[K01808]"%in%KOs)
    )

  #M00166
  #Reductive pentose phosphate cycle, ribulose-5P => glyceraldehyde-3P
  #K00855 (K01601-K01602) K00927 (K05298,K00150,K00134)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00166","Reductive pentose phosphate cycle, ribulose-5P => glyceraldehyde-3P",
      "[K00855]"%in%KOs&("[K01601]"%in%KOs)&"[K00927]"%in%KOs&("[K05298]"%in%KOs|"[K00150]"%in%KOs|"[K00134]"%in%KOs)
    )
  
  #M00167
  #Reductive pentose phosphate cycle, glyceraldehyde-3P => ribulose-5P
  #(K01623,K01624) (K03841,K02446,K11532,K01086) K00615 (K01623,K01624) (K01100,K11532,K01086) K00615 (K01807,K01808)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00167","Reductive pentose phosphate cycle, glyceraldehyde-3P => ribulose-5P",
      ("[K01623]"%in%KOs|"[K01624]"%in%KOs)&("[K03841]"%in%KOs|"[K02446]"%in%KOs|"[K11532]"%in%KOs|"[K01086]"%in%KOs)&"[K00615]"%in%KOs&("[K01623]"%in%KOs|"[K01624]"%in%KOs)&("[K01100]"%in%KOs|"[K11532]"%in%KOs|"[K01086]"%in%KOs)&"[K00615]"%in%KOs&("[K01807]"%in%KOs|"[K01808]"%in%KOs)
    )
  
  #M00168
  #CAM (Crassulacean acid metabolism), dark
  #K01595 (K00025,K00026,K00024)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00168","CAM (Crassulacean acid metabolism), dark",
      "[K01595]"%in%KOs&("[K00025]"%in%KOs|"[K00026]"%in%KOs|"[K00024]"%in%KOs)
    )
  
  #M00169
  #CAM (Crassulacean acid metabolism), light
  #K00029 K01006
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00169","CAM (Crassulacean acid metabolism), light",
      "[K00029]"%in%KOs&"[K01006]"%in%KOs
    )
  
  #M00172
  #C4-dicarboxylic acid cycle, NADP - malic enzyme type
  #K01595 K00051 K00029 K01006
  Res_DF[nrow(Res_DF)+1,] <- 
  c("Carbon fixation","M00172","C4-dicarboxylic acid cycle, NADP - malic enzyme type",
    "[K01595]"%in%KOs&"[K00051]"%in%KOs&"[K00029]"%in%KOs&"[K01006]"%in%KOs
  )
  
  #M00171
  #C4-dicarboxylic acid cycle, NAD - malic enzyme type
  #K01595 K14454 K14455 (K00025,K00026) K00028 (K00814,K14272) K01006
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00171","C4-dicarboxylic acid cycle, NAD - malic enzyme type",
      "[K01595]"%in%KOs&"[K14454]"%in%KOs&"[K14455]"%in%KOs&("[K00025]"%in%KOs|"[K00026]"%in%KOs)&"[K00028]"%in%KOs&("[K00814]"%in%KOs|"[K14272]"%in%KOs)&"[K01006]"%in%KOs
    )
  
  #M00170
  #C4-dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type
  #K01595 K14454 K14455 K01610
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00170","C4-dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type",
      "[K01595]"%in%KOs&"[K14454]"%in%KOs&"[K14455]"%in%KOs&"[K01610]"%in%KOs
    )
  
  #M00173
  #Reductive citrate cycle (Arnon-Buchanan cycle)
  #(K00169+K00170+K00171+K00172,K03737) ((K01007,K01006) K01595,K01959+K01960,K01958) K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K00244+K00245-K00246-K00247,K18556+K18557+K18558+K18559+K18560) (K01902+K01903) (K00174+K00175-K00177-K00176) K00031 (K01681,K01682) (K15230+K15231,K15232+K15233 K15234)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00173","Reductive citrate cycle (Arnon-Buchanan cycle)",
      (("[K00169]"%in%KOs&"[K00170]"%in%KOs&"[K00171]"%in%KOs&"[K00172]"%in%KOs)|"[K03737]"%in%KOs)&((("[K01007]"%in%KOs|"[K01006]"%in%KOs)&"[K01595]"%in%KOs)|("[K01959]"%in%KOs&"[K01960]"%in%KOs)|"[K01958]"%in%KOs)&"[K00024]"%in%KOs&("[K01676]"%in%KOs|"[K01679]"%in%KOs|("[K01677]"%in%KOs&"[K01678]"%in%KOs))&(("[K00239]"%in%KOs&"[K00240]"%in%KOs)|("[K00244]"%in%KOs&"[K00245]"%in%KOs)|("[K18556]"%in%KOs&"[K18557]"%in%KOs&"[K18558]"%in%KOs&"[K18559]"%in%KOs&"[K18560]"%in%KOs))&("[K01902]"%in%KOs&"[K01903]"%in%KOs)&("[K00174]"%in%KOs&"[K00175]"%in%KOs)&"[K00031]"%in%KOs&("[K01681]"%in%KOs|"[K01682]"%in%KOs)&(("[K15230]"%in%KOs&"[K15231]"%in%KOs)|(("[K15232]"%in%KOs&"[K15233]"%in%KOs)&"[K15234]"%in%KOs))
    )
  
  #M00376
  #3-Hydroxypropionate bi-cycle
  #(K02160+K01961+K01962+K01963) K14468 K14469 K15052 K05606 (K01847,K01848+K01849) (K14471+K14472) (K00239+K00240+K00241) K01679
  #K08691 K14449 K14470 K09709
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00376","3-Hydroxypropionate bi-cycle",
      (("[K02160]"%in%KOs&"[K01961]"%in%KOs&"[K01962]"%in%KOs&"[K01963]"%in%KOs)&"[K14468]"%in%KOs&"[K14469]"%in%KOs&"[K15052]"%in%KOs&"[K05606]"%in%KOs&("[K01847]"%in%KOs|("[K01848]"%in%KOs&"[K01849]"%in%KOs))&("[K14471]"%in%KOs&"[K14472]"%in%KOs)&("[K00239]"%in%KOs&"[K00240]"%in%KOs&"[K00241]"%in%KOs)&"[K01679]"%in%KOs)|("[K08691]"%in%KOs&"[K14449]"%in%KOs&"[K14470]"%in%KOs&"[K09709]"%in%KOs)
    )
  
  #M00375 
  #Hydroxypropionate-hydroxybutylate cycle
  #(K01964+K15037+K15036) K15017 K15039 K15018 K15019 K15020 K05606 (K01848+K01849) (K15038,K15017) K14465 (K14466,K18861) K14534 K15016 K00626
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00375","Hydroxypropionate-hydroxybutylate cycle",
      ("[K01964]"%in%KOs&"[K15037]"%in%KOs&"[K15036]"%in%KOs)&"[K15017]"%in%KOs&"[K15039]"%in%KOs&"[K15018]"%in%KOs&"[K15019]"%in%KOs&"[K15020]"%in%KOs&"[K05606]"%in%KOs&("[K01848]"%in%KOs&"[K01849]"%in%KOs)&("[K15038]"%in%KOs|"[K15017]"%in%KOs)&"[K14465]"%in%KOs&("[K14466]"%in%KOs|"[K18861]"%in%KOs)&"[K14534]"%in%KOs&"[K15016]"%in%KOs&"[K00626]"%in%KOs
    )
  
  #M00374
  #Dicarboxylate-hydroxybutyrate cycle
  #(K00169+K00170+K00171+K00172) K01007 K01595 K00024 (K01677+K01678) (K00239+K00240-K00241-K18860) (K01902+K01903) (K15038,K15017) K14465 (K14467,K18861) K14534 K15016 K00626
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00374","Dicarboxylate-hydroxybutyrate cycle",
      ("[K00169]"%in%KOs&"[K00170]"%in%KOs&"[K00171]"%in%KOs&"[K00172]"%in%KOs)&"[K01007]"%in%KOs&"[K01595]"%in%KOs&"[K00024]"%in%KOs&("[K01677]"%in%KOs&"[K01678]"%in%KOs)&("[K00239]"%in%KOs&"[K00240]"%in%KOs)&("[K01902]"%in%KOs&"[K01903]"%in%KOs)&("[K15038]"%in%KOs|"[K15017]"%in%KOs)&"[K14465]"%in%KOs&("[K14467]"%in%KOs|"[K18861]"%in%KOs)&"[K14534]"%in%KOs&"[K15016]"%in%KOs&"[K00626]"%in%KOs
    )
  
  #M00377
  #Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)
  #K00198 (K05299-K15022,K22015+K25123+K25124) K01938 K01491-K01500 K00297-K25007-K25008 K15023 K14138+K00197+K00194
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00377","Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)",
      "[K00198]"%in%KOs&(("[K05299]"%in%KOs)|("[K22015]"%in%KOs&"[K25123]"%in%KOs&"[K25124]"%in%KOs))&"[K01938]"%in%KOs&"[K01491]"%in%KOs&"[K00297]"%in%KOs&"[K15023]"%in%KOs&("[K14138]"%in%KOs&"[K00197]"%in%KOs&"[K00194]"%in%KOs)
    )
  
  #M00579
  #Phosphate acetyltransferase-acetate kinase pathway, acetyl-CoA => acetate
  #(K00625,K13788,K15024) K00925
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00579","Phosphate acetyltransferase-acetate kinase pathway, acetyl-CoA => acetate",
      ("[K00625]"%in%KOs|"[K13788]"%in%KOs|"[K15024]"%in%KOs)&"[K00925]"%in%KOs
    )
  
  #M00620
  #Incomplete reductive citrate cycle, acetyl-CoA => oxoglutarate
  #(K00169+K00170+K00171+K00172) (K01959+K01960) K00024 (K01677+K01678) (K18209+K18210) (K01902+K01903) (K00174+K00175+K00176+K00177)
  Res_DF[nrow(Res_DF)+1,] <- 
    c("Carbon fixation","M00620","Incomplete reductive citrate cycle, acetyl-CoA => oxoglutarate",
      ("[K00169]"%in%KOs&"[K00170]"%in%KOs&"[K00171]"%in%KOs&"[K00172]"%in%KOs)&("[K01959]"%in%KOs&"[K01960]"%in%KOs)&"[K00024]"%in%KOs&("[K01677]"%in%KOs&"[K01678]"%in%KOs)&("[K18209]"%in%KOs&"[K18210]"%in%KOs)&("[K01902]"%in%KOs&"[K01903]"%in%KOs)&("[K00174]"%in%KOs&"[K00175]"%in%KOs&"[K00176]"%in%KOs&"[K00177]"%in%KOs)
    )
  
  return(Res_DF)
}

#####LOAD DATA#####

#Counts Transcriptome
Read_Counting_RNA <- featureCounts(files = c("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH02_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH03_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH07_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH08_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH15_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH16_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH19_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH20_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH21_NO_rRNA.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MT_FOMES_MAPPED_READS/VBRNAH22_NO_rRNA.sorted.bam"),
                                    annot.ext = "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_metatrans_Trinity_contigs_fgs.gff",
                                    isGTFAnnotationFile = TRUE,GTF.featureType = "CDS",isPairedEnd = FALSE, nthreads = 6)

rawCounts_RNA <- as.data.frame(Read_Counting_RNA$counts)
#remove the unnecessary part of names and simplify
colnames(rawCounts_RNA) <- gsub("_NO_rRNA.sorted.bam","",colnames(rawCounts_RNA))
colnames(rawCounts_RNA) <- gsub("VBRNA","RNA_",colnames(rawCounts_RNA))

#Import the sample metadata
sampleData_RNA <- read.csv("~/PostDoc/02_Projects/05_Fomes_transcription/06_Miscellaneous/RNA_Metadata.csv", header = T)
sampleData_RNA$Condition <- as.factor(sampleData_RNA$Condition)
sampleData_RNA$Sample <- gsub("VBRNA","RNA_",sampleData_RNA$Sample)

#Import gene taxonomy
GeneTaxonomy_RNA <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_MG_MT_TABLES/FOMES_METATRANS_TAX_CAZy.tab",header = TRUE, sep = "\t", comment.char = "", quote = "")
GeneTaxonomy_RNA <- GeneTaxonomy_RNA[,c("ctg_name","ctg_length","bitscore","TAX_BEST","e.val","CAZy")]
colnames(GeneTaxonomy_RNA) <- c("Gene","Transcript Length","bitscore","TAX_Key","e.val","CAZ_model")

#Counts Genomes
Read_Counting_DNA <- featureCounts(files = c("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH02.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH03.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH07.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH08.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH15.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH16.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH19.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH20.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH21.sorted.bam",
                                              "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_FOMES_MAPPED_READS/VBDNAH22.sorted.bam"),
                                    annot.ext = "~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_metagenome_contigs_fgs.gff",
                                    isGTFAnnotationFile = TRUE,GTF.featureType = "CDS",isPairedEnd = FALSE, nthreads = 6)

rawCounts_DNA <- as.data.frame(Read_Counting_DNA$counts)
#remove the unnecessary part of names and simplify
colnames(rawCounts_DNA) <- gsub(".sorted.bam","",colnames(rawCounts_DNA))
colnames(rawCounts_DNA) <- gsub("VBDNA","DNA_",colnames(rawCounts_DNA))

#Import the sample metadata
sampleData_DNA <- sampleData_RNA
sampleData_DNA$Sample <- gsub("RNA_","DNA_",sampleData_DNA$Sample)

#Import gene taxonomy
GeneTaxonomy_DNA <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_Annotation/FOMES_METAGENOME_TAX_CAZy_KOFAM.tab",sep = "\t", header = T, quote = "", comment.char = "")
GeneTaxonomy_DNA <- GeneTaxonomy_DNA[,c("ctg_name","ctg_length","bitscore","TAX_BEST","e.val","CAZy","KOFAM")]
colnames(GeneTaxonomy_DNA) <- c("Gene","Transcript Length","bitscore","TAX_Key","e.val","CAZ_model","KOFAM")
#Remove the non-CAZyme entries from CAZ_model
GeneTaxonomy_DNA$CAZ_model <- gsub("SLH","-",GeneTaxonomy_DNA$CAZ_model)
GeneTaxonomy_DNA$CAZ_model <- gsub("cohesin","-",GeneTaxonomy_DNA$CAZ_model)

#Load taxonomies of the trees, CAZymes and KOFAMs
Taxonomy_Organisms <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_MG_MT_TABLES/MG_MT_taxonomy_tree.tab",sep = "\t", header = T, quote = "", comment.char = "")
Taxonomy_CAZymes <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_MG_MT_TABLES/MG_MT_CAZy_tree.tab",sep = "\t", header = T)
Taxonomy_KOFAM <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/MG_Annotation/KOFAM_MT_MG_tree.tab",sep = "\t", header = T, quote = "", comment.char = "")
#Clean them up for better linking
colnames(Taxonomy_Organisms) <- c("TAX_Domain", "TAX_Phylum", "TAX_Class", "TAX_Order", "TAX_Family", "TAX_Genus", "TAX_Key")
colnames(Taxonomy_CAZymes) <- c("CAZ_Family", "CAZ_Subfamily", "CAZ_model")
colnames(Taxonomy_KOFAM) <- c("KO_Level_1","KO_Level_2","KO_Level_3","KO_Function","KOFAM")
#Add entries for unidentified genes
Taxonomy_Organisms[nrow(Taxonomy_Organisms)+1,] <- c( c("Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "-"))
Taxonomy_CAZymes[nrow(Taxonomy_CAZymes)+1,] <- c("-", "-", "-")
Taxonomy_KOFAM[nrow(Taxonomy_KOFAM)+1,] <- c("-", "-", "-", "-", "-")

#Load Ruben's CAZyme table
Ruben_CAZymes <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/06_Miscellaneous/Classification_per_target_substrate.csv",header = TRUE,sep = ",")
#Correct old label
Ruben_CAZymes$activity <- gsub("glycoprotein/glucosaminglycan","Glycoconjugates",Ruben_CAZymes$activity)
Ruben_CAZymes$activity <- gsub("beta glucans","beta-glucans",Ruben_CAZymes$activity)
Ruben_CAZymes$activity <- gsub("alpha glucans","alpha-glucans",Ruben_CAZymes$activity)
#Change all to use capitals
Ruben_CAZymes$activity <- str_to_sentence(Ruben_CAZymes$activity)

#####Analysis#####

###Set colours###

condition_colours <- c("#66c2a5","#fc8d62")
names(condition_colours) <- c("Fresh","Rotten")

###Pre-Analysis Quality Check###

#Create a dataframe of read counts
rawCounts_RNA_DF <- as.data.frame(colSums(rawCounts_RNA))
rawCounts_DNA_DF <- as.data.frame(colSums(rawCounts_DNA))
colnames(rawCounts_RNA_DF)[1] <- "Reads"
colnames(rawCounts_DNA_DF)[1] <- "Reads"

#Combine with sample data
rawCounts_RNA_DF <- cbind(rawCounts_RNA_DF, sampleData_RNA)
rawCounts_DNA_DF <- cbind(rawCounts_DNA_DF, sampleData_DNA)

#Plot read numbers per sample
rawCounts_RNA_Plot <- ggplot(rawCounts_RNA_DF, aes(x = Sample, y = Reads, fill = Condition)) +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(n.breaks = 7, labels = function(x) format(x, scientific = TRUE)) +
  labs(x=NULL, y = "Reads", title = "Read Counts: Metatranscriptome") + 
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
        legend.position = "none", plot.title = element_text(hjust = 0.5))
rawCounts_DNA_Plot <- ggplot(rawCounts_DNA_DF, aes(x = Sample, y = Reads, fill = Condition)) +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(n.breaks = 7, labels = function(x) format(x, scientific = TRUE)) +
  labs(x=NULL, y = "Reads", title = "Read Counts: Metagenome") + 
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
        legend.position = "none", plot.title = element_text(hjust = 0.5))

#Save plots
rawCounts_Plot <- plot_grid(rawCounts_RNA_Plot + labs(tag = "A"), rawCounts_DNA_Plot + labs(tag = "B"), ncol = 2)
ggsave("rawCounts_Plot.svg",rawCounts_Plot,width = 12, height = 4)

##Calculate TPM##
#https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/

#Gene lengths
GFF_RNA <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_metatrans_Trinity_contigs_fgs.gff", header = F, sep="\t")
GFF_RNA$Gene <- gsub("gene_id=","",GFF_RNA$V9)
GFF_RNA$Length <- GFF_RNA$V5 - GFF_RNA$V4
GFF_DNA <- read.table("~/PostDoc/02_Projects/05_Fomes_transcription/02_Cleaned_Data/Fomes_metagenome_contigs_fgs.gff", header = F, sep="\t")
GFF_DNA$Gene <- gsub("gene_id=","",GFF_DNA$V9)
GFF_DNA$Length <- GFF_DNA$V5 - GFF_DNA$V4

#Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
RPK_RNA <- rawCounts_RNA
for (COL in colnames(RPK_RNA)) {
  RPK_RNA[,COL] <- RPK_RNA[,COL]/(GFF_RNA$Length/1000)
}
RPK_DNA <- rawCounts_DNA
for (COL in colnames(RPK_DNA)) {
  RPK_DNA[,COL] <- RPK_DNA[,COL]/(GFF_DNA$Length/1000)
}

#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
scaling_factor_RNA <- colSums(RPK_RNA)/1000000
scaling_factor_DNA <- colSums(RPK_DNA)/1000000

#Divide the RPK values by the “per million” scaling factor. This gives you TPM.
TPM_RNA <- RPK_RNA
for (n in 1:ncol(TPM_RNA)) {
  TPM_RNA[,n] <- TPM_RNA[,n]/scaling_factor_RNA[n]
}
TPM_DNA <- RPK_DNA
for (n in 1:ncol(TPM_DNA)) {
  TPM_DNA[,n] <- TPM_DNA[,n]/scaling_factor_DNA[n]
}

##Test sample clustering##

#Determine euclidean distances of samples
sampleDists_RNA <- dist(t(TPM_RNA))
sampleDists_DNA <- dist(t(TPM_DNA))

#Make a labelled matrix of distances
sampleDistsMatrix_RNA <- as.matrix(sampleDists_RNA)
rownames(sampleDistsMatrix_RNA) <- paste(sampleData_RNA$Condition, sampleData_RNA$tree, sep = "_" )
colnames(sampleDistsMatrix_RNA) <- sampleData_RNA$Sample
sampleDistsMatrix_DNA <- as.matrix(sampleDists_DNA)
rownames(sampleDistsMatrix_DNA) <- paste(sampleData_DNA$Condition, sampleData_DNA$tree, sep = "_" )
colnames(sampleDistsMatrix_DNA) <- sampleData_DNA$Sample

#Show the distances on a heatmap
heatmap_colours <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
Heatmap_TPM_Euclidian_RNA <- 
  as.ggplot(pheatmap(sampleDistsMatrix_RNA,
                                      clustering_distance_rows = sampleDists_RNA,
                                      clustering_distance_cols = sampleDists_RNA,
                                      col = heatmap_colours)) +
  labs(title = "Metatranscriptome heatmap") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))
Heatmap_TPM_Euclidian_DNA <- 
  as.ggplot(pheatmap(sampleDistsMatrix_DNA,
                                      clustering_distance_rows = sampleDists_DNA,
                                      clustering_distance_cols = sampleDists_DNA,
                                      col = heatmap_colours)) +
  labs(title = "Metagenome heatmap") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))

#Save plots
Heatmap_Plots <- plot_grid(as.ggplot(Heatmap_TPM_Euclidian_RNA) + labs(tag = "A"), as.ggplot(Heatmap_TPM_Euclidian_DNA) + labs(tag = "B"), ncol = 2)
ggsave("Heatmap_Plots.svg",Heatmap_Plots,width = 12, height = 4)

#Show samples on a PCA
pca_res_RNA <- prcomp(TPM_RNA, scale. = TRUE)
dtp_RNA <- data.frame('Condition' = as.character(sampleData_RNA$Condition),'Tree' = as.character(sampleData_RNA$tree), pca_res_RNA$rotation[,1:2]) 
PCAloadings_RNA <- data.frame(Variables = rownames(pca_res_RNA$rotation), pca_res_RNA$rotation)
percentVar_RNA <- round(100 * summary(pca_res_RNA)$importance[2,1:2],digits = 2)
pca_res_DNA <- prcomp(TPM_DNA, scale. = TRUE)
dtp_DNA <- data.frame('Condition' = as.character(sampleData_DNA$Condition),'Tree' = as.character(sampleData_DNA$tree), pca_res_DNA$rotation[,1:2]) 
PCAloadings_DNA <- data.frame(Variables = rownames(pca_res_DNA$rotation), pca_res_DNA$rotation)
percentVar_DNA <- round(100 * summary(pca_res_DNA)$importance[2,1:2],digits = 2)

PCA_RNA <- 
  ggplot(data = dtp_RNA) + 
  geom_point(aes(x = PC1, y = PC2, col = Condition)) + 
  xlab(paste0("PC1: ", percentVar_RNA[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_RNA[2], "% variance")) +
  scale_y_continuous(n.breaks = 7) + 
  scale_x_continuous(n.breaks = 7) +
  geom_text_repel(aes(x = PC1, y = PC2, label = row.names(dtp_RNA)), max.overlaps = 15, nudge_y = -0.1, force = 2) +
  ggtitle("PCA of TPM-transformed Metatranscriptomes") +
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))
PCA_DNA <- 
  ggplot(data = dtp_DNA) + 
  geom_point(aes(x = PC1, y = PC2, col = Condition)) + 
  xlab(paste0("PC1: ", percentVar_DNA[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_DNA[2], "% variance")) +
  scale_y_continuous(n.breaks = 7) + 
  scale_x_continuous(n.breaks = 7) +
  geom_text_repel(aes(x = PC1, y = PC2, label = row.names(dtp_DNA)), max.overlaps = 15, force = 2) +
  ggtitle("PCA of TPM-transformed Metagenomes") +
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))

#Save plots
PCA_Plots <- plot_grid(PCA_RNA + labs(tag = "A"), PCA_DNA + labs(tag = "B"), ncol = 2)
ggsave("PCA_Plots.svg",PCA_Plots,width = 12, height = 4)

#Perform clustering analysis

#Same distances calculated previously for heatmap

#Find the best Linkage Method to use for clustering 
#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
#function to compute agglomerative coefficient
ac_RNA <- function(x) {
  agnes(x = sampleDists_RNA, diss = TRUE, method = x)$ac
}
ac_DNA <- function(x) {
  agnes(x = sampleDists_DNA, diss = TRUE, method = x)$ac
}
#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac_RNA)
sapply(m, ac_DNA)
#Ward has highest ac in both cases, so best method to use

#Cluster samples
final_clust_RNA <- hclust(sampleDists_RNA, method = "ward.D2")
final_clust_DNA <- hclust(sampleDists_DNA, method = "ward.D2")

#Perform bootstrapping of clusters
set.seed(123456)
Bootstrap_clust_RNA <- Bclust(t(TPM_RNA), method.d = "euclidean", method.c = "ward.D2", iter=1000, mc.cores = 7)
set.seed(123456)
#Needs at least 19 GB free ram. (~4.2  GB per core at max)
Bootstrap_clust_DNA <- Bclust(t(TPM_DNA), method.d = "euclidean", method.c = "ward.D2", iter=1000, mc.cores = 7)

#Colour clusters by vector is not officially supported but seems to work for now.
clust_colours_RNA <- condition_colours[sampleData_RNA$Condition][final_clust_RNA$order]
clust_colours_DNA <- condition_colours[sampleData_DNA$Condition][final_clust_DNA$order]

#Get labels and co-ordinates 
#Not how this is meant to be used but it seems to work
plot(Bootstrap_clust_RNA)
Bootstrap_labels_RNA <- Bclabels(hcl = final_clust_RNA,values = Bootstrap_clust_RNA$values)
plot(Bootstrap_clust_DNA)
Bootstrap_labels_DNA <- Bclabels(hcl = final_clust_DNA,values = Bootstrap_clust_DNA$values)

#Final set up
dendro_clust_RNA <- dendro_data(final_clust_RNA)
dendro_clust_DNA <- dendro_data(final_clust_DNA)
dendro_clust_RNA_labels <- as.data.frame(cbind(Bootstrap_labels_RNA$coords[,"x"],cbind(Bootstrap_labels_RNA$coords[,"y"],Bootstrap_labels_RNA$labels)))
colnames(dendro_clust_RNA_labels) <- c("x","y","label")
dendro_clust_DNA_labels <- as.data.frame(cbind(Bootstrap_labels_DNA$coords[,"x"],cbind(Bootstrap_labels_DNA$coords[,"y"],Bootstrap_labels_DNA$labels)))
colnames(dendro_clust_DNA_labels) <- c("x","y","label")

#Plot clustering
Cluster_plot_RNA <- 
  ggdendrogram(final_clust_RNA) + 
  geom_text(data = dendro_clust_RNA_labels, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 4, vjust = 1.1) +
  labs(y = "", x = "", title = "Metatranscriptomic clustering") +
  scale_y_continuous(breaks = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = clust_colours_RNA),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5))
Cluster_plot_DNA <- 
  ggdendrogram(final_clust_DNA) + 
  geom_text(data = dendro_clust_DNA_labels, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 4, vjust = 1.1) +
  labs(y = "", x = "", title = "Metagenomic clustering") +
  scale_y_continuous(breaks = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = clust_colours_DNA),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5))

#Save plots
Cluster_plots <- plot_grid(Cluster_plot_RNA + labs(tag = "A"), Cluster_plot_DNA + labs(tag = "B"), ncol=2)
ggsave(filename = paste("Cluster_plots.svg",sep=""), Cluster_plots, width = 20, height = 6)



#The data does not match what is expected and has quality problems.

##Fixing problems with data##

#Sample VBRNAH22 will be excluded from the RNA data due to low read count.
rawCounts_RNA_FIXED <- subset(rawCounts_RNA,select = -RNA_H22)
TPM_RNA_FIXED <- subset(TPM_RNA,select = -RNA_H22)
sampleData_RNA_FIXED <- sampleData_RNA[!sampleData_RNA$Sample=="RNA_H22",]
rawCounts_RNA_DF_FIXED <- rawCounts_RNA_DF[!rownames(rawCounts_RNA_DF)=="RNA_H22",]

#Sample VBRNAH08 will be recoded as a fresh sample.
sampleData_RNA_FIXED[sampleData_RNA_FIXED$Sample=="RNA_H08","Condition"] <- "Fresh"
rawCounts_RNA_DF_FIXED[rawCounts_RNA_DF_FIXED$Sample=="RNA_H08","Condition"] <- "Fresh"

#Recode 8 and 22 as fresh for DNA
sampleData_DNA_FIXED <- sampleData_DNA
sampleData_DNA_FIXED[sampleData_DNA_FIXED$Sample=="DNA_H08","Condition"] <- "Fresh"
sampleData_DNA_FIXED[sampleData_DNA_FIXED$Sample=="DNA_H22","Condition"] <- "Fresh"
rawCounts_DNA_DF_FIXED <- rawCounts_DNA_DF
rawCounts_DNA_DF_FIXED$Condition <- sampleData_DNA_FIXED$Condition

#Figure with the updated classifications

#Plot read numbers per fixed sample
rawCounts_RNA_Plot_FIXED <- ggplot(rawCounts_RNA_DF_FIXED, aes(x = Sample, y = Reads, fill = Condition)) +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(n.breaks = 7, labels = function(x) format(x, scientific = TRUE)) +
  labs(x=NULL, y = "Reads", title = "Read Counts: Metatranscriptome") + 
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
        legend.position = "none", plot.title = element_text(hjust = 0.5))
rawCounts_DNA_Plot_FIXED <- ggplot(rawCounts_DNA_DF_FIXED, aes(x = Sample, y = Reads, fill = Condition)) +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(n.breaks = 7, labels = function(x) format(x, scientific = TRUE)) +
  labs(x=NULL, y = "Reads", title = "Read Counts: Metagenome") + 
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
        legend.position = "none", plot.title = element_text(hjust = 0.5))

##Fixed sample clustering##

#Determine euclidean distances of samples
sampleDists_RNA_FIXED <- dist(t(TPM_RNA_FIXED))
sampleDists_DNA_FIXED <- dist(t(TPM_DNA))

#Make a labelled matrix of distances
sampleDistsMatrix_RNA_FIXED <- as.matrix(sampleDists_RNA_FIXED)
rownames(sampleDistsMatrix_RNA_FIXED) <- paste(sampleData_RNA_FIXED$Condition, sampleData_RNA_FIXED$tree, sep = "_" )
colnames(sampleDistsMatrix_RNA_FIXED) <- sampleData_RNA_FIXED$Sample
sampleDistsMatrix_DNA_FIXED <- as.matrix(sampleDists_DNA_FIXED)
rownames(sampleDistsMatrix_DNA_FIXED) <- paste(sampleData_DNA_FIXED$Condition, sampleData_DNA_FIXED$tree, sep = "_" )
colnames(sampleDistsMatrix_DNA_FIXED) <- sampleData_DNA_FIXED$Sample

#Show the distances on a heatmap
Heatmap_TPM_Euclidian_RNA_FIXED <- 
  as.ggplot(pheatmap(sampleDistsMatrix_RNA_FIXED,
                     clustering_distance_rows = sampleDists_RNA_FIXED,
                     clustering_distance_cols = sampleDists_RNA_FIXED,
                     col = heatmap_colours)) +
  labs(title = "Metatranscriptome heatmap") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))
Heatmap_TPM_Euclidian_DNA_FIXED <- 
  as.ggplot(pheatmap(sampleDistsMatrix_DNA_FIXED,
                     clustering_distance_rows = sampleDists_DNA_FIXED,
                     clustering_distance_cols = sampleDists_DNA_FIXED,
                     col = heatmap_colours)) +
  labs(title = "Metagenome heatmap") +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))

#Show fixed samples on a PCA
pca_res_RNA_FIXED <- prcomp(TPM_RNA_FIXED, scale. = TRUE)
dtp_RNA_FIXED <- data.frame('Condition' = as.character(sampleData_RNA_FIXED$Condition),'Tree' = as.character(sampleData_RNA_FIXED$tree), pca_res_RNA_FIXED$rotation[,1:2]) 
PCAloadings_RNA_FIXED <- data.frame(Variables = rownames(pca_res_RNA_FIXED$rotation), pca_res_RNA_FIXED$rotation)
percentVar_RNA_FIXED <- round(100 * summary(pca_res_RNA_FIXED)$importance[2,1:2],digits = 2)
pca_res_DNA_FIXED <- prcomp(TPM_DNA, scale. = TRUE)
dtp_DNA_FIXED <- data.frame('Condition' = as.character(sampleData_DNA_FIXED$Condition),'Tree' = as.character(sampleData_DNA_FIXED$tree), pca_res_DNA_FIXED$rotation[,1:2]) 
PCAloadings_DNA_FIXED <- data.frame(Variables = rownames(pca_res_DNA_FIXED$rotation), pca_res_DNA_FIXED$rotation)
percentVar_DNA_FIXED <- round(100 * summary(pca_res_DNA_FIXED)$importance[2,1:2],digits = 2)

PCA_RNA_FIXED <- 
  ggplot(data = dtp_RNA_FIXED) + 
  geom_point(aes(x = PC1, y = PC2, col = Condition)) + 
  xlab(paste0("PC1: ", percentVar_RNA[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_RNA[2], "% variance")) +
  scale_y_continuous(n.breaks = 7) + 
  scale_x_continuous(n.breaks = 7) +
  geom_text_repel(aes(x = PC1, y = PC2, label = row.names(dtp_RNA_FIXED)), max.overlaps = 15, nudge_y = -0.1, force = 2) +
  ggtitle("PCA of TPM-transformed Metatranscriptomes") +
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))
PCA_DNA_FIXED <- 
  ggplot(data = dtp_DNA_FIXED) + 
  geom_point(aes(x = PC1, y = PC2, col = Condition)) + 
  xlab(paste0("PC1: ", percentVar_DNA[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_DNA[2], "% variance")) +
  scale_y_continuous(n.breaks = 7) + 
  scale_x_continuous(n.breaks = 7) +
  geom_text_repel(aes(x = PC1, y = PC2, label = row.names(dtp_DNA_FIXED)), max.overlaps = 15, force = 2) +
  ggtitle("PCA of TPM-transformed Metagenomes") +
  scale_fill_manual(values = condition_colours) +
  theme_classic() +
  theme(text=element_text(size=16), plot.title = element_text(hjust = 0.5))

#Perform clustering analysis

#Find the best Linkage Method to use for clustering 
#Only the RNA distances changed
#function to compute agglomerative coefficient
ac_RNA_FIXED <- function(x) {
  agnes(x = sampleDists_RNA_FIXED, diss = TRUE, method = x)$ac
}
#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac_RNA_FIXED)
#Ward is still the best method

#Cluster samples
final_clust_RNA_FIXED <- hclust(sampleDists_RNA_FIXED, method = "ward.D2")
final_clust_DNA_FIXED <- hclust(sampleDists_DNA_FIXED, method = "ward.D2")

#Perform bootstrapping of clusters
set.seed(123456)
Bootstrap_clust_RNA_FIXED <- Bclust(t(TPM_RNA_FIXED), method.d = "euclidean", method.c = "ward.D2", iter=1000, mc.cores = 7)
set.seed(123456)
#Needs ~19 GB of free RAM
Bootstrap_clust_DNA_FIXED <- Bclust(t(TPM_DNA), method.d = "euclidean", method.c = "ward.D2", iter=1000, mc.cores = 7)

#Colour clusters by vector is not officially supported but seems to work for now.
clust_colours_RNA_FIXED <- condition_colours[sampleData_RNA_FIXED$Condition][final_clust_RNA_FIXED$order]
clust_colours_DNA_FIXED <- condition_colours[sampleData_DNA_FIXED$Condition][final_clust_DNA_FIXED$order]

#Get labels and co-ordinates 
#Not how this is meant to be used but it seems to work
plot(Bootstrap_clust_RNA_FIXED)
Bootstrap_labels_RNA_FIXED <- Bclabels(hcl = final_clust_RNA_FIXED,values = Bootstrap_clust_RNA_FIXED$values)
plot(Bootstrap_clust_DNA_FIXED)
Bootstrap_labels_DNA_FIXED <- Bclabels(hcl = final_clust_DNA_FIXED,values = Bootstrap_clust_DNA_FIXED$values)

#Final set up
dendro_clust_RNA_FIXED <- dendro_data(final_clust_RNA_FIXED)
dendro_clust_DNA_FIXED <- dendro_data(final_clust_DNA_FIXED)
dendro_clust_RNA_labels_FIXED <- as.data.frame(cbind(Bootstrap_labels_RNA_FIXED$coords[,"x"],cbind(Bootstrap_labels_RNA_FIXED$coords[,"y"],Bootstrap_labels_RNA_FIXED$labels)))
colnames(dendro_clust_RNA_labels_FIXED) <- c("x","y","label")
dendro_clust_DNA_labels_FIXED <- as.data.frame(cbind(Bootstrap_labels_DNA_FIXED$coords[,"x"],cbind(Bootstrap_labels_DNA_FIXED$coords[,"y"],Bootstrap_labels_DNA_FIXED$labels)))
colnames(dendro_clust_DNA_labels_FIXED) <- c("x","y","label")

#Plot clustering
Cluster_plot_RNA_FIXED <- 
  ggdendrogram(final_clust_RNA_FIXED) + 
  geom_text(data = dendro_clust_RNA_labels_FIXED, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 4, vjust = 1.1) +
  labs(y = "", x = "", title = "Metatranscriptomic clustering") +
  scale_y_continuous(breaks = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = clust_colours_RNA_FIXED),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5))
Cluster_plot_DNA_FIXED <- 
  ggdendrogram(final_clust_DNA_FIXED) + 
  geom_text(data = dendro_clust_DNA_labels_FIXED, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 4, vjust = 1.1) +
  labs(y = "", x = "", title = "Metagenomic clustering") +
  scale_y_continuous(breaks = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = clust_colours_DNA_FIXED),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5))

#Figure S1
Figure_S1 <-
  plot_grid(rawCounts_DNA_Plot_FIXED + labs(tag = "A"), rawCounts_RNA_Plot_FIXED + labs(tag = "B"),
            PCA_DNA_FIXED + labs(tag = "C"), PCA_RNA_FIXED + labs(tag = "D"),
            as.ggplot(Heatmap_TPM_Euclidian_DNA_FIXED) + labs(tag = "E"), as.ggplot(Heatmap_TPM_Euclidian_RNA_FIXED) + labs(tag = "F"),
            Cluster_plot_DNA_FIXED + labs(tag = "G"), Cluster_plot_RNA_FIXED + labs(tag = "H"),
            ncol = 2)
ggsave("Figure_S1.eps",Figure_S1,width = 2*168, height = 2*237, units = "mm")

###Main Analysis###

##Create tidy tables

#RNA
Tidy_Taxonomy_RNA <- data.frame()
for (SAMPLENAME in colnames(TPM_RNA_FIXED)) {
  #Sets up the tidy dataframe and additionally removes any sample with 0 TPM.
  #This is necessary due to the much larger table size than the RNA.
  #SAMPLE, GENE, TPM, TAX_Domain, TAX_Phylum, TAX_Class, TAX_Order, TAX_Family, TAX_Genus, CAZ_Family, CAZ_Subfamily, CAZ_model
  cat(paste ("Performing ", which(colnames(TPM_RNA_FIXED)==SAMPLENAME)," of ",length(colnames(TPM_RNA_FIXED)),"\n",sep = ""))
  WIP <- data.frame(Sample = rep(SAMPLENAME,nrow(TPM_RNA_FIXED)), Condition = rep(sampleData_RNA_FIXED[sampleData_RNA_FIXED$Sample==SAMPLENAME,"Condition"],nrow(TPM_RNA_FIXED)), 
                    Gene = row.names(TPM_RNA_FIXED), TPM = TPM_RNA_FIXED[,SAMPLENAME])
  WIP <- WIP[WIP$TPM > 0,]
  Tidy_Taxonomy_RNA <- rbind(Tidy_Taxonomy_RNA,WIP)
}
#Add in the keys for the predicted taxonomy and CAZymes.
Tidy_Taxonomy_RNA <- merge.data.frame(Tidy_Taxonomy_RNA,GeneTaxonomy_RNA,by="Gene")
#Expands the taxonomy
Tidy_Taxonomy_RNA <- merge.data.frame(Tidy_Taxonomy_RNA,Taxonomy_Organisms,by="TAX_Key")
#Expands the CAZyme classification
Tidy_Taxonomy_RNA <- merge.data.frame(Tidy_Taxonomy_RNA,Taxonomy_CAZymes,by="CAZ_model")
#Add CAZyme target according to Ruben's information
for (CAZy in unique(Tidy_Taxonomy_RNA$CAZ_model)) {
  cat(paste ("Performing ", which(unique(Tidy_Taxonomy_RNA$CAZ_model)==CAZy)," of ",length(unique(Tidy_Taxonomy_RNA$CAZ_model)),"\n",sep = ""))
  if (CAZy == "-") {
    Tidy_Taxonomy_RNA[which(Tidy_Taxonomy_RNA$CAZ_model==CAZy),"CAZ_target"] <- "-"
  }
  else if (CAZy%in%Ruben_CAZymes$cazy) {
    Tidy_Taxonomy_RNA[which(Tidy_Taxonomy_RNA$CAZ_model==CAZy),"CAZ_target"] <-
      paste(Ruben_CAZymes[Ruben_CAZymes$cazy==CAZy,"activity"], collapse = " / ")
  }
  else {
    Tidy_Taxonomy_RNA[which(Tidy_Taxonomy_RNA$CAZ_model==CAZy),"CAZ_target"] <- "Unknown"
  }
}
Tidy_Taxonomy_RNA <- Tidy_Taxonomy_RNA[,c("Gene","Sample","Condition","TPM","Transcript Length","TAX_Key","TAX_Domain","TAX_Phylum",
                                          "TAX_Class","TAX_Order","TAX_Family","TAX_Genus","CAZ_Family","CAZ_Subfamily","CAZ_model","CAZ_target")]
#DNA
Tidy_Taxonomy_DNA <- data.frame()
for (SAMPLENAME in colnames(TPM_DNA)) {
  #Sets up the tidy dataframe and additionally removes any sample with 0 TPM.
  #This is necessary due to the much larger table size than the RNA.
  #SAMPLE, GENE, TPM, TAX_Domain, TAX_Phylum, TAX_Class, TAX_Order, TAX_Family, TAX_Genus, CAZ_Family, CAZ_Subfamily, CAZ_model
  cat(paste ("Performing ", which(colnames(TPM_DNA)==SAMPLENAME)," of ",length(colnames(TPM_DNA)),"\n",sep = ""))
  WIP <- data.frame(Sample = rep(SAMPLENAME,nrow(TPM_DNA)), Condition = rep(sampleData_DNA_FIXED[sampleData_DNA_FIXED$Sample==SAMPLENAME,"Condition"],nrow(TPM_DNA)), 
                    Gene = row.names(TPM_DNA), TPM = TPM_DNA[,SAMPLENAME])
  WIP <- WIP[WIP$TPM > 0,]
  Tidy_Taxonomy_DNA <- rbind(Tidy_Taxonomy_DNA,WIP)
}
#Add in the keys for the predicted taxonomy and CAZymes.
Tidy_Taxonomy_DNA <- merge.data.frame(Tidy_Taxonomy_DNA,GeneTaxonomy_DNA,by="Gene")
#Expands the taxonomy
Tidy_Taxonomy_DNA <- merge.data.frame(Tidy_Taxonomy_DNA,Taxonomy_Organisms,by="TAX_Key")
#Expands the CAZyme classification
Tidy_Taxonomy_DNA <- merge.data.frame(Tidy_Taxonomy_DNA,Taxonomy_CAZymes,by="CAZ_model")
# #Expands the KOFAM classification
# #Some KOs have multiple functions assigned, this disrupts the correct TPM counting and can not be automatically added at this stage
# Tidy_Taxonomy_DNA <- merge.data.frame(Tidy_Taxonomy_DNA,Taxonomy_KOFAM,by="KOFAM")
#Add CAZyme target according to Ruben's information
for (CAZy in unique(Tidy_Taxonomy_DNA$CAZ_model)) {
  cat(paste ("Performing ", which(unique(Tidy_Taxonomy_DNA$CAZ_model)==CAZy)," of ",length(unique(Tidy_Taxonomy_DNA$CAZ_model)),"\n",sep = ""))
  if (CAZy == "-") {
    Tidy_Taxonomy_DNA[which(Tidy_Taxonomy_DNA$CAZ_model==CAZy),"CAZ_target"] <- "-"
  }
  else if (CAZy%in%Ruben_CAZymes$cazy) {
    Tidy_Taxonomy_DNA[which(Tidy_Taxonomy_DNA$CAZ_model==CAZy),"CAZ_target"] <-
      paste(Ruben_CAZymes[Ruben_CAZymes$cazy==CAZy,"activity"], collapse = " / ")
  }
  else {
    Tidy_Taxonomy_DNA[which(Tidy_Taxonomy_DNA$CAZ_model==CAZy),"CAZ_target"] <- "Unknown"
  }
}
Tidy_Taxonomy_DNA <- Tidy_Taxonomy_DNA[,c("Gene","Sample","Condition","TPM","Transcript Length","TAX_Key","TAX_Domain","TAX_Phylum",
                                          "TAX_Class","TAX_Order","TAX_Family","TAX_Genus","CAZ_Family","CAZ_Subfamily","CAZ_model","CAZ_target",
                                          "KOFAM")]

##Percentage Fomes##

#Percent Fomes reads in Fresh vs Rotten
GeneTaxonomy_RNA_Clean <- GeneTaxonomy_RNA
GeneTaxonomy_RNA_Clean$Genus <- gsub("\\[Genome_Fomes\\]","\\[Fomes",GeneTaxonomy_RNA_Clean$TAX_Key)
GeneTaxonomy_RNA_Clean$Genus <- gsub("\\[Transcriptome_Fomes\\]","\\[Fomes",GeneTaxonomy_RNA_Clean$Genus)
GeneTaxonomy_RNA_Clean$Genus <- gsub("Candidatus ","Candidatus_",GeneTaxonomy_RNA_Clean$Genus)
GeneTaxonomy_RNA_Clean$Genus <- gsub("-","\\[Unidentified",GeneTaxonomy_RNA_Clean$Genus)
GeneTaxonomy_RNA_Clean$Genus <- gsub("\\'","",GeneTaxonomy_RNA_Clean$Genus)
GeneTaxonomy_RNA_Clean$Genus <- gsub("\\[","",str_extract(GeneTaxonomy_RNA_Clean$Genus,"\\[[A-z]+"))

rawCounts_RNA_DF_FIXED <- as.data.frame(colSums(rawCounts_RNA_FIXED))
colnames(rawCounts_RNA_DF_FIXED)[1] <- "Reads"
rawCounts_RNA_DF_FIXED$Fomes <- colSums(rawCounts_RNA_FIXED[GeneTaxonomy_RNA_Clean[GeneTaxonomy_RNA_Clean$Genus=="Fomes","Gene"],])
rawCounts_RNA_DF_FIXED$Fomes_TPM <- unlist(aggregate(Tidy_Taxonomy_RNA[Tidy_Taxonomy_RNA$TAX_Genus=="Fomes",c("Sample","TPM")], .~Sample,sum)["TPM"])
rawCounts_RNA_DF_FIXED$Fomes_Percent <- rawCounts_RNA_DF_FIXED$Fomes/rawCounts_RNA_DF_FIXED$Reads*100
rawCounts_RNA_DF_FIXED$Fomes_Percent_TPM <- rawCounts_RNA_DF_FIXED$Fomes_TPM/(1e+06)*100
rawCounts_RNA_DF_FIXED <- cbind(rawCounts_RNA_DF_FIXED, sampleData_RNA_FIXED)

RawCountPlot_RNA_FIXED <- ggplot(rawCounts_RNA_DF_FIXED, aes(x = Sample, y = Reads, fill = Condition)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None") +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(n.breaks = 7) +
  labs(x=NULL, y = "Reads") + 
  scale_fill_manual(values = condition_colours)
RawCountPlot_RNA_FIXED_Fomes <- ggplot(rawCounts_RNA_DF_FIXED, aes(x = Sample, y = Fomes/Reads*100, fill = Condition)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None") +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70)) +
  labs(x=NULL, y = "Percent Fomes Reads") + 
  scale_fill_manual(values = condition_colours)
RawCountPlot_RNA_FIXED <- plot_grid(RawCountPlot_RNA_FIXED + labs(tag = "A"), RawCountPlot_RNA_FIXED_Fomes + labs(tag = "B"), ncol = 2)
ggsave("RawCountPlot_RNA_FIXED.svg",RawCountPlot_RNA_FIXED,width = 12, height = 4)

GeneTaxonomy_DNA_Clean <- GeneTaxonomy_DNA
GeneTaxonomy_DNA_Clean$Genus <- gsub("\\[Genome_Fomes\\]","\\[Fomes",GeneTaxonomy_DNA_Clean$TAX_Key)
GeneTaxonomy_DNA_Clean$Genus <- gsub("\\[Transcriptome_Fomes\\]","\\[Fomes",GeneTaxonomy_DNA_Clean$Genus)
GeneTaxonomy_DNA_Clean$Genus <- gsub("Candidatus ","Candidatus_",GeneTaxonomy_DNA_Clean$Genus)
GeneTaxonomy_DNA_Clean$Genus <- gsub("-","\\[Unidentified",GeneTaxonomy_DNA_Clean$Genus)
GeneTaxonomy_DNA_Clean$Genus <- gsub("\\'","",GeneTaxonomy_DNA_Clean$Genus)
GeneTaxonomy_DNA_Clean$Genus <- gsub("\\[","",str_extract(GeneTaxonomy_DNA_Clean$Genus,"\\[[A-z]+"))

rawCounts_DNA_DF_FIXED <- as.data.frame(colSums(rawCounts_DNA))
colnames(rawCounts_DNA_DF_FIXED)[1] <- "Reads"
rawCounts_DNA_DF_FIXED$Fomes <- colSums(rawCounts_DNA[GeneTaxonomy_DNA_Clean[GeneTaxonomy_DNA_Clean$Genus=="Fomes","Gene"],])
rawCounts_DNA_DF_FIXED$Fomes_TPM <- unlist(aggregate(Tidy_Taxonomy_DNA[Tidy_Taxonomy_DNA$TAX_Genus=="Fomes",c("Sample","TPM")], .~Sample,sum)["TPM"])
rawCounts_DNA_DF_FIXED$Fomes_Percent <- rawCounts_DNA_DF_FIXED$Fomes/rawCounts_DNA_DF_FIXED$Reads*100
rawCounts_DNA_DF_FIXED$Fomes_Percent_TPM <- rawCounts_DNA_DF_FIXED$Fomes_TPM/(1e+06)*100
rawCounts_DNA_DF_FIXED <- cbind(rawCounts_DNA_DF_FIXED, sampleData_DNA_FIXED)

#Bacteria
rawCounts_DNA_DF_FIXED$Bacteria <- colSums(rawCounts_DNA[GeneTaxonomy_DNA_Clean[GeneTaxonomy_DNA_Clean$Genus%in%unique(Tidy_Taxonomy_DNA[Tidy_Taxonomy_DNA$TAX_Domain=="Bacteria","TAX_Genus"]),"Gene"],])
rawCounts_DNA_DF_FIXED$Bacteria_TPM <- unlist(aggregate(Tidy_Taxonomy_DNA[Tidy_Taxonomy_DNA$TAX_Domain=="Bacteria",c("Sample","TPM")], .~Sample,sum)["TPM"])
rawCounts_DNA_DF_FIXED$Bacteria_Percent <- rawCounts_DNA_DF_FIXED$Bacteria/rawCounts_DNA_DF_FIXED$Reads*100
rawCounts_DNA_DF_FIXED$Bacteria_Percent_TPM <- rawCounts_DNA_DF_FIXED$Bacteria_TPM/(1e+06)*100

#Averages
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Fomes_Percent")],mean)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Fomes_Percent")],median)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Bacteria_Percent")],mean)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Bacteria_Percent")],median)

aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Fomes_Percent_TPM")],median)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Bacteria_Percent_TPM")],median)

aggregate(.~Condition, rawCounts_RNA_DF_FIXED[,c("Condition","Fomes_Percent")],median)
aggregate(.~Condition, rawCounts_RNA_DF_FIXED[,c("Condition","Fomes_Percent_TPM")],median)



RawCountPlot_DNA <- ggplot(rawCounts_DNA_DF_FIXED, aes(x = Sample, y = Reads, fill = Condition)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None") +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(n.breaks = 7) +
  labs(x=NULL, y = "Reads") + 
  scale_fill_manual(values = condition_colours)
RawCountPlot_DNA_Fomes <- ggplot(rawCounts_DNA_DF_FIXED, aes(x = Sample, y = Fomes/Reads*100, fill = Condition)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None") +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80)) +
  labs(x=NULL, y = "Percent Fomes Reads") + 
  scale_fill_manual(values = condition_colours)
RawCountPlot_DNA_Bacteria <- ggplot(rawCounts_DNA_DF_FIXED, aes(x = Sample, y = Bacteria/Reads*100, fill = Condition)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None") +
  geom_col() +
  facet_grid(. ~ Condition, scales = "free_x", space = "free") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80)) +
  labs(x=NULL, y = "Percent Bacterial Reads") + 
  scale_fill_manual(values = condition_colours)
RawCountPlot_DNA <- plot_grid(RawCountPlot_DNA + labs(tag = "A"), RawCountPlot_DNA_Fomes + labs(tag = "B"), ncol = 2)
ggsave("RawCountPlot_DNA.svg",RawCountPlot_DNA,width = 12, height = 4)

#Create a vector of Fomes Metagenomic TPM percentage for ordering
Sample_Order_DNA <- rownames(rawCounts_DNA_DF_FIXED[order(rawCounts_DNA_DF_FIXED$Fomes_Percent_TPM,decreasing = T),])
Sample_Order_RNA <- gsub("DNA","RNA",Sample_Order_DNA)

#Count family-level taxa in each sample
for (x in 1:nrow(rawCounts_RNA_DF_FIXED)) {
  rawCounts_RNA_DF_FIXED[x,"Family"] <- length(unique(Tidy_Taxonomy_RNA[Tidy_Taxonomy_RNA$Gene%in%row.names(rawCounts_RNA_FIXED[rawCounts_RNA_FIXED[,rawCounts_RNA_DF_FIXED[x,"Sample"]]>0,]),"TAX_Family"]))
}

aggregate(.~Condition, rawCounts_RNA_DF_FIXED[,c("Condition","Family")],min)
aggregate(.~Condition, rawCounts_RNA_DF_FIXED[,c("Condition","Family")],mean)
aggregate(.~Condition, rawCounts_RNA_DF_FIXED[,c("Condition","Family")],median)
aggregate(.~Condition, rawCounts_RNA_DF_FIXED[,c("Condition","Family")],max)

for (x in 1:nrow(rawCounts_DNA_DF_FIXED)) {
  rawCounts_DNA_DF_FIXED[x,"Family"] <- length(unique(Tidy_Taxonomy_DNA[Tidy_Taxonomy_DNA$Gene%in%row.names(rawCounts_DNA[rawCounts_DNA[,rawCounts_DNA_DF_FIXED[x,"Sample"]]>0,]),"TAX_Family"]))
}

aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Family")],min)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Family")],mean)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Family")],median)
aggregate(.~Condition, rawCounts_DNA_DF_FIXED[,c("Condition","Family")],max)

###What is the overall taxonomy?

#Plot the data
Taxonomy_RNA_TPM_Barplots <- list()
Taxonomy_DNA_TPM_Barplots <- list()

for (TERM in c("TAX_Phylum", "TAX_Order", "TAX_Family","TAX_Genus")) {
  #Set the threshold
  THRESH <- 0.01
  
  #Combine low abundance for better visibility
  Fomes_tidy_shrunk_RNA <- Combine_Low_Abundance(Tidy_Taxonomy_RNA,TERM,"TPM",THRESH)
  Fomes_tidy_shrunk_DNA <- Combine_Low_Abundance(Tidy_Taxonomy_DNA,TERM,"TPM",THRESH)
  
  #Set the order according to Fomes percentage
  Fomes_tidy_shrunk_RNA$Sample <- factor(Fomes_tidy_shrunk_RNA$Sample,levels = Sample_Order_RNA)
  Fomes_tidy_shrunk_DNA$Sample <- factor(Fomes_tidy_shrunk_DNA$Sample,levels = Sample_Order_DNA)
  
  #Set the taxa name order
  Taxa_names <- sort(unique(c(unique(Fomes_tidy_shrunk_RNA[[TERM]]),unique(Fomes_tidy_shrunk_DNA[[TERM]]))))
  TERM_Fomes <- as.character(Fomes_tidy_shrunk_RNA[Fomes_tidy_shrunk_RNA$TAX_Genus=="Fomes",][1,TERM])
  Taxa_names<- Taxa_names[c(which(Taxa_names==TERM_Fomes), which(!Taxa_names%in%c(TERM_Fomes,paste("Below ",THRESH*100,"%",sep = ""),"Unknown")),
               which(Taxa_names==paste("Below ",THRESH*100,"%",sep = "")), which(Taxa_names=="Unknown"))]
  Fomes_tidy_shrunk_DNA[,TERM] <- factor(Fomes_tidy_shrunk_DNA[,TERM], levels = Taxa_names)
  Fomes_tidy_shrunk_RNA[,TERM] <- factor(Fomes_tidy_shrunk_RNA[,TERM], levels = Taxa_names)

  #All same colours in all figures
  Taxa_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(Taxa_names))
  names(Taxa_colours) <- Taxa_names
  
  #Figures
  Fomes_tidy_shrunk_RNA_fig <- 
    ggplot(Fomes_tidy_shrunk_RNA, aes(x = Sample)) +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM), colour = "black") +
    #facet_wrap(~Condition, scales = "free_x") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "TPM", title = "Metatranscriptomic taxonomy", fill = gsub("TAX_","",TERM)) + 
    scale_fill_manual(values = Taxa_colours, limits = names(Taxa_colours), guide = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  Fomes_tidy_shrunk_DNA_fig <- 
    ggplot(Fomes_tidy_shrunk_DNA, aes(x = Sample)) +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM), colour = "black") +
    #facet_wrap(~Condition, scales = "free_x") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "TPM", title = "Metagenomic taxonomy", fill = gsub("TAX_","",TERM)) + 
    scale_fill_manual(values = Taxa_colours, limits = names(Taxa_colours), guide = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  #Add the figures to the lists
  Taxonomy_RNA_TPM_Barplots[TERM] <- list(Fomes_tidy_shrunk_RNA_fig)
  Taxonomy_DNA_TPM_Barplots[TERM] <- list(Fomes_tidy_shrunk_DNA_fig)
}

#Save the figures
for (TERM in c("TAX_Phylum", "TAX_Order", "TAX_Family")) {
  ggsave(paste("Taxonomy_RNA_TPM_Barplots_",TERM,".svg",sep = ""),Taxonomy_RNA_TPM_Barplots[[TERM]], width = 7, height = 7)
  ggsave(paste("Taxonomy_DNA_TPM_Barplots_",TERM,".svg",sep = ""),Taxonomy_DNA_TPM_Barplots[[TERM]], width = 7, height = 7)
}

#Figure 1
Fig_1_legend <- get_legend(Taxonomy_RNA_TPM_Barplots$TAX_Family)
Figure_1 <-
  plot_grid(
    plot_grid(Taxonomy_DNA_TPM_Barplots$TAX_Family + labs(tag = "A") + theme(legend.position = "none"),
              Taxonomy_RNA_TPM_Barplots$TAX_Family + labs(tag = "B") + theme(legend.position = "none"),
              ncol = 1),
    Fig_1_legend,
    ncol = 2, rel_widths = c(1,0.4))
ggsave("Figure_1.eps",Figure_1,width = 2*114, height = 2*100, units = "mm", )

#Figure S2
Fig_S2_legend <- get_legend(Taxonomy_RNA_TPM_Barplots$TAX_Phylum)
Figure_S2 <-
  plot_grid(
    plot_grid(Taxonomy_DNA_TPM_Barplots$TAX_Phylum + labs(tag = "A") + theme(legend.position = "none"),
              Taxonomy_RNA_TPM_Barplots$TAX_Phylum + labs(tag = "B") + theme(legend.position = "none"),
              ncol = 1),
    Fig_S2_legend,
    ncol = 2, rel_widths = c(1,0.3))
ggsave("Figure_S2.eps",Figure_S2,width = 2*114, height = 2*100, units = "mm")

#Clean up figures
#The figures are easy/quick to generate if needed but are exceptional large which puts pressure on memory resources 
#and greatly slows down saving the session.
rm(Taxonomy_DNA_TPM_Barplots,Fomes_tidy_shrunk_DNA,Fomes_tidy_shrunk_DNA_fig)

#Most abundant families
Fomes_tidy_shrunk_DNA <- Combine_Low_Abundance(Tidy_Taxonomy_DNA,"TAX_Family","TPM",0.01)
Rotten_Families <- as.data.frame(aggregate(.~TAX_Family, Fomes_tidy_shrunk_DNA[Fomes_tidy_shrunk_DNA$Condition=="Rotten",c("TAX_Family","TPM")],sum))
Rotten_Families[,"Percent"] <- Rotten_Families$TPM/sum(Rotten_Families$TPM)*100
Fomes_tidy_shrunk_DNA_NoH19 <- Fomes_tidy_shrunk_DNA[!Fomes_tidy_shrunk_DNA$Sample=="DNA_H19",]
Fresh_Families <- as.data.frame(aggregate(.~TAX_Family, Fomes_tidy_shrunk_DNA_NoH19[Fomes_tidy_shrunk_DNA_NoH19$Condition=="Fresh",c("TAX_Family","TPM")],sum))
Fresh_Families[,"Percent"] <- Fresh_Families$TPM/sum(Fresh_Families$TPM)*100

##CAZyme Expression##

#Limit the tidy tables to only CAZymes
Tidy_CAZymes_RNA <- subset(Tidy_Taxonomy_RNA, !Tidy_Taxonomy_RNA$CAZ_model=="-")
Tidy_CAZymes_DNA <- subset(Tidy_Taxonomy_DNA, !Tidy_Taxonomy_DNA$CAZ_model=="-")

#Who is expressing or what is being expressed?

Cazymes_TPM_RNA_Barplots_Counts <- list()
Cazymes_TPM_RNA_Barplots_Props <- list()
Cazymes_TPM_DNA_Barplots_Counts <- list()
Cazymes_TPM_DNA_Barplots_Props <- list()

for (TERM in c("TAX_Phylum", "TAX_Order", "TAX_Genus", "CAZ_Family", "CAZ_Subfamily", "CAZ_target")) {
  #Set threshhold
  THRESH <- 0.01
  
  #Combine low abundance for better visibility
  Tidy_CAZymes_RNA_shrunk <- Combine_Low_Abundance(Tidy_CAZymes_RNA,TERM,"TPM",THRESH)
  Tidy_CAZymes_DNA_shrunk <- Combine_Low_Abundance(Tidy_CAZymes_DNA,TERM,"TPM",THRESH)
  
  #Set the order according to Fomes percentage
  Tidy_CAZymes_DNA_shrunk$Sample <- factor(Tidy_CAZymes_DNA_shrunk$Sample,levels = Sample_Order_DNA)
  Tidy_CAZymes_RNA_shrunk$Sample <- factor(Tidy_CAZymes_RNA_shrunk$Sample,levels = Sample_Order_RNA)
  
  #Set the term order
  Term_names <- sort(unique(c(unique(Tidy_CAZymes_RNA_shrunk[[TERM]]),unique(Tidy_CAZymes_DNA_shrunk[[TERM]]))))
  Term_names<- Term_names[c(which(!Term_names%in%c(paste("Below ",THRESH*100,"%",sep = ""),"Unknown")),
                            which(Term_names==paste("Below ",THRESH*100,"%",sep = "")), which(Term_names=="Unknown"))]
  Tidy_CAZymes_DNA_shrunk[,TERM] <- factor(Tidy_CAZymes_DNA_shrunk[,TERM], levels = Term_names, ordered = T)
  Tidy_CAZymes_RNA_shrunk[,TERM] <- factor(Tidy_CAZymes_RNA_shrunk[,TERM], levels = Term_names, ordered = T)
  
  #All same colours in all figures
  Term_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(Term_names))
  if (length(Term_names) < 9) {
    Term_colours <- colorRampPalette(brewer.pal(8, "Set2"))(length(Term_names))
  }
  names(Term_colours) <- Term_names
  
  #Identify which CAZymes come from Fomes
  Tidy_CAZymes_RNA_shrunk[,"Source"] <- Tidy_CAZymes_RNA_shrunk$TAX_Genus=="Fomes"
  Tidy_CAZymes_RNA_shrunk[,"Source"] <- gsub("TRUE","Fomes",Tidy_CAZymes_RNA_shrunk[,"Source"])
  Tidy_CAZymes_RNA_shrunk[,"Source"] <- gsub("FALSE","Non-Fomes",Tidy_CAZymes_RNA_shrunk[,"Source"])
  Tidy_CAZymes_DNA_shrunk[,"Source"] <- Tidy_CAZymes_DNA_shrunk$TAX_Genus=="Fomes"
  Tidy_CAZymes_DNA_shrunk[,"Source"] <- gsub("TRUE","Fomes",Tidy_CAZymes_DNA_shrunk[,"Source"])
  Tidy_CAZymes_DNA_shrunk[,"Source"] <- gsub("FALSE","Non-Fomes",Tidy_CAZymes_DNA_shrunk[,"Source"])
  Source_colours <- c("#d8b365","#5ab4ac")
  names(Source_colours) <- c("Non-Fomes","Fomes")
  
  #Counts (TPM)
  Tidy_CAZymes_RNA_shrunk_count_fig <- 
    ggplot(Tidy_CAZymes_RNA_shrunk, aes(x = Sample)) +
    geom_bar(aes(fill = Source, weight = TPM), just = 1.55,width = 0.3, colour = "black") +
    scale_fill_manual(values = Source_colours) +
    new_scale_fill() +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM, group =factor(paste(Tidy_CAZymes_RNA_shrunk[,"Source"],Tidy_CAZymes_RNA_shrunk[,TERM]),levels = c(paste("Fomes ",levels(Tidy_CAZymes_RNA_shrunk[,TERM]),sep = ""),paste("Non-Fomes ",levels(Tidy_CAZymes_RNA_shrunk[,TERM]),sep = "")))), just = 0.25,width = 0.6, colour = "black") +
    scale_fill_manual(values = Term_colours, limits = names(Term_colours)) +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "Count (TPM)", title = "Metatranscriptome CAZymes", fill = str_to_sentence(gsub("[A-Z]{3}_","",TERM))) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16),
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  Tidy_CAZymes_DNA_shrunk_count_fig <- 
    ggplot(Tidy_CAZymes_DNA_shrunk, aes(x = Sample)) +
    geom_bar(aes(fill = Source, weight = TPM), just = 1.55,width = 0.3, colour = "black") +
    scale_fill_manual(values = Source_colours) +
    new_scale_fill() +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM, group =factor(paste(Tidy_CAZymes_DNA_shrunk[,"Source"],Tidy_CAZymes_DNA_shrunk[,TERM]),levels = c(paste("Fomes ",levels(Tidy_CAZymes_DNA_shrunk[,TERM]),sep = ""),paste("Non-Fomes ",levels(Tidy_CAZymes_DNA_shrunk[,TERM]),sep = "")))), just = 0.25,width = 0.6, colour = "black") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "Count (TPM)", title = "Metagenome CAZymes", fill = str_to_sentence(gsub("[A-Z]{3}_","",TERM))) + 
    scale_fill_manual(values = Term_colours, limits = names(Term_colours)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16),
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  #Proportion
  Tidy_CAZymes_DNA_shrunk_prop_fig <- 
    ggplot(Tidy_CAZymes_DNA_shrunk, aes(x = Sample)) +
    geom_bar(aes(fill = Source, weight = TPM), position = "fill", just = 1.55,width = 0.3, colour = "black") +
    scale_fill_manual(values = Source_colours) +
    new_scale_fill() +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM, group =factor(paste(Tidy_CAZymes_DNA_shrunk[,"Source"],Tidy_CAZymes_DNA_shrunk[,TERM]),levels = c(paste("Fomes ",levels(Tidy_CAZymes_DNA_shrunk[,TERM]),sep = ""),paste("Non-Fomes ",levels(Tidy_CAZymes_DNA_shrunk[,TERM]),sep = "")))),position = "fill", just = 0.25,width = 0.6, colour = "black") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "Proportion of TPM counts", title = "Metagenome CAZymes", fill = str_to_sentence(gsub("[A-Z]{3}_","",TERM))) + 
    scale_fill_manual(values = Term_colours, limits = names(Term_colours)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16),
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  Tidy_CAZymes_RNA_shrunk_prop_fig <- 
    ggplot(Tidy_CAZymes_RNA_shrunk, aes(x = Sample)) +
    geom_bar(aes(fill = Source, weight = TPM), position = "fill", just = 1.55,width = 0.3, colour = "black") +
    scale_fill_manual(values = Source_colours) +
    new_scale_fill() +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM, group =factor(paste(Tidy_CAZymes_RNA_shrunk[,"Source"],Tidy_CAZymes_RNA_shrunk[,TERM]),levels = c(paste("Fomes ",levels(Tidy_CAZymes_RNA_shrunk[,TERM]),sep = ""),paste("Non-Fomes ",levels(Tidy_CAZymes_RNA_shrunk[,TERM]),sep = "")))),position = "fill", just = 0.25,width = 0.6, colour = "black") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "Proportion of TPM counts", title = "Metatranscriptome CAZymes", fill = str_to_sentence(gsub("[A-Z]{3}_","",TERM))) + 
    scale_fill_manual(values = Term_colours, limits = names(Term_colours)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16),
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  #Add the figures to the lists
  Cazymes_TPM_RNA_Barplots_Counts[TERM] <- list(Tidy_CAZymes_RNA_shrunk_count_fig)
  Cazymes_TPM_RNA_Barplots_Props[TERM] <- list(Tidy_CAZymes_RNA_shrunk_prop_fig)
  Cazymes_TPM_DNA_Barplots_Counts[TERM] <- list(Tidy_CAZymes_DNA_shrunk_count_fig)
  Cazymes_TPM_DNA_Barplots_Props[TERM] <- list(Tidy_CAZymes_DNA_shrunk_prop_fig)
  
  #Save the figures
  #If not done during the loop, the colours will not save correctly.
  ggsave(paste("Cazymes_TPM_RNA_Barplots_Counts_",TERM,".svg",sep = ""),Cazymes_TPM_RNA_Barplots_Counts[[TERM]], width = 10, height = 7)
  ggsave(paste("Cazymes_TPM_RNA_Barplots_Props_",TERM,".svg",sep = ""),Cazymes_TPM_RNA_Barplots_Props[[TERM]], width = 10, height = 7)
  ggsave(paste("Cazymes_TPM_DNA_Barplots_Counts_",TERM,".svg",sep = ""),Cazymes_TPM_DNA_Barplots_Counts[[TERM]], width = 10, height = 7)
  ggsave(paste("Cazymes_TPM_DNA_Barplots_Props_",TERM,".svg",sep = ""),Cazymes_TPM_DNA_Barplots_Props[[TERM]], width = 10, height = 7)
}

#For graphical abstract
ggsave("Graphical_Abstract_panel.eps",Cazymes_TPM_DNA_Barplots_Props$CAZ_target + theme(axis.text.x=element_blank(), text=element_text(size=14)),width = 168.000, height = 2*93.8, units = "mm")

#Construct figure 2
Fig_2_legend <- get_legend(Cazymes_TPM_RNA_Barplots_Counts$CAZ_target)
Figure_2 <- 
  plot_grid(
    plot_grid(Cazymes_TPM_DNA_Barplots_Counts$CAZ_target + labs(tag = "A") + theme(legend.position = "none"),
              Cazymes_TPM_RNA_Barplots_Counts$CAZ_target + labs(tag = "B") + theme(legend.position = "none"),
              Cazymes_TPM_DNA_Barplots_Props$CAZ_target + labs(tag = "C") + theme(legend.position = "none"),
              Cazymes_TPM_RNA_Barplots_Props$CAZ_target + labs(tag = "D") + theme(legend.position = "none"),
              ncol = 2),
    Fig_2_legend,
    ncol = 2, rel_widths = c(1,0.35))
ggsave("Figure_2.eps",Figure_2,width = 2*168, height = 2*168, units = "mm")

##Who is secreting what?

#Set threshold
THRESH <- 0.01

#Combine low abundance for better visibility
Tidy_CAZymes_RNA_shrunk_grouped <- Tidy_CAZymes_RNA[Tidy_CAZymes_RNA$Condition=="Rotten",]
Tidy_CAZymes_RNA_shrunk_grouped <- Combine_Low_Abundance(Tidy_CAZymes_RNA_shrunk_grouped,"TAX_Phylum","TPM",THRESH)
Tidy_CAZymes_RNA_shrunk_grouped <- Combine_Low_Abundance(Tidy_CAZymes_RNA_shrunk_grouped,"CAZ_target","TPM",THRESH)

#Sort the targets
Target_levels <- sort(unique(Tidy_CAZymes_RNA_shrunk_grouped$CAZ_target),decreasing = TRUE)
Target_levels <- Target_levels[c(which(Target_levels=="Unknown"), which(Target_levels==paste("Below ",THRESH*100,"%",sep = "")),
                               which(!Target_levels%in%c(paste("Below ",THRESH*100,"%",sep = ""),"Unknown")))]
Tidy_CAZymes_RNA_shrunk_grouped$CAZ_target <- factor(Tidy_CAZymes_RNA_shrunk_grouped$CAZ_target, levels = Target_levels)

#Sort the phyla
Phylum_levels <- sort(unique(Tidy_CAZymes_RNA_shrunk_grouped$TAX_Phylum))
Phylum_levels <- Phylum_levels[c(which(!Phylum_levels%in%paste("Below ",THRESH*100,"%",sep = "")),
                                 which(Phylum_levels==paste("Below ",THRESH*100,"%",sep = "")))]
Tidy_CAZymes_RNA_shrunk_grouped$TAX_Phylum <- factor(Tidy_CAZymes_RNA_shrunk_grouped$TAX_Phylum, levels = Phylum_levels)

#Safe colours
F3_taxa <- sort(unique(Tidy_CAZymes_RNA_shrunk_grouped$TAX_Phylum))
F3_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(F3_taxa))
names(F3_colours) <- F3_taxa

#Counts (TPM)
Figure_3 <- 
  ggplot(Tidy_CAZymes_RNA_shrunk_grouped, aes(y = CAZ_target)) +
  geom_bar(aes(fill = TAX_Phylum, weight = TPM), colour = "black") +
  facet_wrap(c("Sample"), ncol = 1) +
  labs(y=NULL, x = "Count (TPM)", title = "CAZyme expressing taxa", fill = "Phylum") + 
  scale_fill_manual(values = F3_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1), text=element_text(size=16),
        legend.position = "bottom", plot.title = element_text(hjust = 0.5))
ggsave(paste("Figure_3.eps",sep = ""),Figure_3,width = 2*168, height = 2*115, units = "mm")

#Generate tables of what is expressed by whom

Cazyme_gene_expression_tables <- list()

for (TERM_2 in c("CAZ_Family", "CAZ_Subfamily", "CAZ_target")) {
  
  #Individual tables by condition
  fresh_table <- as.data.frame(table(Tidy_CAZymes_RNA[Tidy_CAZymes_RNA$Condition=="Fresh"&Tidy_CAZymes_RNA$TPM>0,TERM_2],Tidy_CAZymes_RNA[Tidy_CAZymes_RNA$Condition=="Fresh"&Tidy_CAZymes_RNA$TPM>0,"TAX_Family"]))
  colnames(fresh_table) <- c(TERM_2,"TAX_Family","Fresh_Genes")
  rotten_table <- as.data.frame(table(Tidy_CAZymes_RNA[Tidy_CAZymes_RNA$Condition=="Rotten"&Tidy_CAZymes_RNA$TPM>0,TERM_2],Tidy_CAZymes_RNA[Tidy_CAZymes_RNA$Condition=="Rotten"&Tidy_CAZymes_RNA$TPM>0,"TAX_Family"]))
  colnames(rotten_table) <- c(TERM_2,"TAX_Family","Rotten_Genes")
  
  #Merge tables
  merged_table <- merge.data.frame(fresh_table,rotten_table, all = TRUE)
  
  #Convert NAs to 0 for merging
  merged_table[is.na(merged_table$Fresh_Genes),"Fresh_Genes"] <- 0
  merged_table[is.na(merged_table$Rotten_Genes),"Rotten_Genes"] <- 0
  
  #Remove unneeded lines
  merged_table <- merged_table[!(merged_table$Fresh_Genes+merged_table$Rotten_Genes)==0,]
  
  #Remove factors
  merged_table[,"TAX_Family"] <- as.character(merged_table[,"TAX_Family"])
  merged_table[,TERM_2] <- as.character(merged_table[,TERM_2])
  
  #Add in the other taxonomic levels
  merged_table <- merge.data.frame(merged_table,unique(Tidy_CAZymes_RNA[,c("TAX_Phylum","TAX_Class","TAX_Order","TAX_Family")]))
  
  #Sort
  merged_table <- merged_table[,c(TERM_2,"TAX_Phylum","TAX_Class","TAX_Order","TAX_Family","Fresh_Genes","Rotten_Genes")]
  merged_table <- merged_table[order(merged_table[,TERM_2],merged_table[,"TAX_Phylum"],merged_table[,"TAX_Class"],merged_table[,"TAX_Order"],merged_table[,"TAX_Family"]),]
  
  #Add in TPM
  for (ROW in 1:nrow(merged_table)) {
    merged_table[ROW,"Fresh_TPM"] <- round(sum(Tidy_CAZymes_RNA[Tidy_CAZymes_RNA[,TERM_2]==merged_table[ROW,TERM_2] & 
                                                        Tidy_CAZymes_RNA$TAX_Family==merged_table[ROW,"TAX_Family"] &
                                                        Tidy_CAZymes_RNA$Condition=="Fresh","TPM"]),digits = 3)
    merged_table[ROW,"Rotten_TPM"] <- round(sum(Tidy_CAZymes_RNA[Tidy_CAZymes_RNA[,TERM_2]==merged_table[ROW,TERM_2] & 
                                                            Tidy_CAZymes_RNA$TAX_Family==merged_table[ROW,"TAX_Family"] &
                                                            Tidy_CAZymes_RNA$Condition=="Rotten","TPM"]),digits = 3)
  }
  
  #Save table
  Cazyme_gene_expression_tables[paste(TERM_2,sep = "_")] <- list(merged_table)
}

#Save tables
write.table(Cazyme_gene_expression_tables$CAZ_target,"Table_S2.csv",row.names = FALSE,sep = ",")


#Bacterial functional potential

#Focus on bacteria due to difficulty of dealing with introns
#Unknown, archaea and viruses are also excluded. Archaea and viruses are a tiny number of sequences.
#Uknown are sizeable but few KOs and no complete pathways according to EvaluatePathways(). They can all be safely excluded without affecting functional analysis.
Tidy_Taxonomy_DNA_Bacteria <- subset(Tidy_Taxonomy_DNA,Tidy_Taxonomy_DNA$TAX_Domain=="Bacteria")

#Overview of energy metabolism in different conditions
for (COND in unique(Tidy_Taxonomy_DNA_Bacteria$Condition)) {
  #Identify all the families and count their TPM
  MG_Bact_Families <- aggregate(.~TAX_Family, Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$Condition == COND,c("TPM","TAX_Family")],sum)
  MG_Bact_Families <- MG_Bact_Families[order(MG_Bact_Families[,"TPM"],decreasing = T),]
  MG_Bact_Families$Percent <- MG_Bact_Families$TPM/sum(MG_Bact_Families$TPM)*100
  assign(paste("MG_Bact_Families",COND,sep = "_"),MG_Bact_Families)
  
  #Identify how many pathways are present overall and in each family
  MG_Bact_Pathways <- EvaluatePathways(Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$Condition == COND,"KOFAM"])
  colnames(MG_Bact_Pathways) <- c("Energy Metabolism","Module","Name","Bacteria")
  for (TAX in MG_Bact_Families[MG_Bact_Families$Percent>1,"TAX_Family"]) {
    TAX_path <- EvaluatePathways(Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$Condition == COND & Tidy_Taxonomy_DNA_Bacteria$TAX_Family==TAX,"KOFAM"])
    colnames(TAX_path) <- c("Energy Metabolism","Module","Name",TAX)
    MG_Bact_Pathways <- merge.data.frame(MG_Bact_Pathways,TAX_path)
  }
  
  #Count how many pathways are present
  for (ROW in 1:nrow(MG_Bact_Pathways)) {
    MG_Bact_Pathways[ROW,"Number of Families"] <- sum(as.logical(MG_Bact_Pathways[ROW,colnames(MG_Bact_Pathways)%in%MG_Bact_Families$TAX_Family]))
  }
  
  #Save a table
  assign(paste("MG_Bact_Pathways",COND,sep = "_"),MG_Bact_Pathways)
}

#Combine Family Abundances
MG_Bact_Families_Combined <- merge.data.frame(MG_Bact_Families_Fresh,MG_Bact_Families_Rotten,by="TAX_Family",all = TRUE)
colnames(MG_Bact_Families_Combined) <- c("TAX_Family","TPM_Fresh","Percentage_Fresh","TPM_Rotten","Percentage_Rotten")
MG_Bact_Families_Combined <- merge.data.frame(MG_Bact_Families_Combined, unique(Taxonomy_Organisms[,c("TAX_Phylum","TAX_Class","TAX_Order","TAX_Family")]), by = "TAX_Family")
MG_Bact_Families_Combined <- MG_Bact_Families_Combined[,c("TAX_Phylum","TAX_Class","TAX_Order","TAX_Family","TPM_Fresh","Percentage_Fresh","TPM_Rotten","Percentage_Rotten")]
MG_Bact_Families_Combined_Sub <- MG_Bact_Families_Combined[MG_Bact_Families_Combined$Percentage_Fresh>1 | MG_Bact_Families_Combined$Percentage_Rotten>1,]
MG_Bact_Families_Combined_Sub <- MG_Bact_Families_Combined_Sub[!is.na(MG_Bact_Families_Combined_Sub$TAX_Family),]

#Broad abundance of bacteria phyla
MG_Bact_Phylum_Abundance <- merge.data.frame(aggregate(x = MG_Bact_Families_Combined[,c("TAX_Phylum","Percentage_Fresh")],by = .~TAX_Phylum, sum),
                 aggregate(x = MG_Bact_Families_Combined[,c("TAX_Phylum","Percentage_Rotten")],by = .~TAX_Phylum, sum))

#Export table
#Combine manually for display
write.table(MG_Bact_Pathways_Fresh,file = "MG_Bact_Pathways_Fresh.csv",quote = T,row.names = FALSE,sep = ",")
write.table(MG_Bact_Pathways_Rotten,file = "MG_Bact_Pathways_Rotten.csv",quote = T,row.names = FALSE,sep = ",")

#Nitrogen cycling

#NifH missing
"[KO2588]"%in%Tidy_Taxonomy_DNA_Bacteria$KOFAM
"[KO2588]"%in%Tidy_Taxonomy_RNA$KOFAM
#nifD present
"[K02586]"%in%Tidy_Taxonomy_DNA_Bacteria$KOFAM
"[K02586]"%in%Tidy_Taxonomy_RNA$KOFAM
#nifK  missing
"[K02591]"%in%Tidy_Taxonomy_DNA_Bacteria$KOFAM
"[K02591]"%in%Tidy_Taxonomy_RNA$KOFAM
#None in metatranscriptome

#nifD only in Proteobacteria (Family: Xanthobacteraceae) and then only in 4 samples!
Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$KOFAM=="[K02586]",]

#Check for dasABC chitin transporter proteins
#dasA present in metagenome only
"[K17329]"%in%Tidy_Taxonomy_DNA_Bacteria$KOFAM
"[K17329]"%in%Tidy_Taxonomy_RNA$KOFAM
#dasB present in metagenome only
"[K17330]"%in%Tidy_Taxonomy_DNA_Bacteria$KOFAM
"[K17330]"%in%Tidy_Taxonomy_RNA$KOFAM
#dasC present in metagenome only
"[K17331]"%in%Tidy_Taxonomy_DNA_Bacteria$KOFAM
"[K17331]"%in%Tidy_Taxonomy_RNA$KOFAM

#dasA only found in Microbacteriaceae and Micromonosporaceae
Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$KOFAM=="[K17329]",]
#dasB only found in Microbacteriaceae and Micromonosporaceae
Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$KOFAM=="[K17330]",]
#dasC only found in Microbacteriaceae
Tidy_Taxonomy_DNA_Bacteria[Tidy_Taxonomy_DNA_Bacteria$KOFAM=="[K17331]",]


#Get a breakdown of Ascomycota

Ascos <- Tidy_Taxonomy_DNA[Tidy_Taxonomy_DNA$TPM>0&Tidy_Taxonomy_DNA$TAX_Phylum=="Ascomycota",]

sort(table(Ascos$TAX_Class),decreasing = T)
sort(table(Ascos$TAX_Order),decreasing = T)
sort(table(Ascos$TAX_Family),decreasing = T)


#Plot the data
Asco_DNA_TPM_Barplots_Count <- list()
Asco_DNA_TPM_Barplots_Prop <- list()

for (TERM in c("TAX_Class", "TAX_Order", "TAX_Family")) {
  #Set the threshold
  THRESH <- 0.01
  
  #Combine low abundance for better visibility
  Asco_tidy_shrunk_DNA <- Combine_Low_Abundance(Ascos,TERM,"TPM",THRESH)
  
  #Set the order according to Fomes percentage
  Asco_tidy_shrunk_DNA$Sample <- factor(Asco_tidy_shrunk_DNA$Sample,levels = Sample_Order_DNA)
  
  # #Set the taxa name order
  Asco_names <- sort(unique(Asco_tidy_shrunk_DNA[[TERM]]))
  Asco_names<- Asco_names[c(which(!Asco_names%in%paste("Below ",THRESH*100,"%",sep = "")),
                            which(Asco_names%in%paste("Below ",THRESH*100,"%",sep = "")))]
  Asco_tidy_shrunk_DNA[,TERM] <- factor(Asco_tidy_shrunk_DNA[,TERM], levels = Asco_names)
  
  #All same colours in all figures
  Asco_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(Asco_names))
  names(Asco_colours) <- Asco_names
  
  Asco_tidy_shrunk_DNA_Count <- 
    ggplot(Asco_tidy_shrunk_DNA, aes(x = Sample)) +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM), colour = "black") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "Count (TPM)", title = "MG Ascomycota", fill = gsub("TAX_","",TERM)) + 
    scale_fill_manual(values = Asco_colours, limits = names(Asco_colours), guide = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  Asco_tidy_shrunk_DNA_Prop <- 
    ggplot(Asco_tidy_shrunk_DNA, aes(x = Sample)) +
    geom_bar(aes(fill = .data[[TERM]], weight = TPM), position = "fill", colour = "black") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free") +
    labs(x=NULL, y = "Proportion of TPM Counts", title = "MG Ascomycota", fill = gsub("TAX_","",TERM)) + 
    scale_fill_manual(values = Asco_colours, limits = names(Asco_colours), guide = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text=element_text(size=16), 
          legend.position = "right", plot.title = element_text(hjust = 0.5))
  
  #Add the figures to the lists
  Asco_DNA_TPM_Barplots_Count[TERM] <- list(Asco_tidy_shrunk_DNA_Count)
  Asco_DNA_TPM_Barplots_Prop[TERM] <- list(Asco_tidy_shrunk_DNA_Prop)
}

Ascos_graph <- plot_grid(Asco_DNA_TPM_Barplots_Count$TAX_Class + labs(tag = "A"),
                         Asco_DNA_TPM_Barplots_Prop$TAX_Class + labs(tag = "B"),
                         Asco_DNA_TPM_Barplots_Count$TAX_Order + labs(tag = "C"),
                         Asco_DNA_TPM_Barplots_Prop$TAX_Order + labs(tag = "D"),
                         Asco_DNA_TPM_Barplots_Count$TAX_Family + labs(tag = "E"),
                         Asco_DNA_TPM_Barplots_Prop$TAX_Family + labs(tag = "F"),
                         ncol = 2)
ggsave("Figure_S3.eps",Ascos_graph,width = 2*168, height = 2*237, units = "mm")



#####

#Session info

# sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Arch Linux
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblas.so.0.3;  LAPACK version 3.11.0
# 
# locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Prague
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] Rsubread_2.14.2    ggnewscale_0.4.9   ggplotify_0.1.2    stringr_1.5.0      ggdendro_0.1.23    shipunov_1.17.1    NbClust_3.0.1     
# [8] cluster_2.1.4      cowplot_1.1.1      ggrepel_0.9.3      pheatmap_1.0.12    RColorBrewer_1.1-3 ggplot2_3.4.3      groundhog_3.1.1   
# 
# loaded via a namespace (and not attached):
# [1] yulab.utils_0.0.8  utf8_1.2.3         generics_0.1.3     stringi_1.7.12     lattice_0.21-8     magrittr_2.0.3     grid_4.3.2        
# [8] fastmap_1.1.1      Matrix_1.6-1       fansi_1.0.4        scales_1.2.1       textshaping_0.3.7  cli_3.6.1          rlang_1.1.1       
# [15] munsell_0.5.0      withr_2.5.0        cachem_1.0.8       tools_4.3.2        parallel_4.3.2     memoise_2.0.1      dplyr_1.1.3       
# [22] colorspace_2.1-0   vctrs_0.6.3        R6_2.5.1           gridGraphics_0.5-1 lifecycle_1.0.3    MASS_7.3-60        ragg_1.2.6        
# [29] pkgconfig_2.0.3    pillar_1.9.0       gtable_0.3.4       glue_1.6.2         Rcpp_1.0.11        systemfonts_1.0.5  tibble_3.2.1      
# [36] tidyselect_1.2.0   rstudioapi_0.15.0  farver_2.1.1       svglite_2.1.2      labeling_0.4.3     compiler_4.3.2   