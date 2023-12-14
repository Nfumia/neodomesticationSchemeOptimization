library(dplyr)
# POP1 REP1
### List files
setwd("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP1/")
data_files1 <- list.files("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP1/")  # Identify file names
data_files1
### Read files loop
foo1 <- list()
for(i in 1:length(data_files1)) {                              
  assign(paste0("output_", i),                                   
         read.csv(paste0("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP1/",
                          data_files1[i])))
}
foo1 <- lapply(data_files1,read.csv)
all1 <- do.call(rbind,foo1)
all1$rep <- "Rep1"
# POP1 REP2
### List files
setwd("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP2/")
data_files2 <- list.files("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP2/")  # Identify file names
data_files2
### Read files loop
foo2 <- list()
for(i in 1:length(data_files2)) {                              
  assign(paste0("output_", i),                                   
         read.csv(paste0("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP2/",
                         data_files2[i])))
}
foo2 <- lapply(data_files2,read.csv)
all2 <- do.call(rbind,foo2)
all2$rep <- "Rep2"
# POP1 REP2
### List files
setwd("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP3/")
data_files3 <- list.files("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP3/")  # Identify file names
data_files3
### Read files loop
foo3 <- list()
for(i in 1:length(data_files3)) {                              
  assign(paste0("output_", i),                                   
         read.csv(paste0("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/1POP/REP3/",
                         data_files3[i])))
}
foo3 <- lapply(data_files3,read.csv)
all3 <- do.call(rbind,foo3)
all3$rep <- "Rep3"

all <- rbind(all1,all2,all3)

rm(all1,all2,all3)
rm(list = ls(pattern = 'output_'))

# What Phase is each cycle?
all <- all %>% mutate(type=case_when(
                  Cycle >= 0 & Cycle <= 40 ~ "burnin",
                  Cycle >= 41 & Cycle <= 42 ~ "germplasmPRS",
                  Cycle >= 43 & Cycle <= 80 ~ "popdev"
                  ))

# Summarize each scheme across replications 
allSummary <- all %>% group_by(NeStart,Nparents,Nprogeny,si,h2,sel,pop,Cycle,type) %>%
  summarize(meanP1=mean(meanPz1),
            meanP2=mean(meanPz2),
            sdP1=sd(meanPz1),
            sdP2=sd(meanPz2),
            var1=mean(varPz1),
            var2=mean(varPz2),
            sdV1=sd(varPz1),
            sdV2=sd(varPz2))
write.csv(allSummary,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/allSummary.csv")
# Apply costs to each type of crop and the scheme used. Phenotyping/plot ranges from $10-20 in field crops, 
# $30-50 in horticultural crops, $40-60 in forestry crops. Genomic selection ranges from $5-15 per sample with
# phenotyping occurring every other cycle meaning half of the applied phenotyping cost is applied per cycle. 
# Maximum avoidance matches PRS variable costs. 
allCost <- all %>% 
  mutate(pLfield=case_when(type=="burnin" & sel=="PRS" ~ 0,
                      type=="popdev" & sel=="PRS" ~ 10,
                      type=="burnin" & sel=="GS" ~ 0,
                      type=="popdev" & sel=="GS" ~ (5+5),
                      type=="burnin" & sel=="MxAv" ~ 0,
                      type=="popdev" & sel=="MxAv" ~ 10,
                      type=="germplasmPRS" ~ 5)) %>% 
  mutate(pMfield=case_when(type=="burnin" & sel=="PRS" ~ 0,
                      type=="popdev" & sel=="PRS" ~ 15,
                      type=="burnin" & sel=="GS" ~ 0,
                      type=="popdev" & sel=="GS" ~ (10+7),
                      type=="burnin" & sel=="MxAv" ~ 0,
                      type=="popdev" & sel=="MxAv" ~ 15,
                      type=="germplasmPRS" ~ 8)) %>%
  mutate(pHfield=case_when(type=="burnin" & sel=="PRS" ~ 0,
                      type=="popdev" & sel=="PRS" ~ 20,
                      type=="burnin" & sel=="GS" ~ 0,
                      type=="popdev" & sel=="GS" ~ (15+10),
                      type=="burnin" & sel=="MxAv" ~ 0,
                      type=="popdev" & sel=="MxAv" ~ 20,
                      type=="germplasmPRS" ~ 10)) %>% 
  mutate(pLhort=case_when(type=="burnin" & sel=="PRS" ~ 0,
                           type=="popdev" & sel=="PRS" ~ 20,
                           type=="burnin" & sel=="GS" ~ 0,
                           type=="popdev" & sel=="GS" ~ (5+10),
                           type=="burnin" & sel=="MxAv" ~ 0,
                           type=="popdev" & sel=="MxAv" ~ 20,
                           type=="germplasmPRS" ~ 10)) %>% 
  mutate(pMhort=case_when(type=="burnin" & sel=="PRS" ~ 0,
                           type=="popdev" & sel=="PRS" ~ 30,
                           type=="burnin" & sel=="GS" ~ 0,
                           type=="popdev" & sel=="GS" ~ (10+15),
                           type=="burnin" & sel=="MxAv" ~ 0,
                           type=="popdev" & sel=="MxAv" ~ 30,
                           type=="germplasmPRS" ~ 15)) %>%
  mutate(pHhort=case_when(type=="burnin" & sel=="PRS" ~ 0,
                           type=="popdev" & sel=="PRS" ~ 40,
                           type=="burnin" & sel=="GS" ~ 0,
                           type=="popdev" & sel=="GS" ~ (15+20),
                           type=="burnin" & sel=="MxAv" ~ 0,
                           type=="popdev" & sel=="MxAv" ~ 40,
                           type=="germplasmPRS" ~ 20)) %>% 
  mutate(pLforest=case_when(type=="burnin" & sel=="PRS" ~ 0,
                        type=="popdev" & sel=="PRS" ~ 40,
                        type=="burnin" & sel=="GS" ~ 0,
                        type=="popdev" & sel=="GS" ~ (5+20),
                        type=="burnin" & sel=="MxAv" ~ 0,
                        type=="popdev" & sel=="MxAv" ~ 40,
                        type=="germplasmPRS" ~ 20)) %>% 
  mutate(pMforest=case_when(type=="burnin" & sel=="PRS" ~ 0,
                          type=="popdev" & sel=="PRS" ~ 50,
                          type=="burnin" & sel=="GS" ~ 0,
                          type=="popdev" & sel=="GS" ~ (10+25),
                          type=="burnin" & sel=="MxAv" ~ 0,
                          type=="popdev" & sel=="MxAv" ~ 50,
                          type=="germplasmPRS" ~ 25)) %>%
  mutate(pHforest=case_when(type=="burnin" & sel=="PRS" ~ 0,
                          type=="popdev" & sel=="PRS" ~ 60,
                          type=="burnin" & sel=="GS" ~ 0,
                          type=="popdev" & sel=="GS" ~ (15+30),
                          type=="burnin" & sel=="MxAv" ~ 0,
                          type=="popdev" & sel=="MxAv" ~ 60,
                          type=="germplasmPRS" ~ 30))


#all <- all %>% 
#          mutate(pL=recode(sel,"PRS"=10,"GS"=8,"MxAv"=10)) %>% 
#          mutate(pM=recode(sel,"PRS"=15,"GS"=12,"MxAv"=15)) %>% 
#          mutate(pH=recode(sel,"PRS"=20,"GS"=15,"MxAv"=20))

all <- allCost %>% 
          mutate(totalFieldpL=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pLfield) %>% 
          mutate(totalFieldpM=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pMfield) %>% 
          mutate(totalFieldpH=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pHfield) %>% 
          mutate(totalHortpL=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pLhort) %>% 
          mutate(totalHortpM=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pMhort) %>% 
          mutate(totalHortpH=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pHhort) %>% 
          mutate(totalForestpL=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pLforest) %>% 
          mutate(totalForestpM=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pMforest) %>% 
          mutate(totalForestpH=ceiling(factorial(Nparents)/(factorial(2)*factorial(Nparents-2)))*Nprogeny*pHforest)



write.csv(all,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/all.csv")
all <- read.csv("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/all.csv")

wrk <- all %>% group_by(NeStart,pop,Nprogeny,Nparents,h2,sel,si) %>%
  summarize(burnin.GainChangeP1=meanPz1[Cycle==40]-meanPz1[Cycle==0],
            burnin.VarChangeV1=varPz1[Cycle==40]-varPz1[Cycle==0],
            burnin.GainChangeP2=meanPz2[Cycle==40]-meanPz2[Cycle==0],
            burnin.VarChangeV2=varPz2[Cycle==40]-varPz2[Cycle==0],
            popdev.GainChangeP1=meanPz1[Cycle==79]-meanPz1[Cycle==40],
            popdev.VarChangeV1=varPz1[Cycle==79]-varPz1[Cycle==40],
            popdev.GainChangeP2=meanPz2[Cycle==79]-meanPz2[Cycle==40],
            popdev.VarChangeV2=varPz2[Cycle==79]-varPz2[Cycle==40],
            total.GainChangeP1=meanPz1[Cycle==79]-meanPz1[Cycle==0],
            total.VarChangeV1=varPz1[Cycle==79]-varPz1[Cycle==0],
            total.GainChangeP2=meanPz2[Cycle==79]-meanPz2[Cycle==0],
            total.VarChangeV2=varPz2[Cycle==79]-varPz2[Cycle==0],
            burnin.PercChangeP1=((meanPz1[Cycle==40]-meanPz1[Cycle==0])/meanPz1[Cycle==0])*100,
            burnin.PercChangeV1=((varPz1[Cycle==40]-varPz1[Cycle==0])/varPz1[Cycle==0])*100,
            burnin.PercChangeP2=((meanPz2[Cycle==40]-meanPz2[Cycle==0])/meanPz2[Cycle==0])*100,
            burnin.PercChangeV2=((varPz2[Cycle==40]-varPz2[Cycle==0])/varPz2[Cycle==0])*100,
            popdev.PercChangeP1=((meanPz1[Cycle==79]-meanPz1[Cycle==40])/meanPz1[Cycle==40])*100,
            popdev.PercChangeV1=((varPz1[Cycle==79]-varPz1[Cycle==40])/varPz1[Cycle==40])*100,
            popdev.PercChangeP2=((meanPz2[Cycle==79]-meanPz2[Cycle==40])/meanPz2[Cycle==40])*100,
            popdev.PercChangeV2=((varPz2[Cycle==79]-varPz2[Cycle==40])/varPz2[Cycle==40])*100,
            total.PercChangeP1=((meanPz1[Cycle==79]-meanPz1[Cycle==0])/meanPz1[Cycle==0])*100,
            total.PercChangeV1=((varPz1[Cycle==79]-varPz1[Cycle==0])/varPz1[Cycle==0])*100,
            total.PercChangeP2=((meanPz2[Cycle==79]-meanPz2[Cycle==0])/meanPz2[Cycle==0])*100,
            total.PercChangeV2=((varPz2[Cycle==79]-varPz2[Cycle==0])/varPz2[Cycle==0])*100,
            costLowField=sum(totalFieldpL),
            costMedField=sum(totalFieldpM),
            costHighField=sum(totalFieldpH),
            costLowHort=sum(totalHortpL),
            costMedHort=sum(totalHortpM),
            costHighHort=sum(totalHortpH),
            costLowForest=sum(totalForestpL),
            costMedForest=sum(totalForestpM),
            costHighForest=sum(totalForestpH))

write.csv(wrk,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/wrk.csv")


### Scale a single dataframe around pop="orphan",h2="hMM",sel="PRS",si=0.50
# Split the large dataframe into groupings to easily scale
wrkList <- split(wrk, list(wrk$NeStart,wrk$Nprogeny,wrk$Nparents)) # split all dataframe into subsets by scheme parameters
lapply(wrkList, \(x) data.frame(x)) |> # form and name split dataframe schemes
  setNames(paste0('schemeAll_', names(wrkList))) |>
  list2env(envir=.GlobalEnv)

# Specify the identifying values for the reference row
reference_values <- data.frame(pop = "orphan", h2 = "hHL", sel = "PRS", si = 0.50)

# Filter the reference row based on identifying values
reference_row <- schemeAll_100.1.10 %>%
  filter(pop == reference_values$pop, h2 == reference_values$h2, sel == reference_values$sel, si == reference_values$si)

# Columns not to scale
columns_not_to_scale <- colnames(schemeAll_100.1.10)[1:7]

# Dataframe of columns to scale
columns_to_scale <- schemeAll_100.1.10 %>%
  select(-one_of(columns_not_to_scale))%>%
  mutate(across(.fns = list(~ . / reference_row[[cur_column()]]), .names = "scaled_{.col}"))

# Combine the scaled columns with the columns not to scale
scaled_df <- cbind(schemeAll_100.1.10[, columns_not_to_scale], columns_to_scale)

###################################################################
### Scale by dividing by the mean of pop="orphan",h2="hHL",sel="PRS",si=0.50 occurrences for each scheme case
# Split the large dataframe into groupings to easily scale
wrkList <- split(wrk, list(wrk$NeStart,wrk$Nprogeny,wrk$Nparents)) # split all dataframe into subsets by scheme parameters
lapply(wrkList, \(x) data.frame(x)) |> # form and name split dataframe schemes
  setNames(paste0('schemeAll_', names(wrkList))) |>
  list2env(envir=.GlobalEnv)

# List all objects in the global environment
all_objects <- ls()

# Filter for data frame names containing 'schemeAll_'
df_names <- all_objects[grep('^schemeAll_', all_objects)]

# Specify the identifying values for the reference row
reference_values <- data.frame(pop = "orphan", h2 = "hHL", sel = "PRS", si = 0.50)

# Create an empty list to store the mean reference rows
mean_reference_rows <- list()

# Loop through data frame names and calculate the mean reference row for each data frame
for (df_name in df_names) {
  df <- get(df_name)  # Get the data frame by name
  
  # Filter the reference row based on identifying values
  reference_row <- df %>%
    filter(pop == reference_values$pop, h2 == reference_values$h2, sel == reference_values$sel, si == reference_values$si)
  
  if (nrow(reference_row) == 0) {
    cat("No reference row found for", df_name, "\n")
    return(NULL)
  }
  
  # Calculate the mean reference row
  mean_reference_row <- reference_row %>%
    summarise(across(.cols = everything(), .fns = mean))
  
  # Store the mean reference row in the list
  mean_reference_rows[[df_name]] <- mean_reference_row
  
  # Assign the mean reference row to the global environment with a specific name
  mean_reference_name <- paste("mean_reference_", df_name, sep = "")
  assign(mean_reference_name, mean_reference_row, envir = .GlobalEnv)
}

# Specify the identifying values for the reference row
reference_values <- data.frame(pop = "orphan", h2 = "hHL", sel = "PRS", si = 0.50)

# List all objects in the global environment
all_objects <- ls()

# Filter for data frame names containing 'schemeAll_'
df_names <- all_objects[grep('^schemeAll_', all_objects)]

# Columns not to scale
columns_not_to_scale <- colnames(schemeAll_100.1.10)[1:7]

# Create a function to perform the scaling operation for a single data frame
scale_single_data_frame <- function(df_name) {
  df <- get(df_name)  # Get the data frame by name
  
  # Filter the reference row based on identifying values
  reference_row <- df %>%
    filter(pop == reference_values$pop, h2 == reference_values$h2, sel == reference_values$sel, si == reference_values$si)
  
  if (nrow(reference_row) == 0) {
    cat("No reference row found for", df_name, "\n")
    return(NULL)
  }
  
  cat("Scaling", df_name, "\n")
  
  # Calculate the mean reference row for the current data frame
  mean_reference_name <- paste("mean_reference_", df_name, sep = "")
  mean_reference_row <- get(mean_reference_name)
  
  # Divide all rows in the selected columns by the mean reference row
  columns_to_scale <- df %>%
    select(-one_of(columns_not_to_scale)) %>%
    mutate(across(.cols = -one_of(columns_not_to_scale), .fns = list(scaled = ~ . / mean_reference_row[[cur_column()]])))
  
  # Combine the scaled columns with the columns not to scale
  scaled_df <- cbind(df[, columns_not_to_scale], columns_to_scale)
  
  # Assign the scaled data frame back to the global environment
  assign(df_name, scaled_df, envir = .GlobalEnv)
}

# Loop through data frame names and scale each one
for (df_name in df_names) {
  scale_single_data_frame(df_name)
}


###################################################################
### Combine all the scaled dataframes from above
# List all objects in the global environment
all_objects <- ls()

# Filter for data frame names containing 'schemeAll_'
df_names <- all_objects[grep('^schemeAll_', all_objects)]

# Create an empty list to store the scaled data frames
scaled_data_frames <- list()

# Loop through data frame names and collect the scaled data frames
for (df_name in df_names) {
  df <- get(df_name)  # Get the scaled data frame by name
  scaled_data_frames[[df_name]] <- df
}

# Bind the scaled data frames back together
combined_df <- do.call(bind_rows, scaled_data_frames)

### Pull unit change (gain and variance for z1 and z2) per unit cost from scaled data
unitByunit <- combined_df %>% group_by(NeStart,pop,Nprogeny,Nparents,h2,sel,si) %>%
  summarize(P1deltaBYcost_lowField=(popdev.GainChangeP1_scaled/costLowField_scaled),
            P1deltaBYcost_medField=(popdev.GainChangeP1_scaled/costMedField_scaled),
            P1deltaBYcost_highField=(popdev.GainChangeP1_scaled/costHighField_scaled),
            P1deltaBYcost_lowHort=(popdev.GainChangeP1_scaled/costLowHort_scaled),
            P1deltaBYcost_medHort=(popdev.GainChangeP1_scaled/costMedHort_scaled),
            P1deltaBYcost_HighHort=(popdev.GainChangeP1_scaled/costHighHort_scaled),
            P1deltaBYcost_lowForest=(popdev.GainChangeP1_scaled/costLowForest_scaled),
            P1deltaBYcost_medForest=(popdev.GainChangeP1_scaled/costMedForest_scaled),
            P1deltaBYcost_highForest=(popdev.GainChangeP1_scaled/costHighForest_scaled),
            V1deltaBYcost_lowField=(popdev.VarChangeV1_scaled/costLowField_scaled),
            V1deltaBYcost_medField=(popdev.VarChangeV1_scaled/costMedField_scaled),
            V1deltaBYcost_highField=(popdev.VarChangeV1_scaled/costHighField_scaled),
            V1deltaBYcost_lowHort=(popdev.VarChangeV1_scaled/costLowHort_scaled),
            V1deltaBYcost_medHort=(popdev.VarChangeV1_scaled/costMedHort_scaled),
            V1deltaBYcost_HighHort=(popdev.VarChangeV1_scaled/costHighHort_scaled),
            V1deltaBYcost_lowForest=(popdev.VarChangeV1_scaled/costLowForest_scaled),
            V1deltaBYcost_medForest=(popdev.VarChangeV1_scaled/costMedForest_scaled),
            V1deltaBYcost_highForest=(popdev.VarChangeV1_scaled/costHighForest_scaled),
            P2deltaBYcost_lowField=(popdev.GainChangeP2_scaled/costLowField_scaled),
            P2deltaBYcost_medField=(popdev.GainChangeP2_scaled/costMedField_scaled),
            P2deltaBYcost_highField=(popdev.GainChangeP2_scaled/costHighField_scaled),
            P2deltaBYcost_lowHort=(popdev.GainChangeP2_scaled/costLowHort_scaled),
            P2deltaBYcost_medHort=(popdev.GainChangeP2_scaled/costMedHort_scaled),
            P2deltaBYcost_HighHort=(popdev.GainChangeP2_scaled/costHighHort_scaled),
            P2deltaBYcost_lowForest=(popdev.GainChangeP2_scaled/costLowForest_scaled),
            P2deltaBYcost_medForest=(popdev.GainChangeP2_scaled/costMedForest_scaled),
            P2deltaBYcost_highForest=(popdev.GainChangeP2_scaled/costHighForest_scaled),
            V2deltaBYcost_lowField=(popdev.VarChangeV2_scaled/costLowField_scaled),
            V2deltaBYcost_medField=(popdev.VarChangeV2_scaled/costMedField_scaled),
            V2deltaBYcost_highField=(popdev.VarChangeV2_scaled/costHighField_scaled),
            V2deltaBYcost_lowHort=(popdev.VarChangeV2_scaled/costLowHort_scaled),
            V2deltaBYcost_medHort=(popdev.VarChangeV2_scaled/costMedHort_scaled),
            V2deltaBYcost_HighHort=(popdev.VarChangeV2_scaled/costHighHort_scaled),
            V2deltaBYcost_lowForest=(popdev.VarChangeV2_scaled/costLowForest_scaled),
            V2deltaBYcost_medForest=(popdev.VarChangeV2_scaled/costMedForest_scaled),
            V2deltaBYcost_highForest=(popdev.VarChangeV2_scaled/costHighForest_scaled))
            


#####################################################
### Cost per percentage change in gain and variance in population development
costPerPerc <- wrk %>% group_by(NeStart,pop,Nprogeny,Nparents,h2,sel,si) %>% 
              summarize(low.costPerPerc_gain1=(costLow/popdev.PercChangeP1),
                        med.costPerPerc_gain1=(costMed/popdev.PercChangeP1),
                        high.costPerPerc_gain1=(costHigh/popdev.PercChangeP1),
                        low.costPerPerc_gain2=(costLow/popdev.PercChangeP2),
                        med.costPerPerc_gain2=(costLow/popdev.PercChangeP2),
                        high.costPerPerc_gain2=(costLow/popdev.PercChangeP2),
                        low.costPerPerc_var1=(costLow/popdev.PercChangeV1),
                        med.costPerPerc_var1=(costLow/popdev.PercChangeV1),
                        high.costPerPerc_var1=(costLow/popdev.PercChangeV1),
                        low.costPerPerc_var2=(costLow/popdev.PercChangeV2),
                        med.costPerPerc_var2=(costLow/popdev.PercChangeV2),
                        high.costPerPerc_var2=(costLow/popdev.PercChangeV2))


####################################################
### Visualization of all parameter combinations
library(ggplot2)
all$iter <- rep(1:(nrow(all)/81),each=81)
allsub <- subset(all,NeStart==25)
allSummary$iter <- rep(1:(nrow(allSummary)/81),each=81)
allSumsub <- subset(allSummary,NeStart==25)
allsub$si <- as.factor(allsub$si)

# Get means of each cycle by the different parameter grouping
popSum <- all %>% group_by(pop,NeStart,Cycle) %>% 
            summarize(meanP1=mean(meanPz1),
                      meanP2=mean(meanPz2),
                      sdP1=sd(meanPz1),
                      sdP2=sd(meanPz2),
                      var1=mean(varPz1),
                      var2=mean(varPz2),
                      sdV1=sd(varPz1),
                      sdV2=sd(varPz2))
popSumsub <- subset(popSum,NeStart==25)
h2Sum <- all %>% group_by(h2,NeStart,Cycle) %>% 
            summarize(meanP1=mean(meanPz1),
                      meanP2=mean(meanPz2),
                      sdP1=sd(meanPz1),
                      sdP2=sd(meanPz2),
                      var1=mean(varPz1),
                      var2=mean(varPz2),
                      sdV1=sd(varPz1),
                      sdV2=sd(varPz2))
h2Sumsub <- subset(h2Sum,NeStart==25)
selSum <- all %>% group_by(sel,NeStart,Cycle) %>% 
            summarize(meanP1=mean(meanPz1),
                      meanP2=mean(meanPz2),
                      sdP1=sd(meanPz1),
                      sdP2=sd(meanPz2),
                      var1=mean(varPz1),
                      var2=mean(varPz2),
                      sdV1=sd(varPz1),
                      sdV2=sd(varPz2))
selSumsub <- subset(selSum,NeStart==25)
siSum <- all %>% group_by(si,NeStart,Cycle) %>% 
            summarize(meanP1=mean(meanPz1),
                      meanP2=mean(meanPz2),
                      sdP1=sd(meanPz1),
                      sdP2=sd(meanPz2),
                      var1=mean(varPz1),
                      var2=mean(varPz2),
                      sdV1=sd(varPz1),
                      sdV2=sd(varPz2))
siSumsub <- subset(siSum,NeStart==25)
siSumsub$si <- as.factor(siSumsub$si)
popselSum <- all %>% group_by(pop,sel,NeStart,Cycle) %>% 
  summarize(meanP1=mean(meanPz1),
            meanP2=mean(meanPz2),
            sdP1=sd(meanPz1),
            sdP2=sd(meanPz2),
            var1=mean(varPz1),
            var2=mean(varPz2),
            sdV1=sd(varPz1),
            sdV2=sd(varPz2))
popselSumsub <- subset(popselSum,NeStart==25)

pop.labels <- c("Landrace","Orphan","Wild")
pop.colors <- c("#009E73","#D55E00","#CC79A7")
h2.colors <- c("#009E73","#D55E00","#CC79A7","#0072B2")
sel.colors <- c("#009E73","#D55E00","#CC79A7")
si.colors <- c("#009E73","#D55E00","#CC79A7")
popsel.colors <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
                   "#AA4499","#44AA99","#999933","#882255")

pP1 = ggplot() +
  geom_line(allsub,mapping=aes(x = Cycle, y = meanPz1, group = iter, color=pop),linewidth = 0.2, alpha = 0.1) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = meanP1, group = pop,color=pop),linewidth=2) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = meanP1, group = pop),color="black",linewidth=1,alpha=0.25) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 40, linetype = "dashed", color = "black") +
  scale_color_manual(values=pop.colors) +
  labs(x="Cycle",y="Phenotypic Value (z1)") +
  labs(color="Population")

pP2 = ggplot() +
  geom_line(allsub,mapping=aes(x = Cycle, y = meanPz2, group = iter, color=pop),linewidth = 0.2, alpha = 0.1) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = meanP2, group = pop,color=pop),linewidth=2) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = meanP2, group = pop),color="black",linewidth=1,alpha=0.25) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 40, linetype = "dashed", color = "black") +
  scale_color_manual(labels=pop.labels,values=pop.colors) +
  labs(x="Cycle",y="Phenotypic Value (z2)") +
  labs(color="Population")

pV1 = ggplot() +
  geom_line(allsub,mapping=aes(x = Cycle, y = varPz1, group = iter, color=pop),linewidth = 0.2, alpha = 0.1) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = var1, group = pop,color=pop),linewidth=2) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = var1, group = pop),color="black",linewidth=1,alpha=0.25) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 40, linetype = "dashed", color = "black") +
  scale_color_manual(labels=pop.labels,values=pop.colors) +
  labs(x="Cycle",y="Phenotypic Variance (z1)") +
  labs(color="Population")

pV2 = ggplot() +
  geom_line(allsub,mapping=aes(x = Cycle, y = varPz2, group = iter, color=pop),linewidth = 0.2, alpha = 0.1) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = var2, group = pop,color=pop),linewidth=2) +
  geom_line(popSumsub,mapping=aes(x = Cycle, y = var2, group = pop),color="black",linewidth=1,alpha=0.25) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 40, linetype = "dashed", color = "black") +
  scale_color_manual(labels=pop.labels,values=pop.colors) +
  labs(x="Cycle",y="Phenotypic Variance (z2)") +
  labs(color="Population")

####################################################
### Calculate BLUPs and BLUEs
library(lme4)       
out1 <- lmer(P1deltaBYcost_HighHort~0+pop+sel*h2+h2*si++(1|Nparents)+(1|Nprogeny)+(1|NeStart),data=unitByunit)
summary(out1)
anova(out1)

out2 <- lmer(V1deltaBYcost_HighHort~0+pop+sel*h2+h2*si++(1|Nparents)+(1|Nprogeny)+(1|NeStart),data=unitByunit)
summary(out2)
anova(out2)

out3 <- lmer(P2deltaBYcost_HighHort~0+pop+sel*h2+h2*si+(1|Nparents)+(1|Nprogeny)+(1|NeStart),data=unitByunit)
summary(out3)
anova(out3)

out4 <- lmer(V2deltaBYcost_HighHort~0+pop+sel*h2+h2*si++(1|Nparents)+(1|Nprogeny)+(1|NeStart),data=unitByunit)
summary(out4)
anova(out4)

library(lme4)
library(dplyr)
#formula <- paste(response_var, "~ pop + sel + h2 + si + (1|Nparents) + (1|Nprogeny) + (1|NeStart)")
# List of response variable names
response_vars <- c(colnames(unitByunit[8:43]))

# Create empty lists to store fixed and random effects
fixed_effects_list <- list()
random_effects_list <- list()
unitByunit$si <- as.factor(unitByunit$si)

# Loop through each response variable and fit a mixed-effects model
for (response_var in response_vars) {
  formula <- paste(response_var, "~0 + pop + sel*h2 + h2*si + (1|Nparents) + (1|Nprogeny) + (1|NeStart)")
  model_fit <- lmer(formula, data = unitByunit, REML = FALSE)  # You can set REML to TRUE if needed
  
  # Extract fixed and random effects and store them in lists
  fixed_effects_list[[response_var]] <- fixef(model_fit)
  random_effects_list[[response_var]] <- ranef(model_fit)
}

# Save fixed effects and random effects as separate objects
for (response_var in response_vars) {
  fixed_effects <- fixed_effects_list[[response_var]]
  random_effects <- random_effects_list[[response_var]]
  
  # Save fixed effects as a separate object (replace "fixed_effects_response1", etc., with appropriate names)
  assign(paste("fixed_effects_", response_var, sep = ""), fixed_effects)
  
  # Save random effects as a separate object (replace "random_effects_response1", etc., with appropriate names)
  assign(paste("random_effects_", response_var, sep = ""), random_effects)
}

# Create dataframe of Random Effects for each mixed model
ranef_df <- do.call(cbind, lapply(names(random_effects_list), function(response_var) {
  df_list <- random_effects_list[[response_var]]
  df <- do.call(rbind, df_list)
  colnames(df) <- paste(response_var, colnames(df), sep = "_")
  return(df)
}))
# Set the column names as the lists
colnames(ranef_df) <- names(random_effects_list)


# Create an empty dataframe with column names
fixed_effects_df <- data.frame(matrix(NA, nrow = length(fixed_effects_list), ncol = length(unlist(fixed_effects_list[[1]]))))
colnames(fixed_effects_df) <- names(unlist(fixed_effects_list[[1]]))

# Fill in the dataframe with values
for (i in 1:length(fixed_effects_list)) {
  fixed_effects_df[i, ] <- unlist(fixed_effects_list[[i]])
}

# Add row names
rownames(fixed_effects_df) <- names(fixed_effects_list)
fixed_effects_df <- t(fixed_effects_df)
fixef_df <- data.frame(fixed_effects_df)
rownames(fixef_df)[1] <- "poplandrace.selGS.h2hHL.si0.25"
rownames(fixef_df)[2] <- "poporphan.selGS.h2hHL.si0.25"
rownames(fixef_df)[3] <- "popwild.selGS.h2hHL.si0.25"

write.csv(ranef_df,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/ranef_df.csv")
# Get BLUES for all schemes
write.csv(fixef_df,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/fixef_df.csv")
# Done in Excel
blues <- readxl::read_xlsx("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/blues.xlsx")
blues <- blues[-c(4:22),]

# Split 'scheme' column into the 4 identifying criteria
split_strings <- strsplit(blues$scheme, "\\.")
# Create a new dataframe with the split values in separate columns
new_df <- data.frame(do.call(rbind, split_strings))
# Combine the two columns with a period between them
new_df$X6 <- paste(new_df$X4, new_df$X5, sep = ".")
new_df <- new_df[,-c(4,5)]
colnames(new_df) <- c("pop", "sel", "h2", "si")
new_df$pop <- sub("^pop", "", new_df$pop)
new_df$sel <- sub("^sel", "", new_df$sel)
new_df$h2 <- sub("^h2", "", new_df$h2)
new_df$si <- sub("^si", "", new_df$si)

blues <- cbind(new_df,blues)
write.csv(blues,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/bluesProc.csv")

library(tidyr)
long_blues <- pivot_longer(blues, cols = 6:41, names_to = "goalBYcost", values_to = "value")
long_blues$pop.sel <- paste(long_blues$pop,long_blues$sel,sep=".")
split_df <- split(long_blues, sub("_.*", "", long_blues$goalBYcost))
P1 <- split_df[["P1deltaBYcost"]]
P2 <- split_df[["P2deltaBYcost"]]
V1 <- split_df[["V1deltaBYcost"]]
V2 <- split_df[["V2deltaBYcost"]]

P1$pop.sel <- factor(P1$pop.sel,levels=c("landrace.PRS","landrace.MxAv","landrace.GS",
                                         "orphan.PRS","orphan.MxAv","orphan.GS",
                                         "wild.PRS","wild.MxAv","wild.GS"))
V1$pop.sel <- factor(V1$pop.sel,levels=c("landrace.PRS","landrace.MxAv","landrace.GS",
                                         "orphan.PRS","orphan.MxAv","orphan.GS",
                                         "wild.PRS","wild.MxAv","wild.GS"))
P2$pop.sel <- factor(P2$pop.sel,levels=c("landrace.PRS","landrace.MxAv","landrace.GS",
                                         "orphan.PRS","orphan.MxAv","orphan.GS",
                                         "wild.PRS","wild.MxAv","wild.GS"))
V2$pop.sel <- factor(V2$pop.sel,levels=c("landrace.PRS","landrace.MxAv","landrace.GS",
                                         "orphan.PRS","orphan.MxAv","orphan.GS",
                                         "wild.PRS","wild.MxAv","wild.GS"))



P1 <- P1 %>% mutate(h2 = recode(h2,hHL = "h2High",hLH = 'h2Low',hLL = 'h2Low',hMM = 'h2Med'))
V1 <- V1 %>% mutate(h2 = recode(h2,hHL = "h2High",hLH = 'h2Low',hLL = 'h2Low',hMM = 'h2Med'))
P2 <- P2 %>% mutate(h2 = recode(h2,hHL = "h2Low",hLH = 'h2High',hLL = 'h2Low',hMM = 'h2Med'))
V2 <- V2 %>% mutate(h2 = recode(h2,hHL = "h2Low",hLH = 'h2High',hLL = 'h2Low',hMM = 'h2Med'))

P1$si <- as.factor(P1$si)
V1$si <- as.factor(V1$si)
P2$si <- as.factor(P2$si)
V2$si <- as.factor(V2$si)

P1 <- P1 %>% mutate(si = recode(si, '0.25'="siLow",'0.5'="siMed",'0.75'="siHigh"))
V1 <- V1 %>% mutate(si = recode(si, '0.25'="siLow",'0.5'="siMed",'0.75'="siHigh"))
P2 <- P2 %>% mutate(si = recode(si, '0.25'="siHigh",'0.5'="siMed",'0.75'="siLow"))
V2 <- V2 %>% mutate(si = recode(si, '0.25'="siHigh",'0.5'="siMed",'0.75'="siLow"))

P1$h2.si <- paste(P1$h2,P1$si,sep=".")
V1$h2.si <- paste(V1$h2,V1$si,sep=".")
P2$h2.si <- paste(P2$h2,P2$si,sep=".")
V2$h2.si <- paste(V2$h2,V2$si,sep=".")

# Plot the BLUES as a facet wrapped around the trait by cost increments (columns)
library(ggplot2)
ggplot(P1, aes(x=scheme,y=value,fill=interaction(pop,sel))) +
  geom_bar(stat="identity") +
  facet_wrap(~goalBYcost,ncol=2,
             labeller=labeller(goalBYcost=c("P1deltaBYcost_highForest"="Forestry High","P1deltaBYcost_medForest"="Forestry Medium","P1deltaBYcost_lowForest"="Forestry Low",
                                            "P1deltaBYcost_HighHort"="Hort High","P1deltaBYcost_medHort"="Hort Medium","P1deltaBYcost_lowHort"="Hort Low",
                                            "P1deltaBYcost_highField"="Field High","P1deltaBYcost_medField"="Field Medium","P1deltaBYcost_lowField"="Field Low"))) +
  labs(x="Scheme",y="Value",title="Simple Oligogenic Trait (z1) Gain by Cost Unit Change") +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=4),
        strip.text.x=element_text(angle=0,hjust=0.5),
        axis.text.y=element_text(size=5)) +
  labs(fill="Population and Selection") +
  scale_fill_manual(labels=c("Landrace GS","Orphan GS","Wild GS","Landrace MxAv","Orphan MxAv",
                             "Wild MxAv","Landrace PRS","Orphan PRS","Wild PRS"),
                    values=c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
                             "#AA4499","#44AA99","#999933","#882255"))


# Plot the BLUEs as a facet wrapped around the 
library(ggplot2)
ggplot(V2,aes(x=pop.sel,y=value,)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=h2.si),size=1,position=position_jitter(width=.15,height=0.15)) +
  labs(x="Population and Selection",y="Baseline Scaled Unit",title="Complex Oligogenic Trait (z2) Scaled Variance by Cost Unit Change") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.0001,size=12),
        strip.text.x=element_text(angle=0,hjust=0.5),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=15)) +
  scale_x_discrete(labels=c("Landrace PRS","Landrace MxAv","Landrace GS",
                            "Orphan PRS","Orphan MxAv","Orphan GS",
                            "Wild PRS","Wild MxAv","Wild GS")) +
  labs(color="Heritability and Index Weighting") +
  scale_color_manual(labels=c("High and High","High and Low","High and Med","Low and High","Low and Low",
                             "Low and Med","Med and High","Med and Low","Med and Med"),
                    values=c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
                             "#AA4499","#44AA99","#999933","#882255"))


