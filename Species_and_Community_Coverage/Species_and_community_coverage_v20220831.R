#This script was used to calculate the species and community coverage
#of the agricultural soil data set as described in the manuscript:
#Tremblay,J., Schreiber,L., Greer, C.W. 2022. "High resolution shotgun metagenonimcs: the more data the better?"

#The script was executed under the following R enviroment
#R version 4.0.5 (2021-03-31)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS High Sierra 10.13.6

#And with the attached packages:
#[1] ggplot2_3.3.5    entropart_1.6-11 vegan_2.5-7      lattice_0.20-44 
#[5] permute_0.9-5    readr_1.4.0

library("readr")
library("vegan")
library("entropart")
library("ggplot2")

#Load all rpoB count tables (supplied as sample data with this script
agri_ct_all <- readRDS(file = "ct_rpoB_all_20220819.RData")
agri_ct_12M <- readRDS(file = "ct_rpoB_12M_20220819.RData")
agri_ct_8M <- readRDS(file = "ct_rpoB_8M_20220819.RData")
agri_ct_4M <- readRDS(file = "ct_rpoB_4M_20220819.RData")
agri_ct_1M <- readRDS(file = "ct_rpoB_1M_20220819.RData")
agri_ct_05M <- readRDS(file = "ct_rpoB_0.5M_20220819.RData")
agri_ct_01M <- readRDS(file = "ct_rpoB_0.1M_20220819.RData")

#For each rpoB count table create a dataframe that summarizes the read counts, the species counts and the community evenness

#all clusters
{
  #Load the sequencing depth file of the original (not subsampled) dataset to calculate the mean cluster count of the "all" dataset
  all_depth_df <- read_tsv(file = "Agricultural_soil_sequencing_depth.all_clusters.tsv")
  mean_cluster_count <- mean(all_depth_df$Bases/all_depth_df$AvgSpotLen)
  #-> ca 15.6 million clusters
  
  #Calculate richness and community evennes with the "vegan" package
  H_sub <- diversity(t(agri_ct_all))
  S_sub <- specnumber(t(agri_ct_all))
  J_sub <- H_sub/log(S_sub)
  
  #Create summary dataframe
  stats_all <- data.frame("sample" = colnames(agri_ct_all),
                          "reads" = colSums(agri_ct_all),
                          "S_obs" = specnumber(t(agri_ct_all)),
                          "S_coverage_sample" = specnumber(t(agri_ct_all)) / specnumber(t(agri_ct_all)),
                          "S_coverage_meta" = specnumber(t(agri_ct_all)) / specnumber(colSums(t(agri_ct_all))),
                          "effort" = 15.6,
                          "J" = J_sub)
  }


#12M
{
  H_sub <- diversity(t(agri_ct_12M))
  S_sub <- specnumber(t(agri_ct_12M))
  J_sub <- H_sub/log(S_sub)
  
   stats_12M <- data.frame("sample" = colnames(agri_ct_12M),
                          "reads" = colSums(agri_ct_12M),
                          "S_obs" = specnumber(t(agri_ct_12M)),
                          "S_coverage_sample" = specnumber(t(agri_ct_12M)) / specnumber(t(agri_ct_all)),
                          "S_coverage_meta" = specnumber(t(agri_ct_12M)) / specnumber(colSums(t(agri_ct_all))),
                          "effort" = 12,
                          "J" = J_sub)
}

#8M
{
  H_sub <- diversity(t(agri_ct_8M))
  S_sub <- specnumber(t(agri_ct_8M))
  J_sub <- H_sub/log(S_sub)
  
  stats_8M <- data.frame("sample" = colnames(agri_ct_8M),
                          "reads" = colSums(agri_ct_8M),
                          "S_obs" = specnumber(t(agri_ct_8M)),
                          "S_coverage_sample" = specnumber(t(agri_ct_8M)) / specnumber(t(agri_ct_all)),
                          "S_coverage_meta" = specnumber(t(agri_ct_8M)) / specnumber(colSums(t(agri_ct_all))),
                          "effort" = 8,
                         "J" = J_sub)
}


#4M
{
  H_sub <- diversity(t(agri_ct_4M))
  S_sub <- specnumber(t(agri_ct_4M))
  J_sub <- H_sub/log(S_sub)
  
  stats_4M <- data.frame("sample" = colnames(agri_ct_4M),
                         "reads" = colSums(agri_ct_4M),
                         "S_obs" = specnumber(t(agri_ct_4M)),
                         "S_coverage_sample" = specnumber(t(agri_ct_4M)) / specnumber(t(agri_ct_all)),
                         "S_coverage_meta" = specnumber(t(agri_ct_4M)) / specnumber(colSums(t(agri_ct_all))),
                         "effort" = 4,
                         "J" = J_sub)
}


#1M
{
  H_sub <- diversity(t(agri_ct_1M))
  S_sub <- specnumber(t(agri_ct_1M))
  J_sub <- H_sub/log(S_sub)
  
  
  stats_1M <- data.frame("sample" = colnames(agri_ct_1M),
                         "reads" = colSums(agri_ct_1M),
                         "S_obs" = specnumber(t(agri_ct_1M)),
                         "S_coverage_sample" = specnumber(t(agri_ct_1M)) / specnumber(t(agri_ct_all)),
                         "S_coverage_meta" = specnumber(t(agri_ct_1M)) / specnumber(colSums(t(agri_ct_all))),
                         "effort" = 1,
                         "J" = J_sub)
}

#0.5M
{
  H_sub <- diversity(t(agri_ct_05M))
  S_sub <- specnumber(t(agri_ct_05M))
  J_sub <- H_sub/log(S_sub)
  
  stats_05M <- data.frame("sample" = colnames(agri_ct_05M),
                         "reads" = colSums(agri_ct_05M),
                         "S_obs" = specnumber(t(agri_ct_05M)),
                         "S_coverage_sample" = specnumber(t(agri_ct_05M)) / specnumber(t(agri_ct_all)),
                         "S_coverage_meta" = specnumber(t(agri_ct_05M)) / specnumber(colSums(t(agri_ct_all))),
                         "effort" = 0.5,
                         "J" = J_sub)
}

#0.1M
{
  H_sub <- diversity(t(agri_ct_01M))
  S_sub <- specnumber(t(agri_ct_01M))
  J_sub <- H_sub/log(S_sub)
  
  stats_01M <- data.frame("sample" = colnames(agri_ct_01M),
                          "reads" = colSums(agri_ct_01M),
                          "S_obs" = specnumber(t(agri_ct_01M)),
                          "S_coverage_sample" = specnumber(t(agri_ct_01M)) / specnumber(t(agri_ct_all)),
                          "S_coverage_meta" = specnumber(t(agri_ct_01M)) / specnumber(colSums(t(agri_ct_all))),
                          "effort" = 0.1,
                          "J" = J_sub)
}

#Combine all the individual dataframes into one
stats_combined <- rbind.data.frame(stats_all,stats_12M,stats_8M,stats_4M,
                                   stats_1M,stats_05M,stats_01M)


#Estimate the community coverages of the subsets on a per sample basis

#Create results column in combined dataframe
stats_combined$com_coverage_sample <- NA

#Some of the subsets did not contain a single read mapped to an rpoB gene, we'll have to add pseudo-count to not break "Coverage" function ("entropart" package), which does not work for levels <2
stats_combined$reads[stats_combined$reads <=1] <- 2

#Create a list of all samples to iterate through
sample_list <- unique(stats_combined$sample)

#Iterate through samples
for(i in sample_list){
  #Iterate through all total read counts of this sample
  for(j in stats_combined$reads[stats_combined$sample == i]){
    #Interpolate the community coverage using the "all" sample as the reference and the subset read count as the level
    stats_combined$com_coverage_sample[stats_combined$sample == i &
                                         stats_combined$reads == j] <- 
      Coverage(Ns = agri_ct_all[,colnames(agri_ct_all) %in% i], Level = j)
  }
}

#Plot the data as a boxplot using ggplot

#Species coverage
ggplot(data = stats_combined) + 
  geom_boxplot(mapping = aes(x = factor(effort),
                             y = S_coverage_sample*100),
               fill = "#4477AA")  +
  theme_classic() +
  scale_y_continuous(name = "Species coverage [%]",
                     breaks = seq(0,100,by = 10)) +
  scale_x_discrete(name = "Subset") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size =14, color = "black"),
        axis.text.y = element_text(size =14, color = "black"))

#Community coverage
ggplot(data = stats_combined) + 
  geom_boxplot(mapping = aes(x = factor(effort),
                             y = com_coverage_sample*100),
               fill = "#4477AA")  +
  theme_classic() +
  scale_y_continuous(name = "Community coverage [%]",
                     breaks = seq(0,100,by = 10)) +
  scale_x_discrete(name = "Subset") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size =14, color = "black"),
        axis.text.y = element_text(size =14, color = "black"))


#Evenness
ggplot(data = stats_combined) + 
  geom_boxplot(mapping = aes(x = factor(effort),
                             y = J),
               fill = "#4477AA")  +
  theme_classic() +
  scale_y_continuous(name = "Eveness J",
                     breaks = seq(0,1,by = 0.1)) +
  scale_x_discrete(name = "Subset") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size =14, color = "black"),
        axis.text.y = element_text(size =14, color = "black"))
