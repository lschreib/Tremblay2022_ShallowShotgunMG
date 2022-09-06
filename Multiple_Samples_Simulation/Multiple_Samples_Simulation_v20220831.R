#This script was used to calculate the species and community coverage
#of the agricultural soil data set as described in the manuscript:
#Tremblay,J., Schreiber,L., Greer, C.W. 2022. "High resolution shotgun metagenonimcs: the more data the better?"

#The script was executed under the following R enviroment
#R version 4.0.5 (2021-03-31)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS High Sierra 10.13.6

#And with the attached packages:
#[1] ggplot2_3.3.5    vegan_2.5-7      lattice_0.20-44 
#[5] permute_0.9-5    readr_1.4.0

library("vegan")
library("ggplot2")

#Load the rpoB count table (supplied as sample data with this script)
ct_all <- readRDS(file = "Agricultural_Soil_rpoB_ct_20220831.RData")

#Rarefy the samples of the count table to a common depth, i.e. the depth of the sample with the fewest reads
rare_threshold <- min(colSums(ct_all))
ct_all_rare <- t(rrarefy(t(ct_all), sample = rare_threshold))

#Remove all species that now have a 0 count in all samples
ct_all_rare <- ct_all_rare[rowSums(ct_all_rare) > 0,]

#Calculate the overall richness of the meta-community represented by the rpoB count table
overall_richness <- length(rownames(ct_all_rare))

#Initiate the output dataframe of the simulation
out_df <- data.frame("Combination" = NA,
                     "Replicate" = NA,
                     "Depth" = NA,
                     "S" = NA,
                     "S_coverage" = NA,
                     "com_coverage" = NA,
                     "Eveness" = NA)

#Calculate the relative abundance of species' across across all samples, i.e. within the meta-community
met_com <- data.frame("otu" = rownames(ct_all_rare),
                      "rel_abund" = rowSums(ct_all_rare)/sum(rowSums(ct_all_rare)))

#Initiate some values of the simulation loop
number_of_samples <- ncol(ct_all_rare)
closed_sum <- rare_threshold
#The "closed_sum" is represents the constant (ceiling) sequencing depth for all simulated combination datasets

#Outer loop to iterate through combinations of 1,2, ...,nmax samples
for(a in seq(1,number_of_samples,by=1)){
  #Initiate a temporary matrix for the iteration that will store the species counts
  it_out <- matrix(nrow = number_of_samples,ncol=a)
  #Create the list of sample names that we will randomly draw from
  sample_vector <- colnames(ct_all_rare)
  
  #Inner loop 1 to generate "number_of_samples" random sample combinations
  for (y in seq(1,number_of_samples,by=1))
  {
    it_out[y,] <- sample(sample_vector,size = a,replace = FALSE)
  } 
  
  #If we are at a = 1, each sample combinations is represented only by one sample,
  #This makes sure that really every sample is used when a=1
  if(a==1){
    it_out <- matrix(nrow = number_of_samples,ncol=1)
    it_out[,1] <- sample_vector
  }
  
  #Inner loop 2 that iterates through all generated sample combinations and calculates the key metrics for these
  for(b in seq(1,nrow(it_out),by=1)){
    
    #Generate the combined community of the sample combination
    if(a == 1){
      achieved_com <- data.frame("otu" = rownames(ct_all_rare),
                                 "abund" = ct_all_rare[,it_out[b,]]) }
    else {
      achieved_com <- data.frame("otu" = rownames(ct_all_rare),
                                 "abund" = rowSums(ct_all_rare[,it_out[b,]]))}
    
    #Calculate the species richness and coverage of the sample combination with the constraint of the constant sequencing depth
    S_sub <- rarefy(achieved_com$abund,sample = closed_sum)
    S_coverage <- S_sub/overall_richness
    
    #If we draw from this community, how much of the meta-community
    #will be covered?
    
    #Rarefy the community of the sample conbintaion to the ceiling sampling depth 
    achieved_com_rare <- data.frame("otu" = rownames(met_com),
                                    "abund_rare" = t(rrarefy(achieved_com$abund, sample = closed_sum)))
    
    #Merge the sample combination counts information with the information of the meta-community
    temp_df <- merge(met_com,achieved_com_rare,by = "otu",all.x = TRUE)
    #Convert the species counts of the sample combination into presence/absence information
    temp_df$abund_rare[temp_df$abund_rare > 0] <- 1
    #Calculate the community coverage (relative to the meta-community) of the sample community
    com_coverage <- sum(temp_df$abund_rare * temp_df$rel_abund)
    
    #Calculate the Shannon index of the sample combination
    H <- diversity(achieved_com_rare$abund_rare)
    #Calculate the species richness of the sample combination
    S <- specnumber(achieved_com_rare$abund_rare)
    #Calculate the community evenness of the sample combination
    J_rare <- H/log(S)
    
    #Put it all the metrics into one df
    out_df_tmp <- data.frame("Combination" = a,
                             "Replicate" = b,
                             "Depth" = closed_sum,
                             "S" = round(S_sub,digits = 0),
                             "S_coverage" = S_coverage,
                             "com_coverage" = com_coverage,
                             "Eveness" = J_rare)
    
    #Append the data of the current sample combination to the overall output dataframe
    out_df <- rbind.data.frame(out_df,out_df_tmp)
    #Print where we are in the simulation
    print(paste("k =",a,"Replicate =",b))
  }
}
#Remove the NA column of the output dataframe

out_df <- out_df[!(is.na(out_df$Combination)),]

#Finally, let's plot the data with ggplot
#First community coverage
ggplot(data = out_df,
       mapping = aes(x = Combination,
                     y = com_coverage*100)) +
  geom_point(color = "#4477AA") +
  stat_summary(fun = "mean",
               geom = "line",
               size = 1.5) +
  theme_classic() +
  scale_x_continuous(breaks = c(1,25,50,71), 
                     name = "# of samples") +
  scale_y_continuous(name = "Community coverage [%]",
                     breaks = seq(0,100,by = 25)) +
  theme(axis.text.x = element_text(size =14, color = "black"),
        axis.text.y = element_text(size =14, color = "black"),
        axis.ticks = element_line(color = "black")) +
  coord_cartesian(ylim = c(0,100))

#Now richness
ggplot(data = out_df,
       mapping = aes(x = Combination,
                     y = S_coverage*100)) +
  geom_point(color = "#4477AA") +
  stat_summary(fun = "mean",
               geom = "line",
               size = 1.5) +
  theme_classic() +
  scale_x_continuous(breaks = c(1,25,50,71), 
                     name = "# of samples") +
  scale_y_continuous(name = "Species coverage [%]",
                     breaks = seq(0,100,by = 25)) +
  theme(axis.text.x = element_text(size =14, color = "black"),
        axis.text.y = element_text(size =14, color = "black"),
        axis.ticks = element_line(color = "black")) +
  coord_cartesian(ylim = c(0,100))

#And finally evenness J
ggplot(data = out_df,
       mapping = aes(x = Combination,
                     y = Eveness)) +
  geom_point(color = "#4477AA") +
  stat_summary(fun = "mean",
               geom = "line",
               size = 1.5) +
  theme_classic() +
  scale_x_continuous(breaks = seq(1,71,by =5), 
                     name = "# of samples") +
  scale_y_continuous(name = "Community eveness J") +
  theme(axis.text.x = element_text(size =14, color = "black"),
        axis.text.y = element_text(size =14, color = "black"),
        axis.ticks = element_line(color = "black"))