
# install.packages("PerformanceAnalytics") # Turn on to install current version
library(PerformanceAnalytics)

#Set up paths
data_dir <- "G:/Shared drives/IMLS MFA/occurrence_points/outputs/exsitu_coverage"
output_dir <- "G:/Shared drives/IMLS MFA/Genetic diversity value"
	# on a Mac
	#data_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/occurrence_points/outputs/exsitu_coverage"
	#output_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/Genetic diversity value"


#Read in file
eco_geo_results<-read.csv("ExSituCoverage_BufferTable_6_30_21.csv")
eco_geo_results<-read.csv(file.path(data_dir,"ExSituCoverage_BufferTable_6_30_21.csv"))
#Pull species names
sp_names_all<-eco_geo_results[,1]
#Identify species with no collections
no_collections<-sp_names_all[which(is.na(eco_geo_results[,2]))]
sp_names_wcoll<-sp_names_all[-which(is.na(eco_geo_results[,2]))]

sp_names_TH<-sp_names_all[intersect(which(eco_geo_results[,14]!="LC"),which(!is.na(eco_geo_results[,2])))]
sp_names_Q_TH<-sp_names_all[intersect(intersect(which(eco_geo_results[,14]!="LC"),which(!is.na(eco_geo_results[,2]))),which(substr(eco_geo_results[,1],1,1)=="Q"))]
sp_names_usa<-sp_names_all[intersect(which(!is.na(eco_geo_results[,9])),which(!is.na(eco_geo_results[,2])))]

#Remove species with no collections
eco_geo_results<-eco_geo_results[-which(is.na(eco_geo_results[,2])),]

#Option to focus on just RL Threatened or not
eco_geo_results[,ncol(eco_geo_results)]

#Get just the percentage
eco_geo_results[,2:(ncol(eco_geo_results)-2)]<-as.numeric(gsub("\\%.*","",as.matrix(eco_geo_results[,2:(ncol(eco_geo_results)-2)])))


#####################
#	CORRELATIONS	#
#####################
cols_eco_geo<-2:11
pdf(file="eco_geo_corr_plots.pdf")
#Trying to figure out if any of the measures are especially bad/ different
#Look for correlations among geographic measures and identify when correlation is low, with raw and ranked values
#This function calculates correlations, outputs top correlated metrics
eco_geo_corr<-function(eg_matrix){
	print("Correlation among <<measures>> themselves")
	print(sort(rowSums(cor(eg_matrix[,cols_eco_geo],use="complete.obs")>.80),decreasing=T)[1:4])	#COUNT those corr w/ others above 0.80
	print(sort(rowMeans(cor(eg_matrix[,cols_eco_geo],use="complete.obs")),decreasing=T)[1:4])		#just take MEAN corr w/ all others
	print("Correlation among <<ranks>>- when put in order from measures")
	print(sort(rowSums(cor(apply(eg_matrix[,cols_eco_geo],2,rank))>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(apply(eg_matrix[,cols_eco_geo],2,rank))),decreasing=T)[1:4])
	chart.Correlation(eg_matrix[,cols_eco_geo], histogram = TRUE, method = "pearson")
	rbind(rowSums(cor(eg_matrix[,cols_eco_geo],use="complete.obs")>.80),rowSums(cor(apply(eg_matrix[,cols_eco_geo],2,rank))>.80), rowMeans(cor(eg_matrix[,cols_eco_geo],use="complete.obs")), rowMeans(cor(apply(eg_matrix[,cols_eco_geo],2,rank))))
	#To also add the one on ranking?
}
#USA species
eco_geo_results_usa<-eco_geo_results[-which(is.na(eco_geo_results[,9])),]
#within rare and not rare sets
eco_geo_results_LC<-eco_geo_results[which(eco_geo_results[,14]=="LC"),]
eco_geo_results_TH<-eco_geo_results[which(eco_geo_results[,14]!="LC"),]
eco_geo_results_Q_TH<-eco_geo_results[intersect(which(eco_geo_results[,14]!="LC"),which(substr(eco_geo_results[,1],1,1)=="Q")),]

#This applies the function to subsets of:
#				all, 			US, 				LC, 			Threatened
#and saves as matrix with each row being one of those in the above function
list_all_corr<-lapply(list(eco_geo_results,eco_geo_results_usa,eco_geo_results_LC,eco_geo_results_TH),eco_geo_corr)
dev.off()
#Over all species: geo10, eco10 and eco50 lower correlation But not bad... also epa ecoregions when looking at corr among ranks... because in ranking it will deal with the empty data in a funny way
#US only: : geo10, eco10 and eco50 lower correlation But not bad- basically the same as all species, but corr EPA are higher prob because of the ecoregions
#LC: geo10 especially low correlation , the rest are similar
#TH: eco10, eco50 and to a lesser degree geo 10 and geo500 are low correlation
#The best are typically geo50, geo100 and eco 100, sometimes eco50.  When the US EPA is good, the 50 and 100 are always 'better' than the 10

#BUT should probably find highest pairs and drop one of each pair... iteratively
#geo100 and eco100 have highest.. drop.  Maybe also drop geo10


#####################
#	RANKING			#
#####################

#So, we have 10 different measures, some of which are highly correlated, some not. How to reconcile, choose among them for tanking
all_ranks<-apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)
#Could take the mean or majority decision...
#There are two ways to get agreement across all of them
#One way to actually rank species is to identify those that most frequently are ranked in a given bunch, say in the top 10
species_ranked1<-data.frame(sp_names_wcoll,rowSums(all_ranks<10))
#Another way to do actually rank species is to take the mean across rows in the rank order
species_ranked2<-data.frame(sp_names_wcoll,rowMeans(all_ranks))
#species_ranked<-data.frame(sp_names_wcoll,rank(rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank))))
colnames(species_ranked1)<-c("sp","rank"); colnames(species_ranked2)<-c("sp","rank")
#It really doesn't matter which approach :)
cbind(species_ranked1[order(species_ranked1$rank),],species_ranked2[order(species_ranked2$rank,decreasing=T),])

#But let's examine more closely the individual columns and how they differ
#Could look at those that might be most different and figure out why
#Examine them by eye

#Can look at all of them or just subsets such as just for Threatened (well, non LC) by commenting in/out the following
#these_results<-eco_geo_results
#these_results<-eco_geo_results_Q_TH
#these_results<-eco_geo_results_usa
these_results<-eco_geo_results_TH

all_ranks<-apply(these_results[,2:(ncol(these_results)-2)],2,rank)

#these_names<-sp_names_wcoll
#these_names<-sp_names_Q_TH
#these_names<-sp_names_usa
these_names<-sp_names_TH

species_ranked_geo10<-data.frame("sp"=these_names,"rank-geo10"=all_ranks[,1])
species_ranked_geo50<-data.frame("sp"=these_names,"rank-geo50"=all_ranks[,2])
species_ranked_geo100<-data.frame("sp"=these_names,"rank-geo100"=all_ranks[,3])
species_ranked_geo500<-data.frame("sp"=these_names,"rank-geo500"=all_ranks[,4])
species_ranked_eco10<-data.frame("sp"=these_names,"rank-eco10"=all_ranks[,5])
species_ranked_eco50<-data.frame("sp"=these_names,"rank-eco50"=all_ranks[,6])
species_ranked_eco100<-data.frame("sp"=these_names,"rank-eco100"=all_ranks[,7])
species_ranked_ecous10<-data.frame("sp"=these_names,"rank-ecous10"=all_ranks[,8])
species_ranked_ecous50<-data.frame("sp"=these_names,"rank-ecous50"=all_ranks[,9])
species_ranked_ecous100<-data.frame("sp"=these_names,"rank-ecous100"=all_ranks[,10])

#Examine them by eye
cbind(species_ranked_geo50[order(species_ranked_geo50$rank),],
	species_ranked_geo100[order(species_ranked_geo100$rank),],
	species_ranked_eco50[order(species_ranked_eco50$rank),],
	species_ranked_ecous50[order(species_ranked_ecous50$rank),])
write.csv(cbind(species_ranked_geo10[order(species_ranked_geo10$rank),],
	species_ranked_geo50[order(species_ranked_geo50$rank),],
	species_ranked_geo100[order(species_ranked_geo100$rank),],
	species_ranked_geo500[order(species_ranked_geo500$rank),],
	species_ranked_eco10[order(species_ranked_eco10$rank),],
	species_ranked_eco50[order(species_ranked_eco50$rank),],
	species_ranked_eco100[order(species_ranked_eco100$rank),],
	species_ranked_ecous10[order(species_ranked_ecous10$rank),],
	species_ranked_ecous50[order(species_ranked_ecous50$rank),],
	species_ranked_ecous100[order(species_ranked_ecous100$rank),]),file="compare_genetic_ranks.csv")
	
#TO DO ADD IN EMILY CODE FOR LINES

#How to get the difference in new rank from old rank
#match(species_ranked_geo50[order(species_ranked_geo50$rank),1],species_ranked_geo100[order(species_ranked_geo100$rank),1])-1:41

count_changes<-function(list_ranks,base=1){
	base_order<-list_ranks[[base]][order(list_ranks[[base]]$rank),1]
	examine_changes<-data.frame(base_order)
	for (i in 1:length(list_ranks)){
		ranks_diff<-match(base_order,list_ranks[[i]][order(list_ranks[[i]]$rank),1])
		examine_changes<-cbind(examine_changes,match(base_order,list_ranks[[i]][order(list_ranks[[i]]$rank),1])-1:length(base_order))
		colnames(examine_changes)[i+1]<-names(list_ranks[[i]][2])
	}
	#examine_changes<-examine_changes[,-2]
	examine_changes
}

#Running count_changes without specifying base will compare to the first in the list
examine_changes<-count_changes(list(species_ranked_geo50,species_ranked_geo10,species_ranked_geo100,species_ranked_geo500,species_ranked_eco10,species_ranked_eco50,species_ranked_eco100,species_ranked_ecous10, species_ranked_ecous50,species_ranked_ecous100))
colSums(abs(examine_changes[,-1])>4)

#Running examine_changes iteratively through comparing to base list
for (j in 1:10){
	examine_changes<-count_changes(list(species_ranked_geo10,species_ranked_geo50,species_ranked_geo100,species_ranked_geo500,species_ranked_eco10,species_ranked_eco50,species_ranked_eco100,species_ranked_ecous10, species_ranked_ecous50,species_ranked_ecous100),base=j)
	print(colnames(examine_changes)[j+1])
	print(colSums(abs(examine_changes[,-1])>4))
	#print(mean(colSums(abs(examine_changes[,-1])>5)))
}

#look at some examples of species lists
#low difference
cbind(as.character(species_ranked_geo10[order(species_ranked_geo10[,2]),][,1]),as.character(species_ranked_geo50[order(species_ranked_geo50[,2]),][,1]))
#high difference 
cbind(as.character(species_ranked_geo50[order(species_ranked_geo50[,2]),][,1]),as.character(species_ranked_eco50[order(species_ranked_eco50[,2]),][,1]),as.character(species_ranked_ecous50[order(species_ranked_ecous50[,2]),][,1]))




eco_geo_results_usa<-eco_geo_results[-which(is.na(eco_geo_results[,10])),]
sp_names_usa<-eco_geo_results_usa[,1]
cor(eco_geo_results_usa[,2:(ncol(eco_geo_results_usa)-2)],use="complete.obs")<.75
cor(apply(eco_geo_results_usa[,2:(ncol(eco_geo_results_usa)-2)],2,rank))<.75
sp_names_usa[rowSums(apply(eco_geo_results_usa[,2:(ncol(eco_geo_results_usa)-2)],2,rank)<10)]


#So far I've looked at all species together the next step is look at RL Threatened vs. LC separate ranking
