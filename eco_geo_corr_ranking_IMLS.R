
# install.packages("PerformanceAnalytics") # Turn on to install current version
library(PerformanceAnalytics)

#Set up paths
data_dir <- "G:/Shared drives/IMLS MFA/occurrence_points/outputs"
output_dir <- "G:/Shared drives/IMLS MFA/Genetic diversity value"

#Read in file
eco_geo_results<-read.csv(file.path(data_dir,"ExSituCoverage_BufferTable_6_30_21.csv"))
#Pull species names
sp_names_all<-eco_geo_results[,1]
#Identify species with no collections
no_collections<-sp_names_all[which(is.na(eco_geo_results[,2]))]
sp_names_wcoll<-sp_names_all[-which(is.na(eco_geo_results[,2]))]
#Remove species with no collections
eco_geo_results<-eco_geo_results[-which(is.na(eco_geo_results[,2])),]
#Maybe remove these rows? These are the non USA ones
#eco_geo_results<-eco_geo_results[-which(is.na(eco_geo_results[,9])),]

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
eco_geo_results_LC<-eco_geo_results[eco_geo_results[,ncol(eco_geo_results)]=="LC",]
eco_geo_results_TH<-eco_geo_results[eco_geo_results[,ncol(eco_geo_results)]!="LC",]

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



#####################
#	RANKING			#
#####################

#So, we have 10 different measures, some of which are highly correlated, some not. How to reconcile, choose among them for tanking
#Could take the mean or majority decision...
#One way to actually rank species is to identify those that most frequently are ranked in a given bunch, say in the top 10
species_ranked1<-data.frame(sp_names_wcoll,rowSums(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)<10))
#Another way to do actually rank species is to take the mean across rows in the rank order
species_ranked2<-data.frame(sp_names_wcoll,rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)))
#species_ranked<-data.frame(sp_names_wcoll,rank(rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank))))
colnames(species_ranked1)<-c("sp","rank"); colnames(species_ranked2)<-c("sp","rank")
#It really doesn't matter which approach :)
cbind(species_ranked1[order(species_ranked1$rank),],species_ranked2[order(species_ranked2$rank,decreasing=T),])

#Could look at those that might be most different and figure out why
#Geo 50
species_ranked_geo50<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,2])
#Geo 100
species_ranked_geo100<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,3])
#Eco 50
species_ranked_eco50<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,6])
#Eco USA 50
species_ranked_ecous50<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,9])
colnames(species_ranked_geo50)<-c("sp","rank-geo50"); colnames(species_ranked_geo100)<-c("sp","rank-geo100")
colnames(species_ranked_eco50)<-c("sp","rank-eco50"); colnames(species_ranked_ecous50)<-c("sp","rank-ecous50")
cbind(species_ranked_geo50[order(species_ranked_geo50$rank),],species_ranked_geo100[order(species_ranked_geo100$rank),],
	species_ranked_eco50[order(species_ranked_eco50$rank),],species_ranked_ecous50[order(species_ranked_ecous50$rank),])
#Maybe 10 is just too fine scale!!
 match(species_ranked_geo50[order(species_ranked_geo50$rank),1],species_ranked_geo100[order(species_ranked_geo100$rank),1])-1:41
examine_changes<-data.frame(species=as.character(species_ranked_geo50[order(species_ranked_geo50$rank),1]),
	change_eco50=(match(species_ranked_geo50[order(species_ranked_geo50$rank),1],species_ranked_eco50[order(species_ranked_eco50$rank),1])-1:41),
	chage_ecoUS50=(match(species_ranked_geo50[order(species_ranked_geo50$rank),1],species_ranked_ecous50[order(species_ranked_ecous50$rank),1])-1:41))



eco_geo_results_usa<-eco_geo_results[-which(is.na(eco_geo_results[,10])),]
sp_names_usa<-eco_geo_results_usa[,1]
cor(eco_geo_results_usa[,2:(ncol(eco_geo_results_usa)-2)],use="complete.obs")<.75
cor(apply(eco_geo_results_usa[,2:(ncol(eco_geo_results_usa)-2)],2,rank))<.75
sp_names_usa[rowSums(apply(eco_geo_results_usa[,2:(ncol(eco_geo_results_usa)-2)],2,rank)<10)]


#So far I've looked at all species together the next step is look at RL Threatened vs. LC separate ranking
