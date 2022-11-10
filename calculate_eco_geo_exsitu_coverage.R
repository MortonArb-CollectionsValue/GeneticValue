################################################################################

## 5-0_calculate_eco_geo_exsitu_coverage.R
### Author: Emily Beckman ### Date: 05/11/2021

### DESCRIPTION:
  # This script creates circular buffers around in situ points (wild occurrence
	#		records) and ex situ points (wild locations where seeds were collected for
	#		cultivation in botanic gardens or arboreta) to calculate the geographic
	#		(area within buffers) and ecological (number of ecoregions within buffers)
	#		diversity conserved in ex situ living collections. A map showing 50km
	#		buffers around in situ and ex situ points is also created.

### DATA IN:
	## Global country boundaries
	#		UIA World Countries Boundaries, UNIGIS Geospatial Education Resources, via
	#			ArcGIS Hub Shapefile
	#			https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0
	## In situ occurrence points (latitude and longitude in decimal degrees)
	# 	Can use the output from 3-1_refine_occurrence_points.R
	# 		https://github.com/MortonArb-CollectionsValue/OccurrencePoints/tree/master/scripts
	# 		Each file has data for only one species and is named "Genus_species.csv"
	# 		You can read in data for mult. species in one file but need to edit the
	#			code to split after reading in
	## Ex situ wild localities (latitude and longitude in decimal degrees)
	# 	Can use the output from 3-1_refine_occurrence_points.R, which has a
	#			"database" column that has "Ex_situ" to distinguish the ex situ records
	#			from the rest of the in situ records
	## Global Ecoregions
	# 	Terrestrial Ecoregions of the World, via WWF (Olson et al., 2001)
	#			https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
	## U.S. EPA Level III and Level IV Ecoregions
	#		Level III (us_eco_l3) -> https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip
	#		Level IV (us_eco_l4)-> https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l4_state_boundaries.zip

### DATA OUT:
  ## Table of geographic and ecological coverage based on three buffer sizes
	#		(10km, 50km, and 100km) and both ecoregion levels (III and IV)
	## Leaflet map with 50 km buffers around in situ and ex situ points, with
	#		ecoregions in background

################################################################################
# Load libraries
################################################################################

## [code chunk from Shannon M. Still]
my.packages <- c("leaflet","rgdal","raster","rgeos","GeoRange","rnaturalearth",
	"dplyr","sf","ggplot2")
#install.packages (my.packages) # turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Set working directory
################################################################################

# either set manually:
#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/occurrence_points"
	### for gap analysis of nine genera dataset
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/Global Tree Conservation Program/Gap Analyses/Conservation Gap Analysis of Priority Genera/occurrence_points"
	### for gap analysis of US oaks 2020
#main_dir <- "/Users/emily/Desktop/work/ArcGIS Map Layers Quercus_2020"

# or use 0-1_set_workingdirectory.R script from OccurrencePoints repo:
#source("./Documents/GitHub/OccurrencePoints/scripts/0-1_set_workingdirectory.R")

################################################################################
# Load functions
################################################################################

# format text in cell for output table
format.cell <- function(ex_result,in_result,final_result){
	cell <- paste0(round(final_result,2),"%","\n",
								 "(",format(round(ex_result,0),format="d",big.mark=",")," / ",
								 format(round(in_result,0),format="d",big.mark=","),")")
	return(cell)
}

# create buffers around points, using specified projection
create.buffers <- function(df,radius,pt_proj,buff_proj){
	# select coordinate columns
	latlong <- df %>% select(decimalLongitude,decimalLatitude)
	# turn occurrence point data into a SpatialPointsDataFrame
	sp_df <- SpatialPointsDataFrame(latlong, df, proj4string = pt_proj)
	# reproject SpatialPointsDataFrame to specified projection
	proj_df <- spTransform(sp_df,buff_proj)
	# place buffer around each point
	buffers <- buffer(proj_df,width=radius,dissolve=T)
	# return buffer polygons
	return(buffers)
}

# create buffers around in situ and ex situ spatial points, calculate areas,
#		then compare to calculate percent coverage
compare.buff.area <- function(insitu,exsitu,radius,pt_proj,buff_proj){
	# create buffers
	buffer_insitu <- create.buffers(insitu,radius,pt_proj,buff_proj)
	buffer_exsitu <- create.buffers(exsitu,radius,pt_proj,buff_proj)
	# calculate buffer area
	print(paste("Based on ",radius/1000," km radius..."))
	area_exsitu <- buffer_exsitu@polygons[[1]]@area/1000000
	print(paste("Area covered by ex situ buffers:", round(area_exsitu,0),"km²"))
	area_insitu <- buffer_insitu@polygons[[1]]@area/1000000
	print(paste("Area covered by in situ buffers:", round(area_insitu,0),"km²"))
	# calculate difference between in situ and ex situ buffer areas (% coverage)
	area_diff_percent <- (area_exsitu/area_insitu)*100
	print(paste0("Percent geographic coverage: ", round(area_diff_percent,2), "%"))
	txt <- format.cell(area_exsitu,area_insitu,area_diff_percent)
	return(txt)
}

# create data frame with ecoregion data extracted for area covered by buffers
intersect.eco.buff <- function(df,radius,pt_proj,buff_proj,eco){
	# create buffers
	buffers <- create.buffers(df,radius,pt_proj,buff_proj)
	# make sure ecoregions are in same projection as buffers
	eco_proj <- spTransform(eco,buff_proj)
	# intersect buffers with ecoregions
	buff_join_eco <- raster::intersect(buffers,eco_proj)
	return(buff_join_eco)
}

# create data frame with ecoregion data extracted for area covered by buffers,
#		for both in situ and ex situ points, then compare count of ecoregions
compare.eco.count <- function(insitu,exsitu,radius,pt_proj,buff_proj,eco){
	# create data frame of ecoregion-buffer intersection
	eco_insitu <- intersect.eco.buff(insitu,radius,pt_proj,buff_proj,eco)
	eco_exsitu <- intersect.eco.buff(exsitu,radius,pt_proj,buff_proj,eco)
	# count number of ecoregions under buffers
	print(paste("Based on ",radius/1000," km radius..."))
	count_exsitu <- nrow(eco_exsitu@data %>% distinct(ECO_ID))
	print(paste0("Number of ecoregions under ex situ buffers: ",count_exsitu))
	count_insitu <- nrow(eco_insitu@data %>% distinct(ECO_ID))
	print(paste0("Number of ecoregions under in situ buffers: ",count_insitu))
	# calculate difference in number of ecoregions
	eco_diff_percent <- (count_exsitu/count_insitu)*100
	print(paste0("Percent ecological coverage: ", round(eco_diff_percent,2), "%"))
	txt <- format.cell(count_exsitu,count_insitu,eco_diff_percent)
	return(txt)
}

# create data frame with Level III ecoregion data extracted for area covered by
#		buffers, for both in situ and ex situ points, then compare ecoregion counts
#compare.ecol3.count <- function(insitu,exsitu,radius,pt_proj,buff_proj,eco){
#	# create data frame of ecoregion-buffer intersection
#	eco_insitu <- intersect.eco.buff(insitu,radius,pt_proj,buff_proj,eco)
#	eco_exsitu <- intersect.eco.buff(exsitu,radius,pt_proj,buff_proj,eco)
#	# count number of ecoregions under buffers
#	print(paste("Based on ",radius/1000," km radius..."))
#	count_exsitu <- nrow(eco_exsitu@data %>% distinct(US_L3CODE))
#	print(paste0("Number of ecoregions under ex situ buffers: ",count_exsitu))
#	count_insitu <- nrow(eco_insitu@data %>% distinct(US_L3CODE))
#	print(paste0("Number of ecoregions under in situ buffers: ",count_insitu))
#	# calculate difference in number of ecoregions
#	eco_diff_percent <- (count_exsitu/count_insitu)*100
#	print(paste0("Percent ecological coverage: ", round(eco_diff_percent,2), "%"))
#	txt <- format.cell(count_exsitu,count_insitu,eco_diff_percent)
#	return(txt)
#}

# create data frame with Level IV ecoregion data extracted for area covered by
#		buffers, for both in situ and ex situ points, then compare ecoregion counts
compare.ecol4.count <- function(insitu,exsitu,radius,pt_proj,buff_proj,eco){
	# create data frame of ecoregion-buffer intersection
	eco_insitu <- intersect.eco.buff(insitu,radius,pt_proj,buff_proj,eco)
	eco_exsitu <- intersect.eco.buff(exsitu,radius,pt_proj,buff_proj,eco)
	# count number of ecoregions under buffers
	print(paste("Based on ",radius/1000," km radius..."))
	count_exsitu <- nrow(eco_exsitu@data %>% distinct(US_L4CODE))
	print(paste0("Number of US ecoregions under ex situ buffers: ",count_exsitu))
	count_insitu <- nrow(eco_insitu@data %>% distinct(US_L4CODE))
	print(paste0("Number of US ecoregions under in situ buffers: ",count_insitu))
	# calculate difference in number of ecoregions
	eco_diff_percent <- (count_exsitu/count_insitu)*100
	print(paste0("Percent US ecological coverage: ", round(eco_diff_percent,2), "%"))
	txt <- format.cell(count_exsitu,count_insitu,eco_diff_percent)
	return(txt)
}

# clip points by boundary so only in target area
# (helpful if focusing on one country/region)
clip.by.boundary <- function(pts,pt_proj,boundary){
	# select coordinate columns
	latlong <- pts %>% select(decimalLongitude,decimalLatitude)
	# turn occurrence point data into a SpatialPointsDataFrame
	spatial_pts <- SpatialPointsDataFrame(latlong, pts, proj4string = pt_proj)
	# clip by boundary created earlier
	spatial_pts <- spatial_pts[boundary, ]
	# keep just the data (not the spatial info you added)
	pts_new <- spatial_pts@data
	return(pts_new)
}

################################################################################
################################################################################
# Set up workspace
################################################################################

# define projection of points
#		WGS84
wgs.proj <- sp::CRS(SRS_string="EPSG:4326")
# define projection for calculations (meters must be the unit)
#		Equal Earth Projection
#		Warning shouldn't cause any errors; metadata says they're "overly cautious"
aea.proj <- sp::CRS(SRS_string="EPSG:8857")

################################################################################
# Choose buffer sizes to use
################################################################################

# buffer size in kilometers = value/1000
largest_buff <- 500000
large_buff <- 100000
med_buff <- 50000
small_buff <- 20000

################################################################################
# Choose target species
################################################################################

# read in target taxa list
taxon_list <- read.csv(file.path(main_dir, "inputs", "taxa_list",
  "target_species_with_syn.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
head(taxon_list)

# OPTIONAL, depending on your workflow:
#		read in manual edits to target species maps
pt_edits <- read.csv(file.path(main_dir, "inputs", "known_distribution",
  "manual_point_edits.csv"),
  header = T, na.strings = c("","NA"),colClasses = "character")
head(pt_edits)

# OPTIONAL: select threatened (CR, EN, VU) and Near Threatened species
#target_sp <- taxon_list %>% filter(rl_category == "CR" | rl_category == "EN" |
#																	 rl_category == "VU" | rl_category == "NT")
#target_sp <- unique(target_sp$species_name_acc)
# randomly select Least Concern species
# 	out of... 8 Malus ; 83 Quercus ; 16 Tilia ; 15 Ulmus
#lc <- taxon_list %>% filter(is.na(map_flag) & list == "desiderata")
#lc_sp <- lc[sample(nrow(lc), 29), ]
#lc_sp <- lc_sp$species_name_acc
	# add Matt's selection
#lc_sp <- c(lc_sp,"Ulmus bergmaniana")
#sort(lc_sp)
# OPTIONAL: select target species
target_sp <- taxon_list %>% filter(grepl("^MAP",map_flag))
target_sp <- sort(gsub(" ","_",target_sp$species_name_acc))
length(target_sp)
#target_sp <- target_sp[c(1:2,4:11,13:15,17:19,21:53)]
	### for gap analysis of US oaks 2020
	#target_sp <- c("Qacerifolia","Qajoensis","Qarkansana","Qaustrina","Qboyntonii",
	#"Qcarmenensis","Qcedrosensis","Qchapmanii","Qcorneliusmulleri","Qdumosa",
	#"Qengelmannii","Qgeorgiana","Qgraciliformis","Qhavardii","Qhinckleyi",
	#"Qinopina","Qlaceyi","Qlobata","Qoglethorpensis","Qpacifica","Qpalmeri",
	#"Qparvula","Qpumila","Qrobusta","Qsadleriana","Qsimilis","Qtomentella",
	#"Qtoumeyi")
# get native dist information to see if using RL (default) or GTS
native_dist <- read.csv(file.path(main_dir,"inputs","known_distribution",
  "target_taxa_with_native_dist.csv"), header = T, na.strings = c("","NA"),
  colClasses = "character")
native_dist <- native_dist %>% dplyr::select(species_name_acc,rl_native_dist,
	gts_native_dist)

################################################################################
# Read in polygon data
################################################################################

# trialing protected areas layer from Khoury et al. 2019 (raster version of WDPA)
#		https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/WTLNRG
#pa_raster <- raster(file.path(main_dir,"inputs","gis_data","dataverse_files",
#	"wdpa_reclass.tif"))

# haven't gotten the PA layer working yet - it is very large and I think needs
#		to be flatted/dissolved (also called "union") so its a single feature instead
#		of multiple features; this should help with the size also
## read in shapefile of global protected areas
	## CITATION: UNEP-WCMC and IUCN (2022), Protected Planet: The World Database on Protected Areas (WDPA) [Online], February 2022, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.
# As of January 2022, the WDPA includes data on 269,457 protected areas of
#		which 96% were polygons and 4% points.

# since the shapefile is large, the source has it split into three;
#		we will load each then merge into one

	# read in shapefile
pa0 <- sf::st_read(file.path(main_dir,"inputs","gis_data","WDPA_Jan2022_Public_shp",
	"WDPA_Jan2022_Public_shp_0","WDPA_Jan2022_Public_shp-polygons.shp"))

pa0 <- readOGR(file.path(main_dir,"inputs","gis_data","WDPA_Jan2022_Public_shp",
	"WDPA_Jan2022_Public_shp_0","WDPA_Jan2022_Public_shp-polygons.shp"))

	# look at its properties
st_geometry_type(pa0)
unique(st_geometry_type(pa0))
st_crs(pa0)
st_bbox(pa0)
pa0 #Simple feature collection with 85,783 features and 30 fields
unique(pa0$ISO3)
	# select terrestrial PAs that are designated/established (not pending),
	#		and in the U.S.
pa0_terr_usa <- pa0 %>%
	filter(MARINE != 2) %>%
	filter(STATUS == "Inscribed" | STATUS == "Designated" | STATUS == "Established") %>% #metadata here: https://pp-import-production.s3.amazonaws.com/WDPA_Manual_1_5.pdf
	filter(grepl("USA",ISO3)) %>%
	select(WDPAID,NAME,DESIG_ENG,DESIG_TYPE,IUCN_CAT,REP_AREA,STATUS,OWN_TYPE,MANG_AUTH,ISO3)
pa0_terr_usa #Simple feature collection with 1,919 features and 10 fields
	# see if any polygons are invalid
pa0_terr_check <- st_is_valid(pa0_terr_usa, reason = TRUE)
which(pa0_terr_check != "Valid Geometry") # rows that are invalid
pa0_terr_usa[142,] #Yukon Delta National Wildlife Refuge Marine Protected Area
	# can plot if not too big/complicated of a PA; otherwise takes too long
#ggplot() +
#  geom_sf(data = pa0_terr_usa[142,], size = 1, color = "black", fill = "cyan1") +
#  coord_sf()
	# we will just remove that invalid polygon for now
pa0_terr_usa <- pa0_terr_usa[-142,]
	# I've tried 'make valid' but doesn't seem to resolve the union issues
#iv1 <- st_make_valid(pa0_terr[16306,])
	# extract the geometry of the conservation areas to a new variable
	# to get rid of all other attributes
pa0_outline <- st_geometry(pa0_terr_usa) %>% st_as_sf()
pa0_outline #Simple feature collection with 1,918 features and 0 fields
	# see if we can do a union with just the U.S. PAs
pa0_terr_usa_union <- st_union(pa0_outline) # hung; ran all night and finished

st_write(pa0_terr_usa_union, file.path(main_dir,"inputs","gis_data",
	"WDPA_Jan2022_Public_shp","WDPA_Jan2022_Public_shp_0",
	"WDPA_Jan2022_Public_shp-polygons-USA_terrestrial_union.shp"))


#pa1 <- sf::st_read(file.path(main_dir,"inputs","gis_data","WDPA_Oct2021_Public_shp",
#	"WDPA_Oct2021_Public_shp_1","WDPA_Oct2021_Public_shp-polygons.shp"))
#pa1 <- st_cast(pa1$geometry[[4]], "POLYGON")
#pa <- rbind(pa0,pa1)
#plot(pa)
#pa2 <- sf::st_read(file.path(main_dir,"inputs","gis_data","WDPA_Oct2021_Public_shp",
#	"WDPA_Oct2021_Public_shp_2","WDPA_Oct2021_Public_shp-polygons.shp"))
#pa_list <- list(pa0,pa1,pa2)
#pa <- do.call(rbind,pa_list)
#rm(pa0);rm(pa1);rm(pa2)
	# reproject to make sure same crs as buffer layer for later intersect
#st_crs(pa) # look at final ID number
#pa_proj <- st_transform(pa,crs = 4326)
#rm(pa)
#	# dissolve polygons into one layer (multi to single)
#	## !! throws an error due to invalid polygons.. also takes a long time
#pa_diss <- st_union(pa_proj)
#rm(pa_proj)


	### for gap analysis of nine genera dataset
	#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/occurrence_points"

# read in shapefile of global ecoregions
ecoregions <- readOGR(file.path(main_dir,"inputs","gis_data",
	"official","wwf_terr_ecos.shp"))

# read in shapefile of U.S. EPA Level III ecoregions
#ecol3 <- readOGR(file.path(main_dir,"inputs","gis_data",
#	"us_eco_l3/us_eco_l3.shp"))

# read in shapefile of U.S. EPA Level IV ecoregions
ecol4 <- readOGR(file.path(main_dir,"inputs","gis_data",
	"us_eco_l4_state_boundaries/us_eco_l4.shp"))

# read in shapefile of country boundaries
world_countries <- readOGR(file.path(main_dir,"inputs","gis_data",
	"UIA_World_Countries_Boundaries-shp/World_Countries__Generalized_.shp"))
target_iso <- c("US")
target_countries <- world_countries[world_countries@data$ISO %in% target_iso,]
	## create polygon for clipping points later, one in each projection
target_countries.wgs <- spTransform(target_countries,wgs.proj)
boundary.wgs <- aggregate(target_countries.wgs,dissolve = TRUE)

	### for gap analysis of nine genera dataset
	#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/Global Tree Conservation Program/Gap Analyses/Conservation Gap Analysis of Priority Genera/occurrence_points"
	### for gap analysis of US oaks 2020
	#main_dir <- "/Volumes/GoogleDrive-103729429307302508433/My Drive/ArcGIS Map Layers Quercus_2020"

################################################################################
## Calculate geographic and ecological coverage of ex situ collections
################################################################################

### START SUMMARY TABLE

# we add each target species as we go along
summary_tbl <- data.frame(
	species = "start",
	geo_sm = "start", geo_md = "start",	geo_lg = "start", geo_vlg = "start",
	eco_sm = "start", eco_md = "start", eco_lg = "start",
	eco_usl4_sm = "start", eco_usl4_md = "start", eco_usl4_lg = "start",
	EOO = "start",
	dist_filter = "start",
	#pa_coverage = "start",
	stringsAsFactors=F)

### CYCLE THROUGH TARGET SPECIES TO CALCULATE EX SITU COVERAGE

for(sp in 1:length(target_sp)){

## can test with one species first
#sp <- 6

	# print progress
  cat("\tStarting ", target_sp[sp], "\n")

	### READ IN AND PREP POINT DATA

	## read in occurrence points (includes ex situ)
	insitu_raw <- read.csv(file.path(main_dir,"outputs","spp_edited_points",
		paste0(target_sp[sp],".csv")), na.strings=c("","NA"), stringsAsFactors = F)
		nrow(insitu_raw)
	 spp.rl.dist <- native_dist[which(native_dist$species_name_acc == gsub("_"," ",target_sp[sp])),]
	## filter as desired
	 insitu <- insitu_raw %>%
	   filter(database == "Ex_situ" |
	     (.cen & .inst & .con & .outl &
					### for gap analysis of nine genera dataset
				  #.bonapnative & .yr1950 & .yrna &
			 #.urb & .yr1950 & .yr1980 & .yrna &
	     #(.gtsnative | is.na(.gtsnative)) &
	     #(.rlnative  | is.na(.rlnative)) &
	     #(.rlintroduced | is.na(.rlintroduced)) &
	     basisOfRecord != "FOSSIL_SPECIMEN" & basisOfRecord != "LIVING_SPECIMEN" &
					### for gap analysis of nine genera dataset
					#basisOfRecord != "H?" &
	     establishmentMeans != "INTRODUCED" & establishmentMeans != "MANAGED" &
	     establishmentMeans != "INVASIVE"))
	  if(!is.na(spp.rl.dist$rl_native_dist)){
	     insitu <- insitu %>%
	     		filter(.rlnative | is.na(.rlnative))
			 dist_filter_val <- "RL"
	  } else if(!is.na(spp.rl.dist$gts_native_dist)){
    	 insitu <- insitu %>%
	     		filter(.gtsnative | is.na(.gtsnative))
			 dist_filter_val <- "GTS"
	  } else {
			 dist_filter_val <- "N/A"
		}
		nrow(insitu)

	## check document with manual point edits to see if anything needs to be added back or removed
	manual.edit <- pt_edits[which(pt_edits$species_name_acc == gsub("_"," ",target_sp[sp])),]
		# bounding box
	if(!is.na(manual.edit$bounding_box)){
		bounds <- unlist(strsplit(manual.edit$bounding_box,"; "))
		for(i in 1:length(bounds)){
			within <- unlist(strsplit(bounds[i],", "))
			insitu <- insitu %>% filter(!(decimalLongitude > as.numeric(within[1]) &
																	  decimalLongitude < as.numeric(within[3]) &
																	 	decimalLatitude > as.numeric(within[2]) &
																		decimalLatitude < as.numeric(within[4])))
		}
	}; nrow(insitu)
		# remove
	if(!is.na(manual.edit$remove)){
		remove <- unlist(strsplit(manual.edit$remove,"; "))
		insitu <- insitu %>% filter(!(UID %in% remove))
	}; nrow(insitu)
			### for gap analysis of nine genera dataset; select US points only
			#insitu <- clip.by.boundary(insitu,wgs.proj,boundary.wgs)
		# add back
	if(!is.na(manual.edit$keep)){
		keep <- unlist(strsplit(manual.edit$keep,"; "))
		add <- insitu_raw %>% filter(UID %in% keep)
		insitu <- suppressMessages(full_join(insitu,add))
	}; nrow(insitu)

		### for gap analysis of US oaks 2020
		#exsitu_raw <- read.csv(file.path(main_dir,target_sp[sp],
		#	paste0(target_sp[sp],"_exsitu.csv")), na.strings=c("","NA"), stringsAsFactors = F)
		#	nrow(exsitu_raw)
		#exsitu <- exsitu_raw %>%
		#	filter(!is.na(latitude) & !is.na(longitude)) %>%
		#	rename(decimalLatitude = latitude,
		#			 	 decimalLongitude = longitude)
		#	nrow(exsitu)
		#insitu_raw <- read.csv(file.path(main_dir,target_sp[sp],
		#	paste0(target_sp[sp],"_insitu.csv")), na.strings=c("","NA"), stringsAsFactors = F)
		#	nrow(insitu_raw)
		#insitu <- insitu_raw %>%
		#	filter(!is.na(latitude) & !is.na(longitude)) %>%
		#	rename(decimalLatitude = latitude,
		#			 	 decimalLongitude = longitude)
		#insitu <- full_join(insitu,exsitu)
		#nrow(insitu)

	### CALCULATE EOO (convex hull)

		# uses package GeoRange to calculate area of convex hull in km2
	hull_area <- CHullAreaEarth(insitu$decimalLongitude,insitu$decimalLatitude)
	hull_area <- round(hull_area,0)
	print(paste("EOO:",hull_area,"km²"))

	### CALCULATE PROTECTED AREA COVERAGE

		# create 10km buffers around in situ points
	#buff_10_insitu <- create.buffers(insitu,10000,wgs.proj,wgs.proj)
		# make sure buffer layer is in same projection as PA layer
	#buff_10_insitu_sf <- st_transform(st_as_sf(buff_10_insitu),crs = 4326)
	#buff_10_insitu_sf
		# intersect the buffer and PA layers to get overlap
	#buff_pa_inter <- st_intersection(buff_10_insitu_sf, pa0_terr_usa_union)
	#	# compare the buffer area to the buff_pa_inter area to get % coverage
	#.....

	## create df with just ex situ points
	exsitu <- insitu %>% filter(database == "Ex_situ")
	# get number of individuals ex situ with lat-long data
	print(paste("Number of ex situ individuals:",sum(as.numeric(exsitu$establishmentMeans))))

	### TURN THIS ON TO CALCULATE FOR MORTON ACCESSIONS ONLY
	#exsitu <- exsitu %>% filter(datasetName == "MortonArb")

		# check there are ex situ points, if not skip to end (can't do calculations)
#	if(nrow(exsitu) == 0){
#	  # add text results to summary table
#		summary_add <- data.frame(
#			species = gsub("_"," ",target_sp[sp]),
#			geo_sm = NA, geo_md = NA,	geo_lg = NA, geo_vlg = NA,
#			eco_sm = NA, eco_md = NA, eco_lg = NA,
#			eco_usl4_sm = NA, eco_usl4_md = NA, eco_usl4_lg = NA,
#			EOO = hull_area,
#				### turn the next row off for gap analysis of US oaks 2020
#			#dist_filter = dist_filter_val,
#			stringsAsFactors=F)
#		print("No ex situ points; skipping buffer calculations")
#	} else {
#
#		### CALCULATE EX SITU COVERAGE
#
#			## Geographic coverage
#
#		# calculate area based on largest buffers
#		geo_coverage_vlg <- compare.buff.area(insitu,exsitu,largest_buff,wgs.proj,aea.proj)
#		# calculate area based on large buffers
#		geo_coverage_lg <- compare.buff.area(insitu,exsitu,large_buff,wgs.proj,aea.proj)
#		# calculate area based on medium buffers
#		geo_coverage_md <- compare.buff.area(insitu,exsitu,med_buff,wgs.proj,aea.proj)
#		# calculate area based on small buffers
#		geo_coverage_sm <- compare.buff.area(insitu,exsitu,small_buff,wgs.proj,aea.proj)
#
#			## Ecological coverage
#
#		## Global ecoregions
#		# count ecoregions under large buffers
#		eco_coverage_lg <- compare.eco.count(insitu,exsitu,large_buff,wgs.proj,aea.proj,ecoregions)
#		# count ecoregions under medium buffers
#		eco_coverage_md <- compare.eco.count(insitu,exsitu,med_buff,wgs.proj,aea.proj,ecoregions)
#		# count ecoregions under small buffers
#		eco_coverage_sm <- compare.eco.count(insitu,exsitu,small_buff,wgs.proj,aea.proj,ecoregions)
#
#		## U.S. Level 4 (most specific) ecoregions
#		# get just points that are in the U.S.
#		us_insitu <- clip.by.boundary(insitu,wgs.proj,boundary.wgs)
#		us_exsitu <- clip.by.boundary(exsitu,wgs.proj,boundary.wgs)
#		# if there are in situ and ex situ points in the U.S., then calculate coverage
#		if(nrow(us_exsitu) > 0 & nrow(us_insitu) > 0){
#			# count ecoregions under large buffers
#			ecol4_coverage_lg <- compare.ecol4.count(us_insitu,us_exsitu,large_buff,wgs.proj,aea.proj,ecol4)
#			# count ecoregions under medium buffers
#			ecol4_coverage_md <- compare.ecol4.count(us_insitu,us_exsitu,med_buff,wgs.proj,aea.proj,ecol4)
#			# count ecoregions under small buffers
#			ecol4_coverage_sm <- compare.ecol4.count(us_insitu,us_exsitu,small_buff,wgs.proj,aea.proj,ecol4)
#		# if there's distribution in the U.S. but no ex situ points, assign 0%
#		} else if(nrow(us_exsitu) == 0 & nrow(us_insitu) > 0){
#			ecol4_coverage_lg <- "0%"
#			ecol4_coverage_md <- "0%"
#			ecol4_coverage_sm <- "0%"
#		# if not in U.S. then NA
#		} else {
#			ecol4_coverage_lg <- NA
#			ecol4_coverage_md <- NA
#			ecol4_coverage_sm <- NA
#		}
#
#		### SUMMARY TABLE
#
#	  ## Add text results to summary table
#	  summary_add <- data.frame(
#			species = gsub("_"," ",target_sp[sp]),
#			geo_sm = geo_coverage_sm,
#			geo_md = geo_coverage_md,
#			geo_lg = geo_coverage_lg,
#			geo_vlg = geo_coverage_vlg,
#			eco_sm = eco_coverage_sm,
#			eco_md = eco_coverage_md,
#			eco_lg = eco_coverage_lg,
#			eco_usl4_sm = ecol4_coverage_sm,
#			eco_usl4_md = ecol4_coverage_md,
#			eco_usl4_lg = ecol4_coverage_lg,
#			EOO = hull_area,
#				### turn the next row off for gap analysis of US oaks 2020
#			#dist_filter = dist_filter_val,
#			stringsAsFactors=F)
#	}
#	summary_tbl[sp,] <- summary_add
}

## write summary table
summary_tbl
write.csv(summary_tbl, file.path(main_dir,"outputs","exsitu_coverage",
	"ExSituCoverage_BufferTable_1_28_22.csv"), row.names = F)

	### for gap analysis of US oaks 2020
#write.csv(summary_tbl, file.path(main_dir,
#	"Oaks2020_BufferTable.csv"), row.names = F)


################################################################################
# Map

#		I've left this section separate, for experimenting with leaflet and
# 	visualizing data for individual species (you have to manually select your
#		target species); it could easily be part of the loop above if desired in
#		the future, if we want all these maps saved.

################################################################################

# if you'd like to have state borders:
# read in state polygons using rnaturalearth package
state_bound <- ne_states(country="united states of america")

# transform ecoregion polygons to WGS84, for mapping
ecol4_wgs <- spTransform(ecol4,wgs.proj)

## CHOOSE TARGET SPECIES TO MAP
target_sp # view target species
sp <- 4 # select one to work with
cat("\tMapping ", target_sp[sp])

## the next section is simply copied from above; the flow could eventually be
#		improved so its not just repeated, but the map is not a necessary part of
#		the above calculations and is sometimes easier to view individually as
#		desired, so I'm just repeating for now
## --------------------------------------------------------------------
### READ IN AND PREP POINT DATA
	## read in occurrence points (includes ex situ)
	insitu_raw <- read.csv(file.path(main_dir,"outputs","spp_edited_points",
		paste0(target_sp[sp],".csv")), na.strings=c("","NA"), stringsAsFactors = F)
		nrow(insitu_raw)
	 spp.rl.dist <- native_dist[which(native_dist$species_name_acc == gsub("_"," ",target_sp[sp])),]
	## filter as desired
	 insitu <- insitu_raw %>%
	   filter(database == "Ex_situ" |
	     (.cen & .inst & .con & .outl &
					### for gap analysis of nine genera dataset
				  #.bonapnative & .yr1950 & .yrna &
			 #.urb & .yr1950 & .yr1980 & .yrna &
	     #(.gtsnative | is.na(.gtsnative)) &
	     #(.rlnative  | is.na(.rlnative)) &
	     #(.rlintroduced | is.na(.rlintroduced)) &
	     basisOfRecord != "FOSSIL_SPECIMEN" & basisOfRecord != "LIVING_SPECIMEN" &
					### for gap analysis of nine genera dataset
					#basisOfRecord != "H?" &
	     establishmentMeans != "INTRODUCED" & establishmentMeans != "MANAGED" &
	     establishmentMeans != "INVASIVE"))
	  if(!is.na(spp.rl.dist$rl_native_dist)){
	     insitu <- insitu %>%
	     		filter(.rlnative | is.na(.rlnative))
			 dist_filter_val <- "RL"
	  } else if(!is.na(spp.rl.dist$gts_native_dist)){
    	 insitu <- insitu %>%
	     		filter(.gtsnative | is.na(.gtsnative))
			 dist_filter_val <- "GTS"
	  } else {
			 dist_filter_val <- "N/A"
		}
		nrow(insitu)
	## check document with manual point edits to see if anything needs to be added back or removed
	manual.edit <- pt_edits[which(pt_edits$species_name_acc == gsub("_"," ",target_sp[sp])),]
		# bounding box
	if(!is.na(manual.edit$bounding_box)){
		bounds <- unlist(strsplit(manual.edit$bounding_box,"; "))
		for(i in 1:length(bounds)){
			within <- unlist(strsplit(bounds[i],", "))
			insitu <- insitu %>%
				filter(!(decimalLongitude > as.numeric(within[1]) &
				decimalLongitude < as.numeric(within[3]) &
				decimalLatitude > as.numeric(within[2]) &
				decimalLatitude < as.numeric(within[4]))) }
	}; nrow(insitu)
		# remove
	if(!is.na(manual.edit$remove)){
		remove <- unlist(strsplit(manual.edit$remove,"; "))
		insitu <- insitu %>% filter(!(UID %in% remove))
	}; nrow(insitu)
			### for gap analysis of nine genera dataset; select US points only
			#insitu <- clip.by.boundary(insitu,wgs.proj,boundary.wgs)
		# add back
	if(!is.na(manual.edit$keep)){
		keep <- unlist(strsplit(manual.edit$keep,"; "))
		add <- insitu_raw %>% filter(UID %in% keep)
		insitu <- suppressMessages(full_join(insitu,add))
	}; nrow(insitu)
	## create df with just ex situ points
	exsitu <- insitu %>% filter(database == "Ex_situ")
## --------------------------------------------------------------------

### CREATE MAP

# select only ecoregions that are within the buffers; otherwise map takes a
#		long time to load in browser and colors are not distinct enough
inter <- intersect.eco.buff(insitu,50000,wgs.proj,wgs.proj,ecol4)
codes <- unique(inter@data$US_L4CODE)
ecol4_sel <- ecol4_wgs[ecol4_wgs@data$US_L4CODE %in% codes,]

# create ecoregion color palette; can run display.brewer.all() to see options
ecol4_pal <- colorFactor(palette = "Set2",domain = ecol4_sel@data$US_L4NAME,
	reverse = F, na.color = "white")

# create map to visualze species distribution and ecoregions
# 	note that you can comment out any section to remove that element
coverage_map <- leaflet() %>%
	## background
	addProviderTiles("Esri.WorldGrayCanvas",#"Esri.WorldShadedRelief",#"Esri.WorldTerrain",
		options = providerTileOptions(maxZoom = 10)) %>%
	## EPA Level IV ecoregions
	addPolygons(
		data = ecol4_sel,
		fillColor = ~ecol4_pal(ecol4_sel@data$US_L4NAME), fillOpacity = 0.6,
		color = ~ecol4_pal(ecol4_sel@data$US_L4NAME), weight = 0.5, opacity = 1) %>%
	## state boundaries
	addPolygons(data = state_bound,
		fillOpacity = 0, color = "#969696", weight = 1.5, opacity = 1) %>%
	## in situ buffers
	addPolygons(data = create.buffers(insitu,50000,wgs.proj,wgs.proj),
	## ex situ buffers
	addPolygons(data = create.buffers(exsitu,50000,wgs.proj,wgs.proj),
		smoothFactor = 0.5,	weight = 1, opacity = 1, color = "#b0f6f7",
		fillOpacity = 0.3) %>%
	## in situ points
	addCircleMarkers(data = insitu,
		lng = ~decimalLongitude, lat = ~decimalLatitude,
		#popup = ~paste("Source(s):", all_source_databases, UID),
		radius = 4, fillOpacity = 1, stroke = F, color = "#1c1c1b") %>%
	## ex situ points
	addCircleMarkers(data = exsitu,
		lng = ~decimalLongitude, lat = ~decimalLatitude,
		#popup = ~paste("Garden:", datasetName, UID),
		radius = 4, fillOpacity = 1, stroke = F, color = "#b0f6f7") %>%
			# could add another layer that highlights just the Morton points or
			#		could use a color palette based on the "database" column :)
	## title
	addControl(paste0(target_sp[sp],
			" in situ distribution and representation in ex situ collections"),
		position = "topright") %>%
	## legend
	addLegend(labels =
		c("Geographic range (occurrence points and 50 km buffers)",
			"Populations sampled for ex situ (wild collection locations and 50 km buffers)",
			"Background colors show U.S. EPA Level IV Ecoregions"),
		colors = c("#1c1c1b","#b0f6f7","#d6ccb4"), title = "Legend",
		position = "bottomright", opacity = 0.8) %>%
	## scale bar
	addScaleBar(position = "bottomleft",
		options = scaleBarOptions(maxWidth = 150)) %>%
	## pick coords for center of frame and zoom level when map first appears
	setView(-97, 40, zoom = 5)

# view map
coverage_map

# save map
# if you have too many points or layers, the map will be too big and can't save
htmlwidgets::saveWidget(coverage_map,
	file = file.path(main_dir,"outputs","exsitu_coverage",
		paste0(target_sp[sp],"-exsitu_coverage_interactive_map.html")))
