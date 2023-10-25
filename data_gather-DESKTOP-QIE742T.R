source(here::here('packages.R'))
source(here::here('functions.R'))

#Data from the database#####
# con <- dbConnect(PostgreSQL(), dbname = "fishecu_complex-survey_db", user = "fishecuuser", host = "172.21.3.20", password = "f1sh3cuus3r!")

#Connecting to the database to extract info about Lipno
#Latin names
# specs <- data.table(dbGetQuery(conn = con, statement =  "SELECT * FROM fishecu.species;"))#to see the whole table (reservoir). Important to see the reservoir ID number.
# specs <- specs[,c(1,2,10)]
# write.xlsx(specs, here::here('Data', 'specs.xlsx'))
specs <- setDT(readxl::read_xlsx(here::here('Data', 'specs.xlsx')))

# Gillnet deployments####
# all_gill_depl <- data.table(dbGetQuery(conn = con, statement = "SELECT *
# FROM fishecu.gillnet_sampling_merged_wide
# WHERE reservoirid IN (2);"))
# write.xlsx(all_gill_depl, here::here('Data', 'all_gill_depl.xlsx'))
all_gill_depl <- setDT(readxl::read_xlsx(here::here('Data', 'all_gill_depl.xlsx')))

#Vilém's function renaming
setnames(x = all_gill_depl,old = 'gearsize', new = 'Effort')#changing the name of releffortm to Effort to allow the function to run
all_gill_depl[, year := year(date_start)]

#Catches ####
# catches_db <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.catch_merged
#                                                                   WHERE  sa_samplingid IN ('",paste(all_gill_depl$sa_samplingid, collapse = "','"), "')
#                                                                   ;", sep = "")))
# write.xlsx(catches_db, here::here('Data', 'catches_db.xlsx'))
catches_db <- setDT(readxl::read_xlsx(here::here('Data', 'catches_db.xlsx')))

#separating the target species
catches_db <- merge(catches_db, specs[, .(sp_speciesid, sp_scientificname, sp_taxonomicorder)], by='sp_speciesid')
catches_db <- catches_db[!sp_speciesid == 'EA']
catches_carp <- catches_db[sp_scientificname %in% c("Cyprinus carpio", "Silurus glanis")]
catches_carp <- merge(catches_carp, all_gill_depl[, .(sa_samplingid, year)], by='sa_samplingid')
#write.xlsx(catches_carp, here::here('Data', 'catches_carp.xlsx'))

#Vpue
splitfactors <- c("sa_samplingid", "gg_gearid", "dl_layertype", "year")
cpue_carp <- getVPUE(samplings = all_gill_depl, catch = catches_carp, split.factors.catch = c("sp_scientificname"), 
                split.factors.samplings = splitfactors, value.var = "ct_abundancestar", 
                effort.colname = "Effort", id.colname = "sa_samplingid")
bpue_carp <- getVPUE(samplings = all_gill_depl, catch = catches_carp, split.factors.catch = c("sp_scientificname"), 
                split.factors.samplings = splitfactors, value.var = "ct_weightstar", 
                effort.colname = "Effort", id.colname = "sa_samplingid")

vpue_carp <- merge(cpue_carp, bpue_carp, by = c("sa_samplingid", "gg_gearid", "dl_layertype", "sp_scientificname", "year"))

#changing the name of variables
setnames(x = vpue_carp, old = c('ct_weightstar.mean','ct_weightstar.se', 'ct_abundancestar.mean','ct_abundancestar.se'),
         new = c('bpue_mean','bpue_se', 'cpue_mean','cpue_se'))#rename the outputs
#tranforming 1000m? per net
vpue_carp[, ':='(cpue_mean = cpue_mean*1000)]
#write.xlsx(vpue_carp, here::here('Data', 'vpue_carp.xlsx'))

sum_size_carp <- catches_carp[!ct_sl == 0,.(Mean = round(mean(ct_sl, na.rm = T), 2),
                                                     SE = round(plotrix::std.error(ct_sl), 2),
                                                     Max = max(ct_sl),
                                                     Min = min(ct_sl),
                                            Ind = sum(ct_abundance)),
                                       by =.(sp_scientificname, year)]
sum_size_carp$year <- as.factor(sum_size_carp$year)
#write.xlsx(sum_size_carp, here::here('Data', 'sum_size_carp.xlsx'))

ggplot(sum_size_carp, aes(x = year, y = Mean)) + 
  geom_ribbon(aes(ymin = Min, ymax = Max), fill = "grey70", alpha = 0.4) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), colour="black", width=.1, position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1), size=3, shape=21, fill="white") + 
  facet_grid(sp_scientificname ~ ., scales="free_y")+
  xlab("Year") +
  ylab("Size (mm)") +
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))

catches_carp[, ':='(ct_weight_kg = ct_weight/1000)]
sum_weight_carp <- catches_carp[!ct_weight_kg == 0,.(Mean = round(mean(ct_weight_kg, na.rm = T), 2),
                                            SE = round(plotrix::std.error(ct_weight_kg), 2),
                                            Max = round(max(ct_weight_kg), 2),
                                            Min = round(min(ct_weight_kg), 2)),
                              by =.(sp_scientificname, year)]
#write.xlsx(sum_weight_carp, here::here('Data', 'sum_weight_carp.xlsx'))

#####length weight correlation
catches_carp$logL <- log(catches_carp$ct_sl)
catches_carp$logW <- log(catches_carp$ct_weight)
catches_carpo <- catches_carp[sp_scientificname == "Cyprinus carpio"]
models <- lapply(split(catches_carp, catches_carp$sp_scientificname), 'lm', formula = logW ~ logL)
models2 <- lapply(split(catches_carpo, catches_carpo$year), 'lm', formula = logW ~ logL)
#need fix
# catches_carpo[sp_scientificname == "Cyprinus carpio", ct_wg_comp := predict.lm(object = models$`Cyprinus carpio`,
#                                                                    newdata = data.frame(logL = log(catches_carpo[sp_scientificname == "Cyprinus carpio"]$ct_sl),
#                                                                                         sp_scientificname = catches_carpo[sp_scientificname == "Cyprinus carpio"]$sp_scientificname))]
# catches_carpo[sp_scientificname == "Cyprinus carpio"]$ct_wg_comp <- exp(catches_carpo[sp_scientificname == "Cyprinus carpio"]$ct_wg_comp) 
# catches_carpo[sp_scientificname == "Cyprinus carpio"]$ct_wg_comp <- exp((summary(models$`Cyprinus carpio`)$sigma^2)/2) * catches_carpo[sp_scientificname == "Cyprinus carpio"]$ct_wg_comp
