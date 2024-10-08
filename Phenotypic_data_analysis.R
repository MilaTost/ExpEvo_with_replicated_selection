rm(list = ls())
### Load packages --------------------------------------------------------------
library(ggplot2)
library(stringr)
library(data.table)
library(gridExtra)
### Load data ------------------------------------------------------------------
setwd("C:/Users/mtost/Documents/Masterarbeit/data/Shoepeg 2021/")
data <- read.table("Shoepeg_all_years_2016_2017_2018_2020_2021.txt")
data$year <- as.factor(data$year)
### Prepare data ---------------------------------------------------------------
PlantHeight_group <- vector(length=(nrow(data)))
PlantHeight_group[which(data$Population == "Tall_plants_1" | data$Population == "Tall_plants_2")] <- "Selected for tall plant height"
PlantHeight_group[which(data$Population == "Short_plants_1" | data$Population == "Short_plants_2")] <- "Selected for short plant height"
data$PlantHeight_group <- PlantHeight_group
Population_labs <- vector(length=(nrow(data)))
Population_labs[which(data$Population == "Tall_plants_1")] <- "Tall plants 1"
Population_labs[which(data$Population == "Tall_plants_2")] <- "Tall plants 2"
Population_labs[which(data$Population == "Short_plants_1")] <- "Short plants 1"
Population_labs[which(data$Population == "Short_plants_2")] <- "Short plants 2"
data$Population_labs <- Population_labs
colnames(data) <- c("Plant_Nr","Population","Location","PlantHeight",
                    "EarHeight", "EarNum", "Tillers", "year",
                    "PlantHeight_group", "Population_labs")
#### Plotting ------------------------------------------------------------------
data <- as.data.table(data)
Short_1_2016 <- data[Population == "Short_plants_1" & year == "2016",]
Short_1_2017 <- data[Population == "Short_plants_1" & year == "2017",]
Short_1_2018 <- data[Population == "Short_plants_1" & year == "2018",]
Short_1_2020 <- data[Population == "Short_plants_1" & year == "2020",]

Short_2_2016 <- data[Population == "Short_plants_2" & year == "2016",]
Short_2_2017 <- data[Population == "Short_plants_2" & year == "2017",]
Short_2_2018 <- data[Population == "Short_plants_2" & year == "2018",]
Short_2_2020 <- data[Population == "Short_plants_2" & year == "2020",]

Tall_1_2016 <- data[Population == "Tall_plants_1" & year == "2016",]
Tall_1_2017 <- data[Population == "Tall_plants_1" & year == "2017",]
Tall_1_2018 <- data[Population == "Tall_plants_1" & year == "2018",]
Tall_1_2020 <- data[Population == "Tall_plants_1" & year == "2020",]

Tall_2_2016 <- data[Population == "Tall_plants_2" & year == "2016",]
Tall_2_2017 <- data[Population == "Tall_plants_2" & year == "2017",]
Tall_2_2018 <- data[Population == "Tall_plants_2" & year == "2018",]
Tall_2_2020 <- data[Population == "Tall_plants_2" & year == "2020",]
### Test for significance regarding the differences in plant height ------------
Short_plants_2016 <- data[PlantHeight_group == "Selected for short plant height" & year == "2016",]
Tall_plants_2016 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2016",]

Short_plants_2020 <- data[PlantHeight_group == "Selected for short plant height" & year == "2020",]
Tall_plants_2020 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2020",]

t.test(Short_plants_2020$PlantHeight,Short_plants_2016$PlantHeight,
       alternative = "less",
       paired = TRUE)
t.test(Tall_plants_2020$PlantHeight,Tall_plants_2016$PlantHeight,
       alternative = "greater",
       paired = TRUE)
results_t_test <- t.test(Tall_plants_2020$PlantHeight,Short_plants_2020$PlantHeight,
                         alternative = "two.sided",
                         paired = TRUE)
results_t_test$p.value
results_t_test$estimate
### Create plant height distribution plot V2 -----------------------------------
#### Prepare the data for plotting ---------------------------------------------
Shoepeg_2016 <- data[year == "2016",]
Shoepeg_2020 <- data[year == "2020",]
Shoepeg_2020_short_pop <- data[year == "2020" & PlantHeight_group == "Selected for short plant height",]
Shoepeg_2020_tall_pop <- data[year == "2020" & PlantHeight_group == "Selected for tall plant height",]
Shoepeg_2016$Density <- rep(1, nrow(Shoepeg_2016))
Shoepeg_2016$year <- as.factor(Shoepeg_2016$year)
Shoepeg_2020$PlantHeight_group<- as.factor(Shoepeg_2020$PlantHeight_group)
#### Prepare the data ----------------------------------------------------------
data <- as.data.table(data)
data$year <- as.factor(data$year)
data$Population <- as.factor(data$Population)

summarized_dt <- data[, Mean_PlantHeight := mean(PlantHeight, na.rm = TRUE), by = c("year", "Population")]
summarized_dt <- data[, SD_PlantHeight := sd(PlantHeight, na.rm = TRUE), by = c("year", "Population")]
summarized_dt <- data[, SE_PlantHeight := sd(PlantHeight, na.rm = TRUE)/length(!is.na(PlantHeight)), by = c("year", "Population")]
summarized_dt <- data[, CI1_Mean := Mean_PlantHeight - 0.95*(SD_PlantHeight / sqrt(nrow(data))), by = c("year", "Population")]
summarized_dt <- data[, CI2_Mean := Mean_PlantHeight + 0.95*(SD_PlantHeight / sqrt(nrow(data))), by = c("year", "Population")]
Calc_Response_to_Sel <- function(Measurement_gen_before,
                                 Measurement_gen_after){
  Response_to_sel <- Measurement_gen_after - Measurement_gen_before
  return(Response_to_sel)
}
##### Calculate it for Short plants 1 
Response_to_Sel_Gen1_Short1 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Short_plants_1" & year == "2016",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Short_plants_1" & year == "2017",Mean_PlantHeight])
Response_to_Sel_Gen2_Short1 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Short_plants_1" & year == "2017",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Short_plants_1" & year == "2018",Mean_PlantHeight])
Response_to_Sel_Gen3_Short1 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Short_plants_1" & year == "2018",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Short_plants_1" & year == "2020",Mean_PlantHeight])
##### Calculate it for Short plants 2
Response_to_Sel_Gen1_Short2 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Short_plants_2" & year == "2016",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Short_plants_2" & year == "2017",Mean_PlantHeight])
Response_to_Sel_Gen2_Short2 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Short_plants_2" & year == "2017",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Short_plants_2" & year == "2018",Mean_PlantHeight])
Response_to_Sel_Gen3_Short2 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Short_plants_2" & year == "2018",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Short_plants_2" & year == "2020",Mean_PlantHeight])
##### Calculate it for Tall plants 1
Response_to_Sel_Gen1_Tall1 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Tall_plants_1" & year == "2016",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Tall_plants_1" & year == "2017",Mean_PlantHeight])
Response_to_Sel_Gen2_Tall1 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Tall_plants_1" & year == "2017",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Tall_plants_1" & year == "2018",Mean_PlantHeight])
Response_to_Sel_Gen3_Tall1 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Tall_plants_1" & year == "2018",Mean_PlantHeight],
                                                    Measurement_gen_after = summarized_dt[Population == "Tall_plants_1" & year == "2020",Mean_PlantHeight])
##### Calculate it for Tall plants 2
Response_to_Sel_Gen1_Tall2 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Tall_plants_2" & year == "2016",Mean_PlantHeight],
                                                   Measurement_gen_after = summarized_dt[Population == "Tall_plants_2" & year == "2017",Mean_PlantHeight])
Response_to_Sel_Gen2_Tall2 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Tall_plants_2" & year == "2017",Mean_PlantHeight],
                                                   Measurement_gen_after = summarized_dt[Population == "Tall_plants_2" & year == "2018",Mean_PlantHeight])
Response_to_Sel_Gen3_Tall2 <- Calc_Response_to_Sel(Measurement_gen_before = summarized_dt[Population == "Tall_plants_2" & year == "2018",Mean_PlantHeight],
                                                   Measurement_gen_after = summarized_dt[Population == "Tall_plants_2" & year == "2020",Mean_PlantHeight])
### Test for significance regarding the differences in plant height ------------
Short_plants_2016 <- data[PlantHeight_group == "Selected for short plant height" & year == "2016",]
Tall_plants_2016 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2016",]

Short_plants_2017 <- data[PlantHeight_group == "Selected for short plant height" & year == "2017",]
Tall_plants_2017 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2017",]

Short_plants_2018 <- data[PlantHeight_group == "Selected for short plant height" & year == "2018",]
Tall_plants_2018 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2018",]

Short_plants_2020 <- data[PlantHeight_group == "Selected for short plant height" & year == "2020",]
Tall_plants_2020 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2020",]

Results_T_test_2016 <- t.test(Tall_plants_2016$PlantHeight,Short_plants_2016$PlantHeight,
                              alternative = "two.sided",
                              paired = TRUE)
Results_T_test_2017 <- t.test(Tall_plants_2017$PlantHeight,Short_plants_2017$PlantHeight,
                              alternative = "two.sided",
                              paired = TRUE)
Results_T_test_2018 <- t.test(Tall_plants_2018$PlantHeight,Short_plants_2018$PlantHeight,
                              alternative = "two.sided",
                              paired = TRUE)
Results_T_test_2020 <- t.test(Tall_plants_2020$PlantHeight,Short_plants_2020$PlantHeight,
                              alternative = "two.sided",
                              paired = TRUE)
#### Combine all data sets -----------------------------------------------------
summarized_dt_short <- rbind(c("Short_plants_1", "2016", NA, summarized_dt[Population == "Short_plants_1" & year == "2016", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_1" & year == "2016",SD_PlantHeight][1], Results_T_test_2016$estimate, Results_T_test_2016$p.value, Results_T_test_2016$conf.int[1], Results_T_test_2016$conf.int[2]), 
                             c("Short_plants_2", "2016", NA, summarized_dt[Population == "Short_plants_2" & year == "2016", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_2" & year == "2016",SD_PlantHeight][1], Results_T_test_2016$estimate, Results_T_test_2016$p.value, Results_T_test_2016$conf.int[1], Results_T_test_2016$conf.int[2]),
                             c("Tall_plants_1", "2016", NA, summarized_dt[Population == "Tall_plants_1" & year == "2016", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_1" & year == "2016",SD_PlantHeight][1], Results_T_test_2016$estimate, Results_T_test_2016$p.value, Results_T_test_2016$conf.int[1], Results_T_test_2016$conf.int[2]), 
                             c("Tall_plants_2", "2016", NA, summarized_dt[Population == "Tall_plants_2" & year == "2016", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_2" & year == "2016",SD_PlantHeight][1], Results_T_test_2016$estimate, Results_T_test_2016$p.value, Results_T_test_2016$conf.int[1], Results_T_test_2016$conf.int[2]),
                             c("Short_plants_1", "2017", Response_to_Sel_Gen1_Short1[1], summarized_dt[Population == "Short_plants_1" & year == "2017", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_1" & year == "2017",SD_PlantHeight][1], Results_T_test_2017$estimate, Results_T_test_2017$p.value, Results_T_test_2017$conf.int[1], Results_T_test_2017$conf.int[2]), 
                             c("Short_plants_2", "2017", Response_to_Sel_Gen1_Short2[1], summarized_dt[Population == "Short_plants_2" & year == "2017", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_2" & year == "2017",SD_PlantHeight][1], Results_T_test_2017$estimate, Results_T_test_2017$p.value, Results_T_test_2017$conf.int[1], Results_T_test_2017$conf.int[2]),
                             c("Tall_plants_1", "2017", Response_to_Sel_Gen1_Tall1[1], summarized_dt[Population == "Tall_plants_1" & year == "2017", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_1" & year == "2017",SD_PlantHeight][1], Results_T_test_2017$estimate, Results_T_test_2017$p.value, Results_T_test_2017$conf.int[1], Results_T_test_2017$conf.int[2]), 
                             c("Tall_plants_2", "2017", Response_to_Sel_Gen1_Tall2[1], summarized_dt[Population == "Tall_plants_2" & year == "2017", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_2" & year == "2017",SD_PlantHeight][1], Results_T_test_2017$estimate, Results_T_test_2017$p.value, Results_T_test_2017$conf.int[1], Results_T_test_2017$conf.int[2]),
                             c("Short_plants_1", "2018", Response_to_Sel_Gen2_Short1[1], summarized_dt[Population == "Short_plants_1" & year == "2018", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_1" & year == "2018",SD_PlantHeight][1], Results_T_test_2018$estimate, Results_T_test_2018$p.value, Results_T_test_2018$conf.int[1], Results_T_test_2018$conf.int[2]), 
                             c("Short_plants_2", "2018", Response_to_Sel_Gen2_Short2[1], summarized_dt[Population == "Short_plants_2" & year == "2018", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_2" & year == "2018",SD_PlantHeight][1], Results_T_test_2018$estimate, Results_T_test_2018$p.value, Results_T_test_2018$conf.int[1], Results_T_test_2018$conf.int[2]),
                             c("Tall_plants_1", "2018", Response_to_Sel_Gen2_Tall1[1], summarized_dt[Population == "Tall_plants_1" & year == "2018", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_1" & year == "2018",SD_PlantHeight][1], Results_T_test_2018$estimate, Results_T_test_2018$p.value, Results_T_test_2018$conf.int[1], Results_T_test_2018$conf.int[2]), 
                             c("Tall_plants_2", "2018", Response_to_Sel_Gen2_Tall2[1], summarized_dt[Population == "Tall_plants_2" & year == "2018", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_2" & year == "2018",SD_PlantHeight][1], Results_T_test_2018$estimate, Results_T_test_2018$p.value, Results_T_test_2018$conf.int[1], Results_T_test_2018$conf.int[2]),
                             c("Short_plants_1", "2020", Response_to_Sel_Gen3_Short1[1], summarized_dt[Population == "Short_plants_1" & year == "2020", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_1" & year == "2020",SD_PlantHeight][1], Results_T_test_2020$estimate, Results_T_test_2020$p.value, Results_T_test_2020$conf.int[1], Results_T_test_2020$conf.int[2]), 
                             c("Short_plants_2", "2020", Response_to_Sel_Gen3_Short2[1], summarized_dt[Population == "Short_plants_2" & year == "2020", Mean_PlantHeight][1], summarized_dt[Population == "Short_plants_2" & year == "2020",SD_PlantHeight][1], Results_T_test_2020$estimate, Results_T_test_2020$p.value, Results_T_test_2020$conf.int[1], Results_T_test_2020$conf.int[2]),
                             c("Tall_plants_1", "2020", Response_to_Sel_Gen3_Tall1[1], summarized_dt[Population == "Tall_plants_1" & year == "2020", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_1" & year == "2020",SD_PlantHeight][1], Results_T_test_2020$estimate, Results_T_test_2020$p.value, Results_T_test_2020$conf.int[1], Results_T_test_2020$conf.int[2]), 
                             c("Tall_plants_2", "2020", Response_to_Sel_Gen3_Tall2[1], summarized_dt[Population == "Tall_plants_2" & year == "2020", Mean_PlantHeight][1], summarized_dt[Population == "Tall_plants_2" & year == "2020",SD_PlantHeight][1], Results_T_test_2020$estimate, Results_T_test_2020$p.value, Results_T_test_2020$conf.int[1], Results_T_test_2020$conf.int[2]))
summarized_dt_short <- as.data.frame(summarized_dt_short)
colnames(summarized_dt_short) <- c("Population", "Year", "Response_to_selection","MeanPlantHeight","SDPlantHeight", "T_test_results_estimate", "T_test_results_p_value", "T_test_results_CI_1", "T_test_results_CI_2")
summarized_dt_short
summarized_dt_short$Response_to_selection <- as.numeric(summarized_dt_short$Response_to_selection)
summarized_dt_short$MeanPlantHeight <- as.numeric(summarized_dt_short$MeanPlantHeight)
summarized_dt_short$SDPlantHeight <- as.numeric(summarized_dt_short$SDPlantHeight)
summarized_dt_short$T_test_results_estimate <- as.numeric(summarized_dt_short$T_test_results_estimate)
summarized_dt_short$T_test_results_p_value <- as.numeric(summarized_dt_short$T_test_results_p_value)
summarized_dt_short <- as.data.table(summarized_dt_short)
summarized_dt_short$T_test_results_p_value <- round(summarized_dt_short$T_test_results_p_value,10)
summarized_dt_short$Response_to_selection <- round(summarized_dt_short$Response_to_selection,2)
summarized_dt_short$T_test_results_CI_1 <- as.numeric(summarized_dt_short$T_test_results_CI_1)
summarized_dt_short$T_test_results_CI_2 <- as.numeric(summarized_dt_short$T_test_results_CI_2)
#### Load plotting parameters --------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 12
color_short_1 <- "palegreen3"
color_short_2 <- "palegreen4"
color_tall_1 <- "plum3"
color_tall_2<- "plum4"
colors_Shoepeg <- c("Short_plants_1" = color_short_1,
                    "Short_plants_2" = color_short_2,
                    "Tall_plants_1" = color_tall_1,
                    "Tall_plants_2" = color_tall_2)
#### Plotting ------------------------------------------------------------------
Response_to_sel_plot1 <- ggplot()+
  stat_summary(data = data,
               fun = "mean", 
               geom = "line",
               aes(y = PlantHeight, x = year, group = Population,
                   color = Population))+
  stat_summary(data = data,
               fun = "mean", geom = "point",
               aes(y = PlantHeight, x = year, group = Population,
                   color = Population))+
  geom_errorbar(data = summarized_dt, width = 0.1,
                aes(ymin = CI1_Mean,
                    ymax = CI2_Mean,
                    x = year,
                    group = Population,
                    color = Population))+
  annotate("text",
           label = paste("p-value ~", c(round(Results_T_test_2016$p.value,10),
                                        round(Results_T_test_2017$p.value,10),
                                        round(Results_T_test_2018$p.value,10),
                                        round(Results_T_test_2020$p.value,10)),"***"),
           y = mean(data$PlantHeight, na.rm = TRUE),
           x = 1.2:4.2,
           angle = 90,
           colour = "black", family = "my", size = 3)+
  annotate("segment",
           y = c(as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2016", Mean_PlantHeight][1]),
                 as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2017", Mean_PlantHeight][1]),
                 as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2018", Mean_PlantHeight][1]),
                 as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2020", Mean_PlantHeight][1])),
           yend = c(as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2016", Mean_PlantHeight][1] + Results_T_test_2016$estimate),
                    as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2017", Mean_PlantHeight][1] + Results_T_test_2017$estimate),
                    as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2018", Mean_PlantHeight][1] + Results_T_test_2018$estimate),
                    as.numeric(summarized_dt[Population == "Short_plants_1" & year == "2020", Mean_PlantHeight][1] + Results_T_test_2020$estimate)),
           x = 1.1:4.1, xend = 1.1:4.1,
           colour = "black")+
  annotate("text", 
           label = paste(c(summarized_dt_short[Population == "Short_plants_1" & Year == "2017", Response_to_selection][1],
                           summarized_dt_short[Population == "Short_plants_1" & Year == "2018", Response_to_selection][1],
                           summarized_dt_short[Population == "Short_plants_1" & Year == "2020", Response_to_selection][1])),
           y = c(summarized_dt_short[Population == "Short_plants_1" & Year == "2017", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Short_plants_1" & Year == "2018", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Short_plants_1" & Year == "2020", MeanPlantHeight][1]-5),
           x = 1.35:3.35,
           colour = color_short_1, family = "my", size = 3)+
  annotate("text", 
           label = paste(c(summarized_dt_short[Population == "Short_plants_2" & Year == "2017", Response_to_selection][1],
                           summarized_dt_short[Population == "Short_plants_2" & Year == "2018", Response_to_selection][1],
                           summarized_dt_short[Population == "Short_plants_2" & Year == "2020", Response_to_selection][1])),
           y = c(summarized_dt_short[Population == "Short_plants_2" & Year == "2017", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Short_plants_2" & Year == "2018", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Short_plants_2" & Year == "2020", MeanPlantHeight][1]-5),
           x = 1.75:3.75,
           colour = color_short_2, family = "my", size = 3)+
  annotate("text", 
           label = paste(c(summarized_dt_short[Population == "Tall_plants_1" & Year == "2017", Response_to_selection][1],
                           summarized_dt_short[Population == "Tall_plants_1" & Year == "2018", Response_to_selection][1],
                           summarized_dt_short[Population == "Tall_plants_1" & Year == "2020", Response_to_selection][1])),
           y = c(summarized_dt_short[Population == "Tall_plants_1" & Year == "2017", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Tall_plants_1" & Year == "2018", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Tall_plants_1" & Year == "2020", MeanPlantHeight][1]-5),
           x = 1.35:3.35,
           colour = color_tall_1, family = "my", size = 3)+
  annotate("text", 
           label = paste(c(summarized_dt_short[Population == "Tall_plants_2" & Year == "2017", Response_to_selection][1],
                           summarized_dt_short[Population == "Tall_plants_2" & Year == "2018", Response_to_selection][1],
                           summarized_dt_short[Population == "Tall_plants_2" & Year == "2020", Response_to_selection][1])),
           y = c(summarized_dt_short[Population == "Tall_plants_2" & Year == "2017", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Tall_plants_2" & Year == "2018", MeanPlantHeight][1]-5,
                 summarized_dt_short[Population == "Tall_plants_2" & Year == "2020", MeanPlantHeight][1]-5),
           x = 1.75:3.75,
           colour = color_tall_2, family = "my", size = 3)+
  theme(text = element_text(size = font_size, family = "my", colour = "black"),
        title = element_text(size = font_size-2, family = "my", colour = "black",
                             face = "bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = font_size, family = "my", colour = "black", face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black", face = "bold"),
        axis.title.x= element_blank(),
        legend.position = "right",
        legend.text = element_text(size = font_size-2, family = "my", colour = "black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        axis.text = element_text(size = font_size, family = "my", colour = "black"))+
  labs(y = "Plant height (cm)")+
  scale_color_manual(values = colors_Shoepeg,
                     labels = c("Short plants 1", "Short plants 2",
                                "Tall plants 1", "Tall plants 2"))
Response_to_sel_plot1
ggsave("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/2023_06_NEW_V5_Poster_Plant_height_measurements.png",
       Response_to_sel_plot1,
       device = png,
       height = 4, width = 7,
       dpi =1800) 
### S3: plot to show the preocedure of simultaneously flowering male plants ----
#### Prepare the data ----------------------------------------------------------
### Flowering dates ------------------------------------------------------------
## Correlation between PH and flowering dates
## how to translate the following part in a formula
library(stringr)
library(ggplot2)
library(data.table)
library(foreach)
setwd("C:/Users/mtost/Documents/Masterarbeit/Material_and_methods/Field/")
shoepag <- readxl::read_excel("field_data_collection_shoepag_12_09_2020.xlsx", sheet = 3, col_names = TRUE)

shoepag <- as.data.table(shoepag)
shoepag$population <- str_sub(shoepag$Barcode_of_plant,1,1)
Short_1 <- shoepag[population == "1",]
Tall_2 <- shoepag[population == "2",]
Tall_1 <- shoepag[population == "3",]
Short_2 <- shoepag[population == "4",]
### Translate flowering dates for plotting and correlation calc ----------------
##### TASSELS 
get_flowering_days_tassel <- function(data){
  # Translate the tassel dates
  day_Tassel <- str_sub(data$Flowering_date_Tassel,1,2)
  month_Tassel <- str_sub(data$Flowering_date_Tassel,4,5)
  date_Tassel <- paste("2020",month_Tassel,day_Tassel, sep = "-", collapse = NULL, recycle0 = F)
  date_Tassel[which(date_Tassel=="2020-NA-NA")] <- NA
  days_to_flowering_Tassels <- as.Date(date_Tassel) - as.Date("2020-04-28")
  return(days_to_flowering_Tassels)
}
Short_1_days_to_flowering_T <- get_flowering_days_tassel(data = Short_1)
Short_2_days_to_flowering_T <- get_flowering_days_tassel(data = Short_2)
Tall_1_days_to_flowering_T <- get_flowering_days_tassel(data = Tall_1)
Tall_2_days_to_flowering_T <- get_flowering_days_tassel(data = Tall_2)
Short_1$Days_to_Tassel_flowering <- Short_1_days_to_flowering_T
Short_2$Days_to_Tassel_flowering <- Short_2_days_to_flowering_T
Tall_1$Days_to_Tassel_flowering <- Tall_1_days_to_flowering_T
Tall_2$Days_to_Tassel_flowering <- Tall_2_days_to_flowering_T
##### SILKS 
get_flowering_days_silk <- function(data){
  # Translate the silk dates
  day_Silk <- str_sub(data$Flowering_date_Silk,1,2)
  month_Silk <- str_sub(data$Flowering_date_Silk,4,5)
  date_Silk <- paste("2020",month_Silk, day_Silk, sep = "-", collapse = NULL, recycle0 = F)
  date_Silk[which(date_Silk=="2020--NA")] <- NA
  days_to_mature_Silks <- as.Date(date_Silk) - as.Date("2020-04-28")
  return(days_to_mature_Silks)
}
Short_1_days_to_flowering_S <- get_flowering_days_silk(data = Short_1)
Short_2_days_to_flowering_S <- get_flowering_days_silk(data = Short_2)
Tall_1_days_to_flowering_S <- get_flowering_days_silk(data = Tall_1)
Tall_2_days_to_flowering_S <- get_flowering_days_silk(data = Tall_2)
Short_1$Days_to_Silk_maturity <- Short_1_days_to_flowering_S
Short_2$Days_to_Silk_maturity <- Short_2_days_to_flowering_S
Tall_1$Days_to_Silk_maturity <- Tall_1_days_to_flowering_S
Tall_2$Days_to_Silk_maturity <- Tall_2_days_to_flowering_S
### Get the flowering time intervals for Ne calc -------------------------------
get_flowering_time_intervals <- function(data,
                                         days_plants_are_recipient){
  # Translate the tassel dates
  day_Tassel <- str_sub(data$Flowering_date_Tassel,1,2)
  month_Tassel <- str_sub(data$Flowering_date_Tassel,4,5)
  date_Tassel <- paste("2020",month_Tassel,day_Tassel, sep = "-", collapse = NULL, recycle0 = F)
  date_Tassel[which(date_Tassel=="2020-NA-NA")] <- NA
  days_to_flowering_Tassels <- as.Date(date_Tassel) - as.Date("2020-04-28")
  # Translate the silk dates
  day_Silk <- str_sub(data$Flowering_date_Silk,1,2)
  month_Silk <- str_sub(data$Flowering_date_Silk,4,5)
  date_Silk <- paste("2020",month_Silk, day_Silk, sep = "-", collapse = NULL, recycle0 = F)
  date_Silk[which(date_Silk=="2020--NA")] <- NA
  days_to_mature_Silks <- as.Date(date_Silk) - as.Date("2020-04-28")
  # when is the maturity of silks ready
  mean_FT <- mean(days_to_mature_Silks, na.rm = TRUE)
  SD_FT <- sd(days_to_mature_Silks, na.rm = TRUE)
  calc_silk_time_int <- function(perc_measurements){
    lower_quantile <- quantile(days_to_mature_Silks, probs = perc_measurements, na.rm = TRUE)
    upper_quantile <- quantile(days_to_mature_Silks, probs = 1-perc_measurements, na.rm = TRUE)
    silk_time_interval <- c(perc_measurements, 1-perc_measurements, lower_quantile, upper_quantile)
    return(silk_time_interval)
  }
  silk_time_intervals <- foreach(perc_measurements = seq(0.01,0.5,0.01), .combine = rbind) %do% calc_silk_time_int(perc_measurements)
  silk_time_intervals <- as.data.frame(silk_time_intervals)
  colnames(silk_time_intervals) <- c("Lower_quantile_in_percentage", "Upper_quantile_in_percentage","Lower_quantile", "Upper_quantile")
  silk_time_intervals$time_int_in_days <- silk_time_intervals$Upper_quantile - silk_time_intervals$Lower_quantile
  majority_of_silks_mature <- silk_time_intervals[min(which(silk_time_intervals$time_int_in_days <= 5)),]
  # How many tassels flower during that time
  get_flowering_tassels <- function(i){
    tassels_flower_then <- length(which(days_to_flowering_Tassels >= silk_time_intervals$Lower_quantile[i] & days_to_flowering_Tassels <= silk_time_intervals$Upper_quantile[i]))
    return(tassels_flower_then)
  }
  Nr_of_flowering_tassels <- foreach(i = 1:nrow(silk_time_intervals), .combine = rbind) %do% get_flowering_tassels(i)
  silk_time_intervals$Nr_of_flowering_tassels <- Nr_of_flowering_tassels
  return(silk_time_intervals)
}
Short_1_results_flowering_times <- get_flowering_time_intervals(data = Short_1, days_plants_are_recipient = 5)
Short_2_results_flowering_times <- get_flowering_time_intervals(data = Short_2, days_plants_are_recipient = 5)
Tall_1_results_flowering_times <- get_flowering_time_intervals(data = Tall_1, days_plants_are_recipient = 5)
Tall_2_results_flowering_times <- get_flowering_time_intervals(data = Tall_2, days_plants_are_recipient = 5)

Short_1_male_parents <- Short_1_results_flowering_times$Nr_of_flowering_tassels[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))]
Short_2_male_parents <- Short_2_results_flowering_times$Nr_of_flowering_tassels[min(which(Short_2_results_flowering_times$time_int_in_days <= 5))]
Tall_1_male_parents <- Tall_1_results_flowering_times$Nr_of_flowering_tassels[min(which(Tall_1_results_flowering_times$time_int_in_days <= 5))]
Tall_2_male_parents <- Tall_2_results_flowering_times$Nr_of_flowering_tassels[min(which(Tall_2_results_flowering_times$time_int_in_days <= 5))]

simFlowPer_Short1 <- (Short_1_male_parents/96)
simFlowPer_Short2 <- (Short_2_male_parents/96)
simFlowPer_Tall1 <- (Tall_1_male_parents/96)
simFlowPer_Tall2 <- (Tall_2_male_parents/96)
(nr_male_parents_Short1 <- (Short_1_male_parents/96)*5000)
(nr_male_parents_Short2 <- (Short_2_male_parents/96)*5000)
(nr_male_parents_Tall1 <- (Tall_1_male_parents/96)*5000)
(nr_male_parents_Tall2 <- (Tall_2_male_parents/96)*5000)
simFlowPer_Short1;simFlowPer_Short2; simFlowPer_Tall1; simFlowPer_Tall2
5000*mean(simFlowPer_Short1, simFlowPer_Short2, simFlowPer_Tall1, simFlowPer_Tall2)
### Plotting of flowering dates ------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 12

Short_1$Days_to_Silk_maturity <- as.numeric(Short_1$Days_to_Silk_maturity)
Short_2$Days_to_Silk_maturity <- as.numeric(Short_2$Days_to_Silk_maturity)
Tall_1$Days_to_Silk_maturity <- as.numeric(Tall_1$Days_to_Silk_maturity)
Tall_2$Days_to_Silk_maturity <- as.numeric(Tall_2$Days_to_Silk_maturity)

Short_1_demo_plot <- ggplot()+
  geom_histogram(data = Short_1, aes(x = Days_to_Tassel_flowering),
                 fill ="royalblue1", colour = "royalblue4",alpha = 0.6)+
  geom_histogram(data = Short_1, aes(x = Days_to_Silk_maturity), 
                 fill ="indianred1", colour = "indianred4", alpha = 0.6)+
  annotate(geom = "rect",
           xmin = Short_1_results_flowering_times$Lower_quantile[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))],
           xmax = Short_1_results_flowering_times$Upper_quantile[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))],
           ymin = 0, ymax = 30, fill = "red", alpha = 0.2)+
  theme(text = element_text(size = font_size, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(colour ="black", size = font_size, face = "bold"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.x= element_text(colour ="black", size = font_size, face = "bold"),
        axis.title.y= element_text(colour ="black", size = font_size, face = "bold"))+
  geom_vline(xintercept = Short_1_results_flowering_times$Lower_quantile[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  geom_vline(xintercept = Short_1_results_flowering_times$Upper_quantile[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  ylim(0,30)+
  annotate("text",
           label = paste(Short_1_results_flowering_times$Nr_of_flowering_tassels[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))],"tassels out of 96 are flowering", "\n","in this time interval"),
           y = 27,
           x = mean(Short_1$Days_to_Silk_maturity, na.rm = TRUE),
           colour = "darkred", family = "my", size = 3)+
  labs(x = "Days to flowering or maturity", y = "Density", tag = "A")
Short_2_demo_plot <- ggplot()+
  geom_histogram(data = Short_2, aes(x = Days_to_Tassel_flowering),
                 fill ="royalblue1", colour = "royalblue4",alpha = 0.6)+
  geom_histogram(data = Short_2, aes(x = Days_to_Silk_maturity), 
                 fill ="indianred1", colour = "indianred4", alpha = 0.6)+
  annotate(geom = "rect",
           xmin = Short_2_results_flowering_times$Lower_quantile[min(which(Short_2_results_flowering_times$time_int_in_days <= 5))],
           xmax = Short_2_results_flowering_times$Upper_quantile[min(which(Short_2_results_flowering_times$time_int_in_days <= 5))],
           ymin = 0, ymax = 30, fill = "red", alpha = 0.2)+
  theme(text = element_text(size = font_size, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.tag = element_text(colour ="black", size = font_size, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.x= element_text(colour ="black", size = font_size, face = "bold"),
        axis.title.y= element_text(colour ="black", size = font_size, face = "bold"))+
  geom_vline(xintercept = Short_2_results_flowering_times$Lower_quantile[min(which(Short_2_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  geom_vline(xintercept = Short_2_results_flowering_times$Upper_quantile[min(which(Short_2_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  ylim(0,30)+
  annotate("text",
           label = paste(Short_2_results_flowering_times$Nr_of_flowering_tassels[min(which(Short_1_results_flowering_times$time_int_in_days <= 5))],"tassels out of 96 are flowering", "\n","in this time interval"),
           y = 27,
           x = mean(Short_2$Days_to_Silk_maturity, na.rm = TRUE),
           colour = "darkred", family = "my", size = 3)+
  labs(x = "Days to flowering or maturity", y = "Density", tag = "B")
Tall_1_demo_plot <- ggplot()+
  geom_histogram(data = Tall_1, aes(x = Days_to_Tassel_flowering),
                 fill ="royalblue1", colour = "royalblue4",alpha = 0.6)+
  geom_histogram(data = Tall_1, aes(x = Days_to_Silk_maturity), 
                 fill ="indianred1", colour = "indianred4", alpha = 0.6)+
  annotate(geom = "rect",
           xmin = Tall_1_results_flowering_times$Lower_quantile[min(which(Tall_1_results_flowering_times$time_int_in_days <= 5))],
           xmax = Tall_1_results_flowering_times$Upper_quantile[min(which(Tall_1_results_flowering_times$time_int_in_days <= 5))],
           ymin = 0, ymax = 30, fill = "red", alpha = 0.2)+
  theme(text = element_text(size = font_size, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.tag = element_text(colour ="black", size = font_size, face = "bold"),
        axis.title.x= element_text(colour ="black", size = font_size, face = "bold"),
        axis.title.y= element_text(colour ="black", size = font_size, face = "bold"))+
  geom_vline(xintercept = Tall_1_results_flowering_times$Lower_quantile[min(which(Tall_1_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  geom_vline(xintercept = Tall_1_results_flowering_times$Upper_quantile[min(which(Tall_1_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  ylim(0,30)+
  annotate("text",
           label = paste(Tall_1_results_flowering_times$Nr_of_flowering_tassels[min(which(Tall_1_results_flowering_times$time_int_in_days <= 5))],"tassels out of 96 are flowering", "\n","in this time interval"),
           y = 27,
           x = mean(Tall_1$Days_to_Silk_maturity, na.rm = TRUE),
           colour = "darkred", family = "my", size = 3)+
  labs(x = "Days to flowering or maturity", y = "Density", tag = "C")
Tall_1_demo_plot
Tall_2_demo_plot <- ggplot()+
  geom_histogram(data = Tall_2, aes(x = Days_to_Tassel_flowering),
                 fill ="royalblue1", colour = "royalblue4",alpha = 0.6)+
  geom_histogram(data = Tall_2, aes(x = Days_to_Silk_maturity), 
                 fill ="indianred1", colour = "indianred4", alpha = 0.6)+
  annotate(geom = "rect",
           xmin = Tall_2_results_flowering_times$Lower_quantile[min(which(Tall_2_results_flowering_times$time_int_in_days <= 5))],
           xmax = Tall_2_results_flowering_times$Upper_quantile[min(which(Tall_2_results_flowering_times$time_int_in_days <= 5))],
           ymin = 0, ymax = 30, fill = "red", alpha = 0.2)+
  theme(text = element_text(size = font_size, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.tag = element_text(colour ="black", size = font_size, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.x= element_text(colour ="black", size = font_size, face = "bold"),
        axis.title.y= element_text(colour ="black", size = font_size, face = "bold"))+
  geom_vline(xintercept = Tall_2_results_flowering_times$Lower_quantile[min(which(Tall_2_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  geom_vline(xintercept = Tall_2_results_flowering_times$Upper_quantile[min(which(Tall_2_results_flowering_times$time_int_in_days <= 5))], colour = "red3", linetype = "dashed")+
  ylim(0,30)+
  annotate("text",
           label = paste(Tall_2_results_flowering_times$Nr_of_flowering_tassels[min(which(Tall_2_results_flowering_times$time_int_in_days <= 5))],"tassels out of 96 are flowering", "\n","in this time interval"),
           y = 27,
           x = mean(Tall_2$Days_to_Silk_maturity, na.rm = TRUE),
           colour = "darkred", family = "my", size = 3)+
  labs(x = "Days to flowering or maturity", y = "Density", tag = "D")
Tall_2_demo_plot
### this part is just used to get a legend -------------------------------------
library(agridat)
ilri.sheep
str(ilri.sheep)
get_legend_plot <- ggplot()+
  geom_histogram(data = ilri.sheep, aes(x = weanwt, fill = sex, colour = sex, y = ..density..), alpha = 0.6)+
  scale_fill_manual(values = c("indianred1","royalblue1"), name = "Flowering", labels = c("Silk","Tassel"))+
  scale_colour_manual(values = c("indianred4","royalblue4"), name = "Flowering", labels = c("Silk","Tassel"))+
  theme(text = element_text(size =12, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.x= element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.y= element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, colour ="black", size = 12),
        legend.title = element_text(hjust = 0.5, colour ="black", size = 12, face = "bold"),
        legend.box = "vertical")+
  labs(x = "Flowering date")
get_legend_plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend_flower <- get_legend(get_legend_plot)

flowering_tassel_all_pop <- gridExtra::grid.arrange(Short_1_demo_plot,
                                                    Short_2_demo_plot,
                                                    Tall_1_demo_plot,
                                                    Tall_2_demo_plot,
                                                    legend_flower,
                                                    nrow = 3, heights = c(4,4,0.6),
                                                    layout_matrix = rbind(c(1,2), c(3,4),
                                                                          c(5,5)))
flowering_tassel_all_pop
ggsave("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/S3_Flowering_dates_ind_quant.png",flowering_tassel_all_pop,
       height = 6,
       width = 10,
       dpi =1800)
### Calculate the correlation between flowering time and plant height ----------
str(Short_1)
Short_1$Plant_Height_in_cm <- as.numeric(Short_1$Plant_Height_in_cm)
Short_2$Plant_Height_in_cm <- as.numeric(Short_2$Plant_Height_in_cm)
Tall_1$Plant_Height_in_cm <- as.numeric(Tall_1$Plant_Height_in_cm)
Tall_2$Plant_Height_in_cm <- as.numeric(Tall_2$Plant_Height_in_cm)
Short_1$Days_to_Tassel_flowering <- as.numeric(Short_1$Days_to_Tassel_flowering)
Short_2$Days_to_Tassel_flowering <- as.numeric(Short_2$Days_to_Tassel_flowering)
Tall_1$Days_to_Tassel_flowering <- as.numeric(Tall_1$Days_to_Tassel_flowering)
Tall_2$Days_to_Tassel_flowering <- as.numeric(Tall_2$Days_to_Tassel_flowering)
Short_1$Flowering_time <- (Short_1$Days_to_Tassel_flowering + Short_1$Days_to_Silk_maturity)/2
Short_2$Flowering_time <- (Short_2$Days_to_Tassel_flowering + Short_2$Days_to_Silk_maturity)/2
Tall_1$Flowering_time <- (Tall_1$Days_to_Tassel_flowering + Tall_1$Days_to_Silk_maturity)/2
Tall_2$Flowering_time <- (Tall_2$Days_to_Tassel_flowering + Tall_2$Days_to_Silk_maturity)/2

cor_res_Short1 <- cor.test(Short_1$Flowering_time, Short_1$Plant_Height_in_cm)
cor_res_Short2 <- cor.test(Short_2$Flowering_time, Short_2$Plant_Height_in_cm)
cor_res_Tall1 <- cor.test(Tall_1$Flowering_time, Tall_1$Plant_Height_in_cm)
cor_res_Tall2 <- cor.test(Tall_2$Flowering_time, Tall_2$Plant_Height_in_cm)
round(cor_res_Short1$estimate,4);round(cor_res_Short1$p.value, 4)
round(cor_res_Short2$estimate,4);round(cor_res_Short2$p.value, 4)
round(cor_res_Tall1$estimate,4);round(cor_res_Tall1$p.value, 4)
round(cor_res_Tall2$estimate,4);round(cor_res_Tall2$p.value, 4)