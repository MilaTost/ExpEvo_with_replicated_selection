rm(list = ls())
windowsFonts(my = windowsFont('Calibri'))
library(ggplot2)
library(stringr)
library(data.table)

setwd("//winfs-uni.top.gwdg.de/mtost$/Masterarbeit/Data_analysis/phenotypic_data/")
data <- read.table("data_all_years.txt")
windowsFonts(my = windowsFont('Calibri'))
plot_phenotypic_measurements <- function(data,
                                         phenotype,
                                         names_of_subpopulations,
                                         years,
                                         name_of_low_phenotype,
                                         name_of_high_phenotype){
  PlantHeight_group <- vector(length=(nrow(data)))
  PlantHeight_group[which(str_detect(name_of_low_phenotype))] <- name_of_low_phenotype
  PlantHeight_group[which(str_detect(name_of_high_phenotype))] <- name_of_high_phenotype
  data$PlantHeight_group <- PlantHeight_group
  data <- as.data.table(data)
  fill_low_pheno <- "darkseagreen2"
  colour_low_pheno <- "darkseagreen3"
  fill_high_pheno <- "thistle2"
  colour_high_pheno <- "thistle3"
  create_plot_low_pheno <- function(data,
                                    name_of_low_phenotype,
                                    phenotype,
                                    i){
    dt_plot <- data[Population == names_of_subpopulations[i] & year == years[i],]
    rand1 <- quantile(phenotype, 0.05, na.rm = TRUE)
    plot1 <- ggplot(data = dt_plot, aes(x = phenotype))+
      geom_density(color = colour_low_pheno, fill = fill_low_pheno)+
      theme(text = element_text(size =12, family = "my"),
            panel.background = element_rect(fill = "white"),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "grey81"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "grey100"),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face ="bold"),
            axis.title.x= element_blank(),
            axis.text = element_text(size =12, family = "my", colour = "black"),
            axis.title.y = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
      labs(y = names_of_subpopulations[i])+
      ggtitle(year[i])+
      geom_vline(aes(xintercept=quantile(dt_plot$phenotype, 0.05, na.rm = TRUE)), color="darkred")+
      geom_vline(aes(xintercept=mean(dt_plot$phenotype, na.rm = TRUE)), color="navy", linetype = "dashed")+
      xlim(50,320)
    d1 <-  ggplot_build(plot1)$data[[1]]
    plot1 <- plot1 + geom_area(data = subset(d1, x < rand1), aes(x=x, y=y), fill="tomato") 
    rand1 <- quantile(dt_plot$phenotype, 0.05, na.rm = TRUE)
    plot1 <- ggplot(data = dt_plot, aes(x = phenotype))+
      geom_density(color=colour_short,fill =fill_short)+
      theme(text = element_text(size =12, family = "my"),
            panel.background = element_rect(fill = "white"),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "grey81"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "grey100"),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face ="bold"),
            axis.title.x= element_blank(),
            axis.text = element_text(size =12, family = "my", colour = "black"),
            axis.title.y = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
      labs(y = names_of_subpopulations[i])+
      ggtitle(year[i])+
      geom_vline(aes(xintercept=quantile(dt_plot$phenotype, 0.05, na.rm = TRUE)), color="darkred")+
      geom_vline(aes(xintercept=mean(dt_plot$phenotype, na.rm = TRUE)), color="navy", linetype = "dashed")+
      xlim(50,320)
    d1 <-  ggplot_build(plot1)$data[[1]]
    plot <- plot1 + geom_area(data = subset(d1, x < rand1), aes(x=x, y=y), fill="tomato")
    return(plot)
  }
  create_plot_low_pheno(dt_plot = ,
                        name_of_low_phenotype,
                        phenotype,
                        year)
}
names_of_subpopulations <- c("Short_plants_1","Short_plants_2","Tall_plants_1","Tall_plants_2")
years <- c("2016","2017","2018","2020")
str_match(names_of_subpopulations[1],names_of_subpopulations[2])

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

fill_short <- "darkseagreen2"
colour_short <- "darkseagreen3"
fill_tall <- "thistle2"
colour_tall <- "thistle3"

rand1 <- quantile(Short_1_2016$PlantHeight, 0.05, na.rm = TRUE)
plot1 <- ggplot(data = Short_1_2016, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face ="bold"),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"),
        axis.title.y = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
  labs(y = "Short plants 1")+
  ggtitle("2016")+
  geom_vline(aes(xintercept=quantile(Short_1_2016$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_1_2016$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d1 <-  ggplot_build(plot1)$data[[1]]
plot1 <- plot1 + geom_area(data = subset(d1, x < rand1), aes(x=x, y=y), fill="tomato") 

rand1 <- quantile(Short_1_2016$PlantHeight, 0.05, na.rm = TRUE)
plot1 <- ggplot(data = Short_1_2016, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face ="bold"),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"),
        axis.title.y = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
  labs(y = "Short plants 1")+
  ggtitle("2016")+
  geom_vline(aes(xintercept=quantile(Short_1_2016$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_1_2016$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d1 <-  ggplot_build(plot1)$data[[1]]
plot1 <- plot1 + geom_area(data = subset(d1, x < rand1), aes(x=x, y=y), fill="tomato") 
#### plot 2
rand2 <- quantile(Short_2_2016$PlantHeight, 0.05, na.rm = TRUE)
plot2 <- ggplot(data = Short_2_2016, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 12),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"),
        axis.title.y = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
  labs(y = "Short plants 2")+
  geom_vline(aes(xintercept=quantile(Short_2_2016$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_2_2016$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d2 <-  ggplot_build(plot2)$data[[1]]
plot2 <- plot2 + geom_area(data = subset(d2, x < rand2), aes(x=x, y=y), fill="tomato") 
#### plot 3
rand3 <- quantile(Tall_1_2016$PlantHeight, 0.95, na.rm = TRUE)
plot3 <- ggplot(data = Tall_1_2016, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 12),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"),
        axis.title.y = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
  labs(y = "Tall plants 1")+
  geom_vline(aes(xintercept=quantile(Tall_1_2016$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_1_2016$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d3 <- ggplot_build(plot3)$data[[1]]
plot3 <- plot3 + geom_area(data = subset(d3, x > rand3), aes(x=x, y=y), fill="tomato")
#### Plot 4
rand4 <- quantile(Tall_2_2016$PlantHeight, 0.95, na.rm = TRUE)
plot4 <- ggplot(data = Tall_1_2016, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.text = element_text(size =12, family = "my", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
  labs(y = "Tall plants 2", x = "Plant height (cm)")+
  geom_vline(aes(xintercept=quantile(Tall_2_2016$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_2_2016$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d4 <-  ggplot_build(plot4)$data[[1]]
plot4 <- plot4 + geom_area(data = subset(d4, x > rand4), aes(x=x, y=y), fill="tomato")

#### plot 5
rand5 <- quantile(Short_1_2017$PlantHeight, 0.05, na.rm = TRUE)
plot5 <- ggplot(data = Short_1_2017, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  ggtitle("2017")+
  geom_vline(aes(xintercept=quantile(Short_1_2017$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_1_2017$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d5 <-  ggplot_build(plot5)$data[[1]]
plot5 <- plot5 + geom_area(data = subset(d5, x < rand5), aes(x=x, y=y), fill="tomato") 
#### plot 6
rand6 <- quantile(Short_2_2017$PlantHeight, 0.05, na.rm = TRUE)
plot6 <- ggplot(data = Short_2_2017, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  geom_vline(aes(xintercept=quantile(Short_2_2017$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_2_2017$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d6 <-  ggplot_build(plot6)$data[[1]]
plot6 <- plot6 + geom_area(data = subset(d6, x < rand6), aes(x=x, y=y), fill="tomato") 
#### plot 7
rand7 <- quantile(Tall_1_2017$PlantHeight, 0.95, na.rm = TRUE)
plot7 <- ggplot(data = Tall_1_2017, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  geom_vline(aes(xintercept=quantile(Tall_1_2017$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_1_2017$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d7 <-  ggplot_build(plot7)$data[[1]]
plot7 <- plot7 + geom_area(data = subset(d7, x > rand7), aes(x=x, y=y), fill="tomato")
#### Plot 8
rand8 <- quantile(Tall_2_2017$PlantHeight, 0.95, na.rm = TRUE)
plot8 <- ggplot(data = Tall_1_2017, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x= element_text(size =12, family = "my", colour = "black", face = "bold"),
        axis.text = element_text(size =12, family = "my", colour = "black"),
        axis.title = element_text(hjust = 0.5, colour ="black", size = 12,face = "bold"))+
  labs(x = "Plant height (cm)")+
  geom_vline(aes(xintercept=quantile(Tall_2_2017$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_2_2017$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d8 <-  ggplot_build(plot8)$data[[1]]
plot8 <- plot8 + geom_area(data = subset(d8, x > rand8), aes(x=x, y=y), fill="tomato")
#### plot 9
rand9 <- quantile(Short_1_2018$PlantHeight, 0.05, na.rm = TRUE)
plot9 <- ggplot(data = Short_1_2018, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  ggtitle("2018")+
  geom_vline(aes(xintercept=quantile(Short_1_2018$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_1_2018$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d9 <-  ggplot_build(plot9)$data[[1]]
plot9 <- plot9 + geom_area(data = subset(d9, x < rand9), aes(x=x, y=y), fill="tomato") 
#### plot 10
rand10 <- quantile(Short_2_2018$PlantHeight, 0.05, na.rm = TRUE)
plot10 <- ggplot(data = Short_2_2018, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 12),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  geom_vline(aes(xintercept=quantile(Short_2_2018$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_2_2018$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d10 <-  ggplot_build(plot10)$data[[1]]
plot10 <- plot10 + geom_area(data = subset(d10, x < rand10), aes(x=x, y=y), fill="tomato") 
#### plot 11
rand11 <- quantile(Tall_1_2018$PlantHeight, 0.95, na.rm = TRUE)
plot11 <- ggplot(data = Tall_1_2018, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  geom_vline(aes(xintercept=quantile(Tall_1_2018$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_1_2018$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)+
  ylim(0,0.02)
d11 <-  ggplot_build(plot11)$data[[1]]
plot11 <- plot11 + geom_area(data = subset(d11, x > rand11), aes(x=x, y=y), fill="tomato")
#### Plot 12
rand12 <- quantile(Tall_2_2018$PlantHeight, 0.95, na.rm = TRUE)
plot12 <- ggplot(data = Tall_2_2018, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.y= element_blank(),
        axis.title.x= element_text(size =12, family = "my", colour = "black", face = "bold"),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  labs(x = "Plant height (cm)")+
  geom_vline(aes(xintercept=quantile(Tall_2_2018$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_2_2018$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d12 <-  ggplot_build(plot12)$data[[1]]
plot12 <- plot12 + geom_area(data = subset(d12, x > rand12), aes(x=x, y=y), fill="tomato")
#### plot 13
rand13 <- quantile(Short_1_2020$PlantHeight, 0.05, na.rm = TRUE)
plot13 <- ggplot(data = Short_1_2020, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 14, face = "bold"),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  ggtitle("2020")+
  geom_vline(aes(xintercept=quantile(Short_1_2020$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_1_2020$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d13 <-  ggplot_build(plot13)$data[[1]]
plot13 <- plot13 + geom_area(data = subset(d13, x < rand13), aes(x=x, y=y), fill="tomato") 
#### plot 14
rand14 <- quantile(Short_2_2020$PlantHeight, 0.05, na.rm = TRUE)
plot14 <- ggplot(data = Short_2_2020, aes(x = PlantHeight))+
  geom_density(color=colour_short,fill =fill_short)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 12),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  geom_vline(aes(xintercept=quantile(Short_2_2020$PlantHeight, 0.05, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Short_2_2020$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d14 <-  ggplot_build(plot14)$data[[1]]
plot14 <- plot14 + geom_area(data = subset(d14, x < rand14), aes(x=x, y=y), fill="tomato") 
#### plot 15
rand15 <- quantile(Tall_1_2020$PlantHeight, 0.95, na.rm = TRUE)
plot15 <- ggplot(data = Tall_1_2020, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, colour ="black", size = 12),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  geom_vline(aes(xintercept=quantile(Tall_1_2020$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_1_2020$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")+
  xlim(50,320)
d15 <-  ggplot_build(plot15)$data[[1]]
plot15 <- plot15 + geom_area(data = subset(d15, x > rand15), aes(x=x, y=y), fill="tomato")
#### Plot 16
rand16 <- quantile(Tall_2_2020$PlantHeight, 0.95, na.rm = TRUE)
plot16 <- ggplot(data = Tall_2_2020, aes(x = PlantHeight))+
  geom_density(color=colour_tall,fill =fill_tall)+
  theme(text = element_text(size =12, family = "my", colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.y= element_blank(),
        axis.title.x= element_text(size =12, family = "my", colour = "black", face = "bold"),
        axis.text = element_text(size =12, family = "my", colour = "black"))+
  labs(x = "Plant height (cm)")+
  xlim(50,320)+
  geom_vline(aes(xintercept=quantile(Tall_2_2020$PlantHeight, 0.95, na.rm = TRUE)), color="darkred")+
  geom_vline(aes(xintercept=mean(Tall_2_2020$PlantHeight, na.rm = TRUE)), color="navy", linetype = "dashed")
d16 <-  ggplot_build(plot16)$data[[1]]
plot16 <- plot16 + geom_area(data = subset(d16, x > rand16), aes(x=x, y=y), fill="tomato")
density_plot_2 <- gridExtra::grid.arrange(plot1,plot5,plot9,plot13,
                                          plot2,plot6,plot10,plot14,
                                          plot3,plot7,plot11,plot15,
                                          plot4,plot8,plot12,plot16,ncol=4,
                                          heights= c(1.1,1,1,1.1))
ggsave("C:/Users/mtost/Documents/Masterarbeit/Data_analysis/plots_for_paper/phenotypic_measurements_16_12_2021.png",density_plot_2,
       height = 6,
       width = 10,
       dpi = 1200)
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
t.test(Short_plants_2020$PlantHeight,Tall_plants_2020$PlantHeight,
       alternative = "less",
       paired = TRUE)

