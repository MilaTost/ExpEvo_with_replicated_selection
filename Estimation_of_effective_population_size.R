#### Calculate the effective population size -----------------------------------
library(stringr)
library(ggplot2)
library(data.table)
setwd("YOUR/OWN/PATH")
shoepag <- readxl::read_excel("field_data_collection_shoepag_12_09_2020.xlsx", sheet = 3, col_names = TRUE)

### Calculate the effetive population size -------------------------------------
calculate_eff_pop_size <- function(N_males,
                                   N_females){
  N_eff <- (4*N_males*N_females)/(N_males+N_females)
  return(N_eff)
}
calculate_eff_pop_size(N_males = 5000,
                       N_females = 250)
### Tassel flowering dates -----------------------------------------------------
shoepag <- as.data.table(shoepag)
shoepag$population <- str_sub(shoepag$Barcode_of_plant,1,1)
shoepag_1 <- shoepag[population == "1",]
shoepag_2 <- shoepag[population == "2",]
shoepag_3 <- shoepag[population == "3",]
shoepag_4 <- shoepag[population == "4",]
flowering_data_set_T <- function(data){
  day_T <- str_sub(data$Flowering_date_Tassel,1,2)
  day_Tassel <- day_T[!is.na(as.numeric(day_T))]
  index <- which(!is.na(as.numeric(day_T)))
  data <- data[index,]
  month_Tassel <- str_sub(data$Flowering_date_Tassel,4,5)
  date_T <- paste("2020",month_Tassel,day_Tassel, sep = "-", collapse = NULL, recycle0 = F)

  day_S <- str_sub(data$Flowering_date_Silk,1,2)
  day_Silk <- day_S[!is.na(as.numeric(day_S))]
  index <- which(!is.na(as.numeric(day_S)))
  month_S <- str_sub(data$Flowering_date_Silk,4,5)
  month_Silk <- month_S[index]
  date_S <- paste("2020",month_Silk,day_Silk, sep = "-", collapse = NULL, recycle0 = F)

  Flowering_date_T <- as.Date(date_T)
  Flowering_date_T <- as.data.frame(Flowering_date_T)
  colnames(Flowering_date_T) <- "Flowering_date_Tassel"
  Flowering_date_S <- as.Date(date_S)
  Flowering_date_S <- as.data.frame(Flowering_date_S)
  colnames(Flowering_date_S) <- "Flowering_date_Silk"
  return(Flowering_date_T)
}

flowering_data_set_S <- function(data){
  day_T <- str_sub(data$Flowering_date_Tassel,1,2)
  day_Tassel <- day_T[!is.na(as.numeric(day_T))]
  index <- which(!is.na(as.numeric(day_T)))
  data <- data[index,]
  month_Tassel <- str_sub(data$Flowering_date_Tassel,4,5)
  date_T <- paste("2020",month_Tassel,day_Tassel, sep = "-", collapse = NULL, recycle0 = F)

  day_S <- str_sub(data$Flowering_date_Silk,1,2)
  day_Silk <- day_S[!is.na(as.numeric(day_S))]
  index <- which(!is.na(as.numeric(day_S)))
  month_S <- str_sub(data$Flowering_date_Silk,4,5)
  month_Silk <- month_S[index]
  date_S <- paste("2020",month_Silk,day_Silk, sep = "-", collapse = NULL, recycle0 = F)

  Flowering_date_T <- as.Date(date_T)
  Flowering_date_T <- as.data.frame(Flowering_date_T)
  colnames(Flowering_date_T) <- "Flowering_date_Tassel"
  Flowering_date_S <- as.Date(date_S)
  Flowering_date_S <- as.data.frame(Flowering_date_S)
  colnames(Flowering_date_S) <- "Flowering_date_Silk"
  return(Flowering_date_S)
}

flowering_time_intervals <- function(data){
  day_T <- str_sub(data$Flowering_date_Tassel,1,2)
  day_Tassel <- day_T[!is.na(as.numeric(day_T))]
  index <- which(!is.na(as.numeric(day_T)))
  data <- data[index,]
  month_Tassel <- str_sub(data$Flowering_date_Tassel,4,5)
  date_T <- paste("2020",month_Tassel,day_Tassel, sep = "-", collapse = NULL, recycle0 = F)

  day_S <- str_sub(data$Flowering_date_Silk,1,2)
  day_Silk <- day_S[!is.na(as.numeric(day_S))]
  index <- which(!is.na(as.numeric(day_S)))
  month_S <- str_sub(data$Flowering_date_Silk,4,5)
  month_Silk <- month_S[index]
  date_S <- paste("2020",month_Silk,day_Silk, sep = "-", collapse = NULL, recycle0 = F)

  Flowering_date_T <- as.Date(date_T)
  Flowering_date_T <- as.data.frame(Flowering_date_T)
  colnames(Flowering_date_T) <- "Flowering_date_Tassel"
  Flowering_date_S <- as.Date(date_S)
  Flowering_date_S <- as.data.frame(Flowering_date_S)
  colnames(Flowering_date_S) <- "Flowering_date_Silk"

  days_to_flowering_Tassel <- Flowering_date_T$Flowering_date_Tassel-as.Date("2020-04-28")
  as.Date("2020-04-28") + median(days_to_flowering_Tassel)
  mat <- matrix(nrow = 9, ncol= 4)
  mat[,1:2] <- cbind(seq(0.05,0.45,0.05),sort(seq(0.55,0.95,0.05), decreasing = TRUE))
  for (i in 1:9) {
    mat[i,3:4] <- quantile(days_to_flowering_Tassel, probs = c(mat[i,1],mat[i,2]), na.rm = TRUE)
  }
  mat <- as.data.frame(mat)
  mat$time_int <- mat$V4 -mat$V3
  mat$date_1 <- as.Date("2020-04-28") + mat$V3
  mat$date_2 <- as.Date("2020-04-28") + mat$V4
  return(mat)
}
time_intervals_pop1 <- flowering_time_intervals(shoepag_1)# 0.2;0.8
time_intervals_pop2 <- flowering_time_intervals(shoepag_2)# 0.35;0.65
time_intervals_pop3 <- flowering_time_intervals(shoepag_3)# 0.25;0.75
time_intervals_pop4 <- flowering_time_intervals(shoepag_4)# 0.85;0.15

cbind(time_intervals_pop1,time_intervals_pop2,time_intervals_pop3,time_intervals_pop4)

Flowering_date_T_pop1 <- flowering_data_set_T(shoepag_1)
Flowering_date_T_pop2 <- flowering_data_set_T(shoepag_2)
Flowering_date_T_pop3 <- flowering_data_set_T(shoepag_3)
Flowering_date_T_pop4 <- flowering_data_set_T(shoepag_4)

Flowering_date_S_pop1 <- flowering_data_set_S(shoepag_1)
Flowering_date_S_pop2 <- flowering_data_set_S(shoepag_2)
Flowering_date_S_pop3 <- flowering_data_set_S(shoepag_3)
Flowering_date_S_pop4 <- flowering_data_set_S(shoepag_4)

days_to_flowering_Tassel_pop1 <- Flowering_date_T_pop1$Flowering_date_Tassel-as.Date("2020-04-28")
days_to_flowering_Tassel_pop2 <- Flowering_date_T_pop2$Flowering_date_Tassel-as.Date("2020-04-28")
days_to_flowering_Tassel_pop3 <- Flowering_date_T_pop3$Flowering_date_Tassel-as.Date("2020-04-28")
days_to_flowering_Tassel_pop4 <- Flowering_date_T_pop4$Flowering_date_Tassel-as.Date("2020-04-28")

silk_0_pop1 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop1, probs = 0.2, na.rm = TRUE)
silk_1_pop1 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop1, probs = 0.8, na.rm = TRUE)
silk_0_pop2 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop2, probs = 0.35, na.rm = TRUE)
silk_1_pop2 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop2, probs = 0.65, na.rm = TRUE)
silk_0_pop3 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop3, probs = 0.25, na.rm = TRUE)
silk_1_pop3 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop3, probs = 0.75, na.rm = TRUE)
silk_0_pop4 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop4, probs = 0.15, na.rm = TRUE)
silk_1_pop4 <- as.Date("2020-04-28") + quantile(days_to_flowering_Tassel_pop4, probs = 0.85, na.rm = TRUE)

silk_0_num_pop1 <- length(which(silk_0_pop1 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop1)))
silk_1_num_pop1 <- length(which(silk_1_pop1 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop1)))
num_pop_1 <- silk_0_num_pop1 - silk_1_num_pop1
silk_0_num_pop2 <- length(which(silk_0_pop2 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop2)))
silk_1_num_pop2 <- length(which(silk_1_pop2 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop2)))
num_pop_2 <- silk_0_num_pop2 - silk_1_num_pop2
silk_0_num_pop3 <- length(which(silk_0_pop3 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop3)))
silk_1_num_pop3 <- length(which(silk_1_pop3 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop3)))
num_pop_3 <- silk_0_num_pop3 - silk_1_num_pop3
silk_0_num_pop4 <- length(which(silk_0_pop4 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop4)))
silk_1_num_pop4 <- length(which(silk_1_pop4 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop4)))
num_pop_4 <- silk_0_num_pop4 - silk_1_num_pop4

silk_0_num_pop1 <- length(which(silk_0_pop1 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop1)))
silk_1_num_pop1 <- length(which(silk_1_pop1 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop1)))
num_pop_1 <- silk_0_num_pop1 - silk_1_num_pop1
silk_0_num_pop2 <- length(which(silk_0_pop2 < (as.Date("2020-04-28") + days_to_flowering_Silk_pop2)))
silk_1_num_pop2 <- length(which(silk_1_pop2 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop2)))
num_pop_2 <- silk_0_num_pop2 - silk_1_num_pop2
silk_0_num_pop3 <- length(which(silk_0_pop3 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop3)))
silk_1_num_pop3 <- length(which(silk_1_pop3 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop3)))
num_pop_3 <- silk_0_num_pop3 - silk_1_num_pop3
silk_0_num_pop4 <- length(which(silk_0_pop4 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop4)))
silk_1_num_pop4 <- length(which(silk_1_pop4 < (as.Date("2020-04-28") + days_to_flowering_Tassel_pop4)))
num_pop_4 <- silk_0_num_pop4 - silk_1_num_pop4
### Project the numbers of simulatanously flowering plants onto the entire population
sim_flow_males_pop_1 <- (num_pop_1/96)*5000
sim_flow_males_pop_2 <- (num_pop_2/96)*5000
sim_flow_males_pop_3 <- (num_pop_3/96)*5000
sim_flow_males_pop_4 <- (num_pop_4/96)*5000
### New effective population size considering simultanously flowering plants ---
calculate_eff_pop_size(N_males = sim_flow_males_pop_1,
                       N_females = 250)
calculate_eff_pop_size(N_males = sim_flow_males_pop_2,
                       N_females = 250)
calculate_eff_pop_size(N_males = sim_flow_males_pop_3,
                       N_females = 250)
calculate_eff_pop_size(N_males = sim_flow_males_pop_4,
                       N_females = 250)
