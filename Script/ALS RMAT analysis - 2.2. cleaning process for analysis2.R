# Corestem project - ALS AE SAE Analysis.proj
# 20201113 Sehwan Chun at Corestem, Inc.
# 2.1. cleaning process for ALS AE SAE analysis

#### 1. source Loading ####
load("./Data/ALS RMAT analysis files.image")
source("./Script/Functions/ALS RMAT analysis - functions.R")

#### 2. info check ####
nrow(subset(ALSFRSR_file_allo_P12, GROUP == 0)) #27
nrow(subset(ALSFRSR_file_allo_P12, GROUP == 1)) #45

nrow(subset(ALSFRSR_file_allo_P2cont, GROUP == 0)) #27
nrow(subset(ALSFRSR_file_allo_P2cont, GROUP == 1)) #6

ALSFRSR_file_allo_P12$Changes = ALSFRSR_file_allo_P12$M12 - ALSFRSR_file_allo_P12$M0
ALSFRSR_file_allo_P2cont$Changes = ALSFRSR_file_allo_P2cont$M12 - ALSFRSR_file_allo_P2cont$M0
ALSFRSR_file_allo_P12_FCSPMM$Changes = ALSFRSR_file_allo_P12_FCSPMM$M12 - ALSFRSR_file_allo_P12_FCSPMM$M0

PMS_ALSFRS_table = PMS_ALSFRS_table[PMS_ALSFRS_table$SUBJ_ID %in% PMS_ALSFRS_init_table$SUBJ_ID,]
#PMS_ALSFRS_init_table = PMS_ALSFRS_init_table[PMS_ALSFRS_init_table$SUBJ_ID %in% PMS_ALSFRS_table$SUBJ_ID,]

PMS_ALSFRS_table = PMS_ALSFRS_table[order(PMS_ALSFRS_table$SUBJ_ID),]
PMS_ALSFRS_init_table = PMS_ALSFRS_init_table[order(PMS_ALSFRS_init_table$SUBJ_ID),]

PMS_ALSFRS_init_table$Date_1st2nd = (PMS_ALSFRS_init_table$Date_2nd - PMS_ALSFRS_init_table$Date_1st) / 28
PMS_ALSFRS_init_table$Date_1st3rd = (PMS_ALSFRS_init_table$Date_3rd - PMS_ALSFRS_init_table$Date_1st) / 28
PMS_ALSFRS_init_table$Date_1st4th = (PMS_ALSFRS_init_table$Date_4th - PMS_ALSFRS_init_table$Date_1st) / 28
PMS_ALSFRS_init_table$Date_1st5th = (PMS_ALSFRS_init_table$Date_5th - PMS_ALSFRS_init_table$Date_1st) / 28
PMS_ALSFRS_init_table$Date_1st6th = (PMS_ALSFRS_init_table$Date_6th - PMS_ALSFRS_init_table$Date_1st) / 28

PMS_ALSFRS_init_table[,"init_slope"] = as.numeric(1)
PMS_ALSFRS_init_table[,"M6_slope"] = as.numeric(1)
PMS_ALSFRS_init_table[,"M9_slope"] = as.numeric(1)
PMS_ALSFRS_init_table[,"M12_slope"] = as.numeric(1)
PMS_ALSFRS_init_table[,"M6_month"] = as.numeric(1)
PMS_ALSFRS_init_table[,"M9_month"] = as.numeric(1)
PMS_ALSFRS_init_table[,"M12_month"] = as.numeric(1)
PMS_ALSFRS_init_table$tx_Date_1st2nd = (PMS_ALSFRS_init_table$Tx_2nd - PMS_ALSFRS_init_table$Tx_1st) / 28
PMS_ALSFRS_init_table$tx_Date_1st3rd = (PMS_ALSFRS_init_table$Tx_3rd - PMS_ALSFRS_init_table$Tx_1st) / 28

for (i in 1:nrow(PMS_ALSFRS_init_table)){
    PMS_ALSFRS_init_table[i,"init_slope"] =  (PMS_ALSFRS_init_table[i,"ALSFRS_1st"] - PMS_ALSFRS_init_table[i,"ALSFRS_initial"]) / 
        ((PMS_ALSFRS_init_table[i,"Date_1st"] - PMS_ALSFRS_init_table[i,"Date_initial"]) / 28)
    
    tmp_set = PMS_ALSFRS_init_table[i, 22:26]
    if(sum(is.na(tmp_set)) != 5){
        tmp_set_min_6m = which.min(abs(PMS_ALSFRS_init_table[i, 22:26] - 6))
        tmp_set_min_value_6m = tmp_set[tmp_set_min_6m]
        tmp_set_min_9m = which.min(abs(PMS_ALSFRS_init_table[i, 22:26] - 9))
        tmp_set_min_value_9m = tmp_set[tmp_set_min_9m]
        tmp_set_min_12m = which.min(abs(PMS_ALSFRS_init_table[i, 22:26] - 12.5))
        tmp_set_min_value_12m = tmp_set[tmp_set_min_12m]
    }
    
    if(abs(tmp_set_min_value_6m - 6) > threshold){
        PMS_ALSFRS_init_table[i,"M6_slope"] = NA
        PMS_ALSFRS_init_table[i,"M6_month"] = NA
    }else{
        PMS_ALSFRS_init_table[i,"M6_slope"] = (PMS_ALSFRS_init_table[i, 10+(2*tmp_set_min_6m)] -
                                           PMS_ALSFRS_init_table[i,10]) / tmp_set_min_value_6m
        PMS_ALSFRS_init_table[i,"M6_month"] = tmp_set_min_value_6m
    }
    if(abs(tmp_set_min_value_9m - 9) > threshold){
        PMS_ALSFRS_init_table[i,"M9_slope"] = NA
        PMS_ALSFRS_init_table[i,"M9_month"] = NA
    }else{
        PMS_ALSFRS_init_table[i,"M9_slope"] = (PMS_ALSFRS_init_table[i, 10+(2*tmp_set_min_9m)] -
                                           PMS_ALSFRS_init_table[i,10]) / tmp_set_min_value_9m
        PMS_ALSFRS_init_table[i,"M9_month"] = tmp_set_min_value_9m
    }
    if(abs(tmp_set_min_value_12m - 12) > threshold){
        PMS_ALSFRS_init_table[i,"M12_slope"] = NA
        PMS_ALSFRS_init_table[i,"M12_month"] = NA
    }else{
        PMS_ALSFRS_init_table[i,"M12_slope"] = (PMS_ALSFRS_init_table[i, 10+(2*tmp_set_min_12m)] -
                                           PMS_ALSFRS_init_table[i,10]) / tmp_set_min_value_12m
        PMS_ALSFRS_init_table[i,"M12_month"] = tmp_set_min_value_12m
        
    }
    
    PMS_ALSFRS_init_table[i,"Change_of_slope_6M"] = (PMS_ALSFRS_init_table[i,"init_slope"] -
                                                         PMS_ALSFRS_init_table[i,"M6_slope"])
    PMS_ALSFRS_init_table[i,"Change_of_slope_9M"] = (PMS_ALSFRS_init_table[i,"init_slope"] -
                                                         PMS_ALSFRS_init_table[i,"M9_slope"])
    PMS_ALSFRS_init_table[i,"Change_of_slope_12M"] = (PMS_ALSFRS_init_table[i,"init_slope"] -
                                                         PMS_ALSFRS_init_table[i,"M12_slope"])
    
    if(is.na(PMS_ALSFRS_init_table[i,"tx_Date_1st3rd"])){
        PMS_ALSFRS_init_table[i,"M6_Multi"] = "Normal"
        PMS_ALSFRS_init_table[i,"M9_Multi"] = "Normal"
        PMS_ALSFRS_init_table[i,"M12_Multi"] = "Normal"
    }else{
        if(is.na(PMS_ALSFRS_init_table[i,"M6_month"])){
            PMS_ALSFRS_init_table[i,"M6_Multi"] = NA
        }else{
            if(PMS_ALSFRS_init_table[i,"tx_Date_1st3rd"] > PMS_ALSFRS_init_table[i,"M6_month"]){
                PMS_ALSFRS_init_table[i,"M6_Multi"] = "Normal"
            }else{
                PMS_ALSFRS_init_table[i,"M6_Multi"] = "Multi"
            }
        }
        if(is.na(PMS_ALSFRS_init_table[i,"M9_month"])){
            PMS_ALSFRS_init_table[i,"M9_Multi"] = NA
        }else{
            if(PMS_ALSFRS_init_table[i,"tx_Date_1st3rd"] > PMS_ALSFRS_init_table[i,"M9_month"]){
                PMS_ALSFRS_init_table[i,"M9_Multi"] = "Normal"
            }else{
                PMS_ALSFRS_init_table[i,"M9_Multi"] = "Multi"
            }
        }
        if(is.na(PMS_ALSFRS_init_table[i,"M12_month"])){
            PMS_ALSFRS_init_table[i,"M12_Multi"] = NA
        }else{
            if(PMS_ALSFRS_init_table[i,"tx_Date_1st3rd"] > PMS_ALSFRS_init_table[i,"M12_month"]){
                PMS_ALSFRS_init_table[i,"M12_Multi"] = "Normal"
            }else{
                PMS_ALSFRS_init_table[i,"M12_Multi"] = "Multi"
            }
        }
    }
}




t.test(subset(PMS_ALSFRS_init_table, INJ >= 3)$Change_of_slope_12M,
            subset(PMS_ALSFRS_init_table, INJ < 3)$Change_of_slope_12M,
       var.equal = T)

PMS_ALSFRS_init_table = PMS_ALSFRS_init_table

ggp = ggboxplot(data = PMS_ALSFRS_init_table, x = "M6_Multi", y = "Change_of_slope_6M", color = "M6_Multi", add = "jitter",
          xlab = "Group", ylab = "Change of ALSFRS-R slope", add.params = list(size = 2.5)) + stat_compare_means(method = "t.test", 
                                                                                                                 method.args = list(var.equal = TRUE))
ggp

ggsave("./Output/PMS change of ALSFRS-R slope.tiff", width = 6, height = 6, units = "in", dpi = 300)
save.image(file = "./Data/ALS RMAT analysis files.image")






#

if(sum(is.na(PMS_ALSFRS_table[i,14:16])) == 3){
    PMS_ALSFRS_init_table[i,"M12_slope"] = NA
}else{
    PMS_ALSFRS_init_table[i,"M12_slope"] = as.numeric((max(PMS_ALSFRS_table[i,14:16], na.rm = T) - PMS_ALSFRS_table[i,3]) / 12)
}