#Corestem project - ALS AE SAE Analysis.proj
# 20201113 Sehwan Chun at Corestem, Inc.
# function.R

#### 1. Library Loading ####
packs = c("data.table", "readxl", "ggpubr", "writexl", "car")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. simple LOCF and stat making ####
LOCF_table_run = function(ALS_table){
        ALS_table_LOCF = ALS_table
        for (i in 1:nrow(ALS_table_LOCF)){
                if(is.na(ALS_table_LOCF[i,15])){
                       tmp = na.omit(as.numeric(ALS_table_LOCF[i,]))
                       ALS_table_LOCF[i,15] = tmp[length(tmp)]
                }
        }
        
        ALS_table_LOCF$Changes = ALS_table_LOCF$M12 - ALS_table_LOCF$M0 
 return(ALS_table_LOCF)
}
stat_table_run = function(ALS_table){
        tmp_table  = data.frame(test = NA, pvalue = NA)
        tmp_table[1:4,1] = c("wilcox.test",
                              "var.test",
                              "t.test",
                              "ANCOVA")
        
        ALS_table_cont = subset(ALS_table, GROUP == 0)
        ALS_table_treat = subset(ALS_table, GROUP == 1)
        
        tmp_table[1,2] = wilcox.test(ALS_table_cont$Changes, ALS_table_treat$Changes)$p.value
        tmp_table[2,2] = var.test(ALS_table_cont$Changes, ALS_table_treat$Changes)$p.value

        if(tmp_table[2,2] < 0.05){
                tmp_table[3,2] = t.test(ALS_table_cont$Changes,
                                        ALS_table_treat$Changes,
                                        var.equal = FALSE)$p.value
        }else{
                tmp_table[3,2] = t.test(ALS_table_cont$Changes,
                                        ALS_table_treat$Changes,
                                        var.equal = TRUE)$p.value
        }
        
        tmp_table[4,2] = summary(lm(data = ALS_table, Changes ~ GROUP + M0))$coefficients[2,4]
        
        
        
        return(tmp_table)
}
PMS_slope_comparison_run = function(PMS_data_table, PMS_data_source, threshold = 1.5){
        if(PMS_data_source == "DreamCIS"){
                PMS_ALSFRS_table = PMS_ALSFRS_table[PMS_ALSFRS_table$SUBJ_ID %in% PMS_data_table$SUBJ_ID,]
                PMS_data_table = PMS_data_table[PMS_data_table$SUBJ_ID %in% PMS_ALSFRS_table$SUBJ_ID,]
                PMS_ALSFRS_table = PMS_ALSFRS_table[order(PMS_ALSFRS_table$SUBJ_ID),]
                PMS_data_table = PMS_data_table[order(PMS_data_table$SUBJ_ID),]
                
                PMS_data_table[1,22:26] = as.numeric(1)
                
                for (i in 1:nrow(PMS_data_table)){
                        PMS_data_table[i,22] =  (PMS_data_table[i,10] - PMS_data_table[i,8]) / 
                                ((PMS_data_table[i,11] - PMS_data_table[i,9]) / 28)
                        
                        if(sum(is.na(PMS_ALSFRS_table[i,8:10])) == 3){
                                PMS_data_table[i,23] = NA
                        }else{
                                PMS_data_table[i,23] = as.numeric((max(PMS_ALSFRS_table[i,8:10], na.rm = T) - PMS_ALSFRS_table[i,3]) / 6)
                        }
                        
                        if(sum(is.na(PMS_ALSFRS_table[i,11:13])) == 3){
                                PMS_data_table[i,24] = NA
                        }else{
                                PMS_data_table[i,24] = as.numeric((max(PMS_ALSFRS_table[i,11:13], na.rm = T) - PMS_ALSFRS_table[i,3]) / 9)
                        }
                        
                        if(sum(is.na(PMS_ALSFRS_table[i,14:16])) == 3){
                                PMS_data_table[i,25] = NA
                        }else{
                                PMS_data_table[i,25] = as.numeric((max(PMS_ALSFRS_table[i,14:16], na.rm = T) - PMS_ALSFRS_table[i,3]) / 12)
                        }
                        
                        
                        PMS_data_table[i,26] = (PMS_data_table[i,22] - PMS_data_table[i,23])
                        PMS_data_table[i,27] = (PMS_data_table[i,22] - PMS_data_table[i,24])
                        PMS_data_table[i,28] = (PMS_data_table[i,22] - PMS_data_table[i,25])
                        
                }
                
                colnames(PMS_data_table)[22:28] = c("init_slope",
                                                    "M6_slope",
                                                    "M9_slope",
                                                    "M12_slope",
                                                    "Change_of_slope_6M",
                                                    "Change_of_slope_9M",
                                                    "Change_of_slope_12M")
                
                PMS_data_table$Date_1st2nd = (PMS_data_table$Date_2nd - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st3rd = (PMS_data_table$Date_3rd - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st4th = (PMS_data_table$Date_4th - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st5th = (PMS_data_table$Date_5th - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st6th = (PMS_data_table$Date_6th - PMS_data_table$Date_1st) / 28
                PMS_data_table$tx_Date_1st2nd = (PMS_data_table$Tx_2nd - PMS_data_table$Tx_1st) / 28
                PMS_data_table$tx_Date_1st3rd = (PMS_data_table$Tx_3rd - PMS_data_table$Tx_1st) / 28
                
                PMS_data_table[,"M6_month"] = as.numeric(1)
                PMS_data_table[,"M9_month"] = as.numeric(1)
                PMS_data_table[,"M12_month"] = as.numeric(1)
                
                for (i in 1:nrow(PMS_data_table)){
                        
                        tmp_set = PMS_data_table[i, 29:33]
                        if(sum(is.na(tmp_set)) == 5){
                                PMS_data_table[i,"M6_month"] = NA
                                PMS_data_table[i,"M9_month"] = NA
                                PMS_data_table[i,"M12_month"] = NA
                                
                        }else if(sum(is.na(tmp_set)) != 5){
                                tmp_set_min_6m = which.min(abs(PMS_data_table[i, 29:33] - 6.5))
                                tmp_set_min_value_6m = tmp_set[tmp_set_min_6m]
                                tmp_set_min_9m = which.min(abs(PMS_data_table[i, 29:33] - 9.5))
                                tmp_set_min_value_9m = tmp_set[tmp_set_min_9m]
                                tmp_set_min_12m = which.min(abs(PMS_data_table[i, 29:33] - 12.5))
                                tmp_set_min_value_12m = tmp_set[tmp_set_min_12m]
                        }
                        
                        if(abs(tmp_set_min_value_6m - 6) >= threshold){
                                PMS_data_table[i,"M6_month"] = NA
                        }else{
                                PMS_data_table[i,"M6_month"] = tmp_set_min_value_6m
                        }
                        if(abs(tmp_set_min_value_9m - 9) >= threshold){
                                PMS_data_table[i,"M9_month"] = NA
                        }else{
                                PMS_data_table[i,"M9_month"] = tmp_set_min_value_9m
                        }
                        if(abs(tmp_set_min_value_12m - 12) >= threshold){
                                PMS_data_table[i,"M12_month"] = NA
                        }else{
                                PMS_data_table[i,"M12_month"] = tmp_set_min_value_12m
                        }
                        
                        if(is.na(PMS_data_table[i,"tx_Date_1st3rd"])){
                                PMS_data_table[i,"M6_Multi"] = "Normal"
                                PMS_data_table[i,"M9_Multi"] = "Normal"
                                PMS_data_table[i,"M12_Multi"] = "Normal"
                        }else{
                                if(is.na(PMS_data_table[i,"M6_month"])){
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > 6.5){
                                                PMS_data_table[i,"M6_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M6_Multi"] = "Multi"
                                        }
                                }else{
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > PMS_data_table[i,"M6_month"]){
                                                PMS_data_table[i,"M6_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M6_Multi"] = "Multi"
                                        }
                                }
                                if(is.na(PMS_data_table[i,"M9_month"])){
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > 9.5){
                                                PMS_data_table[i,"M9_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M9_Multi"] = "Multi"
                                        }
                                }else{
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > PMS_data_table[i,"M9_month"]){
                                                PMS_data_table[i,"M9_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M9_Multi"] = "Multi"
                                        }
                                }
                                if(is.na(PMS_data_table[i,"M12_month"])){
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > 12.5){
                                                PMS_data_table[i,"M12_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M12_Multi"] = "Multi"
                                        }
                                }else{
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > PMS_data_table[i,"M12_month"]){
                                                PMS_data_table[i,"M12_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M12_Multi"] = "Multi"
                                        }
                                }
                        }
                }
                
                return(PMS_data_table)
                
        }else if(PMS_data_source == "HYCRC"){
                PMS_data_table = PMS_data_table[order(PMS_data_table$SUBJ_ID),]
                
                PMS_data_table$Date_1st2nd = (PMS_data_table$Date_2nd - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st3rd = (PMS_data_table$Date_3rd - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st4th = (PMS_data_table$Date_4th - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st5th = (PMS_data_table$Date_5th - PMS_data_table$Date_1st) / 28
                PMS_data_table$Date_1st6th = (PMS_data_table$Date_6th - PMS_data_table$Date_1st) / 28
                
                PMS_data_table[,"init_slope"] = as.numeric(1)
                PMS_data_table[,"M6_slope"] = as.numeric(1)
                PMS_data_table[,"M9_slope"] = as.numeric(1)
                PMS_data_table[,"M12_slope"] = as.numeric(1)
                PMS_data_table[,"M6_month"] = as.numeric(1)
                PMS_data_table[,"M9_month"] = as.numeric(1)
                PMS_data_table[,"M12_month"] = as.numeric(1)
                
                PMS_data_table$tx_Date_1st2nd = (PMS_data_table$Tx_2nd - PMS_data_table$Tx_1st) / 28
                PMS_data_table$tx_Date_1st3rd = (PMS_data_table$Tx_3rd - PMS_data_table$Tx_1st) / 28
                
                for (i in 1:nrow(PMS_data_table)){
                        PMS_data_table[i,"init_slope"] =  (PMS_data_table[i,"ALSFRS_1st"] - PMS_data_table[i,"ALSFRS_initial"]) / 
                                ((PMS_data_table[i,"Date_1st"] - PMS_data_table[i,"Date_initial"]) / 28)
                        
                        tmp_set = PMS_data_table[i, 22:26]
                        
                        if(sum(is.na(tmp_set)) == 5){
                                PMS_data_table[i,"M6_slope"] = NA
                                PMS_data_table[i,"M6_month"] = NA
                                PMS_data_table[i,"M9_slope"] = NA
                                PMS_data_table[i,"M9_month"] = NA
                                PMS_data_table[i,"M12_slope"] = NA
                                PMS_data_table[i,"M12_month"] = NA
                                PMS_data_table[i,"Change_of_slope_6M"] = NA
                                PMS_data_table[i,"Change_of_slope_9M"] = NA
                                PMS_data_table[i,"Change_of_slope_12M"] = NA
                                
                        }else if(sum(is.na(tmp_set)) != 5){
                                tmp_set_min_6m = which.min(abs(PMS_data_table[i, 22:26] - 6))
                                tmp_set_min_value_6m = tmp_set[tmp_set_min_6m]
                                tmp_set_min_9m = which.min(abs(PMS_data_table[i, 22:26] - 9))
                                tmp_set_min_value_9m = tmp_set[tmp_set_min_9m]
                                tmp_set_min_12m = which.min(abs(PMS_data_table[i, 22:26] - 12))
                                tmp_set_min_value_12m = tmp_set[tmp_set_min_12m]
                                
                                if(abs(tmp_set_min_value_6m - 6) >= threshold){
                                        PMS_data_table[i,"M6_slope"] = NA
                                        PMS_data_table[i,"M6_month"] = NA
                                }else{
                                        PMS_data_table[i,"M6_slope"] = (PMS_data_table[i, 10+(2*tmp_set_min_6m)] -
                                                                                PMS_data_table[i,10]) / tmp_set_min_value_6m
                                        PMS_data_table[i,"M6_month"] = tmp_set_min_value_6m
                                }
                                if(abs(tmp_set_min_value_9m - 9) >= threshold){
                                        PMS_data_table[i,"M9_slope"] = NA
                                        PMS_data_table[i,"M9_month"] = NA
                                }else{
                                        PMS_data_table[i,"M9_slope"] = (PMS_data_table[i, 10+(2*tmp_set_min_9m)] -
                                                                                PMS_data_table[i,10]) / tmp_set_min_value_9m
                                        PMS_data_table[i,"M9_month"] = tmp_set_min_value_9m
                                }
                                if(abs(tmp_set_min_value_12m - 12) >= threshold){
                                        PMS_data_table[i,"M12_slope"] = NA
                                        PMS_data_table[i,"M12_month"] = NA
                                }else{
                                        PMS_data_table[i,"M12_slope"] = (PMS_data_table[i, 10+(2*tmp_set_min_12m)] -
                                                                                 PMS_data_table[i,10]) / tmp_set_min_value_12m
                                        PMS_data_table[i,"M12_month"] = tmp_set_min_value_12m
                                        
                                }
                                
                                PMS_data_table[i,"Change_of_slope_6M"] = (PMS_data_table[i,"init_slope"] -
                                                                                  PMS_data_table[i,"M6_slope"])
                                PMS_data_table[i,"Change_of_slope_9M"] = (PMS_data_table[i,"init_slope"] -
                                                                                  PMS_data_table[i,"M9_slope"])
                                PMS_data_table[i,"Change_of_slope_12M"] = (PMS_data_table[i,"init_slope"] -
                                                                                   PMS_data_table[i,"M12_slope"])
                        }
                }
                
                for (i in 1:nrow(PMS_data_table)){
                        if(is.na(PMS_data_table[i,"tx_Date_1st3rd"])){
                                PMS_data_table[i,"M6_Multi"] = "Normal"
                                PMS_data_table[i,"M9_Multi"] = "Normal"
                                PMS_data_table[i,"M12_Multi"] = "Normal"
                        }else{
                                if(is.na(PMS_data_table[i,"M6_month"])){
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > 6.5){
                                                PMS_data_table[i,"M6_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M6_Multi"] = "Multi"
                                        }
                                }else{
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > PMS_data_table[i,"M6_month"]){
                                                PMS_data_table[i,"M6_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M6_Multi"] = "Multi"
                                        }
                                }
                                if(is.na(PMS_data_table[i,"M9_month"])){
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > 9.5){
                                                PMS_data_table[i,"M9_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M9_Multi"] = "Multi"
                                        }
                                }else{
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > PMS_data_table[i,"M9_month"]){
                                                PMS_data_table[i,"M9_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M9_Multi"] = "Multi"
                                        }
                                }
                                if(is.na(PMS_data_table[i,"M12_month"])){
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > 12.5){
                                                PMS_data_table[i,"M12_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M12_Multi"] = "Multi"
                                        }
                                }else{
                                        if(PMS_data_table[i,"tx_Date_1st3rd"] > PMS_data_table[i,"M12_month"]){
                                                PMS_data_table[i,"M12_Multi"] = "Normal"
                                        }else{
                                                PMS_data_table[i,"M12_Multi"] = "Multi"
                                        }
                                }
                        }
                }
                
                return(PMS_data_table)
        }
        
}
PMS_slope_comparison_plot_run = function(PMS_extended_data_table, month){
        
        if(PMS_extended_data_table == 1){
                PMS_extended_data_table = PMS_data_table_DreamCIS
                table_name = "PMS_data_table_DreamCIS"
        }else{
                PMS_extended_data_table = PMS_data_table_HYCRC
                table_name = "PMS_data_table_HYCRC"
        }
        
        if(month == 6){
                tmp_6m_table = subset(PMS_extended_data_table, is.na(Change_of_slope_6M) == F)
                tmp_6m_table = subset(tmp_6m_table, ALSFRS_1st >=31 & ALSFRS_1st <= 46)
                tmp_6m_table$GROUP = ifelse(tmp_6m_table$INJ >= 3,
                                            "Multi",
                                            "Normal")
                
                tmp_6m_table$GROUP = factor(tmp_6m_table$GROUP, levels = c("Normal","Multi"))
                tmp_6m_table$M6_Multi = factor(tmp_6m_table$M6_Multi, levels = c("Normal","Multi"))
                
                PMS_6M_AD_Wilcox = ggboxplot(data = tmp_6m_table,
                                             x = "GROUP",
                                             y = "Change_of_slope_6M",
                                             color = "GROUP",
                                             add = "jitter",
                                             xlab = "Group (AD)",
                                             ylab = "Change of ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                PMS_6M_TD_Wilcox = ggboxplot(data = tmp_6m_table,
                                             x = "M6_Multi",
                                             y = "Change_of_slope_6M",
                                             color = "M6_Multi",
                                             add = "jitter",
                                             xlab = "Group (TD)",
                                             ylab = "Change of ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                
                if(length(subset(tmp_6m_table, GROUP != "Normal")$Change_of_slope_6M) > 1){
                        PMS_6M_AD_var = var.test(subset(tmp_6m_table, GROUP == "Normal")$Change_of_slope_6M,
                                                 subset(tmp_6m_table, GROUP != "Normal")$Change_of_slope_6M)
                        PMS_6M_AD_var = ifelse(PMS_6M_AD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_6M_AD_TT = ggboxplot(data = tmp_6m_table,
                                                 x = "GROUP",
                                                 y = "Change_of_slope_6M",
                                                 color = "GROUP",
                                                 add = "jitter",
                                                 xlab = "Group (AD)",
                                                 ylab = "Change of ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_6M_AD_var))
                }else{
                        PMS_6M_AD_TT = NULL
                }
                
                if(length(subset(tmp_6m_table, M6_Multi != "Normal")$Change_of_slope_6M) > 1){
                        PMS_6M_TD_var = var.test(subset(tmp_6m_table, M6_Multi == "Normal")$Change_of_slope_6M,
                                                 subset(tmp_6m_table, M6_Multi != "Normal")$Change_of_slope_6M)
                        PMS_6M_TD_var = ifelse(PMS_6M_TD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_6M_TD_TT = ggboxplot(data = tmp_6m_table,
                                                 x = "M6_Multi",
                                                 y = "Change_of_slope_6M",
                                                 color = "M6_Multi",
                                                 add = "jitter",
                                                 xlab = "Group (TD)",
                                                 ylab = "Change of ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_6M_TD_var))
                }else{
                        PMS_6M_TD_TT = NULL
                }
                
                PMS_6M_total = ggarrange(PMS_6M_AD_Wilcox, PMS_6M_TD_Wilcox, PMS_6M_AD_TT, PMS_6M_TD_TT , ncol = 2, nrow = 2)
                ggsave(filename = paste0("./Output/",table_name,"_change_of_slope_",month,"M.tiff"),
                       plot = PMS_6M_total,
                       device = "tiff",
                       width = 12,
                       height = 12,
                       dpi = 300)
        } else if(month == 9){
                tmp_9m_table = subset(PMS_extended_data_table, is.na(Change_of_slope_9M) == F)
                tmp_9m_table = subset(tmp_9m_table, ALSFRS_1st >=31 & ALSFRS_1st <= 46)
                tmp_9m_table$GROUP = ifelse(tmp_9m_table$INJ >= 3,
                                            "Multi",
                                            "Normal")
                
                tmp_9m_table$GROUP = factor(tmp_9m_table$GROUP, levels = c("Normal","Multi"))
                tmp_9m_table$M9_Multi = factor(tmp_9m_table$M9_Multi, levels = c("Normal","Multi"))
                
                PMS_9M_AD_Wilcox = ggboxplot(data = tmp_9m_table,
                                             x = "GROUP",
                                             y = "Change_of_slope_9M",
                                             color = "GROUP",
                                             add = "jitter",
                                             xlab = "Group (AD)",
                                             ylab = "Change of ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                PMS_9M_TD_Wilcox = ggboxplot(data = tmp_9m_table,
                                             x = "M9_Multi",
                                             y = "Change_of_slope_9M",
                                             color = "M9_Multi",
                                             add = "jitter",
                                             xlab = "Group (TD)",
                                             ylab = "Change of ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                
                if(length(subset(tmp_9m_table, GROUP != "Normal")$Change_of_slope_9M) > 1){
                        PMS_9M_AD_var = var.test(subset(tmp_9m_table, GROUP == "Normal")$Change_of_slope_9M,
                                                 subset(tmp_9m_table, GROUP != "Normal")$Change_of_slope_9M)
                        PMS_9M_AD_var = ifelse(PMS_9M_AD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_9M_AD_TT = ggboxplot(data = tmp_9m_table,
                                                 x = "GROUP",
                                                 y = "Change_of_slope_9M",
                                                 color = "GROUP",
                                                 add = "jitter",
                                                 xlab = "Group (AD)",
                                                 ylab = "Change of ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_9M_AD_var))
                }else{
                        PMS_9M_AD_TT = NULL
                }
                
                if(length(subset(tmp_9m_table, M9_Multi != "Normal")$Change_of_slope_9M) > 1){
                        PMS_9M_TD_var = var.test(subset(tmp_9m_table, M9_Multi == "Normal")$Change_of_slope_9M,
                                                 subset(tmp_9m_table, M9_Multi != "Normal")$Change_of_slope_9M)
                        PMS_9M_TD_var = ifelse(PMS_9M_TD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_9M_TD_TT = ggboxplot(data = tmp_9m_table,
                                                 x = "M9_Multi",
                                                 y = "Change_of_slope_9M",
                                                 color = "M9_Multi",
                                                 add = "jitter",
                                                 xlab = "Group (TD)",
                                                 ylab = "Change of ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_9M_TD_var))
                }else{
                        PMS_9M_TD_TT = NULL
                }
                
                PMS_9M_total = ggarrange(PMS_9M_AD_Wilcox, PMS_9M_TD_Wilcox, PMS_9M_AD_TT, PMS_9M_TD_TT , ncol = 2, nrow = 2)
                ggsave(filename = paste0("./Output/",table_name,"_change_of_slope_",month,"M.tiff"),
                       plot = PMS_9M_total,
                       device = "tiff",
                       width = 12,
                       height = 12,
                       dpi = 300)
        }else if(month == 12){
                tmp_12m_table = subset(PMS_extended_data_table, is.na(Change_of_slope_12M) == F)
                tmp_12m_table = subset(tmp_12m_table, ALSFRS_1st >=31 & ALSFRS_1st <= 46)
                tmp_12m_table$GROUP = ifelse(tmp_12m_table$INJ >= 3,
                                             "Multi",
                                             "Normal")
                
                tmp_12m_table$GROUP = factor(tmp_12m_table$GROUP, levels = c("Normal","Multi"))
                tmp_12m_table$M12_Multi = factor(tmp_12m_table$M12_Multi, levels = c("Normal","Multi"))
                
                PMS_12M_AD_Wilcox = ggboxplot(data = tmp_12m_table,
                                              x = "GROUP",
                                              y = "Change_of_slope_12M",
                                              color = "GROUP",
                                              add = "jitter",
                                              xlab = "Group (AD)",
                                              ylab = "Change of ALSFRS-R slope",
                                              add.params = list(size = 2.5),
                                              legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                PMS_12M_TD_Wilcox = ggboxplot(data = tmp_12m_table,
                                              x = "M12_Multi",
                                              y = "Change_of_slope_12M",
                                              color = "M12_Multi",
                                              add = "jitter",
                                              xlab = "Group (TD)",
                                              ylab = "Change of ALSFRS-R slope",
                                              add.params = list(size = 2.5),
                                              legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                
                if(length(subset(tmp_12m_table, GROUP != "Normal")$Change_of_slope_12M) > 1){
                        PMS_12M_AD_var = var.test(subset(tmp_12m_table, GROUP == "Normal")$Change_of_slope_12M,
                                                  subset(tmp_12m_table, GROUP != "Normal")$Change_of_slope_12M)
                        PMS_12M_AD_var = ifelse(PMS_12M_AD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_12M_AD_TT = ggboxplot(data = tmp_12m_table,
                                                  x = "GROUP",
                                                  y = "Change_of_slope_12M",
                                                  color = "GROUP",
                                                  add = "jitter",
                                                  xlab = "Group (AD)",
                                                  ylab = "Change of ALSFRS-R slope",
                                                  add.params = list(size = 2.5),
                                                  legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_12M_AD_var))
                }else{
                        PMS_12M_AD_TT = NULL
                }
                
                if(length(subset(tmp_12m_table, M12_Multi != "Normal")$Change_of_slope_12M) > 1){
                        PMS_12M_TD_var = var.test(subset(tmp_12m_table, M12_Multi == "Normal")$Change_of_slope_12M,
                                                  subset(tmp_12m_table, M12_Multi != "Normal")$Change_of_slope_12M)
                        PMS_12M_TD_var = ifelse(PMS_12M_TD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_12M_TD_TT = ggboxplot(data = tmp_12m_table,
                                                  x = "M12_Multi",
                                                  y = "Change_of_slope_12M",
                                                  color = "M12_Multi",
                                                  add = "jitter",
                                                  xlab = "Group (TD)",
                                                  ylab = "Change of ALSFRS-R slope",
                                                  add.params = list(size = 2.5),
                                                  legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_12M_TD_var))
                }else{
                        PMS_12M_TD_TT = NULL
                }
                
                PMS_12M_total = ggarrange(PMS_12M_AD_Wilcox, PMS_12M_TD_Wilcox, PMS_12M_AD_TT, PMS_12M_TD_TT , ncol = 2, nrow = 2)
                ggsave(filename = paste0("./Output/",table_name,"_change_of_slope_",month,"M.tiff"),
                       plot = PMS_12M_total,
                       device = "tiff",
                       width = 12,
                       height = 12,
                       dpi = 300)
        }
}
PMS_slope_comparison_plot_np_run = function(PMS_extended_data_table, month){
        
        if(PMS_extended_data_table == 1){
                PMS_extended_data_table = PMS_data_table_DreamCIS
                table_name = "PMS_data_table_DreamCIS"
        }else{
                PMS_extended_data_table = PMS_data_table_HYCRC
                table_name = "PMS_data_table_HYCRC"
        }
        
        if(month == 6){
                tmp_6m_table = subset(PMS_extended_data_table, is.na(Change_of_slope_6M) == F)
                #tmp_6m_table = subset(tmp_6m_table, ALSFRS_1st >=31 & ALSFRS_1st <= 46)
                tmp_6m_table$GROUP = ifelse(tmp_6m_table$INJ >= 3,
                                            "Multi",
                                            "Normal")
                
                tmp_6m_table$GROUP = factor(tmp_6m_table$GROUP, levels = c("Normal","Multi"))
                tmp_6m_table$M6_Multi = factor(tmp_6m_table$M6_Multi, levels = c("Normal","Multi"))
                
                PMS_6M_AD_Wilcox = ggboxplot(data = tmp_6m_table,
                                             x = "GROUP",
                                             y = "M6_slope",
                                             color = "GROUP",
                                             add = "jitter",
                                             xlab = "Group (AD)",
                                             ylab = " ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                PMS_6M_TD_Wilcox = ggboxplot(data = tmp_6m_table,
                                             x = "M6_Multi",
                                             y = "M6_slope",
                                             color = "M6_Multi",
                                             add = "jitter",
                                             xlab = "Group (TD)",
                                             ylab = " ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                
                if(length(subset(tmp_6m_table, GROUP != "Normal")$Change_of_slope_6M) > 1){
                        PMS_6M_AD_var = var.test(subset(tmp_6m_table, GROUP == "Normal")$Change_of_slope_6M,
                                                 subset(tmp_6m_table, GROUP != "Normal")$Change_of_slope_6M)
                        PMS_6M_AD_var = ifelse(PMS_6M_AD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_6M_AD_TT = ggboxplot(data = tmp_6m_table,
                                                 x = "GROUP",
                                                 y = "M6_slope",
                                                 color = "GROUP",
                                                 add = "jitter",
                                                 xlab = "Group (AD)",
                                                 ylab = " ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_6M_AD_var))
                }else{
                        PMS_6M_AD_TT = NULL
                }
                
                if(length(subset(tmp_6m_table, M6_Multi != "Normal")$Change_of_slope_6M) > 1){
                        PMS_6M_TD_var = var.test(subset(tmp_6m_table, M6_Multi == "Normal")$Change_of_slope_6M,
                                                 subset(tmp_6m_table, M6_Multi != "Normal")$Change_of_slope_6M)
                        PMS_6M_TD_var = ifelse(PMS_6M_TD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_6M_TD_TT = ggboxplot(data = tmp_6m_table,
                                                 x = "M6_Multi",
                                                 y = "M6_slope",
                                                 color = "M6_Multi",
                                                 add = "jitter",
                                                 xlab = "Group (TD)",
                                                 ylab = " ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_6M_TD_var))
                }else{
                        PMS_6M_TD_TT = NULL
                }
                
                PMS_6M_total = ggarrange(PMS_6M_AD_Wilcox, PMS_6M_TD_Wilcox, PMS_6M_AD_TT, PMS_6M_TD_TT , ncol = 2, nrow = 2)
                ggsave(filename = paste0("./Output/",table_name,"_ALSFRS-R_slope_",month,"M.tiff"),
                       plot = PMS_6M_total,
                       device = "tiff",
                       width = 12,
                       height = 12,
                       dpi = 300)
        } else if(month == 9){
                tmp_9m_table = subset(PMS_extended_data_table, is.na(Change_of_slope_9M) == F)
                #tmp_9m_table = subset(tmp_9m_table, ALSFRS_1st >=31 & ALSFRS_1st <= 46)
                tmp_9m_table$GROUP = ifelse(tmp_9m_table$INJ >= 3,
                                            "Multi",
                                            "Normal")
                
                tmp_9m_table$GROUP = factor(tmp_9m_table$GROUP, levels = c("Normal","Multi"))
                tmp_9m_table$M9_Multi = factor(tmp_9m_table$M9_Multi, levels = c("Normal","Multi"))
                
                PMS_9M_AD_Wilcox = ggboxplot(data = tmp_9m_table,
                                             x = "GROUP",
                                             y = "M9_slope",
                                             color = "GROUP",
                                             add = "jitter",
                                             xlab = "Group (AD)",
                                             ylab = " ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                PMS_9M_TD_Wilcox = ggboxplot(data = tmp_9m_table,
                                             x = "M9_Multi",
                                             y = "M9_slope",
                                             color = "M9_Multi",
                                             add = "jitter",
                                             xlab = "Group (TD)",
                                             ylab = " ALSFRS-R slope",
                                             add.params = list(size = 2.5),
                                             legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                
                if(length(subset(tmp_9m_table, GROUP != "Normal")$Change_of_slope_9M) > 1){
                        PMS_9M_AD_var = var.test(subset(tmp_9m_table, GROUP == "Normal")$Change_of_slope_9M,
                                                 subset(tmp_9m_table, GROUP != "Normal")$Change_of_slope_9M)
                        PMS_9M_AD_var = ifelse(PMS_9M_AD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_9M_AD_TT = ggboxplot(data = tmp_9m_table,
                                                 x = "GROUP",
                                                 y = "M9_slope",
                                                 color = "GROUP",
                                                 add = "jitter",
                                                 xlab = "Group (AD)",
                                                 ylab = " ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_9M_AD_var))
                }else{
                        PMS_9M_AD_TT = NULL
                }
                
                if(length(subset(tmp_9m_table, M9_Multi != "Normal")$Change_of_slope_9M) > 1){
                        PMS_9M_TD_var = var.test(subset(tmp_9m_table, M9_Multi == "Normal")$Change_of_slope_9M,
                                                 subset(tmp_9m_table, M9_Multi != "Normal")$Change_of_slope_9M)
                        PMS_9M_TD_var = ifelse(PMS_9M_TD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_9M_TD_TT = ggboxplot(data = tmp_9m_table,
                                                 x = "M9_Multi",
                                                 y = "M9_slope",
                                                 color = "M9_Multi",
                                                 add = "jitter",
                                                 xlab = "Group (TD)",
                                                 ylab = " ALSFRS-R slope",
                                                 add.params = list(size = 2.5),
                                                 legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_9M_TD_var))
                }else{
                        PMS_9M_TD_TT = NULL
                }
                
                PMS_9M_total = ggarrange(PMS_9M_AD_Wilcox, PMS_9M_TD_Wilcox, PMS_9M_AD_TT, PMS_9M_TD_TT , ncol = 2, nrow = 2)
                ggsave(filename = paste0("./Output/",table_name,"_ALSFRS-R_slope_",month,"M.tiff"),
                       plot = PMS_9M_total,
                       device = "tiff",
                       width = 12,
                       height = 12,
                       dpi = 300)
        }else if(month == 12){
                tmp_12m_table = subset(PMS_extended_data_table, is.na(Change_of_slope_12M) == F)
                #tmp_12m_table = subset(tmp_12m_table, ALSFRS_1st >=31 & ALSFRS_1st <= 46)
                tmp_12m_table$GROUP = ifelse(tmp_12m_table$INJ >= 3,
                                             "Multi",
                                             "Normal")
                
                tmp_12m_table$GROUP = factor(tmp_12m_table$GROUP, levels = c("Normal","Multi"))
                tmp_12m_table$M12_Multi = factor(tmp_12m_table$M12_Multi, levels = c("Normal","Multi"))
                
                PMS_12M_AD_Wilcox = ggboxplot(data = tmp_12m_table,
                                              x = "GROUP",
                                              y = "M12_slope",
                                              color = "GROUP",
                                              add = "jitter",
                                              xlab = "Group (AD)",
                                              ylab = " ALSFRS-R slope",
                                              add.params = list(size = 2.5),
                                              legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                PMS_12M_TD_Wilcox = ggboxplot(data = tmp_12m_table,
                                              x = "M12_Multi",
                                              y = "M12_slope",
                                              color = "M12_Multi",
                                              add = "jitter",
                                              xlab = "Group (TD)",
                                              ylab = " ALSFRS-R slope",
                                              add.params = list(size = 2.5),
                                              legend = "none") +
                        stat_compare_means(method = "wilcox.test")
                
                
                if(length(subset(tmp_12m_table, GROUP != "Normal")$Change_of_slope_12M) > 1){
                        PMS_12M_AD_var = var.test(subset(tmp_12m_table, GROUP == "Normal")$Change_of_slope_12M,
                                                  subset(tmp_12m_table, GROUP != "Normal")$Change_of_slope_12M)
                        PMS_12M_AD_var = ifelse(PMS_12M_AD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_12M_AD_TT = ggboxplot(data = tmp_12m_table,
                                                  x = "GROUP",
                                                  y = "M12_slope",
                                                  color = "GROUP",
                                                  add = "jitter",
                                                  xlab = "Group (AD)",
                                                  ylab = " ALSFRS-R slope",
                                                  add.params = list(size = 2.5),
                                                  legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_12M_AD_var))
                }else{
                        PMS_12M_AD_TT = NULL
                }
                
                if(length(subset(tmp_12m_table, M12_Multi != "Normal")$Change_of_slope_12M) > 1){
                        PMS_12M_TD_var = var.test(subset(tmp_12m_table, M12_Multi == "Normal")$Change_of_slope_12M,
                                                  subset(tmp_12m_table, M12_Multi != "Normal")$Change_of_slope_12M)
                        PMS_12M_TD_var = ifelse(PMS_12M_TD_var$p.value > 0.05, TRUE, FALSE)
                        PMS_12M_TD_TT = ggboxplot(data = tmp_12m_table,
                                                  x = "M12_Multi",
                                                  y = "M12_slope",
                                                  color = "M12_Multi",
                                                  add = "jitter",
                                                  xlab = "Group (TD)",
                                                  ylab = " ALSFRS-R slope",
                                                  add.params = list(size = 2.5),
                                                  legend = "none") +
                                stat_compare_means(method = "t.test",
                                                   method.args = list(var.equal = PMS_12M_TD_var))
                }else{
                        PMS_12M_TD_TT = NULL
                }
                
                PMS_12M_total = ggarrange(PMS_12M_AD_Wilcox, PMS_12M_TD_Wilcox, PMS_12M_AD_TT, PMS_12M_TD_TT , ncol = 2, nrow = 2)
                ggsave(filename = paste0("./Output/",table_name,"_ALSFRS-R_slope_",month,"M.tiff"),
                       plot = PMS_12M_total,
                       device = "tiff",
                       width = 12,
                       height = 12,
                       dpi = 300)
        }
}

PMS_shortage_run = function(PMS_extended_data_table, type){
        if(type == "HYCRC"){
                PMS_extended_data_table = PMS_extended_data_table[,c(1,10:21)]
                PMS_extended_data_table[,"Date_2nd"] =  (PMS_extended_data_table[,"Date_2nd"] - PMS_extended_data_table[,"Date_1st"])
                PMS_extended_data_table[,"Date_3rd"] =  (PMS_extended_data_table[,"Date_3rd"] - PMS_extended_data_table[,"Date_1st"])
                PMS_extended_data_table[,"Date_4th"] =  (PMS_extended_data_table[,"Date_4th"] - PMS_extended_data_table[,"Date_1st"])
                PMS_extended_data_table[,"Date_5th"] =  (PMS_extended_data_table[,"Date_5th"] - PMS_extended_data_table[,"Date_1st"])
                PMS_extended_data_table[,"Date_6th"] =  (PMS_extended_data_table[,"Date_6th"] - PMS_extended_data_table[,"Date_1st"])
                PMS_extended_data_table[,"Date_1st"] = 0
                
                PMS_shortage_data_table = data.frame(NULL)
                for (i in 1:nrow(PMS_extended_data_table)){
                        for (j in 1:6){
                                PMS_shortage_data_table[(6*(i-1))+j, 1] = PMS_extended_data_table[i,1]
                                PMS_shortage_data_table[(6*(i-1))+j, 2] = PMS_extended_data_table[i,(2*j)+1]
                                PMS_shortage_data_table[(6*(i-1))+j, 3] = PMS_extended_data_table[i,(2*j)]
                        }
                }
                
                colnames(PMS_shortage_data_table) = c("subject_id",
                                                      "ALSFRS_Delta",
                                                      "ALSFRS_R_Total")
                
                PMS_shortage_data_table = PMS_shortage_data_table[complete.cases(PMS_shortage_data_table),]
                row.names(PMS_shortage_data_table) = 1:nrow(PMS_shortage_data_table)
                
        }else if(type == "DreamCIS"){
                PMS_extended_data_table = PMS_extended_data_table[PMS_extended_data_table$SUBJ_ID %in% PMS_ALSFRS_table$SUBJ_ID,]
                PMS_ALSFRS_table = PMS_ALSFRS_table[PMS_ALSFRS_table$SUBJ_ID %in% PMS_extended_data_table$SUBJ_ID,]
                
                PMS_extended_data_table = PMS_extended_data_table[order(PMS_extended_data_table$SUBJ_ID),]
                PMS_ALSFRS_table = PMS_ALSFRS_table[order(PMS_ALSFRS_table$SUBJ_ID),]
                
                PMS_shortage_data_table = data.frame(NULL)
                for (i in 1:nrow(PMS_ALSFRS_table)){
                        for (j in 1:15){
                                PMS_shortage_data_table[(15*(i-1))+j, 1] = PMS_ALSFRS_table[i,1]
                                PMS_shortage_data_table[(15*(i-1))+j, 2] = (j-1) * 28
                                PMS_shortage_data_table[(15*(i-1))+j, 3] = PMS_ALSFRS_table[i,(j+2)]
                        }
                }
                
                colnames(PMS_shortage_data_table) = c("subject_id",
                                                      "ALSFRS_Delta",
                                                      "ALSFRS_R_Total")
                
                PMS_shortage_data_table = PMS_shortage_data_table[complete.cases(PMS_shortage_data_table),]
                row.names(PMS_shortage_data_table) = 1:nrow(PMS_shortage_data_table)
                
        }
        return(PMS_shortage_data_table)
}
PROACT_summarized_run = function(ALSFRS, Placebo, cp_type, select_type, ALSFRS_cut = FALSE){
        ALSFRS_subjects = unique(ALSFRS$subject_id)
        PROACT_slope = data.frame(subject_id = ALSFRS_subjects,
                                  M6_slope = NA,
                                  M12_slope = NA)
        
        for (i in 1:nrow(PROACT_slope)){
                tmp = subset(ALSFRS, subject_id == ALSFRS_subjects[i])
                tmp$ALSFRS_Delta = (tmp$ALSFRS_Delta / 28)
                
                if(length(which(tmp$ALSFRS_Delta == 0)) == 1){
                        tmp_base_ALSFRS = tmp[which(tmp$ALSFRS_Delta == 0),3]
                        if(ALSFRS_cut == TRUE){
                                if(tmp_base_ALSFRS < 31 | tmp_base_ALSFRS > 46){
                                        tmp_base_ALSFRS = NA
                                }
                        }
                }else if(length(which(tmp$ALSFRS_Delta == 0)) > 1){
                        tmp_base_ALSFRS = NA
                }
                
                
                if(cp_type == "point_style" & select_type == "max"){
                        tmp_6m = subset(tmp, ALSFRS_Delta >= 5 & ALSFRS_Delta <= 7)
                        if(nrow(tmp_6m) == 1){
                                tmp_6m_interval = tmp_6m[1,2]
                                tmp_6m_ALSFRS = tmp_6m[1,3]
                        }else if(nrow(tmp_6m) == 0){
                                tmp_6m_interval = NA
                                tmp_6m_ALSFRS = NA
                        }else if(nrow(tmp_6m) > 1){
                                tmp_6m_interval = tmp_6m[min(which(tmp_6m$ALSFRS_R_Total == max(tmp_6m[,3]))),2]
                                tmp_6m_ALSFRS = tmp_6m[min(which(tmp_6m$ALSFRS_R_Total == max(tmp_6m[,3]))),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 13)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),2]
                                tmp_12m_ALSFRS = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_6m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_6m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }else if(cp_type == "point_style" & select_type == "most"){
                        tmp_6m = subset(tmp, ALSFRS_Delta >= 5 & ALSFRS_Delta <= 7)
                        if(nrow(tmp_6m) == 1){
                                tmp_6m_interval = tmp_6m[1,2]
                                tmp_6m_ALSFRS = tmp_6m[1,3]
                        }else if(nrow(tmp_6m) == 0){
                                tmp_6m_interval = NA
                                tmp_6m_ALSFRS = NA
                        }else if(nrow(tmp_6m) > 1){
                                tmp_6m_interval = tmp_6m[which.min(abs(tmp_6m$ALSFRS_Delta - 6)),2]
                                tmp_6m_ALSFRS = tmp_6m[which.min(abs(tmp_6m$ALSFRS_Delta - 6)),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 13)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),2]
                                tmp_12m_ALSFRS = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_6m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_6m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }else if(cp_type == "block_style" & select_type == "max"){
                        tmp_6m = subset(tmp, ALSFRS_Delta >= 5 & ALSFRS_Delta <= 8)
                        if(nrow(tmp_6m) == 1){
                                tmp_6m_interval = tmp_6m[1,2]
                                tmp_6m_ALSFRS = tmp_6m[1,3]
                        }else if(nrow(tmp_6m) == 0){
                                tmp_6m_interval = NA
                                tmp_6m_ALSFRS = NA
                        }else if(nrow(tmp_6m) > 1){
                                tmp_6m_interval = tmp_6m[min(which(tmp_6m$ALSFRS_R_Total == max(tmp_6m[,3]))),2]
                                tmp_6m_ALSFRS = tmp_6m[min(which(tmp_6m$ALSFRS_R_Total == max(tmp_6m[,3]))),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 14)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),2]
                                tmp_12m_ALSFRS = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_6m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_6m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }else if(cp_type == "block_style" & select_type == "most"){
                        tmp_6m = subset(tmp, ALSFRS_Delta >= 5 & ALSFRS_Delta <= 8)
                        if(nrow(tmp_6m) == 1){
                                tmp_6m_interval = tmp_6m[1,2]
                                tmp_6m_ALSFRS = tmp_6m[1,3]
                        }else if(nrow(tmp_6m) == 0){
                                tmp_6m_interval = NA
                                tmp_6m_ALSFRS = NA
                        }else if(nrow(tmp_6m) > 1){
                                tmp_6m_interval = tmp_6m[which.min(abs(tmp_6m$ALSFRS_Delta - 6)),2]
                                tmp_6m_ALSFRS = tmp_6m[which.min(abs(tmp_6m$ALSFRS_Delta - 6)),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 14)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),2]
                                tmp_12m_ALSFRS = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_6m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_6m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }
                
        }
        PROACT_slope = merge(PROACT_slope, Placebo[,c(1,5)], by = "subject_id")
        PROACT_slope[,5] = "PMS"
        
        return(PROACT_slope)
}
PROACT_stat_run = function(PROACT_ALSFRS, PROACT_Placebo){
        PROACT_slope_point_max = PROACT_summarized_run(PROACT_ALSFRS,
                                                       PROACT_Placebo,
                                                       cp_type = "point_style",
                                                       select_type = "max")
        
        PROACT_slope_point_most = PROACT_summarized_run(PROACT_ALSFRS,
                                                        PROACT_Placebo,
                                                        cp_type = "point_style",
                                                        select_type = "most")
        
        PROACT_slope_block_max = PROACT_summarized_run(PROACT_ALSFRS,
                                                       PROACT_Placebo,
                                                       cp_type = "block_style",
                                                       select_type = "max")
        
        PROACT_slope_block_most = PROACT_summarized_run(PROACT_ALSFRS,
                                                        PROACT_Placebo,
                                                        cp_type = "block_style",
                                                        select_type = "most")
        
        ####
        PROACT_slope_point_max_cut = PROACT_summarized_run(PROACT_ALSFRS,
                                                           PROACT_Placebo,
                                                           ALSFRS_cut = TRUE,
                                                           cp_type = "point_style",
                                                           select_type = "max")
        
        PROACT_slope_point_most_cut = PROACT_summarized_run(PROACT_ALSFRS,
                                                            PROACT_Placebo,
                                                            ALSFRS_cut = TRUE,
                                                            cp_type = "point_style",
                                                            select_type = "most")
        
        PROACT_slope_block_max_cut = PROACT_summarized_run(PROACT_ALSFRS,
                                                           PROACT_Placebo,
                                                           ALSFRS_cut = TRUE,
                                                           cp_type = "block_style",
                                                           select_type = "max")
        
        PROACT_slope_block_most_cut = PROACT_summarized_run(PROACT_ALSFRS,
                                                            PROACT_Placebo,
                                                            ALSFRS_cut = TRUE,
                                                            cp_type = "block_style",
                                                            select_type = "most")
        
        tmp = data.frame(cp_type = c(rep("point",4),rep("block",4)),
                         select_type = rep(c("max", "most"),4), 
                         ALSFRSR_cut = rep(c("No","No","Yes","Yes"),2))
        
        tmp[1,4] = mean(PROACT_slope_point_max$M6_slope, na.rm = T)
        tmp[1,5] = sd(PROACT_slope_point_max$M6_slope, na.rm = T)
        tmp[1,6] = mean(PROACT_slope_point_max$M12_slope, na.rm = T)
        tmp[1,7] = sd(PROACT_slope_point_max$M12_slope, na.rm = T)
        
        tmp[2,4] = mean(PROACT_slope_point_most$M6_slope, na.rm = T)
        tmp[2,5] = sd(PROACT_slope_point_most$M6_slope, na.rm = T)
        tmp[2,6] = mean(PROACT_slope_point_most$M12_slope, na.rm = T)
        tmp[2,7] = sd(PROACT_slope_point_most$M12_slope, na.rm = T)
        
        tmp[3,4] = mean(PROACT_slope_point_max_cut$M6_slope, na.rm = T)
        tmp[3,5] = sd(PROACT_slope_point_max_cut$M6_slope, na.rm = T)
        tmp[3,6] = mean(PROACT_slope_point_max_cut$M12_slope, na.rm = T)
        tmp[3,7] = sd(PROACT_slope_point_max_cut$M12_slope, na.rm = T)
        
        tmp[4,4] = mean(PROACT_slope_point_most_cut$M6_slope, na.rm = T)
        tmp[4,5] = sd(PROACT_slope_point_most_cut$M6_slope, na.rm = T)
        tmp[4,6] = mean(PROACT_slope_point_most_cut$M12_slope, na.rm = T)
        tmp[4,7] = sd(PROACT_slope_point_most_cut$M12_slope, na.rm = T)
        
        #### 
        
        tmp[5,4] = mean(PROACT_slope_block_max$M6_slope, na.rm = T)
        tmp[5,5] = sd(PROACT_slope_block_max$M6_slope, na.rm = T)
        tmp[5,6] = mean(PROACT_slope_block_max$M12_slope, na.rm = T)
        tmp[5,7] = sd(PROACT_slope_block_max$M12_slope, na.rm = T)
        
        tmp[6,4] = mean(PROACT_slope_block_most$M6_slope, na.rm = T)
        tmp[6,5] = sd(PROACT_slope_block_most$M6_slope, na.rm = T)
        tmp[6,6] = mean(PROACT_slope_block_most$M12_slope, na.rm = T)
        tmp[6,7] = sd(PROACT_slope_block_most$M12_slope, na.rm = T)
        
        tmp[7,4] = mean(PROACT_slope_block_max_cut$M6_slope, na.rm = T)
        tmp[7,5] = sd(PROACT_slope_block_max_cut$M6_slope, na.rm = T)
        tmp[7,6] = mean(PROACT_slope_block_max_cut$M12_slope, na.rm = T)
        tmp[7,7] = sd(PROACT_slope_block_max_cut$M12_slope, na.rm = T)
        
        tmp[8,4] = mean(PROACT_slope_block_most_cut$M6_slope, na.rm = T)
        tmp[8,5] = sd(PROACT_slope_block_most_cut$M6_slope, na.rm = T)
        tmp[8,6] = mean(PROACT_slope_block_most_cut$M12_slope, na.rm = T)
        tmp[8,7] = sd(PROACT_slope_block_most_cut$M12_slope, na.rm = T)
        
        return(tmp)
}

PROACT_summarized_short_run = function(ALSFRS, Placebo, cp_type, select_type, ALSFRS_cut = FALSE){
        ALSFRS_subjects = unique(ALSFRS$subject_id)
        PROACT_slope = data.frame(subject_id = ALSFRS_subjects,
                                  M6_slope = NA,
                                  M12_slope = NA)
        
        for (i in 1:nrow(PROACT_slope)){
                tmp = subset(ALSFRS, subject_id == ALSFRS_subjects[i])
                tmp$ALSFRS_Delta = (tmp$ALSFRS_Delta / 28)
                
                if(length(which(tmp$ALSFRS_Delta == 0)) == 1){
                        tmp_base_ALSFRS = tmp[which(tmp$ALSFRS_Delta == 0),3]
                        if(ALSFRS_cut == TRUE){
                                if(tmp_base_ALSFRS < 31 | tmp_base_ALSFRS > 46){
                                        tmp_base_ALSFRS = NA
                                }
                        }
                }else if(length(which(tmp$ALSFRS_Delta == 0)) > 1){
                        tmp_base_ALSFRS = NA
                }
                
                
                if(cp_type == "point_style" & select_type == "max"){
                        tmp_3m = subset(tmp, ALSFRS_Delta >= 3 & ALSFRS_Delta <= 5)
                        if(nrow(tmp_3m) == 1){
                                tmp_3m_interval = tmp_3m[1,2]
                                tmp_3m_ALSFRS = tmp_3m[1,3]
                        }else if(nrow(tmp_3m) == 0){
                                tmp_3m_interval = NA
                                tmp_3m_ALSFRS = NA
                        }else if(nrow(tmp_3m) > 1){
                                tmp_3m_interval = tmp_3m[min(which(tmp_3m$ALSFRS_R_Total == max(tmp_3m[,3]))),2]
                                tmp_3m_ALSFRS = tmp_3m[min(which(tmp_3m$ALSFRS_R_Total == max(tmp_3m[,3]))),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 13)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),2]
                                tmp_12m_ALSFRS = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_3m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_3m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }else if(cp_type == "point_style" & select_type == "most"){
                        tmp_3m = subset(tmp, ALSFRS_Delta >= 2 & ALSFRS_Delta < 3)
                        if(nrow(tmp_3m) == 1){
                                tmp_3m_interval = tmp_3m[1,2]
                                tmp_3m_ALSFRS = tmp_3m[1,3]
                        }else if(nrow(tmp_3m) == 0){
                                tmp_3m_interval = NA
                                tmp_3m_ALSFRS = NA
                        }else if(nrow(tmp_3m) > 1){
                                tmp_3m_interval = tmp_3m[which.min(abs(tmp_3m$ALSFRS_Delta - 2)),2]
                                tmp_3m_ALSFRS = tmp_3m[which.min(abs(tmp_3m$ALSFRS_Delta - 2)),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 13)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),2]
                                tmp_12m_ALSFRS = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_3m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_3m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }else if(cp_type == "block_style" & select_type == "max"){
                        tmp_3m = subset(tmp, ALSFRS_Delta >= 2 & ALSFRS_Delta <= 5)
                        if(nrow(tmp_3m) == 1){
                                tmp_3m_interval = tmp_3m[1,2]
                                tmp_3m_ALSFRS = tmp_3m[1,3]
                        }else if(nrow(tmp_3m) == 0){
                                tmp_3m_interval = NA
                                tmp_3m_ALSFRS = NA
                        }else if(nrow(tmp_3m) > 1){
                                tmp_3m_interval = tmp_3m[min(which(tmp_3m$ALSFRS_R_Total == max(tmp_3m[,3]))),2]
                                tmp_3m_ALSFRS = tmp_3m[min(which(tmp_3m$ALSFRS_R_Total == max(tmp_3m[,3]))),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 14)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),2]
                                tmp_12m_ALSFRS = tmp_12m[min(which(tmp_12m$ALSFRS_R_Total == max(tmp_12m[,3]))),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_3m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_3m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }else if(cp_type == "block_style" & select_type == "most"){
                        tmp_3m = subset(tmp, ALSFRS_Delta >= 2 & ALSFRS_Delta <= 5)
                        if(nrow(tmp_3m) == 1){
                                tmp_3m_interval = tmp_3m[1,2]
                                tmp_3m_ALSFRS = tmp_3m[1,3]
                        }else if(nrow(tmp_3m) == 0){
                                tmp_3m_interval = NA
                                tmp_3m_ALSFRS = NA
                        }else if(nrow(tmp_3m) > 1){
                                tmp_3m_interval = tmp_3m[which.min(abs(tmp_3m$ALSFRS_Delta - 3)),2]
                                tmp_3m_ALSFRS = tmp_3m[which.min(abs(tmp_3m$ALSFRS_Delta - 3)),3]
                        }
                        
                        tmp_12m = subset(tmp, ALSFRS_Delta >= 11 & ALSFRS_Delta <= 14)
                        if(nrow(tmp_12m) == 1){
                                tmp_12m_interval = tmp_12m[1,2]
                                tmp_12m_ALSFRS = tmp_12m[1,3]
                        }else if(nrow(tmp_12m) == 0){
                                tmp_12m_interval = NA
                                tmp_12m_ALSFRS = NA
                        }else if(nrow(tmp_12m) > 1){
                                tmp_12m_interval = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),2]
                                tmp_12m_ALSFRS = tmp_12m[which.min(abs(tmp_12m$ALSFRS_Delta - 12)),3]
                        }
                        
                        PROACT_slope[i,"M6_slope"] = (tmp_3m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_3m_interval
                        PROACT_slope[i,"M12_slope"] = (tmp_12m_ALSFRS - tmp_base_ALSFRS) /
                                tmp_12m_interval
                }
                
        }
        PROACT_slope = merge(PROACT_slope, Placebo[,c(1,5)], by = "subject_id")
        PROACT_slope[,5] = "PMS"
        return(PROACT_slope)
}
PROACT_stat_short_run = function(PROACT_ALSFRS, PROACT_Placebo){
        PROACT_slope_point_max = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                       PROACT_Placebo,
                                                       cp_type = "point_style",
                                                       select_type = "max")
        
        PROACT_slope_point_most = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                        PROACT_Placebo,
                                                        cp_type = "point_style",
                                                        select_type = "most")
        
        PROACT_slope_block_max = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                       PROACT_Placebo,
                                                       cp_type = "block_style",
                                                       select_type = "max")
        
        PROACT_slope_block_most = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                        PROACT_Placebo,
                                                        cp_type = "block_style",
                                                        select_type = "most")
        
        ####
        PROACT_slope_point_max_cut = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                           PROACT_Placebo,
                                                           ALSFRS_cut = TRUE,
                                                           cp_type = "point_style",
                                                           select_type = "max")
        
        PROACT_slope_point_most_cut = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                            PROACT_Placebo,
                                                            ALSFRS_cut = TRUE,
                                                            cp_type = "point_style",
                                                            select_type = "most")
        
        PROACT_slope_block_max_cut = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                           PROACT_Placebo,
                                                           ALSFRS_cut = TRUE,
                                                           cp_type = "block_style",
                                                           select_type = "max")
        
        PROACT_slope_block_most_cut = PROACT_summarized_short_run(PROACT_ALSFRS,
                                                            PROACT_Placebo,
                                                            ALSFRS_cut = TRUE,
                                                            cp_type = "block_style",
                                                            select_type = "most")
        
        tmp = data.frame(cp_type = c(rep("point",4),rep("block",4)),
                         select_type = rep(c("max", "most"),4), 
                         ALSFRSR_cut = rep(c("No","No","Yes","Yes"),2))
        
        tmp[1,4] = mean(PROACT_slope_point_max$M6_slope, na.rm = T)
        tmp[1,5] = sd(PROACT_slope_point_max$M6_slope, na.rm = T)
        tmp[1,6] = mean(PROACT_slope_point_max$M12_slope, na.rm = T)
        tmp[1,7] = sd(PROACT_slope_point_max$M12_slope, na.rm = T)
        
        tmp[2,4] = mean(PROACT_slope_point_most$M6_slope, na.rm = T)
        tmp[2,5] = sd(PROACT_slope_point_most$M6_slope, na.rm = T)
        tmp[2,6] = mean(PROACT_slope_point_most$M12_slope, na.rm = T)
        tmp[2,7] = sd(PROACT_slope_point_most$M12_slope, na.rm = T)
        
        tmp[3,4] = mean(PROACT_slope_point_max_cut$M6_slope, na.rm = T)
        tmp[3,5] = sd(PROACT_slope_point_max_cut$M6_slope, na.rm = T)
        tmp[3,6] = mean(PROACT_slope_point_max_cut$M12_slope, na.rm = T)
        tmp[3,7] = sd(PROACT_slope_point_max_cut$M12_slope, na.rm = T)
        
        tmp[4,4] = mean(PROACT_slope_point_most_cut$M6_slope, na.rm = T)
        tmp[4,5] = sd(PROACT_slope_point_most_cut$M6_slope, na.rm = T)
        tmp[4,6] = mean(PROACT_slope_point_most_cut$M12_slope, na.rm = T)
        tmp[4,7] = sd(PROACT_slope_point_most_cut$M12_slope, na.rm = T)
        
        #### 
        
        tmp[5,4] = mean(PROACT_slope_block_max$M6_slope, na.rm = T)
        tmp[5,5] = sd(PROACT_slope_block_max$M6_slope, na.rm = T)
        tmp[5,6] = mean(PROACT_slope_block_max$M12_slope, na.rm = T)
        tmp[5,7] = sd(PROACT_slope_block_max$M12_slope, na.rm = T)
        
        tmp[6,4] = mean(PROACT_slope_block_most$M6_slope, na.rm = T)
        tmp[6,5] = sd(PROACT_slope_block_most$M6_slope, na.rm = T)
        tmp[6,6] = mean(PROACT_slope_block_most$M12_slope, na.rm = T)
        tmp[6,7] = sd(PROACT_slope_block_most$M12_slope, na.rm = T)
        
        tmp[7,4] = mean(PROACT_slope_block_max_cut$M6_slope, na.rm = T)
        tmp[7,5] = sd(PROACT_slope_block_max_cut$M6_slope, na.rm = T)
        tmp[7,6] = mean(PROACT_slope_block_max_cut$M12_slope, na.rm = T)
        tmp[7,7] = sd(PROACT_slope_block_max_cut$M12_slope, na.rm = T)
        
        tmp[8,4] = mean(PROACT_slope_block_most_cut$M6_slope, na.rm = T)
        tmp[8,5] = sd(PROACT_slope_block_most_cut$M6_slope, na.rm = T)
        tmp[8,6] = mean(PROACT_slope_block_most_cut$M12_slope, na.rm = T)
        tmp[8,7] = sd(PROACT_slope_block_most_cut$M12_slope, na.rm = T)
        
        return(tmp)
}