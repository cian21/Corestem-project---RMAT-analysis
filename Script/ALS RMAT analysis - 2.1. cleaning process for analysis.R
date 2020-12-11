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
ALSFRSR_file_allo_P2cont = subset(ALSFRSR_file_allo_P2cont, GROUP == 0)
ALSFRSR_file_allo_P2cont$M4_Changes = (ALSFRSR_file_allo_P2cont$M4 - ALSFRSR_file_allo_P2cont$M0) / 4
ALSFRSR_file_allo_P2cont$M6_Changes = (ALSFRSR_file_allo_P2cont$M6 - ALSFRSR_file_allo_P2cont$M0) / 6
ALSFRSR_file_allo_P2cont$M12_Changes = (ALSFRSR_file_allo_P2cont$M12 - ALSFRSR_file_allo_P2cont$M0) / 12
ALSFRSR_file_allo_P12_FCSPMM$Changes = ALSFRSR_file_allo_P12_FCSPMM$M12 - ALSFRSR_file_allo_P12_FCSPMM$M0

#### 3. PMS slope comparison ####
PMS_DreamCIS_stat_run = function(PMS_ALSFRS_init_table){
    tmp_table = data.frame(test = NA, value = NA)
    PMS_data_table_DreamCIS = PMS_slope_comparison_run(PMS_ALSFRS_init_table, "DreamCIS")
    PMS_data_table_DreamCIS_shortage = PMS_shortage_run(PMS_data_table_DreamCIS, "DreamCIS")
    PMS_data_table_DreamCIS_shortage_placebo = data.frame(subject_id = unique(PMS_data_table_DreamCIS_shortage$subject_id),
                                                          Study_Arm = "Placebo")
    
    PMS_data_table_DreamCIS_Multi = PMS_data_table_DreamCIS[c(2,39,41,10)]
    colnames(PMS_data_table_DreamCIS_Multi)[1] = "subject_id"
    PMS_data_table_DreamCIS_shortage_placebo = merge(PMS_data_table_DreamCIS_shortage_placebo,
                                                     PMS_data_table_DreamCIS_Multi,
                                                     by = "subject_id")
    
    PMS_data_table_DreamCIS_shortage_placebo_M6_Normal = subset(PMS_data_table_DreamCIS_shortage_placebo,
                                                                M6_Multi == "Normal")
    PMS_data_table_DreamCIS_shortage_placebo_M6_Multi = subset(PMS_data_table_DreamCIS_shortage_placebo,
                                                               M6_Multi == "Multi")
    
    PMS_data_table_DreamCIS_shortage_placebo_M12_Normal = subset(PMS_data_table_DreamCIS_shortage_placebo,
                                                                 M12_Multi == "Normal")
    PMS_data_table_DreamCIS_shortage_placebo_M12_Multi = subset(PMS_data_table_DreamCIS_shortage_placebo,
                                                                M12_Multi == "Multi")
    
    PMS_DREAM_slope_point_most_cut_short = PROACT_summarized_short_run(PMS_data_table_DreamCIS_shortage,
                                                                       PMS_data_table_DreamCIS_shortage_placebo,
                                                                       ALSFRS_cut = F,
                                                                       cp_type = "point_style",
                                                                       select_type = "most")
    
    PMS_DREAM_slope_point_most_cut = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage,
                                                           PMS_data_table_DreamCIS_shortage_placebo,
                                                           ALSFRS_cut = F,
                                                           cp_type = "point_style",
                                                           select_type = "most")
    
    PMS_data_table_DreamCIS_shortage_M6_Normal = subset(PMS_data_table_DreamCIS_shortage,
                                                        PMS_data_table_DreamCIS_shortage$subject_id %in%
                                                            PMS_data_table_DreamCIS_shortage_placebo_M6_Normal$subject_id)
    
    PMS_data_table_DreamCIS_shortage_M6_Multi = subset(PMS_data_table_DreamCIS_shortage,
                                                       PMS_data_table_DreamCIS_shortage$subject_id %in%
                                                           PMS_data_table_DreamCIS_shortage_placebo_M6_Multi$subject_id)
    
    PMS_data_table_DreamCIS_shortage_M12_Normal = subset(PMS_data_table_DreamCIS_shortage,
                                                         PMS_data_table_DreamCIS_shortage$subject_id %in%
                                                             PMS_data_table_DreamCIS_shortage_placebo_M12_Normal$subject_id)
    
    PMS_data_table_DreamCIS_shortage_M12_Multi = subset(PMS_data_table_DreamCIS_shortage,
                                                        PMS_data_table_DreamCIS_shortage$subject_id %in%
                                                            PMS_data_table_DreamCIS_shortage_placebo_M12_Multi$subject_id)
    
    PMS_DreamCIS_slope_point_most_cut_M6_Normal = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M6_Normal,
                                                                        PMS_data_table_DreamCIS_shortage_placebo_M6_Normal,
                                                                        ALSFRS_cut = F,
                                                                        cp_type = "point_style",
                                                                        select_type = "most")
    
    PMS_DreamCIS_slope_point_most_cut_M6_Multi = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M6_Multi,
                                                                       PMS_data_table_DreamCIS_shortage_placebo_M6_Multi,
                                                                       ALSFRS_cut = F,
                                                                       cp_type = "point_style",
                                                                       select_type = "most")
    
    PMS_DreamCIS_slope_point_most_cut_M12_Normal = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M12_Normal,
                                                                         PMS_data_table_DreamCIS_shortage_placebo_M12_Normal,
                                                                         ALSFRS_cut = F,
                                                                         cp_type = "point_style",
                                                                         select_type = "most")
    
    PMS_DreamCIS_slope_point_most_cut_M12_Multi = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M12_Multi,
                                                                        PMS_data_table_DreamCIS_shortage_placebo_M12_Multi,
                                                                        ALSFRS_cut = F,
                                                                        cp_type = "point_style",
                                                                        select_type = "most")
    tmp_table[1:16,1] = c("wilcox M4 1VS2",
                    "wilcox M6 1VS23",
                    "wilcox M6 1VS2",
                    "wilcox M6 1VS3",
                    "wilcox M6 2VS3",
                    "wilcox M12 1VS23",
                    "wilcox M12 1VS2",
                    "wilcox M12 1VS3",
                    "wilcox M12 2VS3",
                    "ANCOVA M4 1vs23",
                    "ANCOVA M6 1vs23",
                    "ANCOVA M6 1vs2",
                    "ANCOVA M6 1vs3",
                    "ANCOVA M12 1vs23",
                    "ANCOVA M12 1vs2",
                    "ANCOVA M12 1vs3")
    
    tmp_table[1,2] = wilcox.test(PMS_DREAM_slope_point_most_cut_short$M6_slope,
                ALSFRSR_file_allo_P2cont$M4_Changes)$p.value
    
    tmp_table[2,2] = wilcox.test(PMS_DREAM_slope_point_most_cut$M6_slope,
                ALSFRSR_file_allo_P2cont$M6_Changes)$p.value
    
    tmp_table[3,2] = wilcox.test(PMS_DreamCIS_slope_point_most_cut_M6_Normal$M6_slope,
                ALSFRSR_file_allo_P2cont$M6_Changes)$p.value
    
    tmp_table[4,2] = wilcox.test(PMS_DreamCIS_slope_point_most_cut_M6_Multi$M6_slope,
                ALSFRSR_file_allo_P2cont$M6_Changes)$p.value
    
    tmp_table[5,2] = wilcox.test(PMS_DreamCIS_slope_point_most_cut_M6_Normal$M6_slope,
                PMS_DreamCIS_slope_point_most_cut_M6_Multi$M6_slope)$p.value
    
    tmp_table[6,2] = wilcox.test(PMS_DREAM_slope_point_most_cut$M12_slope,
                ALSFRSR_file_allo_P2cont$M12_Changes)$p.value
    
    tmp_table[7,2] = wilcox.test(PMS_DreamCIS_slope_point_most_cut_M12_Normal$M12_slope,
                ALSFRSR_file_allo_P2cont$M12_Changes)$p.value
    
    tmp_table[8,2] = wilcox.test(PMS_DreamCIS_slope_point_most_cut_M12_Multi$M12_slope,
                ALSFRSR_file_allo_P2cont$M12_Changes)$p.value

    tmp_table[9,2] = wilcox.test(PMS_DreamCIS_slope_point_most_cut_M12_Multi$M12_slope,
                PMS_DreamCIS_slope_point_most_cut_M12_Normal$M12_slope)$p.value
    
    
    P12_Cont_tmp = ALSFRSR_file_allo_P2cont[,c(1,17,18,3,2)]
    P12_Cont_tmp_4m = ALSFRSR_file_allo_P2cont[,c(1,16,18,3,2)]
    
    
    colnames(P12_Cont_tmp) = colnames(PMS_DREAM_slope_point_most_cut)
    colnames(P12_Cont_tmp_4m) = colnames(PMS_DREAM_slope_point_most_cut)
    
    
    P12_Cont_tmp_4M_total = rbind(PMS_DREAM_slope_point_most_cut_short,
                                  P12_Cont_tmp_4m)
    
    tmp_table[10,2] = summary(lm(data = P12_Cont_tmp_4M_total, formula = M6_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    
    P12_Cont_tmp_6M_total = rbind(PMS_DREAM_slope_point_most_cut,
                                  P12_Cont_tmp)
    
    tmp_table[11,2] = summary(lm(data = P12_Cont_tmp_6M_total, formula = M6_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    P12_Cont_tmp_6M_normal = rbind(PMS_DreamCIS_slope_point_most_cut_M6_Normal,
                                  P12_Cont_tmp)
    
    tmp_table[12,2] = summary(lm(data = P12_Cont_tmp_6M_normal, formula = M6_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    P12_Cont_tmp_6M_multi = rbind(PMS_DreamCIS_slope_point_most_cut_M6_Multi,
                                   P12_Cont_tmp)
    
    tmp_table[13,2] = summary(lm(data = P12_Cont_tmp_6M_multi, formula = M6_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    P12_Cont_tmp_12M_total = rbind(PMS_DREAM_slope_point_most_cut,
                                  P12_Cont_tmp)
    
    tmp_table[14,2] = summary(lm(data = P12_Cont_tmp_12M_total, formula = M12_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    P12_Cont_tmp_12M_normal = rbind(PMS_DreamCIS_slope_point_most_cut_M12_Normal,
                                   P12_Cont_tmp)
    
    tmp_table[15,2] = summary(lm(data = P12_Cont_tmp_12M_normal, formula = M12_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    P12_Cont_tmp_12M_multi = rbind(PMS_DreamCIS_slope_point_most_cut_M12_Multi,
                                  P12_Cont_tmp)
    
    tmp_table[16,2] = summary(lm(data = P12_Cont_tmp_12M_multi, formula = M12_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    PMS_DreamCIS_slope_point_most_cut_M6_Normal$V5 = "Normal"
    PMS_DreamCIS_slope_point_most_cut_M6_Multi$V5 = "Multi"
    PMS_DreamCIS_slope_point_most_cut_M12_Normal$V5 = "Normal"
    PMS_DreamCIS_slope_point_most_cut_M12_Multi$V5 = "Multi"
    
    M6_Normal_Multi = rbind(PMS_DreamCIS_slope_point_most_cut_M6_Normal,PMS_DreamCIS_slope_point_most_cut_M6_Multi)
    M12_Normal_Multi = rbind(PMS_DreamCIS_slope_point_most_cut_M12_Normal,PMS_DreamCIS_slope_point_most_cut_M12_Multi)
    
    summary(lm(data = M6_Normal_Multi, formula = M6_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    summary(lm(data = M12_Normal_Multi, formula = M12_slope ~ V5 + ALSFRS_1st))$coefficients[2,4]
    
    return(tmp_table)
}

tmp = PMS_DreamCIS_stat_run(PMS_ALSFRS_init_table)
tmp

#### 4. PROACT ####
PROACT_Placebo = subset(PROACT_Placebo, Study_Arm == "Placebo")
PROACT_Placebo = subset(PROACT_Placebo, Treatment_Group_Delta == 0)

PROACT_ALSFRS = subset(PROACT_ALSFRS, subject_id %in% PROACT_Placebo$subject_id)
PROACT_ALSFRS = PROACT_ALSFRS[,c(1,13,15)]
PROACT_ALSFRS = PROACT_ALSFRS[complete.cases(PROACT_ALSFRS),]
PROACT_stat_table = PROACT_stat_run(PROACT_ALSFRS, PROACT_Placebo)
PROACT_stat_short_table = PROACT_stat_short_run(PROACT_ALSFRS, PROACT_Placebo)

#### 4. SLOPE AND CYTOKINE ####
PMS_DreamCIS_stat_run = function(PMS_ALSFRS_init_table){
    PMS_data_table_DreamCIS = PMS_slope_comparison_run(PMS_ALSFRS_init_table, "DreamCIS")
    PMS_data_table_DreamCIS_shortage = PMS_shortage_run(PMS_data_table_DreamCIS, "DreamCIS")
    PMS_data_table_DreamCIS_shortage_placebo = data.frame(subject_id = unique(PMS_data_table_DreamCIS_shortage$subject_id),
                                                          Study_Arm = "Placebo")
    
    PMS_data_table_DreamCIS_Multi = PMS_data_table_DreamCIS[c(2,39,41,10)]
    colnames(PMS_data_table_DreamCIS_Multi)[1] = "subject_id"
    PMS_data_table_DreamCIS_shortage_placebo = merge(PMS_data_table_DreamCIS_shortage_placebo,
                                                     PMS_data_table_DreamCIS_Multi,
                                                     by = "subject_id")
    
    PMS_DREAM_slope_point_most_cut_short = PROACT_summarized_short_run(PMS_data_table_DreamCIS_shortage,
                                                                       PMS_data_table_DreamCIS_shortage_placebo,
                                                                       ALSFRS_cut = F,
                                                                       cp_type = "point_style",
                                                                       select_type = "most")
    
    PMS_DREAM_slope_point_most_cut = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage,
                                                           PMS_data_table_DreamCIS_shortage_placebo,
                                                           ALSFRS_cut = F,
                                                           cp_type = "point_style",
                                                           select_type = "most")
    
    PMS_data_table_DreamCIS_cytokine = PMS_ALSFRS_cytokine_table[c(2,13,14,22,23,31,32,40,41)]
    colnames(PMS_data_table_DreamCIS_cytokine)[1] = "subject_id"
    PMS_DREAM_slope_point_most_cut_short = merge(PMS_DREAM_slope_point_most_cut_short,
                                                 PMS_data_table_DreamCIS_cytokine,
                                                 by = "subject_id")
    
    PMS_DREAM_slope_point_most_cut_short$TGF_change = (PMS_DREAM_slope_point_most_cut_short$TGF_1st -
                                                           PMS_DREAM_slope_point_most_cut_short$TGF_Base) / PMS_DREAM_slope_point_most_cut_short$TGF_Base * 100
    PMS_DREAM_slope_point_most_cut_short$IL10_change = (PMS_DREAM_slope_point_most_cut_short$`IL-10_1st` -
                                                           PMS_DREAM_slope_point_most_cut_short$`IL-10_Base`) / PMS_DREAM_slope_point_most_cut_short$`IL-10_Base` * 100
    PMS_DREAM_slope_point_most_cut_short$IL8_change = (PMS_DREAM_slope_point_most_cut_short$`IL-8_1st` -
                                                           PMS_DREAM_slope_point_most_cut_short$`IL-8_Base`) / PMS_DREAM_slope_point_most_cut_short$`IL-8_Base` * 100
    PMS_DREAM_slope_point_most_cut_short$MCP1_change = (PMS_DREAM_slope_point_most_cut_short$`MCP-1_1st` -
                                                           PMS_DREAM_slope_point_most_cut_short$`MCP-1_Base`) / PMS_DREAM_slope_point_most_cut_short$`MCP-1_Base` * 100
    
    slope_model = lm(PMS_DREAM_slope_point_most_cut_short, formula =  M6_slope ~ TGF_change + IL10_change + MCP1_change)
    vif(slope_model)
    
    PMS_DREAM_slope_point_most_cut_short = subset(PMS_DREAM_slope_point_most_cut_short, is.na(M6_slope) == F)
    PMS_DREAM_slope_point_most_cut_short$TGF_change_subgroup = ifelse(PMS_DREAM_slope_point_most_cut_short$M6_slope >= -0.5, "Good response", "Bad response")
    PMS_DREAM_slope_point_most_cut_short$IL10_change_subgroup = ifelse(PMS_DREAM_slope_point_most_cut_short$M6_slope >= -0.5, "Good response", "Bad response")
    PMS_DREAM_slope_point_most_cut_short$IL8_change_subgroup = ifelse(PMS_DREAM_slope_point_most_cut_short$M6_slope >= -0.5, "Good response", "Bad response")
    PMS_DREAM_slope_point_most_cut_short$MCP1_change_subgroup = ifelse(PMS_DREAM_slope_point_most_cut_short$M6_slope >= -0.5, "Good response", "Bad response")
    
    PMS_DREAM_slope_point_most_cut_short$TGF_change_subgroup = factor(PMS_DREAM_slope_point_most_cut_short$TGF_change_subgroup, levels = c("Bad response","Good response"))
    PMS_DREAM_slope_point_most_cut_short$IL10_change_subgroup = factor(PMS_DREAM_slope_point_most_cut_short$IL10_change_subgroup, levels = c("Bad response","Good response"))
    PMS_DREAM_slope_point_most_cut_short$IL8_change_subgroup = factor(PMS_DREAM_slope_point_most_cut_short$IL8_change_subgroup, levels = c("Bad response","Good response"))
    PMS_DREAM_slope_point_most_cut_short$MCP1_change_subgroup = factor(PMS_DREAM_slope_point_most_cut_short$MCP1_change_subgroup, levels = c("Bad response","Good response"))
    
    
    PMS_DREAM_TGF_Wilcox = ggboxplot(data = PMS_DREAM_slope_point_most_cut_short,
                                 x = "TGF_change_subgroup",
                                 y = "TGF_change",
                                 color = "TGF_change_subgroup",
                                 add = "jitter",
                                 xlab = "Group (by response)",
                                 ylab = "TGF-B Change",
                                 add.params = list(size = 2.5),
                                 legend = "none") +
        stat_compare_means(method = "wilcox.test")
    
    PMS_DREAM_IL10_Wilcox = ggboxplot(data = PMS_DREAM_slope_point_most_cut_short,
                                     x = "IL10_change_subgroup",
                                     y = "IL10_change",
                                     color = "IL10_change_subgroup",
                                     add = "jitter",
                                     xlab = "Group (by response)",
                                     ylab = "IL-10 Change",
                                     add.params = list(size = 2.5),
                                     legend = "none") +
        stat_compare_means(method = "wilcox.test")
    
    PMS_DREAM_MCP1_Wilcox = ggboxplot(data = PMS_DREAM_slope_point_most_cut_short,
                                     x = "MCP1_change_subgroup",
                                     y = "MCP1_change",
                                     color = "MCP1_change_subgroup",
                                     add = "jitter",
                                     xlab = "Group (by response)",
                                     ylab = "MCP-1 Change",
                                     add.params = list(size = 2.5),
                                     legend = "none") +
        stat_compare_means(method = "wilcox.test")
    
    PMS_DREAM_IL8_Wilcox = ggboxplot(data = PMS_DREAM_slope_point_most_cut_short,
                                     x = "IL8_change_subgroup",
                                     y = "IL8_change",
                                     color = "IL8_change_subgroup",
                                     add = "jitter",
                                     xlab = "Group (by response)",
                                     ylab = "IL-8 Change",
                                     add.params = list(size = 2.5),
                                     legend = "none") +
        stat_compare_means(method = "wilcox.test")
    
    PMS_DREAM_total_Wilcox = ggarrange(PMS_DREAM_TGF_Wilcox, PMS_DREAM_IL10_Wilcox, PMS_DREAM_IL8_Wilcox, PMS_DREAM_MCP1_Wilcox)
    
    
    PMS_DREAM_total_Wilcox
    
    
    
    
    
    ggsave(filename = paste0("./Output/PMS_2M_cytokine_subgroup_wilcox.tiff"),
           plot = PMS_DREAM_total_Wilcox,
           device = "tiff",
           width = 12,
           height = 12,
           dpi = 300)
    
    
    return(tmp_table)
}