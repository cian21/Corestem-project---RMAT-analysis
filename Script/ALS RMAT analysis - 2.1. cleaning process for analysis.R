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
PMS_data_table_DreamCIS = PMS_slope_comparison_run(PMS_ALSFRS_init_table, "DreamCIS")
PMS_data_table_DreamCIS_shortage = PMS_shortage_run(PMS_data_table_DreamCIS, "DreamCIS")
PMS_data_table_DreamCIS_shortage_placebo = data.frame(subject_id = unique(PMS_data_table_DreamCIS_shortage$subject_id),
                                                      Study_Arm = "Placebo")

PMS_data_table_DreamCIS_Multi = PMS_data_table_DreamCIS[c(2,39,41)]
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
                                                          ALSFRS_cut = T,
                                                          cp_type = "point_style",
                                                          select_type = "most")

PMS_DREAM_slope_point_most_cut = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage,
                                                             PMS_data_table_DreamCIS_shortage_placebo,
                                                             ALSFRS_cut = T,
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
                                                        ALSFRS_cut = T,
                                                        cp_type = "point_style",
                                                        select_type = "most")

PMS_DreamCIS_slope_point_most_cut_M6_Multi = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M6_Multi,
                                                                   PMS_data_table_DreamCIS_shortage_placebo_M6_Multi,
                                                        ALSFRS_cut = T,
                                                        cp_type = "point_style",
                                                        select_type = "most")

PMS_DreamCIS_slope_point_most_cut_M12_Normal = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M12_Normal,
                                                           PMS_data_table_DreamCIS_shortage_placebo_M12_Normal,
                                                           ALSFRS_cut = T,
                                                           cp_type = "point_style",
                                                           select_type = "most")

PMS_DreamCIS_slope_point_most_cut_M12_Multi = PROACT_summarized_run(PMS_data_table_DreamCIS_shortage_M12_Multi,
                                                           PMS_data_table_DreamCIS_shortage_placebo_M12_Multi,
                                                           ALSFRS_cut = T,
                                                           cp_type = "point_style",
                                                           select_type = "most")

wilcox.test(PMS_DREAM_slope_point_most_cut_short$M6_slope,
            ALSFRSR_file_allo_P2cont$M4_Changes)

wilcox.test(PMS_DREAM_slope_point_most_cut$M6_slope,
            ALSFRSR_file_allo_P2cont$M6_Changes)

wilcox.test(PMS_DREAM_slope_point_most_cut$M12_slope,
            ALSFRSR_file_allo_P2cont$M12_Changes)

wilcox.test(PMS_DreamCIS_slope_point_most_cut_M6_Normal$M6_slope,
            ALSFRSR_file_allo_P2cont$M6_Changes)

wilcox.test(PMS_DreamCIS_slope_point_most_cut_M6_Multi$M6_slope,
            ALSFRSR_file_allo_P2cont$M6_Changes)

wilcox.test(PMS_DreamCIS_slope_point_most_cut_M6_Normal$M6_slope,
            PMS_DreamCIS_slope_point_most_cut_M6_Multi$M6_slope)

wilcox.test(PMS_DreamCIS_slope_point_most_cut_M12_Normal$M12_slope,
            ALSFRSR_file_allo_P2cont$M12_Changes)

wilcox.test(PMS_DreamCIS_slope_point_most_cut_M12_Multi$M12_slope,
            ALSFRSR_file_allo_P2cont$M12_Changes)

wilcox.test(c(PMS_DreamCIS_slope_point_most_cut_M12_Normal$M12_slope,
              PMS_DreamCIS_slope_point_most_cut_M12_Multi$M12_slope),
            ALSFRSR_file_allo_P2cont$M12_Changes)

wilcox.test(PMS_DreamCIS_slope_point_most_cut_M12_Multi$M12_slope,
       PMS_DreamCIS_slope_point_most_cut_M12_Normal$M12_slope)

#### 4. PROACT ####
PROACT_Placebo = subset(PROACT_Placebo, Study_Arm == "Placebo")
PROACT_Placebo = subset(PROACT_Placebo, Treatment_Group_Delta == 0)

PROACT_ALSFRS = subset(PROACT_ALSFRS, subject_id %in% PROACT_Placebo$subject_id)
PROACT_ALSFRS = PROACT_ALSFRS[,c(1,13,15)]
PROACT_ALSFRS = PROACT_ALSFRS[complete.cases(PROACT_ALSFRS),]
PROACT_stat_table = PROACT_stat_run(PROACT_ALSFRS, PROACT_Placebo)
PROACT_stat_short_table = PROACT_stat_short_run(PROACT_ALSFRS, PROACT_Placebo)