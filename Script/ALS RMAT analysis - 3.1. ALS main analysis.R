# Corestem project - FACS analysis.proj
# 20201113 Sehwan Chun at Corestem, Inc.
# 3.1. ALS AE SAE main analysis

#### 1. source Loading ####
load("./Data/ALS RMAT analysis files.image")
source("./Script/Functions/ALS RMAT analysis - functions.R")

#### 2. simple LOCF running ####
ALSFRSR_file_allo_P12_LOCF = LOCF_table_run(ALSFRSR_file_allo_P12)
ALSFRSR_file_allo_P2cont_LOCF = LOCF_table_run(ALSFRSR_file_allo_P2cont)

stat_allo_P12 = stat_table_run(ALSFRSR_file_allo_P12)
stat_allo_P2cont = stat_table_run(ALSFRSR_file_allo_P2cont)
stat_allo_P12_LOCF = stat_table_run(ALSFRSR_file_allo_P12_LOCF)
stat_allo_P2cont_LOCF = stat_table_run(ALSFRSR_file_allo_P2cont_LOCF)
stat_allo_P12_FCSPMM = stat_table_run(ALSFRSR_file_allo_P12_FCSPMM)


stat_bind = cbind(stat_allo_P12,
                  stat_allo_P12_LOCF$pvalue,
                  stat_allo_P2cont$pvalue,
                  stat_allo_P2cont_LOCF$pvalue,
                  stat_allo_P12_FCSPMM$pvalue)
colnames(stat_bind) = c("TEST",
                        "P12_Allo",
                        "P12_Allo_LOCF",
                        "P2CONT_Allo",
                        "P2CONT_Allo_LOCF",
                        "P12_FCSPMM")

