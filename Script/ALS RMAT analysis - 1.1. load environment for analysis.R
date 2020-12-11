# Corestem project - ALS RMAT Analysis.proj
# 20201116 Sehwan Chun at Corestem, Inc.
# 1.1. load environment for analysis

#### 1. Library Loading ####
packs = c("data.table", "readxl", "ggpubr", "writexl")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. Files Loading ####
ALSFRSR_file = "./Data/Rawdata/Allo 및 임상 12상 통합자료.xlsx"
ALSFRSR_file_allo_P2cont = read_xlsx(ALSFRSR_file, sheet = 3)
ALSFRSR_file_allo_P12 = read_xlsx(ALSFRSR_file, sheet = 4)
ALSFRSR_file_allo_P12_FCSPMM = read_xlsx(ALSFRSR_file, sheet = 5)

PMS_file = "./Data/Rawdata/PMS_V2.xlsx"
PMS_Phase12_change_slope = read_xlsx(PMS_file, sheet = 13)
PMS_ALSFRS_init_table = read_xlsx(PMS_file, sheet = 14)
PMS_ALSFRS_table = read_xlsx(PMS_file, sheet = 15)
PMS_ALSFRS_cytokine_table = read_xlsx(PMS_file, sheet = 16)

PROACT_file = "./Data/Rawdata/PROACT.xlsx"
PROACT_ALSFRS = read_xlsx(PROACT_file, sheet = 1)
PROACT_Placebo = read_xlsx(PROACT_file, sheet = 2)

save.image(file = "./Data/ALS RMAT analysis files.image")
