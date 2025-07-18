

# Load packages -------
library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)

# Import data -----



pubmep = read.csv2("./data/raw/BASE_PUBMEP_LONGI_WIDEformat_NIÑAS_Y_NIÑOS_213inds_31_01_2022_EPIC.csv")

pubmep2 = pubmep %>% 
  dplyr::filter(EPIC_ARRAY850K_T1_SI_o_NO == 1) %>% 
  dplyr::mutate(outcome = paste0(Cole_T1, HOMA_0_AUG_T1)) %>% 
  dplyr::filter(outcome != "Normopeso_T11") %>%
  tidyr::drop_na(Tanner2_T1) %>% 
  dplyr::mutate(outcome = case_when(
    outcome == "Normopeso_T10" ~ "1_NW_noIR", 
    outcome == "Sobrepeso_T10"  | outcome == "Obeso_T10" ~ "2_OW/OB_noIR", 
    outcome == "Sobrepeso_T11"  | outcome == "Obeso_T11" ~ "3_OW/OB_IR"))
  


# pubmep2$outcome %>% table(useNA = "ifany")
# 1_NW_noIR 2_OW/OB_noIR   3_OW/OB_IR 
# 29           45           25 


# Baseline ------

pubmep_t1 = pubmep2 %>% 
  dplyr::select(Code, Sex_T2, Tanner_T1, Origen_T1, Age_T1, BMI_T1, BMIz__Orbegozo__T1, 
                Perimetro_de_cintura_T1, Perimetro_de_cadera_T1, 
                Glucose__mg_dl__T1, Insulin__mU_l__T1, HOMA_T1, QUICKI_T1, 
                SBP_T1, DBP_T1, TAG__mg_dl__T1, HDLc__mg_dl__T1, 
                LDLc__mg_dl__T1, hsCRP_mg_l__T1, Adiponectin__mg_l__T1, 
                Leptin__ug_l__T1, Resistin__ug_l__T1, #Selectin__ug_l__T1, 
                TNF_ng_l__T1, IL8_ng_l__T1, MCP1_ng_l__T1, MPO_ug_l__T1, 
                tPAI1__ug_l__T1, sICAM1__mg_l__T1, outcome
  ) %>% 
  dplyr::mutate(Sex_T2 = if_else(Sex_T2 == "Niña_T2", 1, 0)) %>% 
  dplyr::rename(Age = Age_T1, Waist_circumference = Perimetro_de_cintura_T1, 
                SBP_mmHg = SBP_T1, DBP_mmHg = DBP_T1, 
                Triglycerides_mg_dL = TAG__mg_dl__T1, 
                HDLc_mg_dL = HDLc__mg_dl__T1, 
                Glucose_mg_dL = Glucose__mg_dl__T1, 
                HOMA_IR = HOMA_T1, LDLc_mg_dL = LDLc__mg_dl__T1, 
                Insulin_uU_ml = Insulin__mU_l__T1, Sex = Sex_T2, 
                Experimental_group = outcome)

source("./scripts/functions_referencevalues.R")

pubmep_t1 = pubmep_t1 %>% 
  calculate_wc_zscore(ref_wc_f = wc_ferrandez_score_girls, 
                      ref_wc_m = wc_ferrandez_score_boys) %>% 
  sbp_zscore_stravnsbo() %>% dbp_zscore_stravnsbo() %>%
  tag_zscore_stravnsbo() %>% hdl_zscore_stravnsbo() %>%
  ldl_zscore_stravnsbo() %>% glu_zscore_stravnsbo() %>% 
  insulin_zscore_stravnsbo() %>% homa_zscore_stravnsbo()


pubmep_t1 = pubmep_t1 %>% 
  dplyr::select(Code, Sex:Waist_circumference, WC_zscore, Perimetro_de_cadera_T1,  
                Glucose_mg_dL, Glucose_zscore, Insulin_uU_ml, Insulin_zscore, 
                HOMA_IR, HOMA_IR_zscore, QUICKI_T1, SBP_mmHg, SBP_zscore_Stavnsbo, 
                DBP_mmHg, DBP_zscore_Stavnsbo, Triglycerides_mg_dL, TAG_zscore, 
                HDLc_mg_dL, HDLc_zscore, LDLc_mg_dL, LDLc_zscore, 
                hsCRP_mg_l__T1:sICAM1__mg_l__T1, Experimental_group)

fmi = readRDS("/home/usuario/Escritorio/CIBM/Multiomics_VASN/data/processed/EWAS_pubmep_t1.RDS") %>% 
  dplyr::select(Code, FMI)

pubmep_t1 = pubmep_t1 %>% dplyr::left_join(fmi, by = "Code")

## Statistical testing ----

vars_continuas = pubmep_t1 %>% 
  dplyr::select(Age:sICAM1__mg_l__T1, FMI) %>% 
  colnames()


tests_normalidad = sapply(vars_continuas, function(v) {
  x = pubmep_t1[[v]]
  x1 = pubmep_t1[[v]][pubmep_t1$Experimental_group == "1_NW_noIR"] %>% 
    na.omit()
  x2 = pubmep_t1[[v]][pubmep_t1$Experimental_group == "2_OW/OB_noIR"] %>% 
    na.omit()
  x3 = pubmep_t1[[v]][pubmep_t1$Experimental_group == "3_OW/OB_IR"] %>%
    na.omit()
  n1 = length(x1)
  n2 = length(x2)
  n3 = length(x3)
  print(v)
  if(n1 < 50 | n2 < 50 | n3 < 50){
    print("Shapiro-Wilk test")
    p1 = shapiro.test(x1)$p.value
    p2 = shapiro.test(x2)$p.value
    p3 = shapiro.test(x3)$p.value
  } else{
    print("Anderson-Darling test") # Because Kolmogorov-Smirnov Tests gave us problems due to ties
    p1 = nortest::ad.test(x1)$p.value
    p2 = nortest::ad.test(x2)$p.value
    p3 = nortest::ad.test(x3)$p.value
  }
  print(p1)
  print(p2)
  print(p3)
  # h0: normality
  # h0: no normality
  # TRUE no-parametric
  # FALSE parametric 
  any(c(p1, p2, p3) < 0.05)
})


parametrics = vars_continuas[tests_normalidad == FALSE]
non_parametrics = vars_continuas[tests_normalidad == TRUE]

## Descriptive table baseline  -----
pubmep_t1 %>% 
  dplyr::select(Sex:FMI) %>% 
  dplyr::mutate(Sex = as.factor(Sex), 
                Tanner_T1 = as.factor(Tanner_T1), 
                Origen_T1 = as.factor(Origen_T1)
  ) %>% 
  tbl_summary(by = Experimental_group, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{mean} ({sd})"), #, "{median} ({p25}, {p75})"
              digits = list(
                all_continuous() ~ 2,
                all_categorical() ~ c(0,2)
              ), 
              missing_text = "Missing values"
  ) %>%
  modify_table_styling(
    columns = "label",
    rows = row_type == "label",   
    text_format = "bold"
  ) %>% 
  add_p(
    test = list(
      all_of(parametrics) ~ "oneway.test", 
      all_of(non_parametrics) ~ "kruskal.test", 
      Origen_T1 ~ "fisher.test"
    ), 
    test.args = list(oneway.test = list(var.equal = FALSE), 
                     "fisher.test" = list(workspace = 2e7, simulate.p.value = TRUE)
    ),
    pvalue_fun = label_style_pvalue(digits = 2)
  )

# Information about the statistical tests
# help("add_p")
# help("add_p.tbl_summary")
# help("tests")

# For Two Groups
# Continuous Data
# 
# - t-test: For normal data. Compares the means of two groups.
# - Wilcoxon rank-sum test: For non-normal or heterogeneous data. Compares distributions between two groups.
# - Paired t-test: For paired (before/after) data that follows a normal distribution.
# - Paired Wilcoxon test: For paired (before/after) data that does not follow a normal distribution.
# 
# Categorical Data
# 
# - Chi-square test of independence: Compares proportions between two groups.
# - Chi-square test without correction: Same as above but without correction, suitable for large samples.
# - Fisher's exact test: For small contingency tables (frequencies < 5).
# - McNemar's test: For paired or dependent data (e.g., before and after a treatment).
# 
# For More Than Two Groups
# Continuous Data
# 
# - One-way ANOVA (aov): For normal data with homogeneous variances. Compares means across multiple groups.
# - One-way ANOVA (oneway.test): For normal data with heterogeneous variances.
# - Kruskal-Wallis test: For non-normal or heterogeneous data. Compares distributions across multiple groups.
# - Mood's test: Evaluates differences in dispersion across multiple groups.
# 
# Categorical Data
# 
# - Chi-square test of independence: Compares proportions across multiple groups.
# - Fisher's exact test: For small contingency tables.
# - Proportions test (prop.test): Compares proportions across multiple groups.
# 
# Data Adjusted for Covariates
# 
# - ANCOVA: For continuous data, adjusts comparisons for covariates.
# - Random intercept logistic regression (lme4): For categorical data with hierarchical or nested structures.
