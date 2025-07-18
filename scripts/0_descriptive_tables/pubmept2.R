
# Load packages -------
library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)

# Import data -----

pubmep = read.csv2("./data/raw/BASE_PUBMEP_LONGI_WIDEformat_NIÑAS_Y_NIÑOS_213inds_31_01_2022_EPIC.csv")

pubmep2 = pubmep %>% 
  dplyr::filter(EPIC_ARRAY850K_T2_SI_o_NO == 1) %>% 
  dplyr::mutate(outcome = paste0(Cole_T2, HOMA_0_AUG_T2)) %>% 
  dplyr::filter(outcome != "Normopeso_T21") %>%
  dplyr::filter(outcome != "Normopeso_T2NA") %>% 
  dplyr::filter(outcome != "Sobrepeso_T2NA") %>% 
  dplyr::filter(Tanner_T2 != 1) %>% 
  dplyr::mutate(outcome = case_when(
    outcome == "Normopeso_T20" ~ "1_NW_noIR", 
    outcome == "Sobrepeso_T20"  | outcome == "Obeso_T20" ~ "2_OW/OB_noIR", 
    outcome == "Sobrepeso_T21"  | outcome == "Obeso_T21" ~ "3_OW/OB_IR"))

# pubmep2$outcome %>% table(useNA = "ifany")
# 1_NW_noIR 2_OW/OB_noIR   3_OW/OB_IR 
# 47           51           31 


# Follow-up ------

pubmep_t2 = pubmep2 %>% 
  dplyr::select(Code_new_T2, Sex_T2, Tanner_T2, Origen_T1, Age_T2, BMI_T2, BMIz__Orbegozo__T2, 
                Perimetro_de_cintura_T2, Perimetro_de_cadera_T2, 
                Glucose_mgdl_T2, Insulin_mUl_T2, HOMA_IR_T2, QUICKI_T2, 
                SBP_T2, DBP_T2, TAG_mgdl_T2, HDLc_mgdl_T2, 
                LDLc_mgdl_T2, hsCRPmgl_T2, Adiponectin__mg_l__T2, 
                Leptin__ug_l__T2, Resistin__ug_l__T2, #Selectin__ug_l__T2, 
                TNF_ng_l__T2, IL8_ng_l__T2, MCP1_ng_l__T2, MPO_ug_l__T2, 
                tPAI1__ug_l__T2, sICAM1__mg_l__T2, outcome
  ) %>% 
  dplyr::mutate(Sex_T2 = if_else(Sex_T2 == "Niña_T2", 1, 0)
  ) %>% 
  dplyr::rename(Age = Age_T2, Waist_circumference = Perimetro_de_cintura_T2, 
                SBP_mmHg = SBP_T2, DBP_mmHg = DBP_T2, 
                Triglycerides_mg_dL = TAG_mgdl_T2, 
                HDLc_mg_dL = HDLc_mgdl_T2, 
                Glucose_mg_dL = Glucose_mgdl_T2, 
                HOMA_IR = HOMA_IR_T2, LDLc_mg_dL = LDLc_mgdl_T2, 
                Insulin_uU_ml = Insulin_mUl_T2, Sex = Sex_T2, 
                Experimental_group = outcome)

source("./scripts/functions_referencevalues.R")

pubmep_t2 = pubmep_t2 %>% 
  calculate_wc_zscore(ref_wc_f = wc_ferrandez_score_girls, 
                      ref_wc_m = wc_ferrandez_score_boys) %>% 
  sbp_zscore_stravnsbo() %>% dbp_zscore_stravnsbo() %>%
  tag_zscore_stravnsbo() %>% hdl_zscore_stravnsbo() %>%
  ldl_zscore_stravnsbo() %>% glu_zscore_stravnsbo() %>% 
  insulin_zscore_stravnsbo() %>% homa_zscore_stravnsbo()


pubmep_t2 = pubmep_t2 %>% 
  dplyr::select(Code_new_T2, Sex:Waist_circumference, WC_zscore, Perimetro_de_cadera_T2,  
                Glucose_mg_dL, Glucose_zscore, Insulin_uU_ml, Insulin_zscore, 
                HOMA_IR, HOMA_IR_zscore, QUICKI_T2, SBP_mmHg, SBP_zscore_Stavnsbo, 
                DBP_mmHg, DBP_zscore_Stavnsbo, Triglycerides_mg_dL, TAG_zscore, 
                HDLc_mg_dL, HDLc_zscore, LDLc_mg_dL, LDLc_zscore, 
                hsCRPmgl_T2:sICAM1__mg_l__T2, Experimental_group)

fmi = readRDS("/home/usuario/Escritorio/CIBM/Multiomics_VASN/data/processed/EWAS_pubmep_t2.RDS") %>% 
  dplyr::select(Code_new_T2, FMI)

pubmep_t2 = pubmep_t2 %>% dplyr::left_join(fmi, by = "Code_new_T2")

## Statistical testing ----

vars_continuas = pubmep_t2 %>% 
  dplyr::select(Age:sICAM1__mg_l__T2, FMI) %>% 
  colnames()


tests_normalidad = sapply(vars_continuas, function(v) {
  x = pubmep_t2[[v]]
  x1 = pubmep_t2[[v]][pubmep_t2$Experimental_group == "1_NW_noIR"] %>% 
    na.omit()
  x2 = pubmep_t2[[v]][pubmep_t2$Experimental_group == "2_OW/OB_noIR"] %>% 
    na.omit()
  x3 = pubmep_t2[[v]][pubmep_t2$Experimental_group == "3_OW/OB_IR"] %>%
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
pubmep_t2 %>% 
  dplyr::mutate(hsCRPmgl_T2 = if_else(hsCRPmgl_T2 > 100, NA, hsCRPmgl_T2)) %>% 
  dplyr::select(Sex:FMI) %>% 
  dplyr::mutate(Sex = as.factor(Sex), 
                Tanner_T2  = as.factor(Tanner_T2), 
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
