#####
#Libraries loading
library(readr)
library(tidyverse)
library(gtsummary)
#####
#upload data
candidate63 <- read_delim("CopomuS-ranks.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

#paired t.test on differences on average, to assess significant RRI stability changes
set <- candidate63 %>% select(mutation, E_ww, E_wm, E_mw, E_mm)
diff_wm <- set$E_wm - set$E_ww
diff_mw <- set$E_mw - set$E_ww
diff_mm <- set$E_mm - set$E_ww
diff_wm_mm <- set$E_wm - set$E_mm
diff_mw_mm <- set$E_mw - set$E_mm

t.test(diff_wm, mu = 0, alternative = "two.sided")
t.test(diff_mw, mu = 0, alternative = "two.sided")
t.test(diff_mm, mu = 0, alternative = "two.sided")
t.test(diff_wm_mm, mu = 0, alternative = "two.sided")
t.test(diff_mw_mm, mu = 0, alternative = "two.sided")

#quick descriptive on mutations predicted in CoM analysis
table63 <- tbl_summary(
  candidate63,
  by = bpMutated,
  type = all_dichotomous() ~ "categorical"
) %>% 
  modify_header(label = "**bpMutated**") %>%
  bold_labels() %>% 
  italicize_labels()

table63 %>% as_flex_table() %>%
  flextable::save_as_docx(path="table63.docx") #produco un .docx

#####
#upload data
candidate899<- read_delim("cand899_mut.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)

#paired t.test on differences on average, to assess significant RRI stability changes
set <- candidate899 %>% select(mutation, E_ww, E_wm, E_mw, E_mm)
diff_wm <- set$E_wm - set$E_ww
diff_mw <- set$E_mw - set$E_ww
diff_mm <- set$E_mm - set$E_ww
diff_wm_mm <- set$E_wm - set$E_mm
diff_mw_mm <- set$E_mw - set$E_mm

t.test(diff_wm, mu = 0, alternative = "two.sided")
t.test(diff_mw, mu = 0, alternative = "two.sided")
t.test(diff_mm, mu = 0, alternative = "two.sided")
t.test(diff_wm_mm, mu = 0, alternative = "two.sided")
t.test(diff_mw_mm, mu = 0, alternative = "two.sided")

#quick descriptive on mutations predicted in CoM analysis
table899 <- tbl_summary(
  candidate899,
  by = bpMutated,
  type = all_dichotomous() ~ "categorical"
) %>% 
  modify_header(label = "**bpMutated**") %>%
  bold_labels() %>% 
  italicize_labels()

table899 %>% as_flex_table() %>%
  flextable::save_as_docx(path="table899.docx") #produco un .docx

#####
#upload data
candidate1405 <- read_delim("cand1405_Rv0005_mut.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

#paired t.test on differences on average, to assess significant RRI stability changes
set <- candidate1405 %>% select(mutation, E_ww, E_wm, E_mw, E_mm)
diff_wm <- set$E_wm - set$E_ww
diff_mw <- set$E_mw - set$E_ww
diff_mm <- set$E_mm - set$E_ww
diff_wm_mm <- set$E_wm - set$E_mm
diff_mw_mm <- set$E_mw - set$E_mm

t.test(diff_wm, mu = 0, alternative = "two.sided")
t.test(diff_mw, mu = 0, alternative = "two.sided")
t.test(diff_mm, mu = 0, alternative = "two.sided")
t.test(diff_wm_mm, mu = 0, alternative = "two.sided")
t.test(diff_mw_mm, mu = 0, alternative = "two.sided")

#quick descriptive on mutations predicted in CoM analysis
table1405 <- tbl_summary(
  candidate1405,
  by = bpMutated,
  type = all_dichotomous() ~ "categorical"
) %>% 
  modify_header(label = "**bpMutated**") %>%
  bold_labels() %>% 
  italicize_labels()

table1405 %>% as_flex_table() %>%
  flextable::save_as_docx(path="table1405.docx") #produco un .docx

#####
#upload data
candidate2023 <- read_delim("cand2023_Rv1111c_mut.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

#paired t.test on differences on average, to assess significant RRI stability changes
set <- candidate2023 %>% select(mutation, E_ww, E_wm, E_mw, E_mm)
diff_wm <- set$E_wm - set$E_ww
diff_mw <- set$E_mw - set$E_ww
diff_mm <- set$E_mm - set$E_ww
diff_wm_mm <- set$E_wm - set$E_mm
diff_mw_mm <- set$E_mw - set$E_mm

t.test(diff_wm, mu = 0, alternative = "two.sided")
t.test(diff_mw, mu = 0, alternative = "two.sided")
t.test(diff_mm, mu = 0, alternative = "two.sided")
t.test(diff_wm_mm, mu = 0, alternative = "two.sided")
t.test(diff_mw_mm, mu = 0, alternative = "two.sided")

#quick descriptive on mutations predicted in CoM analysis
table2023 <- tbl_summary(
  candidate2023,
  by = bpMutated,
  type = all_dichotomous() ~ "categorical"
) %>% 
  modify_header(label = "**bpMutated**") %>%
  bold_labels() %>% 
  italicize_labels()

table2023 %>% as_flex_table() %>%
  flextable::save_as_docx(path="table2023.docx") #produco un .docx

#####
#upload data
candidate2223 <- read_delim("cand2223_Rv1018c_mut.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

#paired t.test on differences on average, to assess significant RRI stability changes
set <- candidate2223 %>% select(mutation, E_ww, E_wm, E_mw, E_mm)
diff_wm <- set$E_wm - set$E_ww
diff_mw <- set$E_mw - set$E_ww
diff_mm <- set$E_mm - set$E_ww
diff_wm_mm <- set$E_wm - set$E_mm
diff_mw_mm <- set$E_mw - set$E_mm

t.test(diff_wm, mu = 0, alternative = "two.sided")
t.test(diff_mw, mu = 0, alternative = "two.sided")
t.test(diff_mm, mu = 0, alternative = "two.sided")
t.test(diff_wm_mm, mu = 0, alternative = "two.sided")
t.test(diff_mw_mm, mu = 0, alternative = "two.sided")

#quick descriptive on mutations predicted in CoM analysis
table2223 <- tbl_summary(
  candidate2223,
  by = bpMutated,
  type = all_dichotomous() ~ "categorical"
) %>% 
  modify_header(label = "**bpMutated**") %>%
  bold_labels() %>% 
  italicize_labels()

table2223 %>% as_flex_table() %>%
  flextable::save_as_docx(path="table2223.docx") #produco un .docx
