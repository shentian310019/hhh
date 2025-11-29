## 00A_nhanes_T2D_HF_screening.R
## 目的：在 NHANES T2D 队列中，对候选心代谢/炎症/肾功能等指标
##       做单变量 Cox 回归（统一调整年龄/性别/种族），
##       输出 HR 表，并画火山图，用来证明“10 条轴”不是拍脑袋选的。

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(stringr)
  library(purrr)
  library(tibble)
})

#------------------------------------------------------------
# 0. 项目路径与数据读取
#------------------------------------------------------------

project_root <- getwd()
setwd(project_root)

nhanes_path <- file.path(project_root, "data", "nhanes_t2d_hf_cohort.rds")
if (!file.exists(nhanes_path)) {
  stop("找不到数据文件：", nhanes_path,
       "\n请确认 data 目录下已有 nhanes_t2d_hf_cohort.rds")
}

dat0 <- readRDS(nhanes_path)
message("数据维度：", paste(dim(dat0), collapse = " x "))
message("示例列名：", paste(head(names(dat0), 30), collapse = ", "))

#------------------------------------------------------------
# 1. 辅助函数：根据别名匹配变量名
#------------------------------------------------------------

find_first_var <- function(candidates, data_names) {
  idx <- match(tolower(candidates), tolower(data_names))
  hit <- idx[!is.na(idx)][1]
  if (is.na(hit)) return(NA_character_)
  data_names[hit]
}

#------------------------------------------------------------
# 2. 统一关键变量（随访时间、结局、协变量）
#------------------------------------------------------------

key_vars <- list(
  t2d_flag = c("t2d_flag", "t2d", "diabetes_flag", "t2dm"),
  time     = c("follow_up_time", "time", "time_to_event", "surv_time", "followup_time", "person_years"),
  hf_event = c("hf_event", "hf", "hf_incident", "hf_outcome", "heart_failure", "hf_event_flag"),
  age      = c("age", "age_baseline", "baseline_age", "ridageyr"),
  sex      = c("sex", "gender", "riagendr"),
  race     = c("race", "ethnicity", "race_ethnicity", "race_eth")
)

resolved_keys <- vapply(key_vars, find_first_var, character(1), data_names = names(dat0))
if (anyNA(resolved_keys)) {
  missing_keys <- names(resolved_keys)[is.na(resolved_keys)]
  stop("以下关键变量无法在数据中找到，请检查：", paste(missing_keys, collapse = ", "))
}

# 统一命名，方便后续公式书写

dat <- dat0 %>%
  rename(
    t2d_flag = !!resolved_keys["t2d_flag"],
    time     = !!resolved_keys["time"],
    hf_event = !!resolved_keys["hf_event"],
    age      = !!resolved_keys["age"],
    sex      = !!resolved_keys["sex"],
    race     = !!resolved_keys["race"]
  )

# 只保留 T2D 患者，并去掉关键变量缺失

dat_t2d <- dat %>%
  filter(t2d_flag == 1) %>%
  filter(!is.na(time), !is.na(hf_event))

# 将性别、种族转为因子，避免 Cox 将其视为连续变量

if (!is.factor(dat_t2d$sex)) dat_t2d <- dat_t2d %>% mutate(sex = factor(sex))
if (!is.factor(dat_t2d$race)) dat_t2d <- dat_t2d %>% mutate(race = factor(race))

#------------------------------------------------------------
# 3. 定义候选指标及其分组/别名
#------------------------------------------------------------

candidate_vars <- tibble(
  label = c(
    # 代谢
    "BMI", "WC", "WHR", "HbA1c", "FPG",
    # 血压
    "SBP", "DBP", "PP",
    # 血脂
    "LDL", "HDL", "TG", "TC",
    # 肾功能
    "eGFR", "UACR", "urine_protein",
    # 炎症
    "CRP",
    # 血细胞
    "WBC", "NEU", "MONO", "LYM", "PLT",
    # 既往病史
    "prev_MI", "prev_AF", "prev_stroke"
  ),
  group = c(
    rep("Metabolic", 5),
    rep("BloodPressure", 3),
    rep("Lipid", 4),
    rep("Renal", 3),
    rep("Inflammation", 1),
    rep("BloodCell", 5),
    rep("CVD_history", 3)
  ),
  candidates = list(
    # 代谢
    c("BMI", "bmi", "BodyMassIndex"),
    c("WC", "waist_circumference", "waist"),
    c("WHR", "waist_hip_ratio", "whr"),
    c("HbA1c", "hba1c", "A1c"),
    c("FPG", "FBS", "fasting_glucose", "FPG_mmol", "fasting_glucose_mmol"),
    # 血压
    c("SBP", "sbp", "systolic_bp"),
    c("DBP", "dbp", "diastolic_bp"),
    c("PP", "pulse_pressure"),
    # 血脂
    c("LDL", "ldl", "ldl_cholesterol"),
    c("HDL", "hdl", "hdl_cholesterol"),
    c("TG", "triglyceride", "triglycerides"),
    c("TC", "tc", "total_cholesterol"),
    # 肾功能
    c("eGFR", "egfr", "eGFR_epi"),
    c("UACR", "uacr", "albumin_creatinine_ratio"),
    c("urine_protein", "proteinuria", "urine_protein_flag"),
    # 炎症
    c("CRP", "crp", "hsCRP", "hs_crp"),
    # 血细胞
    c("WBC", "wbc", "white_cell_count"),
    c("NEU", "neutrophil", "neutrophils"),
    c("MONO", "monocyte", "monocytes"),
    c("LYM", "lymphocyte", "lymphocytes"),
    c("PLT", "platelet", "platelets", "plt_count"),
    # 既往病史
    c("prev_MI", "prior_mi", "history_mi"),
    c("prev_AF", "prior_af", "history_af"),
    c("prev_stroke", "prior_stroke", "history_stroke")
  )
)

candidate_vars <- candidate_vars %>%
  mutate(
    data_var = map_chr(candidates, find_first_var, data_names = names(dat_t2d))
  )

missing_vars <- candidate_vars %>%
  filter(is.na(data_var))

if (nrow(missing_vars) > 0) {
  message("以下变量在数据中找不到，将在 Cox 中跳过：")
  print(missing_vars %>% select(label, group))
}

candidate_vars <- candidate_vars %>%
  filter(!is.na(data_var))

#------------------------------------------------------------
# 4. 定义单变量 Cox 函数
#------------------------------------------------------------

run_cox_one <- function(label, data_var, group_label, data) {
  if (!data_var %in% names(data)) {
    message("变量 ", label, " 不存在，跳过")
    return(NULL)
  }

  # 尽量把暴露变量转为数值型，便于 scale()
  if (!is.numeric(data[[data_var]])) {
    suppressWarnings({
      data[[data_var]] <- as.numeric(as.character(data[[data_var]]))
    })
  }

  if (!is.numeric(data[[data_var]])) {
    message("变量 ", label, " 类型不合法，无法转为数值，跳过")
    return(NULL)
  }

  data_for_fit <- data %>%
    select(time, hf_event, age, sex, race, all_of(data_var)) %>%
    filter(!is.na(.data[[data_var]]))

  if (nrow(data_for_fit) == 0 || all(is.na(data_for_fit[[data_var]]))) {
    message("变量 ", label, " 全部缺失，跳过")
    return(NULL)
  }

  # 构建公式：Surv(time, hf_event) ~ scale(biomarker) + age + sex + race
  form <- as.formula(
    sprintf("Surv(time, hf_event) ~ scale(`%s`) + age + sex + race", data_var)
  )

  fit <- tryCatch(
    survival::coxph(form, data = data_for_fit),
    error = function(e) {
      message("Cox 拟合失败：", label, "，原因：", e$message)
      return(NULL)
    }
  )
  if (is.null(fit)) return(NULL)

  s <- summary(fit)
  coef_row <- sprintf("scale(%s)", data_var)
  if (!coef_row %in% rownames(s$coefficients)) {
    message("在系数表中没有找到变量行：", coef_row)
    return(NULL)
  }

  ctab <- s$coefficients[coef_row, , drop = FALSE]
  ci    <- s$conf.int[coef_row, , drop = FALSE]

  tibble(
    var_name = label,
    data_var = data_var,
    group    = group_label,
    logHR    = ctab[1, "coef"],
    HR       = exp(ctab[1, "coef"]),
    HR_low   = ci[1, "lower .95"],
    HR_high  = ci[1, "upper .95"],
    se       = ctab[1, "se(coef)"],
    z        = ctab[1, "z"],
    p        = ctab[1, "Pr(>|z|)"],
    n        = nrow(data_for_fit)
  )
}

#------------------------------------------------------------
# 5. 循环跑 Cox
#------------------------------------------------------------

group_lookup <- setNames(candidate_vars$group, candidate_vars$label)

res_list <- map2(
  candidate_vars$label,
  candidate_vars$data_var,
  ~ run_cox_one(
    label = .x,
    data_var = .y,
    group_label = group_lookup[[.x]],
    data = dat_t2d
  )
)

res_list <- res_list[!vapply(res_list, is.null, logical(1))]
res_df <- bind_rows(res_list)

if (nrow(res_df) < 5) {
  warning("成功拟合的变量少于 5 个，请检查输入数据和变量名。")
}

# 保存结果表
out_dir <- file.path(project_root, "results", "nhanes_screening")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(res_df, file.path(out_dir, "nhanes_t2d_hf_single_cox_results.csv"))

#------------------------------------------------------------
# 6. 火山图
#------------------------------------------------------------

res_df <- res_df %>%
  mutate(
    label = var_name,
    neg_log10_p = -log10(p),
    group = fct_inorder(group)
  )

top_n <- 15

top_vars <- res_df %>%
  arrange(p) %>%
  slice_head(n = top_n) %>%
  pull(var_name)

p_volcano <- ggplot(res_df, aes(x = logHR, y = neg_log10_p, color = group)) +
  geom_point(alpha = 0.85, size = 2.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  ggrepel::geom_text_repel(
    data = subset(res_df, var_name %in% top_vars),
    aes(label = label),
    size = 3,
    max.overlaps = 50
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "log(HR) for incident HF per 1 SD increase (adjusted for age/sex/race)",
    y = "-log10(p value)",
    color = "Group",
    title = "Single-variable Cox screening in T2D (NHANES)",
    subtitle = "Higher log(HR) and -log10(p) indicate stronger association with incident HF"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(out_dir, "volcano_T2D_HF_single_Cox.png"),
  plot = p_volcano,
  width = 8, height = 6, dpi = 300
)

#------------------------------------------------------------
# 7. 每个 group 选 1-2 个代表变量，绘制森林图
#------------------------------------------------------------

top_by_group <- res_df %>%
  group_by(group) %>%
  arrange(p, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

forest_vars <- unique(top_by_group$data_var)

if (length(forest_vars) == 0) {
  message("没有可用于森林图的变量，跳过森林图绘制。")
} else {
  forest_formula <- as.formula(
    paste(
      "Surv(time, hf_event) ~",
      paste(sprintf("scale(`%s`)", forest_vars), collapse = " + "),
      "+ age + sex + race"
    )
  )

  forest_data <- dat_t2d %>%
    select(time, hf_event, age, sex, race, all_of(forest_vars)) %>%
    tidyr::drop_na()

  forest_fit <- tryCatch(
    coxph(forest_formula, data = forest_data),
    error = function(e) {
      message("森林图的多变量 Cox 拟合失败：", e$message)
      return(NULL)
    }
  )

  if (!is.null(forest_fit)) {
    forest_summary <- summary(forest_fit)
    coef_df <- as.data.frame(forest_summary$coefficients)
    ci_df   <- as.data.frame(forest_summary$conf.int)
    coef_df$term <- rownames(coef_df)
    ci_df$term   <- rownames(ci_df)

    forest_res <- coef_df %>%
      left_join(ci_df %>% select(term, `lower .95`, `upper .95`), by = "term") %>%
      filter(str_detect(term, "^scale\\(")) %>%
      mutate(
        data_var = str_remove(term, "^scale\\((.*)\\)$"),
        logHR    = coef,
        HR       = exp(coef),
        HR_low   = `lower .95`,
        HR_high  = `upper .95`
      ) %>%
      left_join(top_by_group %>% select(var_name, data_var, group), by = "data_var")

    forest_res <- forest_res %>%
      mutate(var_display = fct_reorder(var_name, HR))

    p_forest <- ggplot(forest_res, aes(x = HR, y = var_display, color = group)) +
      geom_point(size = 2.5) +
      geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.2) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      scale_x_log10() +
      labs(
        x = "Hazard Ratio (log scale)",
        y = NULL,
        color = "Group",
        title = "Representative axes for HF risk in T2D (NHANES)"
      ) +
      theme_bw(base_size = 12)

    ggsave(
      filename = file.path(out_dir, "forest_T2D_HF_main_axes.png"),
      plot = p_forest,
      width = 8, height = 6, dpi = 300
    )
  }
}

#------------------------------------------------------------
# 8. 打印前 20 个变量
#------------------------------------------------------------

res_df %>%
  arrange(p) %>%
  slice_head(n = 20) %>%
  select(var_name, group, HR, HR_low, HR_high, p) %>%
  print(n = 20)
