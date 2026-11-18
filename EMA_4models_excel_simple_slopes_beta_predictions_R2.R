library(readxl)
library(dplyr)
library(lmerTest)
library(emmeans)
library(openxlsx)
library(performance)

folder <- "data"
same <- read_excel(file.path(folder, "datasame.xlsx"))
lag  <- read_excel(file.path(folder, "datalag.xlsx"))

same$ID <- factor(same$ID)
lag$ID  <- factor(lag$ID, levels = levels(same$ID))

means <- same %>%
  group_by(ID) %>%
  summarise(PSI_mean = mean(PSI, na.rm = TRUE),
            ASI_mean = mean(ASI, na.rm = TRUE),
            POS = first(SRMposi),
            DANGER = first(SRMdenger), .groups = "drop") %>%
  mutate(POS_c = POS - mean(POS, na.rm = TRUE),
         DANGER_c = DANGER - mean(DANGER, na.rm = TRUE)) %>%
  select(ID, PSI_mean, ASI_mean, POS_c, DANGER_c)

same <- same %>%
  left_join(means, by = "ID") %>%
  mutate(PSI_wp = PSI - PSI_mean, ASI_wp = ASI - ASI_mean)
lag <- lag %>%
  left_join(means, by = "ID") %>%
  mutate(PSI_wp = PSI - PSI_mean, ASI_wp = ASI - ASI_mean)

m1 <- lmer(ASI ~ PSI_wp * POS_c + (1 + PSI_wp | ID),
           data = same, REML = FALSE)
m2 <- lmer(ASIt1 ~ PSI_wp * POS_c + ASI + (1 + PSI_wp | ID),
           data = lag, REML = FALSE)
m3 <- lmer(DEP ~ ASI_wp * DANGER_c + (1 + ASI_wp | ID),
           data = same, REML = FALSE)
m4 <- lmer(DEPt1 ~ ASI_wp * DANGER_c + DEP + (1 + ASI_wp | ID),
           data = lag, REML = FALSE)

models <- list(ASI_same = m1, ASI_lag = m2, DEP_same = m3, DEP_lag = m4)


r2_table <- bind_rows(lapply(names(models), function(nm) {
  r2 <- performance::r2_nakagawa(models[[nm]])
  data.frame(Model = nm,
             Marginal_R2 = r2$R2_marginal)
}))


coef_table <- bind_rows(lapply(names(models), function(nm) {
  model <- models[[nm]]
  cf <- coef(summary(model))
  X <- model.matrix(model)
  sd_y <- sd(model.response(model.frame(model)))
  beta <- sapply(rownames(cf), function(term) {
    if (term == "(Intercept)") return(NA_real_)
    cf[term, "Estimate"] * sd(X[, term]) / sd_y
  })
  data.frame(Model = nm, Variable = rownames(cf),
             B = cf[, "Estimate"], beta = unname(beta),
             SE = cf[, "Std. Error"], df = cf[, "df"],
             t = cf[, "t value"], p = cf[, "Pr(>|t|)"],
             row.names = NULL)
}))

simple_slope <- function(model, name, x, moderator, moderator_sd) {
  trend <- emtrends(model, specs = moderator, var = x,
                    at = setNames(list(c(-moderator_sd, moderator_sd)), moderator),
                    lmer.df = "satterthwaite")
  s <- as.data.frame(summary(trend, infer = c(TRUE, TRUE)))
  slope_col <- grep("\\.trend$", names(s), value = TRUE)
  sd_x <- sd(model.matrix(model)[, x])
  sd_y <- sd(model.response(model.frame(model)))
  data.frame(Model = name, Moderator = moderator,
             Level = c("-1SD (低)", "+1SD (高)"),
             Moderator_centered = s[[moderator]],
             Variable = x, B = s[[slope_col]],
             beta = s[[slope_col]] * sd_x / sd_y,
             SE = s$SE, df = s$df, t = s$t.ratio, p = s$p.value)
}

pos_sd    <- sd(means$POS_c, na.rm = TRUE)
danger_sd <- sd(means$DANGER_c, na.rm = TRUE)

slopes_table <- bind_rows(
  simple_slope(m1, "ASI_same", "PSI_wp", "POS_c", pos_sd),
  simple_slope(m2, "ASI_lag",  "PSI_wp", "POS_c", pos_sd),
  simple_slope(m3, "DEP_same", "ASI_wp", "DANGER_c", danger_sd),
  simple_slope(m4, "DEP_lag",  "ASI_wp", "DANGER_c", danger_sd)
)

# グラフ用：独立変数と調整変数がそれぞれ ±1SD の場合の予測値
plot_values <- function(model, name, x, moderator, moderator_sd, covariate = NULL) {
  dat <- model.frame(model)
  x_sd <- sd(dat[[x]], na.rm = TRUE)
  at_values <- setNames(list(c(-moderator_sd, moderator_sd), c(-x_sd, x_sd)),
                        c(moderator, x))
  if (!is.null(covariate)) {
    at_values[[covariate]] <- mean(dat[[covariate]], na.rm = TRUE)
  }
  pred <- as.data.frame(emmeans(model, specs = c(moderator, x),
                               at = at_values))
  value <- function(m, v) {
    i <- which.min(abs(pred[[moderator]] - m) + abs(pred[[x]] - v))
    pred$emmean[i]
  }
  data.frame(Model = name, Moderator = paste0(moderator, c("_-1SD", "_+1SD")),
             Predictor = x,
             Predictor_minus_1SD = c(value(-moderator_sd, -x_sd), value(moderator_sd, -x_sd)),
             Predictor_plus_1SD = c(value(-moderator_sd, x_sd), value(moderator_sd, x_sd)),
             Covariate_fixed = if (is.null(covariate)) NA_real_
                               else mean(dat[[covariate]], na.rm = TRUE))
}

graph_table <- bind_rows(
  plot_values(m1, "ASI_same", "PSI_wp", "POS_c", pos_sd),
  plot_values(m2, "ASI_lag",  "PSI_wp", "POS_c", pos_sd, "ASI"),
  plot_values(m3, "DEP_same", "ASI_wp", "DANGER_c", danger_sd),
  plot_values(m4, "DEP_lag",  "ASI_wp", "DANGER_c", danger_sd, "DEP")
)

wb <- createWorkbook()
addWorksheet(wb, "回帰係数")
addWorksheet(wb, "単純主効果")
addWorksheet(wb, "グラフ用予測値")
addWorksheet(wb, "R2")
writeData(wb, "回帰係数", coef_table)
writeData(wb, "単純主効果", slopes_table)
writeData(wb, "グラフ用予測値", graph_table)
writeData(wb, "R2", r2_table)
setColWidths(wb, "回帰係数", cols = 1:ncol(coef_table), widths = "auto")
setColWidths(wb, "単純主効果", cols = 1:ncol(slopes_table), widths = "auto")
setColWidths(wb, "グラフ用予測値", cols = 1:ncol(graph_table), widths = "auto")
setColWidths(wb, "R2", cols = 1:ncol(r2_table), widths = "auto")
freezePane(wb, "回帰係数", firstRow = TRUE)
freezePane(wb, "単純主効果", firstRow = TRUE)
freezePane(wb, "グラフ用予測値", firstRow = TRUE)
freezePane(wb, "R2", firstRow = TRUE)
saveWorkbook(wb, file.path(folder, "EMA_4models_results.xlsx"), overwrite = TRUE)
