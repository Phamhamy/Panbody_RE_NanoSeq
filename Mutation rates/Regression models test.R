
library(tidyverse)
library(lme4)
library(lmerTest)
library(nlme)        
library(broom)
library(broom.mixed)
library(AICcmodavg)

set.seed(123)


burden_model <- burdens_new[burdens_new$tissue != "Skin",]

### linear vs quadratic vs exponential testing

mm_linear <- lmer(burden_wg ~ age + (1|PDid), data = burden_model)  # ML, not REML, so AIC comparisons across

mm_quad <- lmer(burden_wg ~ age + I(age^2) + (1|PDid), data = burden_model)  # ML, not REML, so AIC comparisons across

log_start <- lm(log(burden_wg) ~ age, data = burden_model)
a_start <- exp(coef(log_start)[1])
b_start <- coef(log_start)[2]
c(a_start, b_start)

library(nlme)
mm_exp <- nlme(
  burden_wg ~ a * exp(b * age),
  data = burden_model,
  fixed = a + b ~ 1,
  random = a ~ 1 | PDid,
  start = c(a = a_start, b = b_start)
)

aicc_pooled <- tibble(
  model = c("Linear (mixed)", "Quadratic (mixed)", "Exponential (nlme)"),
  AIC  = c(AIC(mm_linear), AIC(mm_quad), AIC(mm_exp))
) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)

print(aicc_pooled)

anova(mm_exp, mm_linear)

### plotting

age_grid <- seq(min(burden_model$age), max(burden_model$age), length.out = 200)

fe_lin   <- fixef(mm_linear)
pred_lin <- fe_lin["(Intercept)"] + fe_lin["age"] * age_grid

fe_quad  <- fixef(mm_quad)
pred_quad <- fe_quad["(Intercept)"] +
  fe_quad["age"]         * age_grid +
  fe_quad["I(age^2)"]    * age_grid^2

fe_exp   <- fixed.effects(mm_exp)          # nlme uses fixed.effects()
pred_exp <- fe_exp["a"] * exp(fe_exp["b"] * age_grid)

pred_df <- tibble(
  age   = rep(age_grid, 3),
  y     = c(pred_lin, pred_quad, pred_exp),
  model = rep(c("Linear (mixed)",
                "Quadratic (mixed)",
                "Exponential (nlme)"),
              each = length(age_grid))
)


models_compare <- ggplot() +
  geom_point(data  = burden_model,
             aes(x = age, y = burden_wg),
             alpha = 0.3, size = 2, colour = "grey40") +
  geom_line(data  = pred_df,
            aes(x = age, y = y, colour = model, linetype = model),
            linewidth = 1) +
  scale_colour_manual(
    values = c("Linear (mixed)"     = "#4472C4",
               "Quadratic (mixed)"  = "#ED7D31",
               "Exponential (nlme)" = "#70AD47")) +
  scale_linetype_manual(
    values = c("Linear (mixed)"     = "solid",
               "Quadratic (mixed)"  = "dashed",
               "Exponential (nlme)" = "dotdash")) +
  coord_cartesian(ylim = c(0, max(burden_model$burden_wg) * 1.1)) +
  labs(x      = "Age",
       y      = "Number of SBS mutations per diploid cell",
       colour = "Model",
       linetype = "Model") +
  theme_pubr() +
  theme(legend.position = "bottom")


ggsave(models_compare ,filename = "models_compare.pdf",dpi=300, units="in", width = 10, height = 7)

## linear models testing 

model1 <- lmer(burden_wg ~ age + (0 + age | structure1),
               data = burden_model, REML = FALSE,
               control = lmerControl(optimizer = "bobyqa"))

model2 <- lmer(burden_wg ~ age+ (0 + age | structure1) + (1 | PDid),
               data = burden_model, REML = FALSE,
               control = lmerControl(optimizer = "bobyqa"))

model3 <- lmer(burden_wg ~ age+ (1 + age | structure1),
               data = burden_model, REML = FALSE,
               control = lmerControl(optimizer = "bobyqa"))

model4 <- lmer(burden_wg ~ age + (1 + age | structure1) + (1 | PDid),
               data = burden_model, REML = FALSE,
               control = lmerControl(optimizer = "bobyqa"))

anova(model1, model2)

anova(model2, model3)

anova(model3, model4)

isSingular(model2)
isSingular(model4)


