
#calculating RSE and p-value for each cell type

# Define the sequence of ages you want to predict for
age_seq <- data.frame(age = 0:100)

unique_structures <- unique(burdens_new$structure1)

lm_equations_all   <- data.frame()

lm_equations_all <- NULL

for (struct in unique_structures) {
  
  filtered_mut_burden <- burdens_new %>% filter(structure1 == struct)
  
  if(nrow(filtered_mut_burden) < 2) next 
  
  structure_model <- lm(round(burden_wg) ~ age, data = filtered_mut_burden)
  
  model_summary <- summary(structure_model)
  R2 <- model_summary$r.squared
  intercept <- round(coef(structure_model)["(Intercept)"], 2)
  slope     <- round(coef(structure_model)["age"], 2)
  
  p_val_age <- if("age" %in% rownames(model_summary$coefficients)) {
    model_summary$coefficients["age", "Pr(>|t|)"]
  } else {
    NA
  }
  
  equation_row <- data.frame(
    structure1 = struct,
    intercept = intercept,
    R2 = R2,
    slope = slope,
    p_value = p_val_age,
    # This matches the "variation per structure" you wanted
    res_std_error = round(model_summary$sigma), 
    model_eq = paste0("y = ", intercept, " + ", slope, "*Age"),
    y_val = max(filtered_mut_burden$burden_wg, na.rm = TRUE) * 1.1
  )
  
  lm_equations_all <- rbind(lm_equations_all, equation_row)
}

lm_equations_all$p_value_adj <- p.adjust(
  lm_equations_all$p_value, 
  method = "BH"
)


lm_equations_all <- lm_equations_all %>%
  arrange(slope) 

#calculating slope for each cell type 
model_lmer <- lmer(burden_wg ~ age+ (0 + age | structure1) + (1 | PDid),
                   data = burden_model,
                   control = lmerControl(optimizer = "bobyqa"))

summary(model_lmer)

coef <- coefficients(model_lmer)$structure1

coef <- coef %>% rownames_to_column(var = "structure1")


# calculating 95%CI for each cell type

all_predictions <- expand.grid(
  age = seq(0, 100, length.out = 50),
  structure1 = unique(burden_model$structure1)
)

# Prediction function (include random effects for structure1, i.e. re.form = ~(0+age||structure1))
pred_fun <- function(fit) {
  predict(fit, newdata = all_predictions, re.form = ~(0 + age | structure1))
}

# Bootstrap (nsim = 200-1000; higher = smoother CI but slower)
boot_res <- bootMer(model_lmer, FUN = pred_fun, nsim = 1000,
                    use.u = TRUE, type = "parametric",
                    parallel = "multicore", ncpus = 4)
all_predictions$predicted <- predict(model_lmer, newdata = all_predictions, re.form = ~(0 + age | structure1))
all_predictions$conf.low  <- apply(boot_res$t, 2, quantile, probs = 0.025, na.rm = TRUE)
all_predictions$conf.high <- apply(boot_res$t, 2, quantile, probs = 0.975, na.rm = TRUE)


burden_model$pred_marginal <- predict(model_lmer, re.form = NA)
burden_model$pred_conditional <- predict(model_lmer)

final_model_r2 <- burden_model %>%
  group_by(structure1) %>%
  summarize(
    Local_R2_Marginal = cor(burden_wg, pred_marginal)^2,
    Local_R2_Conditional = cor(burden_wg, pred_conditional)^2,
    .groups = "drop"
  )

print(final_model_r2)

lmer_model_per_Sample <- data.frame(structure1=unique(burden_model$structure1))
rownames(lmer_model_per_Sample) <- lmer_model_per_Sample$structure1


for(i in lmer_model_per_Sample$structure1){
  slope <- round(coef$age[coef$structure1 == i], digits = 1)
  slope_p <- format(lm_equations_all$p_value_adj[lm_equations_all$structure1 == i], scientific = TRUE, digits = 1)
  r2 <- round(final_model_r2$Local_R2_Conditional[final_model_r2$structure1 == i], digits = 1) 
  rse <- round(lm_equations_all$res_std_error[lm_equations_all$structure1 == i])
  lmer_model_per_Sample[i, "model"] <- paste("slope = ", slope," , ", "p = ", slope_p, sep = "")
  lmer_model_per_Sample[i, "stats"] <- paste("R\u00B2 = ", r2, " , ", "RSE = ", rse, sep = "")
  lmer_model_per_Sample[i, "slope"] <- slope
  lmer_model_per_Sample[i, "R2"] <- r2
  lmer_model_per_Sample[i, "rse"] <- rse
  lmer_model_per_Sample[i, "slope_p"] <- slope_p
}

#Order structure1 based on increasing slope
lmer_model_per_Sample <- lmer_model_per_Sample[order(lmer_model_per_Sample$slope),]

# Define the order explicitly
new_struct_order_sbs <- unique(as.character(lmer_model_per_Sample[order(lmer_model_per_Sample$slope), ]$structure1))

burden_model$structure1 <- factor(burden_model$structure1, levels = new_struct_order_sbs)
all_predictions$structure1 <- factor(all_predictions$structure1, levels = new_struct_order_sbs)
lmer_model_per_Sample$structure1 <- factor(lmer_model_per_Sample$structure1, levels = new_struct_order_sbs)

p_lmer_all <- ggplot() + 
  geom_point(data = burden_model,aes(x = age, y = burden_wg, fill = structure1),size = 4.5, pch = 21, stroke = 0.25) + 
  xlab("Age (years)") +
  ylab("Number of SBSs per diploid cell") +
  geom_line(data = all_predictions,aes(x = age, y = predicted, color = structure1), size = 0.7, linetype = "dashed") +
  geom_ribbon(data = all_predictions,aes(ymin = conf.low, ymax = conf.high, x = age, fill = structure1), alpha = 0.25) +
  geom_text(data = lmer_model_per_Sample, aes(label = as.character(model), x = 5, y = 9500),fontface = "bold",colour = "black", size = 3.5, parse = FALSE, hjust = 0, show.legend = FALSE) +
  geom_text(data = lmer_model_per_Sample, aes(label = as.character(stats), x = 5, y = 8500),fontface = "bold",colour = "black", size = 3.5, parse = FALSE, hjust = 0, show.legend = FALSE) +
  facet_wrap(~structure1, ncol = 6, axes = "all", axis.labels = "margins") + theme_pubr() + scale_fill_manual(values = my_palette) + scale_color_manual(values = my_palette) +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),strip.text.x.top = element_text(face = "bold", size = 10),panel.spacing = unit(0.2, "inches"),axis.text.y = element_text(face = "bold", size = 12),axis.text.x = element_text(face = "bold", size = 12),axis.title.x = element_text(face = "bold", size = 15),axis.title.y = element_text(face = "bold", size = 15),strip.placement = "outside") + ylim(0,10000)

p_lmer_all


## calculating 95%CI for each slope

slope_fun <- function(fit) {
  fixef(fit)["age"] + ranef(fit)$structure1[, "age"]
}

boot_slopes <- bootMer(model_lmer, FUN = slope_fun, nsim = 1000,
                       use.u = TRUE, type = "parametric",
                       parallel = "multicore", ncpus = 4)

# Structure names, in the same order ranef() returns them
structure_names <- rownames(ranef(model_lmer)$structure1)

sbs_boot_ci <- data.frame(
  structure1 = structure_names,
  slope      = fixef(model_lmer)["age"] + ranef(model_lmer)$structure1[, "age"],
  conf.low   = apply(boot_slopes$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  conf.high  = apply(boot_slopes$t, 2, quantile, probs = 0.975, na.rm = TRUE)
)

sbs_boot_ci

