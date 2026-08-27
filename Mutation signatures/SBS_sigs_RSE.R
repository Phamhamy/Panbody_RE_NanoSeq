
library(tidyverse)
library(MuMIn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(forcats)
library(scales)


sbs_sigs_residual <- final_sigs %>% pivot_longer(cols = matches("SBS"), names_to = "signature", values_to = "sigCount") 


df_agg_sig <- sbs_sigs_residual %>%
  group_by(structure1, signature, PDid, age) %>%
  summarise(mutations = sigCount,
            .groups = "drop")

df_agg_total <- df_agg_sig %>%
  group_by(structure1, PDid, age) %>%
  summarise(total_muts = sum(mutations),
            .groups = "drop")

total_RSE <- df_agg_total %>%
  group_by(structure1) %>%
  group_modify(~{
    fit <- lm(total_muts ~ age, data = .x)
    s   <- summary(fit)
    
    # p-value for the age slope
    p_value <- s$coefficients["age", "Pr(>|t|)"]
    
    tibble(
      RSE   = s$sigma,
      RSE2  = s$sigma^2,
      R2    = s$r.squared,
      slope = s$coefficients["age", "Estimate"],
      p     = p_value,
      n     = nrow(.x)
    )
  })


total_RSE <- total_RSE %>%
  mutate(p_adj = p.adjust(p, method = "BH"))


resid_list <- df_agg_sig %>%
  group_by(structure1, signature) %>%
  group_modify(~{
    fit <- lm(mutations ~ age, data = .x)
    tibble(
      donor    = .x$PDid,
      age      = .x$age,
      residual = residuals(fit)
    )
  }) %>%
  ungroup()

# Verify: one row per structure1 × donor × signature
resid_list %>%
  summarise(n = n(), .by = c(structure1, donor, signature)) %>%
  filter(n > 1)   # should return 0 rows


var_results <- resid_list %>%
  dplyr::select(-age) %>%                          # drop age before pivoting
  pivot_wider(names_from  = signature,
              values_from = residual) %>%
  group_by(structure1) %>%
  group_modify(~{
    
    R <- .x %>%
      dplyr::select(-donor) %>%
      dplyr::select(where(is.numeric))
    
    if (any(sapply(R, is.list))) {
      warning("List columns in: ", .y$structure1, " — duplicates remain")
      return(tibble())
    }
    if (nrow(R) < 3) {
      warning("Too few donors in: ", .y$structure1, " — skipping")
      return(tibble())
    }
    
    cov_mat     <- cov(R, use = "pairwise.complete.obs")
    total_var   <- sum(cov_mat)
    diag_var    <- diag(cov_mat)
    diag_sum <- sum(diag_var)
    cov_contrib <- total_var - sum(diag_var)
    pct         <- diag_var / sum(diag_var) * 100
    bind_cols(
      tibble(
        total_residual_var = total_var,
        RSE_decomp          = sqrt(total_var),
        diag_sum = diag_sum,
        covariance_contrib = cov_contrib,
        cov_pct_of_total   = cov_contrib / total_var * 100
      ),
      as_tibble(t(diag_var)) %>% rename_with(~paste0("var_", .)),
      as_tibble(t(pct))      %>% rename_with(~paste0("pct_", .))
    )
    
  }) %>%
  arrange(desc(total_residual_var))


new_struct_order_rse <- unique(as.character(
  total_RSE[order(total_RSE$RSE), ]$structure1
))

total_RSE$structure1   <- factor(total_RSE$structure1,
                                 levels = new_struct_order_rse)
var_results$structure1 <- factor(var_results$structure1,
                                 levels = new_struct_order_rse)
resid_list$structure1  <- factor(resid_list$structure1,
                                 levels = new_struct_order_rse)


var_results  <- total_RSE %>%
  dplyr::select(structure1, RSE_total = RSE) %>%
  left_join(
    var_results,
    by = "structure1"
  ) %>%
  mutate(
    var_diff = round(diag_sum - abs(covariance_contrib), 2),
    difference     = round(RSE_total - RSE_decomp, 2),
    pct_difference = round(abs(RSE_total - RSE_decomp) / RSE_total * 100, 1)
  )



var_results <- var_results %>%
  ungroup() %>%
  mutate(across(starts_with("pct_"),
                ~ . * RSE_total/100,
                .names = "{str_replace(.col, 'pct_', 'rse_')}"))

sigOrder <- c("SBS1","SBS5","SBS4","SBS7a","SBS7b","SBS7d","SBS9","SBS16",
              "SBS18","SBS19","SBS40a","SBS40b","SBS40c","SBS88","SBS92",
              "SBSA","SBSB","SBSC")

var_long <- var_results %>%
  dplyr::select(structure1, matches("^rse_[A-Z]"), -RSE_total, -RSE_decomp,
                -rse_difference, -rse_pred_conditional_sbs1, -rse_pred_conditional_sbs5,
                -rse_pred_marginal_sbs1, -rse_pred_marginal_sbs5) %>%
  pivot_longer(
    cols         = starts_with("rse_"),
    names_to     = "signature",
    values_to    = "RSE",
    names_prefix = "rse_"
  ) %>%
  mutate(signature = factor(signature, levels = sigOrder)) %>%
  group_by(structure1) %>%
  mutate(pct = RSE / sum(RSE, na.rm = TRUE) * 100) %>%
  ungroup()

var_long$signature <- factor(var_long$signature, levels = rev(sigOrder))

outlier_structs <- c("Skin: Epidermis", "Skin: Hair follicles")

var_long_main <- var_long %>% filter(RSE > 0, !structure1 %in% outlier_structs)
var_long_main$structure1 <- factor(var_long_main$structure1, levels = rev(df$cell_type))

var_long_outlier <- var_long %>% filter(RSE > 0, structure1 %in% outlier_structs)

spacers <- df %>%
  filter(!cell_type %in% outlier_structs) %>%
  arrange(quadrant, desc(residual_var)) %>%
  group_by(quadrant) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  filter(quadrant != last(levels(quadrant))) %>%
  mutate(
    cell_type = paste0(".spacer_", quadrant),
    RSE       = NA
  ) %>%
  dplyr::select(cell_type, quadrant)

ordered_levels <- df %>%
  filter(!cell_type %in% outlier_structs) %>%
  arrange(quadrant, desc(residual_var)) %>%
  group_by(quadrant) %>%
  group_modify(~{
    quad <- unique(.y$quadrant)
    if (quad != last(levels(df$quadrant))) {
      bind_rows(.x, tibble(cell_type = paste0(".spacer_", quad), RSE = NA))
    } else {
      .x
    }
  }) %>%
  ungroup() %>%
  pull(cell_type)

spacer_rows <- tibble(
  structure1 = paste0(".spacer_", levels(df$quadrant)[
    levels(df$quadrant) != last(levels(df$quadrant))
  ]),
  signature  = first(levels(var_long_main$signature)),
  RSE        = NA,
  pct        = NA
)

var_long_main_spaced <- var_long_main %>%
  bind_rows(spacer_rows) %>%
  mutate(structure1 = factor(structure1, levels = rev(ordered_levels))) %>%
  mutate(signature = factor(signature, levels = rev(sigOrder)))

sig_colours <- c(
  "SBSC"   = "salmon",
  "SBSB"   = "gold3",
  "SBSA"   = "firebrick",
  "SBS92"  = "springgreen2",
  "SBS88"  = "chocolate4",
  "SBS40c" = "pink",
  "SBS40b" = "#984ea3",
  "SBS40a" = "darkorchid4",
  "SBS19"  = "maroon2",
  "SBS18"  = "darkturquoise",
  "SBS16"  = "yellow",
  "SBS9"   = "blue",
  "SBS7d"  = "wheat4",
  "SBS7b"  = "wheat3",
  "SBS7a"  = "wheat",
  "SBS4"   = "springgreen4",
  "SBS5"   = "#80b1d3",
  "SBS1"   = "#fdb462"
)

#plotting
p_main <- ggplot(var_long_main_spaced, aes(x = RSE, y = structure1, fill = signature)) +
  geom_col(na.rm = TRUE) +
  labs(x = "RSE values", y = NULL, fill = "Signature") +
  scale_fill_manual(values = sig_colours, na.value = "transparent") +
  scale_y_discrete(labels = function(x) ifelse(startsWith(x, ".spacer_"), "", x)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p_outlier <- ggplot(var_long_outlier, aes(x = RSE, y = structure1, fill = signature)) +
  geom_col() +
  labs(x = "RSE values", y = NULL, fill = "Signature", title = NULL) +
  scale_fill_manual(values = sig_colours, na.value = "transparent") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

var_sigs_summary <- (p_main / p_outlier) +
  plot_layout(heights = c(1.5, 0.08))

var_sigs_summary

