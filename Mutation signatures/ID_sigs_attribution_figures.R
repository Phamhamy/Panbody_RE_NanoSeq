
## ## sinatures attribution by cell type

final_sigs_id <- read.table("ID_sigs_attribution.txt", header=T, check.names=F, sep="\t", quote="") %>%
  rowwise(sequenceid) %>% 
  mutate(totalVars = sum(c_across(starts_with("ID")))) %>% # Use c_across for cleaner code
  ungroup() %>% 
  mutate(across(starts_with("ID"), ~ .x / totalVars)) %>% 
  left_join(panbody_metadata %>% 
              dplyr::select(label, PDid, sequenceid, age, structure, 
                            structure1, smoking, drinking, Ethnicity), 
            by = "sequenceid")

final_sigs_id <- na.omit(final_sigs_id)

sigOrder_id <- c("ID1","ID2","ID3","ID5","ID8","ID9","ID11","ID13","ID21")

final_sigs_id <- left_join(final_sigs_id, sbs_burdens[,c(1,13)], by = "sequenceid")

final_sigs_id <- final_sigs_id[final_sigs_id$structure1 != "Testis: Seminiferous tubules",]


final_sigs_id <- final_sigs_id %>%
  mutate(across(starts_with("ID"), ~ .x * burden_indel))

final_sigs_id <- final_sigs_id %>%
  mutate(across(starts_with("ID"), round))

plot_data_id <- final_sigs_id %>% 
  mutate(sequenceid = fct_reorder(sequenceid, desc(age))) %>%
  pivot_longer(
    cols = starts_with("ID"), 
    names_to = "sig", 
    values_to = "sigCount"
  ) %>% 
  mutate(sigCount = replace_na(sigCount, 0)) %>%
  mutate(sig = fct_relevel(sig, rev(sigOrder_id)))

plot_data_id <- plot_data_id %>%
  filter(sigCount > 0) %>% 
  mutate(across(where(is.factor), droplevels))

common_layers <- list(
  scale_fill_manual(values = c("firebrick","darkblue","pink","darkseagreen2","chocolate3","#984ea3","springgreen4",'#80b1d3','#fdb462')),
  theme_pubr(),
  
  facet_grid(structure1 ~ ., scales = "free", space = "free", switch = "y"),
  theme(
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size=11, face = "bold"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )
)

plot_data_id <- plot_data_id %>%
  filter(sigCount > 0) %>% 
  mutate(across(where(is.factor), droplevels))

p_stack_id <- ggplot(plot_data_id, aes(y = sequenceid, x = sigCount, fill = sig)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  labs(x = "Absolute Count") +
  common_layers

p_fill_id <- ggplot(plot_data_id, aes(y = sequenceid, x = sigCount, fill = sig)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  labs(x = "Proportion") +
  common_layers + 
  theme(strip.text.y.left = element_blank())

p_final_id <- (p_stack_id | p_fill_id) + 
  plot_layout(widths = c(3, 2), guides = "collect") & 
  theme(legend.position = "right")

p_final_id


## Bubble and stacked bar plots


df_avg_id <- final_sigs_id %>% 
  pivot_longer(
    cols = starts_with("ID"), 
    names_to = "sig", 
    values_to = "sigCount"
  ) %>% 
  group_by(sequenceid) %>% # Use your sample ID column here
  mutate(sample_prop = sigCount / sum(sigCount)) %>%
  group_by(structure1, sig) %>%
  summarise(
    avg_proportion = mean(sample_prop, na.rm = TRUE),
    perc_presence = sum(sigCount > 0) / n(),.groups = 'drop') %>%
  mutate(sig = fct_relevel(sig, rev(sigOrder_id))) %>%
  mutate(structure1 = fct_rev(factor(structure1)))

# The Vertical Bar Plot (Summary)
p_bar_vert_id <- ggplot(df_avg_id, aes(x = avg_proportion, y = structure1, fill = sig)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("firebrick","darkblue","pink","darkseagreen2","chocolate4","#984ea3","springgreen4",'#80b1d3','#fdb462')) +
  theme_pubr() +
  labs(x = "Mean Contribution") +
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    legend.position = "right"
  )

# The Vertical Bubble Plot
p_bubble_vert_id <- df_avg_id %>%
  filter(perc_presence > 0.001) %>%
  # Reverse the signature order here
  mutate(sig = fct_relevel(sig, sigOrder_id)) %>%
  ggplot(aes(x = sig, y = structure1, size = perc_presence, color = sig)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(1, 12), name = "% of samples", labels = scales::percent) +
  scale_colour_manual(values = c("#fdb462", "#80b1d3", "springgreen4", "#984ea3", "chocolate4", "darkseagreen2", "pink", "darkblue", "firebrick")) +
  theme_pubr() +
  labs(x = "Signatures", y = "Structure") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    
    axis.text.y = element_text(hjust = 1, vjust = 0.5, face = "bold"),
    legend.position = "none" # patchwork will handle this later
  )

library(patchwork)

combined_plot_vert_id <- p_bubble_vert_id | p_bar_vert_id + 
  plot_layout(widths = c(3, 1), guides = 'collect') & 
  theme(legend.position = "right")

combined_plot_vert_id