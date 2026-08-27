## ## sinatures attribution by cell type

sbs_burdens  <- burdens_new


final_sigs <- read.table("SBS_sig_attribution.txt", header=T,check.names =F, sep="\t",quote = "") %>%
  rowwise(sequenceid) %>% mutate(totalVars = sum(c(SBS1,SBS4,SBS5,SBS7a,SBS7b,SBS7d,SBS9,SBS16,SBS18,SBS19,SBS40a,SBS40b,SBS40c,SBS88,SBS92,SBSA, SBSB, SBSC))) %>% ungroup() %>% mutate(across(matches("SBS"), ~ .x/totalVars)) %>%
  left_join(panbody_metadata %>% dplyr::select(label, PDid, sequenceid, age, structure, structure1, smoking, drinking, Ethnicity), by = "sequenceid") 

final_sigs <- na.omit(final_sigs)

sigOrder <- c("SBS1","SBS5","SBS4","SBS7a","SBS7b","SBS7d","SBS9","SBS16","SBS18","SBS19","SBS40a","SBS40b","SBS40c","SBS88","SBS92","SBSA","SBSB","SBSC")

final_sigs <- left_join(final_sigs, sbs_burdens[,c(1,10)], by = "sequenceid")

left_join(sbs_burdens[,c(1,10)],final_sigs, by = "sequenceid")


final_sigs <- final_sigs %>%
  mutate(across(starts_with("SBS"), ~ .x * burden_wg))

final_sigs <- final_sigs %>%
  mutate(across(starts_with("SBS"), round))

plot_data <- final_sigs[!final_sigs$structure %in% c("Epidermis", "Hair follicles"), ] %>% 
  mutate(sequenceid = fct_reorder(sequenceid, desc(age))) %>%
  pivot_longer(cols = matches("SBS"), names_to = "sig", values_to = "sigCount") %>% 
  mutate(sigCount = replace_na(sigCount, 0)) %>%
  mutate(sig = fct_relevel(sig, rev(sigOrder)))

plot_data <- plot_data %>%
  filter(sigCount > 0) %>% 
  mutate(across(where(is.factor), droplevels))

common_layers <- list(
  scale_fill_manual(values = c("salmon","gold3","firebrick","springgreen2","chocolate4","pink","#984ea3","darkorchid4","maroon2","turquoise","yellow","blue","wheat3","wheat","springgreen4",'#80b1d3','#fdb462')),
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

p_stack <- ggplot(plot_data, aes(y = sequenceid, x = sigCount, fill = sig)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  labs(x = "Absolute Count") +
  common_layers

p_fill <- ggplot(plot_data, aes(y = sequenceid, x = sigCount, fill = sig)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  labs(x = "Proportion") +
  common_layers +
  theme(strip.text.y.left = element_blank()) # Hide duplicate facet labels

p_final <- (p_stack | p_fill) + 
  plot_layout(widths = c(3, 2), guides = "collect") & 
  theme(legend.position = "right")

p_final


## plot Skin samples 
plot_data_skin <- final_sigs[final_sigs$structure %in% c("Epidermis", "Hair follicles"), ] %>% 
  # Use desc(age) so the youngest is at the top of the Y-axis
  mutate(sequenceid = fct_reorder(sequenceid, desc(age))) %>%
  pivot_longer(cols = matches("SBS"), names_to = "sig", values_to = "sigCount") %>% 
  mutate(sigCount = replace_na(sigCount, 0)) %>%
  mutate(sig = fct_relevel(sig, rev(sigOrder)))


common_layers <- list(
  scale_fill_manual(values = c("salmon","gold3","firebrick","springgreen2","chocolate4","pink","#984ea3","darkorchid4","maroon2","turquoise","yellow","blue","wheat4","wheat3","wheat","springgreen4",'#80b1d3','#fdb462')),
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

p_stack_skin <- ggplot(plot_data_skin, aes(y = sequenceid, x = sigCount, fill = sig)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  labs(x = "Absolute Count") +
  common_layers

p_fill_skin <- ggplot(plot_data_skin, aes(y = sequenceid, x = sigCount, fill = sig)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  labs(x = "Proportion") +
  common_layers +
  theme(strip.text.y.left = element_blank()) # Hide duplicate facet labels

p_final_skin <- (p_stack_skin | p_fill_skin) + 
  plot_layout(widths = c(3, 2), guides = "collect") & 
  theme(legend.position = "none")

p_final_skin

## combine figures
p_sigs_all <- (p_final / p_final_skin) /
  wrap_elements(legend) +
  plot_layout(heights = c(1.5,0.06,0.08))   # outlier small, main large


## bubble and stacked bar plot


df_avg <- final_sigs[,-c(30,31,32,33)] %>% 
  pivot_longer(cols = matches("SBS"), names_to = "sig", values_to = "sigCount") %>% 
  group_by(sequenceid) %>% # Use your sample ID column here
  mutate(sample_prop = sigCount / sum(sigCount)) %>%
  group_by(structure1, sig) %>%
  summarise(
    avg_proportion = mean(sample_prop, na.rm = TRUE),
    perc_presence = sum(sigCount > 0) / n(),.groups = 'drop') %>%
  mutate(sig = fct_relevel(sig, rev(sigOrder))) %>%
  mutate(structure1 = fct_rev(factor(structure1)))

# The Vertical Bar Plot (Summary)
p_bar_vert <- ggplot(df_avg, aes(x = avg_proportion, y = structure1, fill = sig)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("salmon","gold3","firebrick","springgreen2","chocolate4","pink","#984ea3","darkorchid4","maroon2","darkturquoise","yellow","blue","wheat4","wheat3","wheat","springgreen4",'#80b1d3','#fdb462')) +
  theme_pubr() +
  labs(x = "Mean Contribution") +
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    legend.position = "right"
  )

# The Vertical Bubble Plot
p_bubble_vert <- df_avg %>%
  filter(perc_presence > 0.001) %>%
  # Reverse the signature order here
  mutate(sig = fct_relevel(sig, sigOrder)) %>%
  ggplot(aes(x = sig, y = structure1, size = perc_presence, color = sig)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(1, 12), name = "% of samples", labels = scales::percent) +
  scale_colour_manual(values = c("#fdb462", "#80b1d3", "springgreen4", "wheat", "wheat3", "wheat4", "blue", "yellow", "darkturquoise", "maroon2", "darkorchid4", "#984ea3", "pink", "chocolate4", "springgreen2", "firebrick", "gold3", "salmon")) +
  theme_pubr() +
  labs(x = "Signatures", y = "Structure") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    
    axis.text.y = element_text(hjust = 1, vjust = 0.5, face = "bold"),
    legend.position = "none" # patchwork will handle this later
  )

library(patchwork)

combined_plot_vert <- p_bubble_vert | p_bar_vert + 
  plot_layout(widths = c(3, 1), guides = 'collect') & 
  theme(legend.position = "right")

combined_plot_vert
