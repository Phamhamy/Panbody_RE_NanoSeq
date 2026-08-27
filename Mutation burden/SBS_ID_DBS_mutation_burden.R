
### order structure by median of burdens

library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggplotify)
library(RColorBrewer)
library(readxl)

burdens_new <- read.table("/Panbody_burden_all_donor.txt",header = TRUE,sep = "\t",quote = "",stringsAsFactors = FALSE)  

structures = burdens_new %>% dplyr::select(structure1) %>% unique()
structures = as.list(structures$structure1)

struct_order  <- data.frame(
  structure=character(),
  burden_mega=numeric(),
  median=numeric(), 
  stringsAsFactors=FALSE) 

for(i in structures){
  
  burden <- burdens_new %>% filter(structure1 == i) %>% dplyr::select(structure1, burden_wg)
  burden$median = median(burden$burden_wg)
  
  struct_order <- rbind(struct_order, burden)
  
}

struct_order = struct_order[,-2]
struct_order = distinct(struct_order, structure1, .keep_all = TRUE)

struct_order <- struct_order[order(struct_order$median),]
rownames(struct_order) <- NULL

burdens_new$structure1 <- factor(burdens_new$structure1 , levels = struct_order$structure1)

#order PDid based on burden_wg
burden_order <- as.data.frame(unique(burdens_new$sequenceid))
colnames(burden_order) <- "sequenceid"
for(i in burden_order$sequenceid){
  burden_order[which(burden_order$sequenceid == i), "Burden"] <- unique(burdens_new[which(burdens_new$sequenceid == i), "burden_wg"])
}
burden_order <- burden_order[order(burden_order$Burden),]
rownames(burden_order) <- NULL

burdens_new$sequenceid <- factor(burdens_new$sequenceid, levels = burden_order$sequenceid, labels = paste0(burden_order$sequenceid,"(",burden_order$Burden,")"))

burdens_new <- burdens_new[order(burdens_new$structure1,burdens_new$burden_wg),]



shading_data <- burdens_new %>%
  distinct(structure1) %>%
  mutate(facet_id = as.numeric(as.factor(structure1)),
         fill_color = ifelse(facet_id %% 2 == 0, "grey95", "white"))

SBS_burden <- ggplot(burdens_new, aes(x = structure1)) + 
  geom_rect(data = shading_data, 
            aes(fill = fill_color), 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_point(aes(y = burden_wg, group = sequenceid), colour = "springgreen4",
             size = 1.5, position = position_dodge2(width = 0.5)) +
  stat_summary(aes(y = burden_wg), fun = median, position = position_dodge2(0.5), 
               colour = "firebrick3", geom = "point", shape = 95, size = 13, show.legend = FALSE) + 
  facet_wrap(~structure1, scales = "free_x", nrow = 1) +
  scale_y_break(c(15000, 50000), scales = 0.2) + 
  labs(x = NULL, y = "SBSs per diploid cell") +
  theme_pubr() + 
  theme(legend.position = "none",
        strip.text = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14, face = "bold"), 
        axis.title.y = element_text (size= 15, face = "bold")) 

ggsave("SBS_burden_by_structure.pdf", plot = SBS_burden, height=5, width=20)


ID_burden <- ggplot(burdens_new, aes(x = structure1)) + 
  geom_rect(data = shading_data, 
            aes(fill = fill_color), 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_point(aes(y = burden_indel, group = sequenceid), colour = "#3f007d",
             size = 1.5, position = position_dodge2(width = 0.5)) +
  stat_summary(aes(y = burden_indel), fun = median, position = position_dodge2(0.5), 
               colour = "gold", geom = "point", shape = 95, size = 13, show.legend = FALSE) + 
  facet_wrap(~structure1, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "IDNELs per diploid cell") +
  theme_pubr() + 
  theme(legend.position = "none",
        strip.text = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14, face = "bold"), 
        axis.title.y = element_text (size= 15, face = "bold"))

ggsave("ID_burden_by_structure.pdf", plot = ID_burden, height=5, width=20)


DBS_burden <- ggplot(burdens_new, aes(x = structure1)) + 
  geom_rect(data = shading_data, 
            aes(fill = fill_color), 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_point(aes(y = burden_dbs, group = sequenceid), colour = "orange2",
             size = 1.5, position = position_dodge2(width = 0.5)) +
  stat_summary(aes(y = burden_dbs), fun = median, position = position_dodge2(0.5), 
               colour = "darkblue", geom = "point", shape = 95, size = 13, show.legend = FALSE) +
  scale_y_break(c(400, 600), scales = 0.2) + 
  
  facet_wrap(~structure1, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "DBSs per diploid cell") +
  theme_pubr() + 
  theme(legend.position = "none",
        strip.text = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size=13, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"), 
        axis.title.y = element_text (size= 15, face = "bold"))

ggsave("DBS_burden_by_structure.pdf", plot = DBS_burden, height=8, width=20)


burden_combined_plot <- (wrap_elements(SBS_burden) / 
                           wrap_elements(ID_burden) / 
                           wrap_elements(DBS_burden)) + 
  plot_layout(ncol = 1, guides = 'collect')

ggsave("SBS_ID_DBS_burden_by_structure.pdf", plot = burden_combined_plot, height=15, width=20)
