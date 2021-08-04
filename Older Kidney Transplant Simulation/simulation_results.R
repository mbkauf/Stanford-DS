library(survminer)
library(survival)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)

##### Functions ######
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


##### Clean data #####
obs_df <- obs_df %>%
  rename(
    gl_time = months_to_gl,
    gs_death_time = months_to_gs_mort,
    gl_death_time = months_to_gl_mort
  ) %>%
  mutate(sim = rep(0, 4417)) %>%
  mutate(sim = factor(sim, levels = c(0, 1), 
                      labels = c("Actual", "Simulated"))) %>%
  mutate(tx_year_grp = ifelse(tx_year_c <= 5, 0, 1)) %>%
  mutate(tx_year_grp = factor(tx_year_grp, levels = c(0, 1),
                              labels = c("2010-2014", "2015+")))

sim_df <- sim_df %>% 
  mutate(tx_year_grp = ifelse(tx_year_c <= 5, 0, 1)) %>%
  mutate(tx_year_grp = factor(tx_year_grp, levels = c(0, 1),
                              labels = c("2010-2014", "2015+")))

subgrp_df <- rbind(sim_df, obs_df)
subgrp_df <- subgrp_df %>%
  mutate(kdpi_cat = ifelse(kdpi < 20, 0, 
                           ifelse(kdpi >= 20 & kdpi < 85, 1, 2))) %>%
  mutate(kdpi_cat = factor(kdpi_cat, levels = c(0, 1, 2),
                           labels = c("KDPI: < 20", "KDPI: 20-84", "KDPI: 85+"))) %>%
  mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
                           labels = c("Diabetes: No", "Diabetes: Yes")))
  

##### For SMDM #####
### Survfit Objects ###
fit_gl <- survfit(Surv(time = gl_time, graft_loss) ~ sim, 
                  data = subgrp_df)  # Graft loss

fit_gs_death <- survfit(Surv(time = gs_death_time, gs_death) ~ sim, 
                        data = subgrp_df)  # Death w/ function

fit_gl_death <- survfit(Surv(time = gl_death_time, gl_death) ~ sim, 
                        data = subgrp_df)  # Death after graft loss

fit_gl_txyr <- survfit(Surv(time = gl_time, graft_loss) ~ sim + tx_year_grp, 
                       data = subgrp_df)  # Graft loss

fit_gs_death_txyr <- survfit(Surv(time = gs_death_time, gs_death) ~ sim + tx_year_grp, 
                             data = subgrp_df)  # Death w/ function

fit_gl_death_txyr <- survfit(Surv(time = gl_death_time, gl_death) ~ sim + tx_year_grp, 
                             data = subgrp_df)  # Death after graft loss

### Survival Plots ###
gl_all <- ggsurvplot(fit_gl,
                      censor = FALSE,
                      title = "Graft Loss",
                      conf.int = TRUE,
                      short.panel.labs = TRUE,
                      legend.title = "",
                      ylim = c(.65, 1),
                      palette = "grey",
                      xlab = "Months")
gl_all$plot <- gl_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position='none')+ 
  theme(strip.text.x = element_text(size = 8))

gs_death_all <- ggsurvplot(fit_gs_death,
                            censor = FALSE,
                            conf.int = TRUE,
                            title = "Death with Function",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            palette = "grey",
                            xlab = "Months")
gs_death_all$plot <- gs_death_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position='none')+ 
  theme(strip.text.x = element_text(size = 8))

gl_death_all <- ggsurvplot(fit_gl_death,
                            censor = FALSE,
                            conf.int = TRUE,
                            title = "Death after Graft Loss",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            xlim = c(0, 36),
                            break.x.by = 6,
                            palette = "grey",
                            xlab = "Months")
gl_death_all$plot <- gl_death_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position='none') + 
  theme(strip.text.x = element_text(size = 8))

gl_year <- ggsurvplot(fit_gl_txyr,
                      censor = FALSE,
                      conf.int = TRUE,
                      title = "",
                      short.panel.labs = TRUE,
                      facet.by = "tx_year_grp",
                      legend.title = "",
                      ylim = c(.65, 1),
                      palette = "grey",
                      xlab = "Months")
gl_year <- gl_year +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position='none')+ 
  theme(strip.text.x = element_text(size = 8))

gs_death_year <- ggsurvplot(fit_gs_death_txyr,
                            censor = FALSE,
                            conf.int = TRUE,
                            title = "",
                            short.panel.labs = TRUE,
                            facet.by = "tx_year_grp",
                            legend.title = "",
                            ylim = c(0, 1),
                            palette = "grey",
                            xlab = "Months")
gs_death_year <- gs_death_year +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position='none')+ 
  theme(strip.text.x = element_text(size = 8))

gl_death_year <- ggsurvplot(fit_gl_death_txyr,
                            censor = FALSE,
                            conf.int = TRUE,
                            title = "",
                            short.panel.labs = TRUE,
                            facet.by = "tx_year_grp",
                            legend.title = "",
                            legend.labs = c("Observed", 
                                            "Simulated"),
                            ylim = c(0, 1),
                            xlim = c(0, 36),
                            break.x.by = 6,
                            palette = "grey",
                            xlab = "Months")
gl_death_year <- gl_death_year +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position='right') + 
  theme(strip.text.x = element_text(size = 8))

mylegend<-g_legend(gl_death_year)

smdm <- grid.arrange(arrangeGrob(gs_death_all$plot, gl_all$plot, gl_death_all$plot,
                         gs_death_year, gl_year, gl_death_year + 
                           theme(legend.position="none"),
                         nrow = 2),
             mylegend, nrow=2,heights=c(10, 1))

ggsave(filename = paste0(out_path, "smdm_bw.png"), plot = smdm , width = 7, height = 5)



##### Kaplan-Meier Plots (5-Year Analysis) #####
### All Recipients
## Graft Loss
fit_gl <- survfit(Surv(time = gl_time, graft_loss) ~ sim, data = df_test_long)

ggsurvplot(fit_gl, data = df_test_long,
           censor = FALSE,
           title = "Graft Loss",
           legend.title = "",
           legend.labs = c("Actual", "Simulated"),
           xlab = "Months",
           ylim = c(0.8, 1),
           xlim = c(0, 60),
           break.x.by = 12)
ggsave(filename = paste0(out_path, "gl_comp5.png"))

## Death w/ functioning graft
fit_gs_death <- survfit(Surv(time = gs_death_time, gs_death) ~ sim, 
                        data = df_test_long)

ggsurvplot(fit_gs_death, data = df_test_long,
           censor = FALSE,                    
           title = "Death with Functioning Graft",
           legend.title = "",
           xlab = "Months",
           legend.labs = c("Actual", "Simulated"),
           ylim = c(0.7, 1),
           break.x.by = 12)
ggsave(filename = paste0(out_path, "gs_death_comp5.png"))

## Death after graft loss
fit_gl_death <- survfit(Surv(time = gl_death_time, gl_death) ~ sim, 
                        data = df_test_long)

ggsurvplot(fit_gl_death, data = df_test_long,
           censor = FALSE,
           title = "Death after Graft Loss",
           legend.title = "",
           legend.labs = c("Actual", "Simulated"),
           xlab = "Months",
           ylim = c(0, 1),
           xlim = c(0, 40),
           break.x.by = 10)
ggsave(filename = paste0(out_path, "gl_death_comp5.png"))

## Combined graft loss and death w/ functioning graft
comb_fit <- list(Graft_Loss = fit_gl, GS_Death = fit_gs_death)

comb_plot <- ggsurvplot(comb_fit, data = df_test_long,
                        combine = TRUE,
                        censor = FALSE,                    
                        title = "Actual vs. Simulated",
                        legend.title = "",
                        legend.labs = c("Actual - Graft Loss", 
                                        "Simulated - Graft Loss",
                                        "Actual - Death w/ Graft", 
                                        "Simulated - Death w/ Graft"),
                        ylim = c(0, 1),
                        linetype = c("solid", "dotdash", "solid", "dotdash"))
comb_plot$plot <- comb_plot$plot + 
  theme(legend.text = element_text(size = 7)) +
  scale_color_grey()
comb_plot
ggsave(filename = paste0(out_path, "combined_gl_gsdeath_bw.png"), scale = 0.8)

### KM by Race
## Graft Loss
fit_gl_race <- survfit(Surv(time = gl_time, graft_loss) ~ sim + race, 
                       data = df_test_long)

ggsurvplot(fit_gl_race, data = df_test_long,
           censor = FALSE,                     # Remove censor points
           facet.by = "race",
           short.panel.labs = TRUE,
           title = "Graft Loss by Race/Ethnicity",
           legend.title = "",
           ylim = c(0.7, 1))
ggsave(filename = paste0(out_path, "gl_race_comp10.png"), scale = 0.8)

## Death w/ functioning graft
fit_gs_death_race <- survfit(Surv(time = gs_death_time, gs_death) ~ sim + race, 
                             data = df_test_long)

ggsurvplot(fit_gs_death_race, data = df_test_long,
           censor = FALSE,                     # Remove censor points,
           facet.by = "race",
           short.panel.labs = TRUE,
           title = "Death with Functioning Graft by Race/Ethnicity",
           legend.title = "",
           ylim = c(0, 1))
ggsave(filename = paste0(out_path, "gs_death_race_comp10.png"), scale = 0.8)

## Death w/ functioning graft
fit_gl_death_race <- survfit(Surv(time = gl_death_time, gl_death) ~ sim + race, 
                             data = df_test_long)

ggsurvplot(fit_gl_death_race, data = df_test_long,
           censor = FALSE,                     # Remove censor points,
           facet.by = "race",
           short.panel.labs = TRUE,
           title = "Death after Graft Loss by Race/Ethnicity",
           legend.title = "",
           ylim = c(0, 1))
ggsave(filename = paste0(out_path, "gl_death_race_comp10.png"), scale = 0.8)

### KM by Tx Year
## Graft Loss
fit_gl_txyr <- survfit(Surv(time = gl_time, graft_loss) ~ sim + tx_year_grp, 
                       data = subgrp_df)

ggsurvplot(fit_gl_txyr, data = subgrp_df,
           censor = FALSE,                     # Remove censor points
           facet.by = "tx_year_grp",
           short.panel.labs = TRUE,
           title = "Graft Loss by Transplant Year",
           legend.title = "",
           xlab = "Months",
           ylim = c(0.7, 1),
           xlim = c(0, 60),
           break.x.by = 12)
ggsave(filename = paste0(out_path, "gl_txyr_comp5.png"), scale = 0.8)

## Death w/ functioning graft
fit_gs_death_txyr <- survfit(Surv(time = gs_death_time, gs_death) ~ sim + tx_year_grp, 
                             data = subgrp_df)

ggsurvplot(fit_gs_death_txyr, data = subgrp_df,
           censor = FALSE,                     # Remove censor points,
           facet.by = "tx_year_grp",
           short.panel.labs = TRUE,
           title = "Death with Functioning Graft by Transplant Year",
           legend.title = "",
           ylim = c(0, 1),
           xlim = c(0, 60),
           break.x.by = 12,
           xlab = "Months")
ggsave(filename = paste0(out_path, "gs_death_txyr_comp5.png"), scale = 0.8)




##### Kaplan-Meier Plots (with CIs) #####
### Whole Sample ###
fit_gl <- survfit(Surv(time = gl_time, graft_loss) ~ sim, 
                  data = subgrp_df)  # Graft loss

fit_gs_death <- survfit(Surv(time = gs_death_time, gs_death) ~ sim, 
                        data = subgrp_df)  # Death w/ function

fit_gl_death <- survfit(Surv(time = gl_death_time, gl_death) ~ sim, 
                        data = subgrp_df)  # Death after graft loss

### Survival Plots ###
gl_all <- ggsurvplot(fit_gl,
                     censor = FALSE,
                     conf.int = TRUE,
                     title = "Observed vs. Simulated: Graft Loss",
                     legend.title = "",
                     legend.labs = c("Actual", "Simulated"),
                     xlab = "Months",
                     ylim = c(0.6, 1))
gl_all$plot <- gl_all$plot +
  theme(legend.text = element_text(size = 7)) +
  theme(axis.text = element_text(size = 7))
ggsave(filename = paste0(out_path, "gl_surv.png"))

gs_death_all <- ggsurvplot(fit_gs_death,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Death with Function",
                           legend.title = "",
                           legend.labs = c("Actual", "Simulated"),
                           xlab = "Months",
                           ylim = c(0, 1))
gs_death_all$plot <- gs_death_all$plot +
  theme(legend.text = element_text(size = 7)) +
  theme(axis.text = element_text(size = 7))
ggsave(filename = paste0(out_path, "gs_death_surv.png"))

gl_death_all <- ggsurvplot(fit_gl_death,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Death after Graft Loss",
                           legend.title = "",
                           legend.labs = c("Actual", "Simulated"),
                           xlab = "Months",
                           ylim = c(0, 1))
gl_death_all$plot <- gl_death_all$plot +
  theme(legend.text = element_text(size = 7)) +
  theme(axis.text = element_text(size = 7))
ggsave(filename = paste0(out_path, "gl_death_surv.png"))


### By Sub Group ###
## Race
fit_gl_race <- survfit(Surv(time = gl_time, graft_loss) ~ race + sim,
                       data = subgrp_df)

fit_gs_death_race <- survfit(Surv(time = gs_death_time, gs_death) ~ race + sim, 
                             data = subgrp_df)

fit_gl_death_race <- survfit(Surv(time = gl_death_time, gl_death) ~ race + sim,
                             data = subgrp_df)

## Tx Year
fit_gl_txyr <- survfit(Surv(time = gl_time, graft_loss) ~ tx_year_grp + sim,
                       data = subgrp_df)

fit_gs_death_txyr <- survfit(Surv(time = gs_death_time, gs_death) ~ tx_year_grp + sim, 
                             data = subgrp_df)

fit_gl_death_txyr <- survfit(Surv(time = gl_death_time, gl_death) ~ tx_year_grp + sim,
                             data = subgrp_df)

## KDPI category
fit_gl_kdpi <- survfit(Surv(time = gl_time, graft_loss) ~ kdpi_cat + sim,
                       data = subgrp_df)

fit_gs_death_kdpi <- survfit(Surv(time = gs_death_time, gs_death) ~ kdpi_cat + sim, 
                             data = subgrp_df)

fit_gl_death_kdpi <- survfit(Surv(time = gl_death_time, gl_death) ~ kdpi_cat + sim,
                             data = subgrp_df)

## Diabetes Status
fit_gl_diab <- survfit(Surv(time = gl_time, graft_loss) ~ diab_stat + sim,
                       data = subgrp_df)

fit_gs_death_diab <- survfit(Surv(time = gs_death_time, gs_death) ~ diab_stat + sim, 
                             data = subgrp_df)

fit_gl_death_diab <- survfit(Surv(time = gl_death_time, gl_death) ~ diab_stat + sim,
                             data = subgrp_df)

### By Race
## graft loss
gl_race_plot <- ggsurvplot(fit_gl_race,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Graft Loss",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           ylim = c(0.5, 1),
                           facet.by = "race")
gl_race_plot$plot <- gl_race_plot$plot + 
  theme(legend.text = element_text(size = 7)) +
  theme(axis.text = element_text(size = 7))
gl_race_plot
ggsave(filename = paste0(out_path, "gl_race.png"))

## death w/ functioning graft
gs_death_race_plot <- ggsurvplot(fit_gs_death_race,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Death w/ Functioning Graft",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           facet.by = "race")
gs_death_race_plot$plot <- gs_death_race_plot$plot + 
  theme(legend.text = element_text(size = 7))
gs_death_race_plot
ggsave(filename = paste0(out_path, "gs_death_race.png"))
## death after graft loss
gl_death_race_plot <- ggsurvplot(fit_gl_death_race,
                                 censor = FALSE,
                                 conf.int = TRUE,
                                 title = "Observed vs. Simulated: Death after Graft Loss",
                                 legend.title = "",
                                 short.panel.labs = TRUE,
                                 xlab = "Months",
                                 facet.by = "race",
                                 xlim = c(0, 40),
                                 break.x.by = 8)
gl_death_race_plot$plot <- gl_death_race_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_death_race_plot
ggsave(filename = paste0(out_path, "gl_death_race.png"))

### By Tx Year
## graft loss
gl_txyr_plot <- ggsurvplot(fit_gl_txyr,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Graft Loss",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           ylim = c(0.7, 1),
                           xlim = c(0, 60),
                           break.x.by = 12,
                           facet.by = "tx_year_grp")
gl_txyr_plot$plot <- gl_txyr_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_txyr_plot
ggsave(filename = paste0(out_path, "gl_txyr.png"))

## death w/ functioning graft
gs_death_txyr_plot <- ggsurvplot(fit_gs_death_txyr,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Death with Function",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           ylim = c(0.25, 1),
                           xlim = c(0, 60),
                           break.x.by = 12,
                           facet.by = "tx_year_grp")
gs_death_txyr_plot$plot <- gs_death_txyr_plot$plot + 
  theme(legend.text = element_text(size = 7))
gs_death_txyr_plot
ggsave(filename = paste0(out_path, "gs_death_txyr.png")) 

## death after graft loss
gl_death_txyr_plot <- ggsurvplot(fit_gl_death_txyr,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Death after Graft Loss",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           ylim = c(0, 1),
                           xlim = c(0, 36),
                           break.x.by = 6,
                           facet.by = "tx_year_grp")
gl_death_txyr_plot$plot <- gl_death_txyr_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_death_txyr_plot
ggsave(filename = paste0(out_path, "gl_death_txyr.png")) 

### By KDPI category
## graft loss
gl_kdpi_plot <- ggsurvplot(fit_gl_kdpi,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Graft Loss",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           ylim = c(0.5, 1),
                           xlim = c(0, 120),
                           break.x.by = 30,
                           facet.by = "kdpi_cat")
gl_kdpi_plot$plot <- gl_kdpi_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_kdpi_plot
ggsave(filename = paste0(out_path, "gl_kdpi.png"))

## death w/ functioning graft
gs_death_kdpi_plot <- ggsurvplot(fit_gs_death_kdpi,
                                 censor = FALSE,
                                 conf.int = TRUE,
                                 title = "Observed vs. Simulated: Death w/ Functioning Graft",
                                 legend.title = "",
                                 short.panel.labs = TRUE,
                                 xlab = "Months",
                                 ylim = c(0, 1),
                                 xlim = c(0, 120),
                                 break.x.by = 30,
                                 facet.by = "kdpi_cat")
gs_death_kdpi_plot$plot <- gs_death_kdpi_plot$plot + 
  theme(legend.text = element_text(size = 7))
gs_death_kdpi_plot
ggsave(filename = paste0(out_path, "gs_death_kdpi.png"))

## death after graft loss
gl_death_kdpi_plot <- ggsurvplot(fit_gl_death_kdpi,
                                 censor = FALSE,
                                 conf.int = TRUE,
                                 title = "Observed vs. Simulated: Death after Graft Loss",
                                 legend.title = "",
                                 short.panel.labs = TRUE,
                                 xlab = "Months",
                                 ylim = c(0, 1),
                                 xlim = c(0, 36),
                                 break.x.by = 6,
                                 facet.by = "kdpi_cat")
gl_death_kdpi_plot$plot <- gl_death_kdpi_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_death_kdpi_plot
ggsave(filename = paste0(out_path, "gl_death_kdpi.png"))


### By Diabetes Status
## graft loss
gl_diab_plot <- ggsurvplot(fit_gl_diab,
                           censor = FALSE,
                           conf.int = TRUE,
                           title = "Observed vs. Simulated: Graft Loss",
                           legend.title = "",
                           short.panel.labs = TRUE,
                           xlab = "Months",
                           ylim = c(0.5, 1),
                           xlim = c(0, 120),
                           break.x.by = 30,
                           facet.by = "diab_stat")
gl_diab_plot$plot <- gl_diab_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_diab_plot
ggsave(filename = paste0(out_path, "gl_diab.png"))

## death w/ functioning graft
gs_death_diab_plot <- ggsurvplot(fit_gs_death_diab,
                                 censor = FALSE,
                                 conf.int = TRUE,
                                 title = "Observed vs. Simulated: Death w/ Functioning Graft",
                                 legend.title = "",
                                 short.panel.labs = TRUE,
                                 xlab = "Months",
                                 ylim = c(0, 1),
                                 xlim = c(0, 120),
                                 break.x.by = 30,
                                 facet.by = "diab_stat")
gs_death_diab_plot$plot <- gs_death_diab_plot$plot + 
  theme(legend.text = element_text(size = 7))
gs_death_diab_plot
ggsave(filename = paste0(out_path, "gs_death_diab.png"))

## death after graft loss
gl_death_diab_plot <- ggsurvplot(fit_gl_death_diab,
                                 censor = FALSE,
                                 conf.int = TRUE,
                                 title = "Observed vs. Simulated: Death after Graft Loss",
                                 legend.title = "",
                                 short.panel.labs = TRUE,
                                 xlab = "Months",
                                 ylim = c(0, 1),
                                 xlim = c(0, 36),
                                 break.x.by = 6,
                                 facet.by = "diab_stat")
gl_death_diab_plot$plot <- gl_death_diab_plot$plot + 
  theme(legend.text = element_text(size = 7))
gl_death_diab_plot
ggsave(filename = paste0(out_path, "gl_death_diab.png"))

##### Examine Results #####
summary(fit_gl, times = c(12,24,36,48,60,72,84,96,108,120))
summary(fit_gs_death, times = c(60,72,84,96,108,120))
summary(fit_gl_death, times = c(0,6,12,18,24,30,36))

summary(fit_gl_race, times = c(12,24,36,48,60,120))
summary(fit_gs_death_race, times = c(12,24,36,48,60,120))
summary(fit_gl_death_race, times = c(6,12,18,24,30,36))


##### Scraps #####
# For SMDM
splots <- list()
splots[[1]] <- ggsurvplot(comb_gl,
                          combine = TRUE,
                          censor = FALSE,
                          conf.int = TRUE,
                          conf.int.style = "ribbon",
                          title = "",
                          legend.title = "",
                          legend.labs = c("Observed - Graft Loss", 
                                          "Simulated - Graft Loss",
                                          "Observed - Death w/ Function", 
                                          "Simulated - Death w/ Function"),
                          ylim = c(0, 1),
                          linetype = c("solid", "dotdash",
                                       "solid", "dotdash"),
                          palette = "grey",
                          xlab = "Months")
splots[[1]]$plot <- splots[[1]]$plot + 
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2))

splots[[2]] <- ggsurvplot(comb_gl_death,
                          combine = TRUE,
                          censor = FALSE,
                          conf.int = TRUE,
                          conf.int.style = "ribbon",
                          title = "",
                          legend.title = "",
                          legend.labs = c("Observed - Death after Graft Loss", 
                                          "Simulated - Death after Graft Loss"),
                          ylim = c(0, 1),
                          xlim = c(0, 36),
                          break.x.by = 6,
                          linetype = c("solid", "dotdash"),
                          palette = "grey",
                          xlab = "Months",
                          ylab = "")
splots[[2]]$plot <- splots[[2]]$plot + 
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2))

lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,NA),
             c(3,3,3,NA))
layout_matrix = lay

