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
  mutate(ddtx = ifelse(event_cd==0, 1, 0)) %>%
  mutate(ldtx = ifelse(event_cd==1, 1, 0)) %>%
  mutate(mort = ifelse(event_cd==2, 1, 0)) %>%
  mutate(remove = ifelse(event_cd==3, 1, 0)) %>%
  select(!event_cd) %>%
  mutate(sim = rep(0, 10766)) %>%
  mutate(sim = factor(sim, levels = c(0, 1), 
                      labels = c("Actual", "Simulated"))) %>%
  mutate(list_year_grp = ifelse(list_year_c <= 5, 0, 1)) %>%
  mutate(list_year_grp = factor(list_year_grp, levels = c(0, 1),
                              labels = c("2010-2014", "2015+")))
obs_df$ddtx[is.na(obs_df$ddtx)] <- 0
obs_df$ldtx[is.na(obs_df$ldtx)] <- 0
obs_df$mort[is.na(obs_df$mort)] <- 0
obs_df$remove[is.na(obs_df$remove)] <- 0

sim_df <- sim_df %>%
  mutate(sim = rep(1, 100000)) %>%
  mutate(sim = factor(sim, levels = c(0, 1), 
                      labels = c("Actual", "Simulated"))) %>%
  mutate(list_year_grp = ifelse(list_year_c <= 5, 0, 1)) %>%
  mutate(list_year_grp = factor(list_year_grp, levels = c(0, 1),
                              labels = c("2010-2014", "2015+")))

subgrp_df <- rbind(sim_df, obs_df)
subgrp_df <- subgrp_df %>%
  mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
                           labels = c("Diabetes: No", "Diabetes: Yes")))
  

##### Multinomial Logit Method #####
### Survfit Objects ###
fit_event   <- survfit(Surv(time = months_to_event, event) ~ sim,
                       data = subgrp_df)  # Any Event

fit_ddtx    <- survfit(Surv(time = months_to_event, ddtx) ~ sim,
                       data = subgrp_df)  # Deceased Donor Tx

fit_ldtx    <- survfit(Surv(time = months_to_event, ldtx) ~ sim,
                       data = subgrp_df)  # Living Donor Tx

fit_mort    <- survfit(Surv(time = months_to_event, mort) ~ sim,
                       data = subgrp_df)  # Waitlist Mortality

fit_remove  <- survfit(Surv(time = months_to_event, remove) ~ sim,
                       data = subgrp_df)  # Other Removal

### Survival Plots ###
event_all <- ggsurvplot(fit_event,
                      censor = FALSE,
                      title = "Any Event",
                      conf.int = TRUE,
                      short.panel.labs = TRUE,
                      legend.title = "",
                      ylim = c(0, 1),
                      palette = "grey",
                      xlab = "Months")
event_all$plot <- event_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8))

ddtx_all <- ggsurvplot(fit_ddtx,
                            censor = FALSE,
                            conf.int = TRUE,
                            title = "Deceased Donor Tx",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            palette = "grey",
                            xlab = "Months")
ddtx_all$plot <- ddtx_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8))

ldtx_all <- ggsurvplot(fit_ldtx,
                            censor = FALSE,
                            conf.int = TRUE,
                            title = "Living Donor Tx",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0.4, 1),
                            palette = "grey",
                            xlab = "Months")
ldtx_all$plot <- ldtx_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8))

mort_all <- ggsurvplot(fit_mort,
                       censor = FALSE,
                       conf.int = TRUE,
                       title = "Wait List Mortality",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0.4, 1),
                       palette = "grey",
                       xlab = "Months")
mort_all$plot <- mort_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8))

remove_all <- ggsurvplot(fit_remove,
                       censor = FALSE,
                       conf.int = TRUE,
                       title = "Wait List Removal",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0, 1),
                       palette = "grey",
                       xlab = "Months")
remove_all$plot <- remove_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8))

plots_all <- grid.arrange(arrangeGrob(ddtx_all$plot, ldtx_all$plot, 
                                      mort_all$plot, remove_all$plot + 
                                        theme(legend.position="none"), 
                                      nrow = 2),
                          nrow=2,heights=c(10, 1))

ggsave(filename = paste0(out_path, "list_mlogit.png"), plot = plots_all , width = 7, height = 5)

##### Parametric Method #####
### Survfit Objects ###
fit_event   <- survfit(Surv(time = months_to_event, event) ~ sim,
                       data = subgrp_df)  # Any Event

fit_ddtx    <- survfit(Surv(time = months_to_event, ddtx) ~ sim,
                       data = subgrp_df)  # Deceased Donor Tx

fit_ldtx    <- survfit(Surv(time = months_to_event, ldtx) ~ sim,
                       data = subgrp_df)  # Living Donor Tx

fit_mort    <- survfit(Surv(time = months_to_event, mort) ~ sim,
                       data = subgrp_df)  # Waitlist Mortality

fit_remove  <- survfit(Surv(time = months_to_event, remove) ~ sim,
                       data = subgrp_df)  # Other Removal

### Survival Plots ###
event_all <- ggsurvplot(fit_event,
                        censor = FALSE,
                        title = "Any Event",
                        conf.int = TRUE,
                        short.panel.labs = TRUE,
                        legend.title = "",
                        ylim = c(0, 1),
                        palette = "grey",
                        xlab = "Months")
event_all$plot <- event_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8))

ddtx_all <- ggsurvplot(fit_ddtx,
                       censor = FALSE,
                       conf.int = TRUE,
                       title = "Deceased Donor Tx",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0, 1),
                       palette = "grey",
                       xlab = "Months")
ddtx_all$plot <- ddtx_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8))

ldtx_all <- ggsurvplot(fit_ldtx,
                       censor = FALSE,
                       conf.int = TRUE,
                       title = "Living Donor Tx",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0.6, 1),
                       palette = "grey",
                       xlab = "Months")
ldtx_all$plot <- ldtx_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8))

mort_all <- ggsurvplot(fit_mort,
                       censor = FALSE,
                       conf.int = TRUE,
                       title = "Wait List Mortality",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0.4, 1),
                       palette = "grey",
                       xlab = "Months")
mort_all$plot <- mort_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8))

remove_all <- ggsurvplot(fit_remove,
                         censor = FALSE,
                         conf.int = TRUE,
                         title = "Wait List Removal",
                         short.panel.labs = TRUE,
                         legend.title = "",
                         ylim = c(0, 1),
                         palette = "grey",
                         xlab = "Months")
remove_all$plot <- remove_all$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8))

plots_all <- grid.arrange(arrangeGrob(ddtx_all$plot, ldtx_all$plot, 
                                      mort_all$plot, remove_all$plot + 
                                        theme(legend.position="none"), 
                                      nrow = 2),
                          nrow=2,heights=c(10, 1))

ggsave(filename = paste0(out_path, "list_parametric.png"), plot = plots_all , width = 7, height = 5)

