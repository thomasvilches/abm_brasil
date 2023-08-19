

library(dplyr)
library(ggplot2)
library(rcartocolor)

setwd("/data/thomas-covid/Brazil/")

herd = c("5", "10", "20", "30")
beta = c("0468", "04695", "053", "0607")

Vaccines = c("No Vaccine", "CoronaVac", "Astrazeneca", "Pfizer")
Vac_efs = c("0", "5038", "7042", "94") # efficacy against symptomatic

names(Vac_efs) = Vaccines

redp  = c("0.0", "0.5", "1.0") #reduction compared to symptomatic
ef_sev = c("0.0", "1.0")

risk_red = c("0.0", "1.0")
vac_rate = c(300, 600)

cores = rcartocolor::carto_pal(12, name = "Bold")

read_file_incidence <- function(folder, type, beta, herdi, efsymp = "0", redp = "0.0", risk_red = "0.0", ef_sev = "0.0", vac_rate = 300, ag = "all"){
  
  folder_sim = paste0(folder, "/results_prob_0_", beta,"_vac_0_", efsymp, "_herd_immu_", herdi,
                      "_redp_", redp, "_fmild_1.0_rp_", risk_red, "_sev_", ef_sev, "_saopaulo_", vac_rate)
  
  data.cases1 = read.table(paste0(folder_sim, "/simlevel_",type,"_inc_",ag,".dat"),',',h = T) 
  data.cases1 = as.matrix(data.cases1[,-1])
  
  return(data.cases1)
}


read_file_incidence2 <- function(folder, type, beta, herdi, efsymp = "0", redp = "0.0", risk_red = "0.0", ef_sev = "0.0", vac_rate = 300, ag = "all", initial = 1){
  
  folder_sim = paste0(folder, "/results_prob_0_", beta,"_vac_0_", efsymp, "_herd_immu_", herdi,
                      "_redp_", redp, "_fmild_1.0_rp_", risk_red, "_sev_", ef_sev, "_saopaulo_", vac_rate, "_", initial)
  
  # print(folder_sim)
  data.cases1 = read.table(paste0(folder_sim, "/simlevel_",type,"_inc_",ag,".dat"),',',h = T) 
  data.cases1 = as.matrix(data.cases1[,-1])
  
  return(data.cases1)
}

fc <- function(x, i) mean(x[i])

return_boot_average <- function(M, n = 1000){
  
  total <- as.vector(colSums(M))
  medias <- boot::boot(total, fc, n)
  medias <- as.vector(medias$t)
  return(medias)
}

create_df <- function(n, folder, type, beta, herdi, Vaccine = "No Vaccine", redp = "0.0", risk_red = "0.0", ef_sev = "0.0", vac_rate = 300, ag = "all", initial = 1){
  
  if(initial == 1){
    
    sims = read_file_incidence(folder, type, beta, herdi, Vac_efs[Vaccine], redp, risk_red, ef_sev, vac_rate, ag)
  
    }else{
      
      sims = read_file_incidence2(folder, type, beta, herdi, Vac_efs[Vaccine], redp, risk_red, ef_sev, vac_rate, ag, initial)
      
  }
  medias <- return_boot_average(sims, n)
  df <- data.frame(result = medias)
  df$type <- type
  df$beta <- beta
  df$herdi <- herdi
  df$Vaccine <- Vaccine
  df$redp <- redp
  df$risk_red <- risk_red
  df$ef_sev <- ef_sev
  df$vac_rate <- vac_rate
  df$ag <- ag
  return(df)
  
}

variable_names <- list(
  "lat" = "Infections" ,
  "hos" = "Hospitalizations",
  "icu" = "ICU",
  "ded" = "Deaths"
)


variable_labeller <- function(variable,value){
  return(variable_names[value])
}


create_df_initial <- function(n, folder, type, beta, herdi, Vaccine = "No Vaccine", redp = "0.0", risk_red = "0.0", ef_sev = "0.0", vac_rate = 300, ag = "all", initial = 1){
  
  if(initial == 1){
    
    sims = read_file_incidence(folder, type, beta, herdi, Vac_efs[Vaccine], redp, risk_red, ef_sev, vac_rate, ag)
    
  }else{
    # print("entrou aqui")
    sims = read_file_incidence2(folder, type, beta, herdi, Vac_efs[Vaccine], redp, risk_red, ef_sev, vac_rate, ag, initial)
    
  }
  medias <- sum(sims[,-1])/ncol(sims[,-1])
  df <- data.frame(result = medias)
  df$type <- type
  df$beta <- beta
  df$herdi <- herdi
  df$Vaccine <- Vaccine
  df$redp <- redp
  df$risk_red <- risk_red
  df$ef_sev <- ef_sev
  df$vac_rate <- vac_rate
  df$ag <- ag
  df$initial <- initial
  return(df)
}



return_plot <- function(df_plot, type, df_l = NULL){
  
  if(is.null(df_l))
    p_plot <- ggplot(df_plot)+
      geom_jitter(aes(x = scenario, y = result, color = scenario), alpha = 0.5, width = 0.3)+
      geom_boxplot(aes(x = scenario, y = result), fill = NA, size = 1.2, width = 0.3)+
      #facet_wrap(.~type, labeller = variable_labeller, scales = "free_y")+
      scale_color_carto_d(palette = "Bold", name = "Scenarios")+
      labs(y = "Number of outcomes occurences", title = variable_names[[type]])+
      guides(color = guide_legend(override.aes = list(size=3)))+
      theme_bw()+
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")
        
      )
  else
    p_plot <- ggplot(df_plot)+
      geom_jitter(aes(x = scenario, y = result, color = scenario), alpha = 0.5, width = 0.3)+
      geom_boxplot(aes(x = scenario, y = result), fill = NA, size = 1.2, width = 0.3)+
      geom_blank(data = df_l, aes(x = scenario, y = lim_max))+
      geom_blank(data = df_l, aes(x = scenario, y = lim_min))+
      #facet_wrap(.~type, labeller = variable_labeller, scales = "free_y")+
      scale_color_carto_d(palette = "Bold", name = "Scenarios")+
      labs(y = "Number of outcomes occurences", title = variable_names[[type]])+
      guides(color = guide_legend(override.aes = list(size=3)))+
      theme_bw()+
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")
        
      )
  return(p_plot)
}


set.seed(10023)
n = 1000
i = 3 # herd immunity

# Simulations mean --------------------------------------------------------

folder <- "normal_sims/"

# no vac
dfn <- create_df(1000, folder, "lat", beta[i], herd[i])
dfn <- bind_rows(dfn,create_df(1000, folder, "hos", beta[i], herd[i]))
dfn <- bind_rows(dfn,create_df(1000, folder, "icu", beta[i], herd[i]))
dfn <- bind_rows(dfn,create_df(1000, folder, "ded", beta[i], herd[i]))



#CoronaVac

df <- dfn
Vaccine = "CoronaVac"

types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0", "1.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
          
          df1 <- create_df(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 300)
          df <- bind_rows(df, df1)
        
      }
    }
  }
}


types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
       
          
          df1 <- create_df(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 600)
          df <- bind_rows(df, df1)
        
      }
    }
  }
}


write.csv(df, "output/save_mean_sims_corona.csv", row.names = F)



####



#Astrazeneca


df <- dfn
Vaccine = "Astrazeneca"



types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0", "1.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){

          
          df1 <- create_df(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 300)
          df <- bind_rows(df, df1)
        
      }
    }
  }
}


types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
          
          df1 <- create_df(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 600)
          df <- bind_rows(df, df1)
        
      }
    }
  }
}



write.csv(df, "output/save_mean_sims_astra.csv", row.names = F)



#Pfizer


df <- dfn
Vaccine = "Pfizer"

types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0", "1.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
          
          df1 <- create_df(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 300)
          df <- bind_rows(df, df1)
        
      }
    }
  }
}


types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
          
          df1 <- create_df(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 600)
          df <- bind_rows(df, df1)
        
      }
    }
  }
}


write.csv(df, "output/save_mean_sims_pfi.csv", row.names = F)



# Making plots ------------------------------------------------------------

df <- read.csv("output/save_mean_sims_corona.csv", colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))

glimpse(df)

## Corona
df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))


df_plot <- df %>% 
  filter(risk_red == "0.0", vac_rate == 300) %>% 
  mutate(
    scenario = case_when(
      Vaccine == "No Vaccine" ~ "No Vaccine",
      TRUE ~ paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%")
    )
  )

df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                        "50% - 100%", "100% - 0%", "100% - 100%")
)

df_l_inf <- data.frame(scenario = unique(df_plot$scenario))
df_l_inf$lim_min <- 16000
df_l_inf$lim_max <- 26000

df_l_hos <- data.frame(scenario = unique(df_plot$scenario))
df_l_hos$lim_min <- 125
df_l_hos$lim_max <- 350

df_l_icu <- data.frame(scenario = unique(df_plot$scenario))
df_l_icu$lim_min <- 50
df_l_icu$lim_max <- 160

df_l_ded <- data.frame(scenario = unique(df_plot$scenario))
df_l_ded$lim_min <- 40
df_l_ded$lim_max <- 140

df_p <- df_plot %>% filter(type == "lat")
p1 <- return_plot(df_p, "lat", df_l_inf)+scale_y_continuous(n.breaks = 3)
df_p <- df_plot %>% filter(type == "hos")
p2 <- return_plot(df_p, "hos", df_l_hos)
df_p <- df_plot %>% filter(type == "icu")
p3 <- return_plot(df_p, "icu", df_l_icu)
df_p <- df_plot %>% filter(type == "ded")
p4 <- return_plot(df_p, "ded", df_l_ded)

p_plot <- ggpubr::ggarrange(p1,p2+theme(axis.title.y = element_blank()),p3,
                            p4+theme(axis.title.y = element_blank()), common.legend = TRUE, ncol = 2, nrow = 2,
                            legend = "right")



ggsave(plot = p_plot, "output/newFigures/outcomes_corona_vac.png", device = "png", dpi = 300, width = 6.5, height = 5)  


## Astra

df <- read.csv("output/save_mean_sims_astra.csv", colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))

glimpse(df)

df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))


df_plot <- df %>% 
  filter(risk_red == "0.0", vac_rate == 300) %>% 
  mutate(
    scenario = case_when(
      Vaccine == "No Vaccine" ~ "No Vaccine",
      TRUE ~ paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%")
    )
  )

df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                        "50% - 100%", "100% - 0%", "100% - 100%")
)


df_p <- df_plot %>% filter(type == "lat")
p1 <- return_plot(df_p, "lat", df_l_inf)+scale_y_continuous(n.breaks = 3)
df_p <- df_plot %>% filter(type == "hos")
p2 <- return_plot(df_p, "hos", df_l_hos)
df_p <- df_plot %>% filter(type == "icu")
p3 <- return_plot(df_p, "icu", df_l_icu)
df_p <- df_plot %>% filter(type == "ded")
p4 <- return_plot(df_p, "ded", df_l_ded)

p_plot <- ggpubr::ggarrange(p1,p2+theme(axis.title.y = element_blank()),p3,
                            p4+theme(axis.title.y = element_blank()), common.legend = TRUE, ncol = 2, nrow = 2,
                            legend = "right")


ggsave(plot = p_plot, "output/newFigures/outcomes_astra_vac.png", device = "png", dpi = 300, width = 6.5, height = 5)  



## Pfizer


df <- read.csv("output/save_mean_sims_pfi.csv", colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))


glimpse(df)

df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))


df_plot <- df %>% 
  filter(risk_red == "0.0", vac_rate == 300) %>% 
  mutate(
    scenario = case_when(
      Vaccine == "No Vaccine" ~ "No Vaccine",
      TRUE ~ paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%")
    )
  )

df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                        "50% - 100%", "100% - 0%", "100% - 100%")
)


df_p <- df_plot %>% filter(type == "lat")
p1 <- return_plot(df_p, "lat", df_l_inf)+scale_y_continuous(n.breaks = 3)
df_p <- df_plot %>% filter(type == "hos")
p2 <- return_plot(df_p, "hos", df_l_hos)
df_p <- df_plot %>% filter(type == "icu")
p3 <- return_plot(df_p, "icu", df_l_icu)
df_p <- df_plot %>% filter(type == "ded")
p4 <- return_plot(df_p, "ded", df_l_ded)

p_plot <- ggpubr::ggarrange(p1,p2+theme(axis.title.y = element_blank()),p3,
                            p4+theme(axis.title.y = element_blank()), common.legend = TRUE, ncol = 2, nrow = 2,
                            legend = "right")

ggsave(plot = p_plot, "output/newFigures/outcomes_pfi_vac.png", device = "png", dpi = 300, width = 6.5, height = 5)  



# Plot reduction ----------------------------------------------------------

vac <- "corona"

df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))

glimpse(df)

df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))

dfnovac <- df %>% filter(Vaccine == "No Vaccine")

df_plot <- df %>% filter(Vaccine != "No Vaccine")


df_plot <- df_plot %>% 
  mutate(
    scenario = paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%"),
    scenario_y = case_when(
      vac_rate == "600" ~ "Double vaccination rate",
      risk_red == "1.0" ~ "Low risk perception",
      TRUE ~ "Baseline"
    )
  ) %>% 
  group_by(scenario_y, type, scenario) %>% 
  mutate(ID = row_number()) %>% 
  ungroup()


df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                        "50% - 100%", "100% - 0%", "100% - 100%")
)



df_plot <- dfnovac %>%
  select(result_nvac = result, type) %>% 
  group_by(type) %>% 
  mutate(ID = row_number()) %>%
  ungroup() %>% 
  right_join(
    df_plot, by = c("ID", "type")
  ) %>% 
  mutate(
    Reduction = (result_nvac - result)/result_nvac
  ) %>% select(-result_nvac)





## plot for Baseline


limits_min <- c("lat" = -0.15, "hos" = 0.0, "icu" = 0.0, "ded" = 0.10)
limits_max <- c("lat" = 0.35, "hos" = 0.6, "icu" = 0.65, "ded" = 0.65)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Baseline")

df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

 
  p_plot <- ggplot(df_plot1)+
    geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
    geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
    scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
    facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
    geom_blank(aes(y = lim_min))+
    geom_blank(aes(y = lim_max))+
    scale_color_manual(values = cores[2:7], name = "Scenarios")+
    guides(color = guide_legend(override.aes = list(size=3)))+
    #scale_y_continuous(expand = c(0,0))+
    labs(y = "Reduction")+
    theme_bw()+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      # legend.key.size = unit(2.0, 'cm'),
      legend.position = "bottom",
      legend.box.spacing = unit(1.5, "cm")
      
    )
  





ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_base.png"), device = "png", dpi = 300, width = 8, height = 5.5)  



## plot for Double

limits_min <- c("lat" = -0.15, "hos" = 0.2, "icu" = 0.2, "ded" = 0.2)
limits_max <- c("lat" = 0.6, "hos" = 0.8, "icu" = 0.8, "ded" = 0.85)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Double vaccination rate")

df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom",
    legend.box.spacing = unit(1.5, "cm")
    
  )



ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_double.png"), device = "png", dpi = 300, width = 8, height = 5.5)  





## plot for Risk perception

limits_min <- c("lat" = -0.2, "hos" = -0.05, "icu" = 0.0, "ded" = 0.05)
limits_max <- c("lat" = 0.35, "hos" = 0.6, "icu" = 0.6, "ded" = 0.65)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Low risk perception")


df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_low.png"), device = "png", dpi = 300, width = 8, height = 5)  


####################################################################################

vac = "astra"

df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))

glimpse(df)

df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))

dfnovac <- df %>% filter(Vaccine == "No Vaccine")

df_plot <- df %>% filter(Vaccine != "No Vaccine")


df_plot <- df_plot %>% 
  mutate(
    scenario = paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%"),
    scenario_y = case_when(
      vac_rate == "600" ~ "Double vaccination rate",
      risk_red == "1.0" ~ "Low risk perception",
      TRUE ~ "Baseline"
    )
  ) %>% 
  group_by(scenario_y, type, scenario) %>% 
  mutate(ID = row_number()) %>% 
  ungroup()


df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                        "50% - 100%", "100% - 0%", "100% - 100%")
)



df_plot <- dfnovac %>%
  select(result_nvac = result, type) %>% 
  group_by(type) %>% 
  mutate(ID = row_number()) %>%
  ungroup() %>% 
  right_join(
    df_plot, by = c("ID", "type")
  ) %>% 
  mutate(
    Reduction = (result_nvac - result)/result_nvac
  ) %>% select(-result_nvac)





## plot for Baseline


limits_min <- c("lat" = -0.15, "hos" = 0.0, "icu" = 0.0, "ded" = 0.10)
limits_max <- c("lat" = 0.35, "hos" = 0.6, "icu" = 0.65, "ded" = 0.65)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Baseline")

df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom",
    legend.box.spacing = unit(1.5, "cm")
    
  )


ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_base.png"), device = "png", dpi = 300, width = 8, height = 5.5)  




## plot for Double

limits_min <- c("lat" = -0.15, "hos" = 0.2, "icu" = 0.2, "ded" = 0.2)
limits_max <- c("lat" = 0.6, "hos" = 0.8, "icu" = 0.8, "ded" = 0.85)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Double vaccination rate")

df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom",
    legend.box.spacing = unit(1.5, "cm")
    
  )
ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_double.png"), device = "png", dpi = 300, width = 8, height = 5.5)


## plot for Risk perception

limits_min <- c("lat" = -0.2, "hos" = -0.05, "icu" = 0.0, "ded" = 0.05)
limits_max <- c("lat" = 0.35, "hos" = 0.6, "icu" = 0.6, "ded" = 0.65)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Low risk perception")


df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_low.png"), device = "png", dpi = 300, width = 8, height = 5)  



####################################################################################

vac = "pfi"

df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))

glimpse(df)

df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))

dfnovac <- df %>% filter(Vaccine == "No Vaccine")

df_plot <- df %>% filter(Vaccine != "No Vaccine")


df_plot <- df_plot %>% 
  mutate(
    scenario = paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%"),
    scenario_y = case_when(
      vac_rate == "600" ~ "Double vaccination rate",
      risk_red == "1.0" ~ "Low risk perception",
      TRUE ~ "Baseline"
    )
  ) %>% 
  group_by(scenario_y, type, scenario) %>% 
  mutate(ID = row_number()) %>% 
  ungroup()


df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                        "50% - 100%", "100% - 0%", "100% - 100%")
)



df_plot <- dfnovac %>%
  select(result_nvac = result, type) %>% 
  group_by(type) %>% 
  mutate(ID = row_number()) %>%
  ungroup() %>% 
  right_join(
    df_plot, by = c("ID", "type")
  ) %>% 
  mutate(
    Reduction = (result_nvac - result)/result_nvac
  ) %>% select(-result_nvac)





## plot for Baseline

limits_min <- c("lat" = -0.15, "hos" = 0.0, "icu" = 0.0, "ded" = 0.10)
limits_max <- c("lat" = 0.35, "hos" = 0.6, "icu" = 0.65, "ded" = 0.65)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Baseline")

df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom",
    legend.box.spacing = unit(1.5, "cm")
    
  )




ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_base.png"), device = "png", dpi = 300, width = 8, height = 5.5)  




## plot for Double

limits_min <- c("lat" = -0.15, "hos" = 0.2, "icu" = 0.2, "ded" = 0.2)
limits_max <- c("lat" = 0.6, "hos" = 0.8, "icu" = 0.8, "ded" = 0.85)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Double vaccination rate")

df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom",
    legend.box.spacing = unit(1.5, "cm")
    #plot.margin = unit(c(1,1,1.5,1), "cm")
    
  )

ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_double.png"), device = "png", dpi = 300, width = 8, height = 5.5)  




## plot for Risk perception

limits_min <- c("lat" = -0.2, "hos" = -0.05, "icu" = 0.0, "ded" = 0.05)
limits_max <- c("lat" = 0.35, "hos" = 0.6, "icu" = 0.6, "ded" = 0.65)


df_plot1 <- df_plot %>% 
  filter(scenario_y == "Low risk perception")


df_plot1$lim_min = limits_min[df_plot1$type]
df_plot1$lim_max = limits_max[df_plot1$type]

# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())


p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1), expand = c(0,0))+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  geom_blank(aes(y = lim_min))+
  geom_blank(aes(y = lim_max))+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(expand = c(0,0))+
  labs(y = "Reduction")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    # legend.key.size = unit(2.0, 'cm'),
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/newFigures/red_", vac,"_vac_low.png"), device = "png", dpi = 300, width = 8, height = 5)  




# Initial -----------------------------------------------------------------

set.seed(10023)
n = 1000
i = 3 # herd immunity

# Simulations mean --------------------------------------------------------

folder <- "change_init/"


create_table_red_initial <- function(df){
  
  df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))
  
  dfnovac <- df %>% filter(Vaccine == "No Vaccine")
  
  df_plot <- df %>% filter(Vaccine != "No Vaccine")
  
  
  df_plot <- df_plot %>% 
    mutate(
      scenario = paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%"),
      scenario_y = case_when(
        vac_rate == "600" ~ "Double vaccination rate",
        risk_red == "1.0" ~ "Low risk perception",
        TRUE ~ "Baseline"
      )
    ) %>% 
    group_by(scenario_y, type, scenario) %>% 
    mutate(ID = row_number()) %>% 
    ungroup()
  
  
  df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                          "50% - 100%", "100% - 0%", "100% - 100%")
  )
  
  
  
  df_plot <- dfnovac %>%
    select(result_nvac = result, type, initial) %>% 
    # group_by(type) %>% 
    # mutate(ID = row_number()) %>%
    # ungroup() %>% 
    right_join(
      df_plot, by = c("initial", "type")
    ) %>% 
    mutate(
      Reduction = (result_nvac - result)/result_nvac
    ) %>% select(-result_nvac)
  return(df_plot)
}

create_plot_red_initial <- function(df_plot, limits_min, limits_max){

  
  df_plot$lim_min = limits_min[df_plot$type]
  df_plot$lim_max = limits_max[df_plot$type]
  
  
 df_plot %>% 
    filter(scenario_y == "Baseline") %>% 
    ggplot()+
    #geom_jitter(aes(x = initial, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
    geom_point(aes(x = initial, y = Reduction, color = scenario), size = 1.8)+
    geom_line(aes(x = initial, y = Reduction, color = scenario), size = 0.8, linetype = 'dashed')+
     scale_y_continuous(labels = scales::percent, expand = c(0,0))+
     facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
     geom_blank(aes(y = lim_min))+
     geom_blank(aes(y = lim_max))+
    scale_color_manual(values = cores[2:7], name = "Scenarios")+
    guides(color = guide_legend(override.aes = list(size=3)))+
    #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
    labs(y = "Reduction", x = "Initial number of infections")+
    theme_bw()+
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
      # axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "bottom"
      
    )
  
}

ii = 5
# no vac
dfn <- create_df_initial(1000, folder, "lat", beta[i], herd[i], initial = ii)
dfn <- bind_rows(dfn,create_df_initial(1000, folder, "hos", beta[i], herd[i], initial = ii))
dfn <- bind_rows(dfn,create_df_initial(1000, folder, "icu", beta[i], herd[i], initial = ii))
dfn <- bind_rows(dfn,create_df_initial(1000, folder, "ded", beta[i], herd[i], initial = ii))

dff_n <- dfn

initials <- seq(5,45,10)


for(ii in initials[-1]){
  
  dfn <- create_df_initial(1000, folder, "lat", beta[i], herd[i], initial = ii)
  dfn <- bind_rows(dfn,create_df_initial(1000, folder, "hos", beta[i], herd[i], initial = ii))
  dfn <- bind_rows(dfn,create_df_initial(1000, folder, "icu", beta[i], herd[i], initial = ii))
  dfn <- bind_rows(dfn,create_df_initial(1000, folder, "ded", beta[i], herd[i], initial = ii))
  
  dff_n <- bind_rows(dfn, dff_n)
}

dff_n


########
df <- dff_n
Vaccine = "CoronaVac"

types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
        
        for(ii in initials[-1]){
          
          df1 <- create_df_initial(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 300, initial = ii)
          df <- bind_rows(df, df1)
        }
        
      }
    }
  }
}


limits_min <- c("lat" = 0.01, "hos" = 0.05, "icu" = 0.05, "ded" = 0.1)
limits_max <- c("lat" = 0.17, "hos" = 0.4, "icu" = 0.45, "ded" = 0.55)



df_plot <- create_table_red_initial(df)
p_plot <- create_plot_red_initial(df_plot, limits_min, limits_max)

vac = "corona"
ggsave(plot = p_plot, paste0("output/newFigures/red_initial_", vac,"_vac.png"), device = "png", dpi = 300, width = 8, height = 5)  



########
df <- dff_n
Vaccine = "Astrazeneca"

types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
        
        for(ii in initials[-1]){
          
          df1 <- create_df_initial(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 300, initial = ii)
          df <- bind_rows(df, df1)
        }
        
      }
    }
  }
}

df_plot <- create_table_red_initial(df)
p_plot <- create_plot_red_initial(df_plot, limits_min, limits_max)


vac = "astra"
ggsave(plot = p_plot, paste0("output/newFigures/red_initial_", vac,"_vac.png"), device = "png", dpi = 300, width = 8, height = 5)  


########
df <- dff_n
Vaccine = "Pfizer"

types = c("lat", "hos", "icu", "ded")
redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")

for(t in types){
  for(red_p in redp){
    for(redr in red_risk){
      for(sev in sev_p){
        
        
        for(ii in initials[-1]){
          
          df1 <- create_df_initial(n, folder, t, beta[i], herd[i], Vaccine, red_p, redr, sev, 300, initial = ii)
          df <- bind_rows(df, df1)
        }
        
      }
    }
  }
}

df_plot <- create_table_red_initial(df)
p_plot <- create_plot_red_initial(df_plot, limits_min, limits_max)


vac = "pfi"
ggsave(plot = p_plot, paste0("output/newFigures/red_initial_", vac,"_vac.png"), device = "png", dpi = 300, width = 8, height = 5)  


# Create tables -----------------------------------------------------------



create_table_raw <- function(df){
  df  %>% 
    mutate(
      scenario = paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%"),
      scenario_y = case_when(
        vac_rate == "600" ~ "Double vaccination rate",
        risk_red == "1.0" ~ "Low risk perception",
        TRUE ~ "Baseline"
      )
    ) %>% 
    group_by(Vaccine, scenario, scenario_y, type) %>% 
    summarise(
      media = mean(result),
      ICmin = unname(quantile(result, 0.025))
      ,ICmax = unname(quantile(result, 0.975))
    )
}

create_table_red <- function(df){
  
  df$type <- factor(df$type, levels = c("lat", "hos", "icu", "ded"))
  
  dfnovac <- df %>% filter(Vaccine == "No Vaccine")
  
  df_plot <- df %>% filter(Vaccine != "No Vaccine")
  
  
  df_plot <- df_plot %>% 
    mutate(
      scenario = paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%"),
      scenario_y = case_when(
        vac_rate == "600" ~ "Double vaccination rate",
        risk_red == "1.0" ~ "Low risk perception",
        TRUE ~ "Baseline"
      )
    ) %>% 
    group_by(scenario_y, type, scenario) %>% 
    mutate(ID = row_number()) %>% 
    ungroup()
  
  
  df_plot$scenario <- factor(df_plot$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                          "50% - 100%", "100% - 0%", "100% - 100%")
  )
  
  
  
  df_plot <- dfnovac %>%
    select(result_nvac = result, type) %>% 
    group_by(type) %>% 
    mutate(ID = row_number()) %>%
    ungroup() %>% 
    right_join(
      df_plot, by = c("ID", "type")
    ) %>% 
    mutate(
      Reduction = (result_nvac - result)/result_nvac
    ) %>% select(-result_nvac) %>% 
    group_by(Vaccine, scenario, scenario_y, type) %>% 
    summarise(
      media = mean(Reduction),
      ICmin = unname(quantile(Reduction, 0.025))
      ,ICmax = unname(quantile(Reduction, 0.975))
    )
  return(df_plot)
}

vac = "pfi"
df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))
head(df)
df_out <- create_table_raw(df)

openxlsx::write.xlsx(df_out, paste0("output/results_raw_", vac, ".xlsx"))


vac = "corona"
df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))
head(df)
df_out <- create_table_raw(df)

openxlsx::write.xlsx(df_out, paste0("output/results_raw_", vac, ".xlsx"))


vac = "astra"
df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))
head(df)
df_out <- create_table_raw(df)

openxlsx::write.xlsx(df_out, paste0("output/results_raw_", vac, ".xlsx"))




vac = "astra"
df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))
head(df)
df_out <- create_table_red(df)
openxlsx::write.xlsx(df_out, paste0("output/results_red_", vac, ".xlsx"))


vac = "pfi"
df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))
head(df)
df_out <- create_table_red(df)
openxlsx::write.xlsx(df_out, paste0("output/results_red_", vac, ".xlsx"))


vac = "corona"
df <- read.csv(paste0("output/save_mean_sims_", vac,".csv"), colClasses = c("numeric", "character", "character", "character","character", "character", "character", "character", "character", "character"))
head(df)
df_out <- create_table_red(df)
openxlsx::write.xlsx(df_out, paste0("output/results_red_", vac, ".xlsx"))






# YLL ---------------------------------------------------------------------



read_data_exp <- function(folder, betai, herdi, efsymp = "0", redp = "0.0", risk_red = "0.0", ef_sev = "0.0", vac_rate = 300, ag = "all"){
  
  folder_sim = paste0(folder, "/results_prob_0_", betai,"_vac_0_", efsymp, "_herd_immu_", herdi,
                      "_redp_", redp, "_fmild_1.0_rp_", risk_red, "_sev_", ef_sev, "_saopaulo_", vac_rate)
  
  data.cases1 = read.table(paste0(folder_sim, "/year_of_death.dat")) 
  data.cases1 = as.matrix(data.cases1)
  
  
  resultados <- purrr::map(1:ncol(data.cases1), function(x) life_exp$Expectancy[1:101]*data.cases1[, x])
  
  resultados <- Reduce(cbind, resultados)
  resultados <- colSums(resultados)
  
  resultados_b <- boot::boot(resultados, function(x, i) mean(x[i]), R = 500)
  
  df <- data.frame(result  = resultados_b$t)
  
  
  df$beta <- betai
  df$herdi <- herdi
  df$Vaccine <- Vaccine
  df$redp <- redp
  df$risk_red <- risk_red
  df$ef_sev <- ef_sev
  df$vac_rate <- vac_rate
  df$idx <- seq(1, nrow(df))
  
  df <- df %>% 
    #filter(risk_red == "0.0", vac_rate == 300) %>% 
    mutate(
      scenario = case_when(
        Vaccine == "No Vaccine" ~ "No Vaccine",
        TRUE ~ paste0(100-as.numeric(redp)*100, "% - ", as.numeric(ef_sev)*100, "%")
      )
    )
  
  df$scenario <- factor(df$scenario, levels = c("No Vaccine", "0% - 0%", "0% - 100%", "50% - 0%", 
                                                          "50% - 100%", "100% - 0%", "100% - 100%")
  )
 df
}


herd = c("5", "10", "20", "30")
beta = c("0468", "04695", "053", "0607")

#Vaccines = c("No Vaccine", "CoronaVac", "Astrazeneca", "Pfizer")
#Vac_efs = c("0", "5038", "7042", "94") # efficacy against symptomatic

# names(Vac_efs) = Vaccines

redp  = c("0.0", "0.5", "1.0") #reduction compared to symptomatic
ef_sev = c("0.0", "1.0")

risk_red = c("0.0", "1.0")
vac_rate = c(300, 600)


gg <- guides(color = guide_legend(override.aes = list(size=3)))

# no vaccine
Vaccine = "No Vaccine"


df_novac <- read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine])


redp <- c("0.0", "0.5", "1.0")
red_risk <- c("0.0")
sev_p <- c("0.0", "1.0")


## CoronaVac
Vaccine = "CoronaVac"

df = NULL
for(red_p in redp){
    for(sev in sev_p){
      
        if(is.null(df))
          df <- read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine], red_p, "0.0", sev, 300)
        else
          df <- bind_rows(df, read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine], red_p, "0.0", sev, 300))
      
    }
  
}


lista <- list(df_novac, df)

df_plot <- Reduce(rbind, lista)

plot_cor <- ggplot(df_plot)+
  geom_jitter(aes(x = scenario, y = result, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = result), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(limits = c(1500, 4000))+
  #facet_wrap(.~type, labeller = variable_labeller, scales = "free_y")+
  scale_color_carto_d(palette = "Bold", name = "Scenarios")+
  labs(y = "Years of life lost", title = "B")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "plain"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
    
  )+gg


## Astrazeneca
Vaccine = "Astrazeneca"

df = NULL
for(red_p in redp){
  for(sev in sev_p){
    
    if(is.null(df))
      df <- read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine], red_p, "0.0", sev, 300)
    else
      df <- bind_rows(df, read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine], red_p, "0.0", sev, 300))
    
  }
  
}


lista <- list(df_novac, df)

df_plot <- Reduce(rbind, lista)

plot_ast <- ggplot(df_plot)+
  geom_jitter(aes(x = scenario, y = result, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = result), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(limits = c(1500, 4000))+
  #facet_wrap(.~type, labeller = variable_labeller, scales = "free_y")+
  scale_color_carto_d(palette = "Bold", name = "Scenarios")+
  labs(y = "Years of life lost", title = "C")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "plain"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
    
  )+gg


## Pfizer
Vaccine = "Pfizer"

df = NULL
for(red_p in redp){
  for(sev in sev_p){
    
    if(is.null(df))
      df <- read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine], red_p, "0.0", sev, 300)
    else
      df <- bind_rows(df, read_data_exp("./", beta[3], herd[3], Vac_efs[Vaccine], red_p, "0.0", sev, 300))
    
  }
  
}


lista <- list(df_novac, df)

df_plot <- Reduce(rbind, lista)

plot_pfi <- ggplot(df_plot)+
  geom_jitter(aes(x = scenario, y = result, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = result), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(limits = c(1500, 4000))+
  #facet_wrap(.~type, labeller = variable_labeller, scales = "free_y")+
  scale_color_carto_d(palette = "Bold", name = "Scenarios")+
  labs(y = "Years of life lost", title = "D")+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "plain"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
    
  )+gg


## juntando



life_exp <- read.table("life_exp.csv", sep = ";", h = FALSE)

names(life_exp) <- c("Age", "Expectancy")
life_exp %>% head()


p.data <- ggplot(life_exp)+
  geom_col(aes(x = Age, y = Expectancy), color = "grey60", fill = "grey60")+
  labs(title = "A", y = "Life expectancy (years)", x = "Age (years)")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "plain"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
    
  )

tt <- theme(axis.title.y = element_blank())


pyll <- ggpubr::ggarrange(ggpubr::ggarrange(p.data), ggpubr::ggarrange(plot_cor+gg, plot_ast+tt+gg,
                                                                       plot_pfi+tt+gg, nrow = 1, widths = c(0.38, .31, .31), common.legend = TRUE, legend = "bottom"), widths = c(0.3, 0.7)
                          ) + ggpubr::bgcolor("white")



ggsave(plot = pyll, paste0("output/newFigures/pyll.png"), device = "png", dpi = 300, width = 10, height = 4)  



# Create table ------------------------------------------------------------


get_prop <- function(Vaccine){
  age_groups <- c("0-4", "5-19", "20-49", "50-64", "65-79", "80-100")
  ags <- paste0("ag", 1:6)
  ags <- c(ags, "all")
  resultados <- lapply(ags, function(x) read_file_incidence("normal_sims/", "lat", beta[3], "20", Vac_efs[Vaccine], "0.0", "0.0", ag = x))
  props <- lapply(resultados, function(x) sum(x[, -1]))
  props <- Reduce(c, props)
  props <- props/props[7]
  
  df <- data.frame(proportion = props[-length(props)], `Age groups` = age_groups)
  df$Vaccine = Vaccine
  df
}

props <- lapply(Vaccines, get_prop)

props <- Reduce(rbind, props)

props <- props %>% 
  tidyr::pivot_wider(names_from = "Age.groups", values_from = "proportion")

xtable::xtable(props)

#[4;19;49;64;79;999]

names(props) <- c("Proportion", "Age Group", "Vaccine")

props$Vaccine <- factor(props$Vaccine, levels = rev(Vaccines))
props$`Age Group` <- factor(props$`Age Group`, levels = c("0-4", "5-19", "20-49", "50-64", "65-79", "80-100"))

ggplot(props)+
  geom_col(aes(x = Proportion, y = Vaccine, color = `Age Group`, fill = `Age Group`))+
  rcartocolor::scale_fill_carto_d(palette = "Bold")+
  rcartocolor::scale_color_carto_d(palette = "Bold")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),#, angle = 90, hjust = 0.5),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    #strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "plain"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
    
  )
ggsave("output/prop_vaccines.png", device = "png", dpi = 300, width = 6.4, height = 3)



