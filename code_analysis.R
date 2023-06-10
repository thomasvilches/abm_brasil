

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



return_plot <- function(df_plot, type){
  
  p_plot <- ggplot(df_plot)+
    geom_jitter(aes(x = scenario, y = result, color = scenario), alpha = 0.5, width = 0.3)+
    geom_boxplot(aes(x = scenario, y = result), fill = NA, size = 1.2, width = 0.3)+
    #facet_wrap(.~type, labeller = variable_labeller, scales = "free_y")+
    scale_color_carto_d(palette = "Bold", name = "Scenarios")+
    labs(y = "Number of outcomes occurences", title = variable_names[[type]])+
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


df_p <- df_plot %>% filter(type == "lat")
p1 <- return_plot(df_p, "lat")+scale_y_continuous(n.breaks = 3)
df_p <- df_plot %>% filter(type == "hos")
p2 <- return_plot(df_p, "hos")
df_p <- df_plot %>% filter(type == "icu")
p3 <- return_plot(df_p, "icu")
df_p <- df_plot %>% filter(type == "ded")
p4 <- return_plot(df_p, "ded")

p_plot <- ggpubr::ggarrange(p1,p2+theme(axis.title.y = element_blank()),p3,
                            p4+theme(axis.title.y = element_blank()), common.legend = TRUE, ncol = 2, nrow = 2,
                            legend = "right")



ggsave(plot = p_plot, "output/outcomes_corona_vac.png", device = "png", dpi = 300, width = 6.5, height = 5)  


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
p1 <- return_plot(df_p, "lat")+scale_y_continuous(n.breaks = 3)
df_p <- df_plot %>% filter(type == "hos")
p2 <- return_plot(df_p, "hos")
df_p <- df_plot %>% filter(type == "icu")
p3 <- return_plot(df_p, "icu")
df_p <- df_plot %>% filter(type == "ded")
p4 <- return_plot(df_p, "ded")

p_plot <- ggpubr::ggarrange(p1,p2+theme(axis.title.y = element_blank()),p3,
                            p4+theme(axis.title.y = element_blank()), common.legend = TRUE, ncol = 2, nrow = 2,
                            legend = "right")


ggsave(plot = p_plot, "output/outcomes_astra_vac.png", device = "png", dpi = 300, width = 6.5, height = 5)  



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
p1 <- return_plot(df_p, "lat")+scale_y_continuous(n.breaks = 3)
df_p <- df_plot %>% filter(type == "hos")
p2 <- return_plot(df_p, "hos")
df_p <- df_plot %>% filter(type == "icu")
p3 <- return_plot(df_p, "icu")
df_p <- df_plot %>% filter(type == "ded")
p4 <- return_plot(df_p, "ded")

p_plot <- ggpubr::ggarrange(p1,p2+theme(axis.title.y = element_blank()),p3,
                            p4+theme(axis.title.y = element_blank()), common.legend = TRUE, ncol = 2, nrow = 2,
                            legend = "right")

ggsave(plot = p_plot, "output/outcomes_pfi_vac.png", device = "png", dpi = 300, width = 6.5, height = 5)  



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

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Baseline")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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


ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_base.png"), device = "png", dpi = 300, width = 8, height = 5)  



## plot for Double

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Double vaccination rate")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7],name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_double.png"), device = "png", dpi = 300, width = 8, height = 5)  




## plot for Risk perception

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Low risk perception")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_low.png"), device = "png", dpi = 300, width = 8, height = 5)  


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

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Baseline")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )


ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_base.png"), device = "png", dpi = 300, width = 8, height = 5)  



## plot for Double

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Double vaccination rate")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_double.png"), device = "png", dpi = 300, width = 8, height = 5)  




## plot for Risk perception

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Low risk perception")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_low.png"), device = "png", dpi = 300, width = 8, height = 5)  



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

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Baseline")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )


ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_base.png"), device = "png", dpi = 300, width = 8, height = 5)  



## plot for Double

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Double vaccination rate")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_double.png"), device = "png", dpi = 300, width = 8, height = 5)  




## plot for Risk perception

df_plot1 <- df_plot %>% 
  filter(scenario_y == "Low risk perception")
# 
# df_plot1 %>% 
#   group_by(scenario, type) %>% 
#   summarise(n = n())

p_plot <- ggplot(df_plot1)+
  geom_jitter(aes(x = scenario, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
  geom_boxplot(aes(x = scenario, y = Reduction), fill = NA, size = 1.2, width = 0.3)+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
  scale_color_manual(values = cores[2:7], name = "Scenarios")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  #scale_y_continuous(breaks = c(18e3, 22e3, 26e3))+
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
    legend.position = "bottom"
    
  )

ggsave(plot = p_plot, paste0("output/red_", vac,"_vac_low.png"), device = "png", dpi = 300, width = 8, height = 5)  




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

create_plot_red_initial <- function(df_plot){
  
 df_plot %>% 
    filter(scenario_y == "Baseline") %>% 
    ggplot()+
    #geom_jitter(aes(x = initial, y = Reduction, color = scenario), alpha = 0.5, width = 0.3)+
    geom_point(aes(x = initial, y = Reduction, color = scenario), size = 1.8)+
    geom_line(aes(x = initial, y = Reduction, color = scenario), size = 0.8, linetype = 'dashed')+
    scale_y_continuous(labels = scales::percent)+
    facet_wrap(.~type, labeller = variable_labeller, scales = "free_y", ncol = 4)+
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

df_plot <- create_table_red_initial(df)
p_plot <- create_plot_red_initial(df_plot)

vac = "corona"
ggsave(plot = p_plot, paste0("output/red_initial_", vac,"_vac.png"), device = "png", dpi = 300, width = 8, height = 5)  



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
p_plot <- create_plot_red_initial(df_plot)


vac = "astra"
ggsave(plot = p_plot, paste0("output/red_initial_", vac,"_vac.png"), device = "png", dpi = 300, width = 8, height = 5)  


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
p_plot <- create_plot_red_initial(df_plot)


vac = "pfi"
ggsave(plot = p_plot, paste0("output/red_initial_", vac,"_vac.png"), device = "png", dpi = 300, width = 8, height = 5)  


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




