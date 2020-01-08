library(knitr)
library(tidyverse)
library(ggthemes)
library(cowplot)
install.packages("cowplot")

rbs_testing <- read.csv(file = "../ChaseWeaver/Desktop/Projects/phage-translation/Data/virus_dataframe.tsv", header = TRUE, sep = "\t")
host_rbs_df <- read.csv(file = "../ChaseWeaver/Desktop/Projects/phage-translation/Data/host_dataframe.tsv", header = TRUE, sep = "\t")

#Ranked sums test on RBS binding
wilcox <- wilcox.test(y = rbs_testing$energy_binding, x = host_rbs_df$energy_binding, exact = FALSE, correct = FALSE)

#fill = "#E69F00"
#fill = ""

#virus and e.coli histogram
ggplot(data = rbs_testing) + 
  geom_histogram(aes(x = energy_binding, y = stat(count/sum(count)),fill="red"),alpha = 0.6, bins = 9) + 
  geom_histogram(data = host_rbs_df, aes(x = energy_binding, y = stat(count/sum(count)),fill="blue"),alpha = 0.6, bins = 9) + 
  xlab("Energy Binding Strength") +
  ylab("Proportion") +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name = "",labels=c("E. coli","Virus"), values=c("#E69F00","#56B4E9")) +
  theme_cowplot(12) + 
  theme(text = element_text(size = 25), legend.text=element_text(size=20), axis.text=element_text(size=20)) 

#ecoli histogram
host_rbs_df %>% 
  ggplot() +
  geom_histogram(aes(x = energy_binding, fill = ""),alpha = 0.6, bins = 9) + 
  xlab("Energy binding strength") +
  scale_fill_manual(name = "", labels = c("E. coli"), values = c("#E69F00")) +
  theme_cowplot(12) + 
  theme(text = element_text(size = 25), legend.text=element_text(size=20), axis.text=element_text(size=20)) +
  scale_y_continuous(expand = c(0,0))



rbs_testing %>% 
  ggplot() +
  geom_histogram(aes(x = energy_binding, fill = "red"),alpha = 0.6, bins = 9) + 
  xlab("Energy binding strength") +
  scale_fill_manual(name= "",labels = c("Virus"), values = c("#56B4E9")) +
  theme_cowplot(12) + 
  theme(text = element_text(size = 25), legend.text=element_text(size=20), axis.text=element_text(size=20)) +
  scale_y_continuous(expand = c(0,0))
# ggplot(data = combined_df) + 
#   geom_density(aes(x = energy_binding, y = ..scaled.., fill=cond), bins = 9,alpha = 0.6) + 
#   xlab("Energy Binding Strength") +
#   ylab("Proportion") +
#   scale_y_continuous() 


bleh_df <- 
  rbs_testing %>% select(energy_binding) %>% mutate(cond="Virus")

host <- 
  host_rbs_df %>% select(energy_binding) %>% mutate(cond="Host")

combined_df <- rbind(bleh_df,host)

stats_df <- read.csv(file = "../ChaseWeaver/Desktop/Projects/phage-translation/Data/stats_df.tsv", header = TRUE, sep = ",")

head(stats_df)

ggplot(data = stats_df) + 
  geom_point(aes(x = mean_difs, y = p_values), color = "#56B4E9") +
  scale_y_continuous(trans = "log10") +
  geom_hline(yintercept = 0.01, color = "red") + 
  geom_vline(xintercept = 0, color = "red") + 
  theme_cowplot(12) + 
  xlab("Mean difference") +
  ylab("P-values") +
  theme(text = element_text(size = 25), axis.text=element_text(size=20))
  

ggplot(data = stats_df) + 
  geom_histogram(aes(x = p_values), fill = "#56B4E9") + 
  geom_vline(xintercept = 0.01, color = "red") +
  xlab("P-values") +
  theme_cowplot(12) +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0)) + 
  theme(text = element_text(size = 25), axis.text=element_text(size=20))

ggplot(data = stats_df) + 
  geom_histogram(aes(x = mean_difs), fill = "#56B4E9") + 
  geom_vline(xintercept = 0.01, color = "red") + 
  xlab("Mean differences in energy binding strength") +
  theme_cowplot(12) +
  scale_y_continuous(expand = c(0,0)) +
  theme(text = element_text(size = 25), axis.text=element_text(size=20))
