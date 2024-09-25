#################################################################
##                  Loading packages and data                  ##
#################################################################
set.seed(4242)

# Loading required packages
library("parSim")
library("bootnet")
library("dplyr")
library("tidyr")
library("ggplot2")
library("psychonetrics")

# Loading data
simres <- read.csv("./data/simres.csv")

# Convert to long format:
longer <- simres %>% pivot_longer(logl:false_neg, names_to = "metric")
longer$n_replication_factor <- factor(longer$n_replication,levels=sort(unique(longer$n_replication)),
                                      labels = paste0("n = ",sort(unique(longer$n_replication))))

# To uppercase:
longer$metric <- toupper(longer$metric)

# Label type:
longer$type <- factor(longer$type, levels = "semi", labels = "semi confirmatory")

# Label reorder:
longer$p_reorder <- factor(longer$p_reorder, levels = sort(unique(longer$p_reorder)), labels = paste0("P(reorder) = ",sort(unique(longer$p_reorder))))

sub_TLI <- longer %>% filter(metric=="TLI")
sub_CFI <- longer %>% filter(metric=="CFI")
sub_RMSEA <- longer %>% filter(metric=="RMSEA")

#getting mean values to report in paper
mean_TLI <- sub_TLI %>%
  filter(((n_replication == 100 | n_replication == 200)) & p_reorder == "P(reorder) = 0") %>% 
  group_by(n_replication) %>% 
  summarise(mean = round(mean(value), digits = 2))
mean_TLI

mean_CFI <- sub_CFI %>%
  filter(((n_replication == 100 | n_replication == 200)) & p_reorder == "P(reorder) = 0") %>% 
  group_by(n_replication) %>% 
  summarise(mean = round(mean(value), digits = 2))
mean_CFI

mean_RMSEA <- sub_RMSEA %>%
  filter(((n_replication == 100 | n_replication == 200)) & p_reorder == "P(reorder) = 0") %>% 
  group_by(n_replication) %>% 
  summarise(mean = round(mean(value), digits = 2))
mean_RMSEA

#################################################################
##                            Plots                            ##
#################################################################
# TLI:
p1 <- ggplot(sub_TLI, aes(x = n_replication_factor, y = value)) + 
  facet_grid(cols = vars(p_reorder)) + geom_boxplot() + theme_bw() + 
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 1.0), minor_breaks =  seq(-10,10,by=0.05)) + 
  geom_hline(yintercept=1) + 
  geom_hline(yintercept=0.97, lty = 2) + geom_hline(yintercept=0.95, lty = 3) + 
  ylab("") + xlab("Replication set sample size") + 
  ggtitle("Tucker-Lewis Index (TLI)","Typical SEM interpretation: > 0.97 = good fit and 0.95 to 0.97 = acceptable fit") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave("./plots/tli.png", width = 6.66, height = 4.63, dpi = 300, scale = 1.2)

# CFI:
p2 <- ggplot(sub_CFI, aes(x = n_replication_factor, y = value)) + 
  facet_grid(cols = vars(p_reorder)) + geom_boxplot() + theme_bw() + 
  scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 0.95, 0.97, 1.0), minor_breaks =  seq(-10,10,by=0.05)) + 
  geom_hline(yintercept=1) + 
  geom_hline(yintercept=0.97, lty = 2) + geom_hline(yintercept=0.95, lty = 3) + 
  ylab("") + xlab("Replication set sample size") + 
  ggtitle("Bentler's Comparative Fit Index (CFI)","Typical SEM interpretation: > 0.97 = good fit and 0.95 to 0.97 = acceptable fit") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave("./plots/cfi.png", width = 6.66, height = 4.63, dpi = 300, scale = 1.2)

# RMSEA:
p3 <- ggplot(sub_RMSEA, aes(x = n_replication_factor, y = value)) + 
  facet_grid(cols = vars(p_reorder)) + geom_boxplot() + theme_bw() + 
  scale_y_continuous(breaks = c(0,0.05,0.08,0.10), minor_breaks =  seq(-10,10,by=0.05)) + 
  geom_hline(yintercept=0) + 
  geom_hline(yintercept=0.05, lty = 2) + 
  geom_hline(yintercept=0.08, lty = 2) + 
  geom_hline(yintercept=0.10, lty = 2) +
  ylab("") + xlab("Replication set sample size") + 
  ggtitle("Root Mean Square Error of Approximation (RMSEA)","Typical SEM interpretation: 
  ≤ 0.05 = good fit, 0.05 to 0.08 = adequate fit, 0.08 to 0.10 = mediocre fit, and > 0.10 = unacceptable")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave("./plots/rmsea.png", width = 6.66, height = 4.63, dpi = 300, scale = 1.2)

# # Print to PDF:
# pdf("./plots/plots.pdf",width=1920/1080 * 7,height=7)
# print(p1)
# print(p2)
# print(p3)
# 
# dev.off()