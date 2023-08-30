# Reproduct Supplementary Figures 3a & 3b
library(lemon)
library(ggplot2)

#---------- data -----------
df <- read.csv(paste(PATH, "results/TableS2plot.csv", sep=''), header=T)

#------ Censored Plot ------
df_cen <- df[which(df$Censor.or.Missing == "Censored" | df$Censor.or.Missing == "Non-censored"),]
df_cen$Censor.or.Missing <- factor(df_cen$Censor.or.Missing, levels = c("Non-censored", "Censored"))

windowsFonts(A = windowsFont("Arial"))
label_both <- c("IDV-RTV=0"="\n IDV-RTV=0\n", "IDV-RTV=1"="\n IDV-RTV=1\n")

#--- Censored Plot Output ---
postscript(paste(PATH, "results/fig3a.eps", sep=""), width = 18, height = 12)
ggplot(df_cen, aes(x = factor(Days), y = Percent, fill = Censor.or.Missing)) + 
  geom_bar(stat = "identity") + 
  geom_text(mapping = aes(label = paste0(round(Percent,1),"%")),
            size = 6, colour = 'black', hjust = .5, position = position_stack(vjust = 0.5)) +
  xlab("Days") + ylab("Censored Percentage") + 
  guides(fill = guide_legend(title = "", reverse = TRUE)) + 
  theme(axis.text = element_text(size = 17, family = "A"), axis.title = element_text(size = 20, face = "bold", family = "A"),
        legend.title = element_text(size = 18, family = "A"), legend.text = element_text(size = 17, family = "A"),
        legend.position = "top", strip.text = element_text(size = 17, face = "bold", family = "A")) +
  scale_fill_manual(values = c("#BAE2F3", "#157CE2")) + 
  facet_rep_grid(IDV ~ ., labeller = as_labeller(label_both), repeat.tick.labels = TRUE) 
dev.off()

#------ Missing Plot ------
df_mis <- df[which(df$Censor.or.Missing == "Missing" | df$Censor.or.Missing == "Observed"),]
df_mis$Censor.or.Missing <- factor(df_mis$Censor.or.Missing, levels = c("Observed", "Missing"))

windowsFonts(A = windowsFont("Arial"))
label_both <- c("IDV-RTV=0"="\n IDV-RTV=0\n", "IDV-RTV=1"="\n IDV-RTV=1\n")

#--- Missing Plot Output ---
postscript(paste(PATH, "results/fig3b.eps", sep=""), width = 18, height = 12)
ggplot(df_mis, aes(x = factor(Days), y = Percent, fill = Censor.or.Missing)) + 
  geom_bar(stat = "identity") + 
  geom_text(mapping = aes(label = paste0(round(Percent,1),"%")),
            size = 6, colour = 'black', hjust = .5, position = position_stack(vjust = 0.5)) +
  xlab("Days") + ylab("Missing Percentage") + 
  guides(fill = guide_legend(title = "", reverse = TRUE)) + 
  theme(axis.text = element_text(size = 17, family = "A"), axis.title = element_text(size = 20, face = "bold", family = "A"),
        legend.title = element_text(size = 18, family = "A"), legend.text = element_text(size = 17, family = "A"),
        legend.position = "top", strip.text = element_text(size = 17, face = "bold", family = "A")) +
  scale_fill_manual(values = c("#F3E3EC","#FF90B1")) + 
  facet_rep_grid(IDV ~., labeller = as_labeller(label_both), repeat.tick.labels = TRUE) 
dev.off()
