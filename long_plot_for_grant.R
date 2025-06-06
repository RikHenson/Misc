
# Plotting data and two-line fits to example grant data, rik.henson@mrc-cbu.cam.ac.uk

library(ggplot2)
library(lme4)
library(lmerTest)

plot <- function(data, x_var, y_var, x_label, y_label, ymin, ymax) {
  ggp = ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(size = 2, stroke = 1, alpha = 0.3) +
    geom_line(aes(group = CCID), lty = 2, colour = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5, lty = 1, aes(group = CCID)) +
    labs(x = x_label, y = y_label) +
    ylim(ymin, ymax) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) 
  return(ggp)
}

### COG

cog_df = read.csv('long_cog_data.csv')
ylabel = "Fluid Intelligence"

(ggp_cog = plot(cog_df, "Age", "Cattell", "Age", ylabel, min(cog_df$Cattell), max(cog_df$Cattell))) 

fit = lmer(Cattell ~ poly(Age,2) + (1|CCID), data=cog_df)
summary(fit)

params = as.data.frame(effects::effect(term = "poly(Age, 2)", mod = fit))

ggp_cog + geom_line(data=params, aes(x=Age, y=fit), linewidth=2, color="blue") 

ggsave(paste("long_plot_","Cog",".png",sep=""), width=2500, height=1500, units="px", dpi=300)


## TGM

tgm_df = read.csv('long_tgm_data.csv')
ylabel = "TotalGM"
tgm_df$TGM = tgm_df$TGM/1e6 # convert to litres

(ggp_tgm = plot(tgm_df, "Age", "TGM", "Age", ylabel, min(tgm_df$TGM), max(tgm_df$TGM))) 

fit = lmer(TGM ~ poly(Age,2) + (1|CCID), data=tgm_df)
summary(fit)

params = as.data.frame(effects::effect(term = "poly(Age, 2)", mod = fit))

ggp_tgm + geom_line(data=params, aes(x=Age, y=fit), linewidth=2, color="blue") 

ggsave(paste("long_plot_","TGM",".png",sep=""), width=2500, height=1500, units="px", dpi=300)



