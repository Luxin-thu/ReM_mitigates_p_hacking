# Clear the workspace and load necessary libraries
rm(list=ls())

library(tidyverse)    # For data manipulation and visualization
library(xtable)       # For creating LaTeX tables
library(latex2exp)    # For LaTeX expression parsing
library(scales)       # For scaling functions in plots

################################### Figure 1 #############################################

# Function to calculate the required acceptance probability based on R2, alpha, and K
c_value_q <- function(R2,alpha,K){
  c_value <- qnorm(1-alpha/2)^2*(R2/K)/(sqrt(1-R2)+sqrt(1-R2+R2/K))^2
  return(pchisq(c_value,df=K))
}

# Set up parameter combinations to compute the required acceptance probabilities
setups <- expand_grid(
  K  = 1:6, # Number of covariates (from 1 to 6)
  R2 = seq(0.1,1,by = 0.01)[-1], # R-squared values from 0.11 to 1
  alpha = c(0.05,0.1) # Significance levels
)

# Calculate required acceptance probabilities for each setup
res_table <- map_dfr(1:nrow(setups),~{
  K <- setups$K[.x]; R2 <- setups$R2[.x]; alpha <- setups$alpha[.x]
  tibble(a = (c_value_q(R2,alpha,K)))
}) %>% bind_cols(setups)


# Define a custom log10 transformation for plotting 
trans_log10_1p <- function() {
  trans_new("log10_1p", function(x) log10(1 + x), function(x) 10^x - 1)
}



# Function to create a LaTeX formatted string for alpha values
appender_alpha <- function(string){
  TeX(paste("$\\alpha = $", string))  
}

# Plot Figure 1
res_table %>%mutate(K = as.factor(K)) %>% ggplot(aes(y=a, x = R2, linetype = K)) + 
  facet_wrap(~alpha,  labeller=  labeller(alpha = as_labeller(appender_alpha, default = label_parsed))) + 
  geom_line() + ylab('Required acceptance probability of ReM') + xlab(TeX("$\\R^2_x$")) + theme(text = element_text(size=25)) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^seq(-7, 0, by = 1)),  # 设置对数刻度
    labels = trans_format("log10", math_format(10^.x))  # 使用科学计数法格式化标签
  )+scale_x_continuous(
    breaks = seq(0.1, 1, by = 0.1)  # 更密集的 x 轴刻度
  )

# Save the plot as a PDF file
ggsave(filename = paste0('Prob_R2.pdf'),units='in',height = 8.5, width = 14)


################################### Figure 2 #############################################

appender_rho <- function(string){
  TeX(paste("$\\rho = $", string))  
}


# Load results for equal coefficients
load("beta_equal_res_all.Rdata")
df_equal <- df_0.05
df_equal$Beta <- 'Equal'

# Load results for unequal coefficients
load("beta_unequal_res_all.Rdata")
df_unequal <- df_0.05
df_unequal$Beta <- 'Unequal'

# Combine both data frames into one for plotting
df_plot <- rbind(df_equal, df_unequal)

# Dataframe of the error bounds
df_e <- data.frame(Bound = e_mat %>% as.vector(), Design = c('ReM', 'ReM', 'ReP', 'ReP'), rho = c(0, 0.8, 0, 0.8))

# plot Figure 2
df_plot %>% mutate(R2 = as.numeric(R2)) %>% mutate(Design = fct_relevel(Design, 'ReM', 'ReP', 'CRE')) %>% ggplot(aes(x = R2, y = ER, color = Design)) + geom_line()  +
  geom_hline(yintercept = 0.05, linetype = 'dashed') + geom_hline(data = df_e, aes(yintercept = Bound, color = Design), linetype = "dashed") + theme(text = element_text(size=15))+
  ylab('Type I error rate') + xlab(TeX("$\\R^2$")) + facet_grid(Beta~rho, labeller =labeller(rho = as_labeller(appender_rho, default = label_parsed)) ) + ylim(c(0,NA))
ggsave(filename = paste0('type_I_error_R2_plot_cre','.pdf'),units='in',width = 8, height = 4.5)

################################### table 1 #############################################


# print table 1
df_plot %>% filter(Beta=='Equal') %>% mutate(ER = paste(sprintf("%.2f", ER*100),' (',sprintf("%.2f", ci*100),')', sep = ''))  %>% dplyr::select(c('ER','R2','rho','Design')) %>% pivot_wider(values_from = ER, names_from = c(R2)) %>% arrange(rho) %>% xtable()

df_plot %>% filter(Beta=='Unequal') %>% mutate(ER = paste(sprintf("%.2f", ER*100),' (',sprintf("%.2f", ci*100),')', sep = ''))  %>% dplyr::select(c('ER','R2','rho','Design')) %>% pivot_wider(values_from = ER, names_from = c(R2)) %>% arrange(rho) %>% xtable()

