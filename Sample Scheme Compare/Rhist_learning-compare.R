

# Load the ggplot2 library
library(ggplot2)


# hist var  hist_var_c1_s1_e1
#####
basevar_c1=output_tab1[,2]
varRatio_c1=rep(1,length(basevar_c1))
varRatio_s1=output_tab_s1[,2]/basevar_c1
varRatio_e1=output_tab_e1[,2]/basevar_c1

replicates <- Jvec

# Create a data frame in long format
my_data <- data.frame(
  replicate = rep(replicates, times = 3),
  dataset = rep(c("Case1", "Starting Group", "Ending Group"), each = length(replicates)),
  variance = c(varRatio_c1, varRatio_s1, varRatio_e1)
)

hist_var_c1_s1_e1=ggplot(my_data, aes(x = factor(replicate), y = variance, fill = dataset)) +
  geom_col(position = position_dodge(width = .4), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Replicate", y = "Variance", title = "Ratio of Variance Compared Replicate and Dataset") +
  coord_cartesian(ylim = c(0.5, 1.2)) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  ) +
  scale_fill_manual(
    values = c("Case1"="orange" ,"Starting Group" = "blue", "Ending Group" = "purple")  # Customize colors
  ) 

hist_var_c1_s1_e1


#####

# hist var  hist_var_c2_s2_e2
#####
basevar_c2=output_tab2[,2]
varRatio_c2=rep(1,length(basevar_c2))
varRatio_s2=output_tab_s2[,2]/basevar_c2
varRatio_e2=output_tab_e2[,2]/basevar_c2

replicates <- Jvec

# Create a data frame in long format
my_data <- data.frame(
  replicate = rep(replicates, times = 3),
  dataset = rep(c("Case2", "Starting Group", "Ending Group"), each = length(replicates)),
  variance = c(varRatio_c2, varRatio_s2, varRatio_e2)
)

hist_var_c2_s2_e2=ggplot(my_data, aes(x = factor(replicate), y = variance, fill = dataset)) +
  geom_col(position = position_dodge(width = .4), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Replicate", y = "Variance", title = "Ratio of Variance Compared Replicate and Dataset") +
  coord_cartesian(ylim = c(0.5, 1.2)) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  ) + 
  scale_fill_manual(
    values = c("Case2"="orange" ,"Starting Group" = "blue", "Ending Group" = "purple")  # Customize colors
  ) 

hist_var_c2_s2_e2
#####

# hist var  hist_var_c3_s3_e3
#####
basevar_c3=output_tab3[,2]
varRatio_c3=rep(1,length(basevar_c3))
varRatio_s3=output_tab_s3[,2]/basevar_c3
varRatio_e3=output_tab_e3[,2]/basevar_c3

replicates <- Jvec

# Create a data frame in long format
my_data <- data.frame(
  replicate = rep(replicates, times = 3),
  dataset = rep(c("Case3", "Starting Group", "Ending Group"), each = length(replicates)),
  variance = c(varRatio_c3, varRatio_s3, varRatio_e3)
)

hist_var_c3_s3_e3=ggplot(my_data, aes(x = factor(replicate), y = variance, fill = dataset)) +
  geom_col(position = position_dodge(width = .4), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Replicate", y = "Variance", title = "Ratio of Variance Compared Replicate and Dataset") +
  coord_cartesian(ylim = c(0.5, 1.2)) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  ) +
  scale_fill_manual(
    values = c("Case3"="orange" ,"Starting Group" = "blue", "Ending Group" = "purple")  # Customize colors
  ) 

hist_var_c3_s3_e3
#####

# hist var  hist_bvar_c1_s1_e1
#####
basevar_c1=output_tab1[,6]
varRatio_c1=rep(1,length(basevar_c1))
varRatio_s1=output_tab_s1[,6]/basevar_c1
varRatio_e1=output_tab_e1[,6]/basevar_c1

replicates <- Jvec

# Create a data frame in long format
my_data <- data.frame(
  replicate = rep(replicates, times = 3),
  dataset = rep(c("Case1", "Starting Group", "Ending Group"), each = length(replicates)),
  variance = c(varRatio_c1, varRatio_s1, varRatio_e1)
)

hist_bvar_c1_s1_e1=ggplot(my_data, aes(x = factor(replicate), y = variance, fill = dataset)) +
  geom_col(position = position_dodge(width = .4), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Replicate", y = "Variance", title = "Ratio of Bootstrap Variance Compared Replicate and Dataset") +
  coord_cartesian(ylim = c(0.5, 1.2)) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

hist_bvar_c1_s1_e1


#####

# hist var  hist_bvar_c2_s2_e2
#####
basevar_c2=output_tab2[,6]
varRatio_c2=rep(1,length(basevar_c2))
varRatio_s2=output_tab_s2[,6]/basevar_c2
varRatio_e2=output_tab_e2[,6]/basevar_c2

replicates <- Jvec

# Create a data frame in long format
my_data <- data.frame(
  replicate = rep(replicates, times = 3),
  dataset = rep(c("Case2", "Starting Group", "Ending Group"), each = length(replicates)),
  variance = c(varRatio_c2, varRatio_s2, varRatio_e2)
)

hist_bvar_c2_s2_e2=ggplot(my_data, aes(x = factor(replicate), y = variance, fill = dataset)) +
  geom_col(position = position_dodge(width = .4), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Replicate", y = "Variance", title = "Ratio of Bootstrap Variance Compared Replicate and Dataset") +
  coord_cartesian(ylim = c(0.5, 1.2)) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

hist_bvar_c2_s2_e2
#####

# hist var  hist_bvar_c3_s3_e3
#####
basevar_c3=output_tab3[,6]
varRatio_c3=rep(1,length(basevar_c3))
varRatio_s3=output_tab_s3[,6]/basevar_c3
varRatio_e3=output_tab_e3[,6]/basevar_c3

replicates <- Jvec

# Create a data frame in long format
my_data <- data.frame(
  replicate = rep(replicates, times = 3),
  dataset = rep(c("Case3", "Starting Group", "Ending Group"), each = length(replicates)),
  variance = c(varRatio_c3, varRatio_s3, varRatio_e3)
)

hist_bvar_c3_s3_e3=ggplot(my_data, aes(x = factor(replicate), y = variance, fill = dataset)) +
  geom_col(position = position_dodge(width = .4), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Replicate", y = "Variance", title = "Ratio of Bootstrap Variance Compared Replicate and Dataset") +
  coord_cartesian(ylim = c(0.5, 1.2)) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

hist_bvar_c3_s3_e3
#####


ggsave(filename = "hist_var_c1_s1_e1.pdf",
       plot = hist_var_c1_s1_e1,
       width = 7, height = 5)

ggsave(filename = "hist_var_c2_s2_e2.pdf",
       plot = hist_var_c2_s2_e2,
       width = 7, height = 5)

ggsave(filename = "hist_var_c3_s3_e3.pdf",
       plot = hist_var_c3_s3_e3,
       width = 7, height = 5)


ggsave(filename = "hist_bvar_c1_s1_e1.pdf",
       plot = hist_bvar_c1_s1_e1,
       width = 7, height = 5)

ggsave(filename = "hist_bvar_c2_s2_e2.pdf",
       plot = hist_bvar_c2_s2_e2,
       width = 7, height = 5)

ggsave(filename = "hist_bvar_c3_s3_e3.pdf",
       plot = hist_bvar_c3_s3_e3,
       width = 7, height = 5)


