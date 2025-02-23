# Load the ggplot2 library
library(ggplot2)


############# joint compare variance s1e1
#####

# Create two separate datasets
var_s1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s1[,2],
  Group = "Starting Group"
)

var_e1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e1[,2],
  Group = "Ending Group"
)

# Combine the datasets
var_s1_e1 <- rbind(var_s1, var_e1)

# Create the plot
plot_var_s1_e1=ggplot(var_s1_e1, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance for Learning Method for Case (i)",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_s1_e1

#####


############# joint compare CR s1e1
#####

# Create two separate datasets
t_s1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s1[,4],
  Group = "Starting Group"
)

t_e1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e1[,4],
  Group = "Ending Group"
)

# Combine the datasets
t_s1_e1 <- rbind(t_s1, t_e1)

# Create the plot
plot_t_s1_e1=ggplot(t_s1_e1 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("Starting Group" = "blue", "Ending Group" = "red")  # Customize colors
  ) +
  labs(
    title = "CR of CI based on t dist. for Learning Method for Case (i)",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_t_s1_e1

#####


############# joint compare variance s2e2
#####

# Create two separate datasets
var_s2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s2[,2],
  Group = "Starting Group"
)

var_e2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e2[,2],
  Group = "Ending Group"
)

# Combine the datasets
var_s2_e2 <- rbind(var_s2, var_e2)

# Create the plot
plot_var_s2_e2=ggplot(var_s2_e2, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance for Learning Method for Case (iii)",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_s2_e2

#####


############# joint compare CR s2e2
#####

# Create two separate datasets
t_s2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s2[,4],
  Group = "Starting Group"
)

t_e2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e2[,4],
  Group = "Ending Group"
)

# Combine the datasets
t_s2_e2 <- rbind(t_s2, t_e2)

# Create the plot
plot_t_s2_e2=ggplot(t_s2_e2 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("Starting Group" = "blue", "Ending Group" = "red")  # Customize colors
  ) +
  labs(
    title = "CR of CI based on t dist. for Learning Method for Case (iii)",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_t_s2_e2

#####


############# joint compare variance s3e3
#####

# Create two separate datasets
var_s3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s3[,2],
  Group = "Starting Group"
)

var_e3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e3[,2],
  Group = "Ending Group"
)

# Combine the datasets
var_s3_e3 <- rbind(var_s3, var_e3)

# Create the plot
plot_var_s3_e3=ggplot(var_s3_e3, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance for Learning Method for Case (ii)",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_s3_e3

#####


############# joint compare CR s3e3
#####

# Create two separate datasets
t_s3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s3[,4],
  Group = "Starting Group"
)

t_e3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e3[,4],
  Group = "Ending Group"
)

# Combine the datasets
t_s3_e3 <- rbind(t_s3, t_e3)

# Create the plot
plot_t_s3_e3=ggplot(t_s3_e3 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("Starting Group" = "blue", "Ending Group" = "red")  # Customize colors
  ) +
  labs(
    title = "CR of CI based on t dist. for Learning Method for Case (ii)",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_t_s3_e3

#####


############# save pdf
#####
ggsave(filename = "plot_var_s1_e1.pdf",
       plot = plot_var_s1_e1,
       width = 7, height = 5)
ggsave(filename = "plot_t_s1_e1.pdf",
       plot = plot_t_s1_e1,
       width = 7, height = 5)
ggsave(filename = "plot_var_s2_e2.pdf",
       plot = plot_var_s2_e2,
       width = 7, height = 5)
ggsave(filename = "plot_t_s2_e2.pdf",
       plot = plot_t_s2_e2,
       width = 7, height = 5)
ggsave(filename = "plot_var_s3_e3.pdf",
       plot = plot_var_s3_e3,
       width = 7, height = 5)
ggsave(filename = "plot_t_s3_e3.pdf",
       plot = plot_t_s3_e3,
       width = 7, height = 5)
#####



