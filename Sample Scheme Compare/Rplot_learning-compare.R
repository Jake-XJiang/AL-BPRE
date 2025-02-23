# Load the ggplot2 library
library(ggplot2)


############# Learning-compare case 1 var
#####

# Create two separate datasets
var_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,2],
  Group = "case1"
)

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
var_c1_s1_e1 <- rbind(var_c1, var_s1, var_e1)

# Create the plot
plot_var_c1_s1_e1=ggplot(var_c1_s1_e1, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance comparison",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_c1_s1_e1

#####

############# Learning-compare case 1 t CR
#####
# Create two separate datasets
t_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,4],
  Group = "case1"
)

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
t_c1_s1_e1  <- rbind(t_c1, t_s1,t_e1)

# Create the plot
plot_t_c1_s1_e1 =ggplot(t_c1_s1_e1 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "t CR comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_t_c1_s1_e1

#####

############# Learning-compare case 1 boots var
#####
# Create two separate datasets
bvar_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,6],
  Group = "case1"
)

bvar_s1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s1[,6],
  Group = "Starting Group"
)

bvar_e1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e1[,6],
  Group = "Ending Group"
)

# Combine the datasets
bvar_c1_s1_e1  <- rbind(bvar_c1, bvar_s1,bvar_e1)

# Create the plot
plot_bvar_c1_s1_e1 =ggplot(bvar_c1_s1_e1 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots Var comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bvar_c1_s1_e1

#####

############# Learning-compare case 1 boots CR
#####
# Create two separate datasets
bCR_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,7],
  Group = "case1"
)

bCR_s1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s1[,7],
  Group = "Starting Group"
)

bCR_e1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e1[,7],
  Group = "Ending Group"
)

# Combine the datasets
bCR_c1_s1_e1  <- rbind(bCR_c1, bCR_s1,bCR_e1)

# Create the plot
plot_bCR_c1_s1_e1 =ggplot(bCR_c1_s1_e1 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots CR comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bCR_c1_s1_e1

#####


############# Learning-compare case 2 var
#####

# Create two separate datasets
var_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,2],
  Group = "case2"
)

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
var_c2_s2_e2 <- rbind(var_c2, var_s2, var_e2)

# Create the plot
plot_var_c2_s2_e2=ggplot(var_c2_s2_e2, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance comparison",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_c2_s2_e2

#####

############# Learning-compare case 2 t CR
#####
# Create two separate datasets
t_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,4],
  Group = "case2"
)

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
t_c2_s2_e2  <- rbind(t_c2, t_s2,t_e2)

# Create the plot
plot_t_c2_s2_e2 =ggplot(t_c2_s2_e2 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "t CR comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_t_c2_s2_e2

#####

############# Learning-compare case 2 boots var
#####
# Create two separate datasets
bvar_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,6],
  Group = "case2"
)

bvar_s2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s2[,6],
  Group = "Starting Group"
)

bvar_e2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e2[,6],
  Group = "Ending Group"
)

# Combine the datasets
bvar_c2_s2_e2  <- rbind(bvar_c2, bvar_s2,bvar_e2)

# Create the plot
plot_bvar_c2_s2_e2 =ggplot(bvar_c2_s2_e2 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots Var comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bvar_c2_s2_e2

#####

############# Learning-compare case 2 boots CR
#####
# Create two separate datasets
bCR_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,7],
  Group = "case2"
)

bCR_s2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s2[,7],
  Group = "Starting Group"
)

bCR_e2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e2[,7],
  Group = "Ending Group"
)

# Combine the datasets
bCR_c2_s2_e2  <- rbind(bCR_c2, bCR_s2,bCR_e2)

# Create the plot
plot_bCR_c2_s2_e2 =ggplot(bCR_c2_s2_e2 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots CR comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bCR_c2_s2_e2

#####


############# Learning-compare case 3 var
#####

# Create two separate datasets
var_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,2],
  Group = "case3"
)

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
var_c3_s3_e3 <- rbind(var_c3, var_s3, var_e3)

# Create the plot
plot_var_c3_s3_e3=ggplot(var_c3_s3_e3, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance comparison",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_c3_s3_e3

#####

############# Learning-compare case 3 t CR
#####
# Create two separate datasets
t_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,4],
  Group = "case3"
)

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
t_c3_s3_e3  <- rbind(t_c3, t_s3,t_e3)

# Create the plot
plot_t_c3_s3_e3 =ggplot(t_c3_s3_e3 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "t CR comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_t_c3_s3_e3

#####

############# Learning-compare case 3 boots var
#####
# Create two separate datasets
bvar_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,6],
  Group = "case3"
)

bvar_s3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s3[,6],
  Group = "Starting Group"
)

bvar_e3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e3[,6],
  Group = "Ending Group"
)

# Combine the datasets
bvar_c3_s3_e3  <- rbind(bvar_c3, bvar_s3,bvar_e3)

# Create the plot
plot_bvar_c3_s3_e3 =ggplot(bvar_c3_s3_e3 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots Var comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bvar_c3_s3_e3

#####

############# Learning-compare case 3 boots CR
#####
# Create two separate datasets
bCR_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,7],
  Group = "case3"
)

bCR_s3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_s3[,7],
  Group = "Starting Group"
)

bCR_e3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab_e3[,7],
  Group = "Ending Group"
)

# Combine the datasets
bCR_c3_s3_e3  <- rbind(bCR_c3, bCR_s3,bCR_e3)

# Create the plot
plot_bCR_c3_s3_e3 =ggplot(bCR_c3_s3_e3 , aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots CR comparison",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bCR_c3_s3_e3

#####



############# save pdf
#####
ggsave(filename = "plot_var_c1_s1_e1.pdf",
       plot = plot_var_c1_s1_e1,
       width = 7, height = 5)
ggsave(filename = "plot_t_c1_s1_e1.pdf",
       plot = plot_t_c1_s1_e1,
       width = 7, height = 5)
ggsave(filename = "plot_bvar_c1_s1_e1.pdf",
       plot = plot_bvar_c1_s1_e1,
       width = 7, height = 5)
ggsave(filename = "plot_bCR_c1_s1_e1.pdf",
       plot = plot_bCR_c1_s1_e1,
       width = 7, height = 5)

ggsave(filename = "plot_var_c2_s2_e2.pdf",
       plot = plot_var_c2_s2_e2,
       width = 7, height = 5)
ggsave(filename = "plot_t_c2_s2_e2.pdf",
       plot = plot_t_c2_s2_e2,
       width = 7, height = 5)
ggsave(filename = "plot_bvar_c2_s2_e2.pdf",
       plot = plot_bvar_c2_s2_e2,
       width = 7, height = 5)
ggsave(filename = "plot_bCR_c2_s2_e2.pdf",
       plot = plot_bCR_c2_s2_e2,
       width = 7, height = 5)

ggsave(filename = "plot_var_c3_s3_e3.pdf",
       plot = plot_var_c3_s3_e3,
       width = 7, height = 5)
ggsave(filename = "plot_t_c3_s3_e3.pdf",
       plot = plot_t_c3_s3_e3,
       width = 7, height = 5)
ggsave(filename = "plot_bvar_c3_s3_e3.pdf",
       plot = plot_bvar_c3_s3_e3,
       width = 7, height = 5)
ggsave(filename = "plot_bCR_c3_s3_e3.pdf",
       plot = plot_bCR_c3_s3_e3,
       width = 7, height = 5)
#####





