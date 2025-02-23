# Load the ggplot2 library
library(ggplot2)


############# T21 Var
#####

# Create two separate datasets
var_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,2],
  Group = "case1"
)

var_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,2],
  Group = "case2"
)

var_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,2],
  Group = "case3"
)

# Combine the datasets
var_c <- rbind(var_c1, var_c2,var_c3)

# Create the plot
plot_var_c=ggplot(var_c, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Variance for different cases",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_var_c

#####


############# T21 phi CR
#####

# Create two separate datasets
phi_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,3],
  Group = "case1"
)

phi_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,3],
  Group = "case2"
)

phi_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,3],
  Group = "case3"
)

# Combine the datasets
phi_c <- rbind(phi_c1, phi_c2,phi_c3)

# Create the plot
plot_phi_c=ggplot(phi_c, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Phi CR for different cases",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_phi_c

#####


############# T21 boots var
#####

# Create two separate datasets
bvar_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,6],
  Group = "case1"
)

bvar_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,6],
  Group = "case2"
)

bvar_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,6],
  Group = "case3"
)

# Combine the datasets
bvar_c <- rbind(bvar_c1, bvar_c2,bvar_c3)

# Create the plot
plot_bvar_c=ggplot(bvar_c, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "bootstrap Variance for different cases",
    x = "Replicates",
    y = "Variance",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bvar_c

#####


############# T21 boots CR
#####

# Create two separate datasets
bCR_c1 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab1[,7],
  Group = "case1"
)

bCR_c3 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab3[,7],
  Group = "case2"
)

bCR_c2 <- data.frame(
  Replicates = Jvec,
  Variance = output_tab2[,7],
  Group = "case3"
)

# Combine the datasets
bCR_c <- rbind(bCR_c1, bCR_c2,bCR_c3)

# Create the plot
plot_bCR_c=ggplot(bCR_c, aes(x = Replicates, y = Variance, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "boots CR for different cases",
    x = "Replicates",
    y = "CR",
    color = "Group"
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom"  # Move the legend to below the plot
  )

plot_bCR_c

#####


############# save pdf
#####
ggsave(filename = "plot_var_c.pdf",
       plot = plot_var_c,
       width = 7, height = 5)
ggsave(filename = "plot_phi_c.pdf",
       plot = plot_phi_c,
       width = 7, height = 5)
ggsave(filename = "plot_bvar_c.pdf",
       plot = plot_bvar_c,
       width = 7, height = 5)
ggsave(filename = "plot_bCR_c.pdf",
       plot = plot_bCR_c,
       width = 7, height = 5)
#####





