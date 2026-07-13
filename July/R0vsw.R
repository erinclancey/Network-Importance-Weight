
# True Parameter Values
w=seq(0,1, length.out = 1000)
p=seq(0,1, length.out = 1000)
beta_par=0.000025
gamma=1/10
S0=5000

R0_random <- S0 * 5 * beta_par / gamma

R0_total <- vector()

#######Calculate R0
for (i in 1:length(w)){
  A_mat <- matrix(rbinom(25,1,p[i]), nrow=5)
  diag(A_mat) <- 1
  A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
  A_vec <- as.vector(A_mat)
  # Build Cmat
  Cmat <- A_mat * w[i] + (1 - w[i])   # off-diagonal formula applied everywhere
  diag(Cmat) <- 1 
  # Compute eigen decomposition
  eig <- eigen(Cmat)
  # Eigenvalues
  dom_eig <- max(eig$values)
  R0_total[i] <- S0*beta_par / gamma * dom_eig
}

plot(w, R0_total)
plot(p, R0_total)
R0_random
min(R0_total)


# Parameter Values
N_grid <- 100 # Reduced resolution for speed
w_vals <- seq(0, 1, length.out = N_grid)
p_vals <- seq(0, 1, length.out = N_grid)
beta_par <- 0.000025
gamma <- 1/10
S0 <- 5000

# Initialize matrix to store R0 results
R0_grid <- matrix(0, nrow = N_grid, ncol = N_grid)

# Calculate R0 over the grid of w and p
for (i in 1:N_grid) {
  for (j in 1:N_grid) {
    w <- w_vals[i]
    p <- p_vals[j]
    
    # Build Adjacency matrix (5x5)
    A_mat <- matrix(rbinom(25, 1, p), nrow=5)
    diag(A_mat) <- 1
    A_mat[lower.tri(A_mat)] <- t(A_mat)[lower.tri(A_mat)]
    
    # Build Cmat
    Cmat <- A_mat * w + (1 - w)
    diag(Cmat) <- 1
    
    # Compute Eigenvalues
    eig <- eigen(Cmat)
    dom_eig <- max(Re(eig$values)) # Using Real part to ensure numeric
    
    R0_grid[i, j] <- S0 * beta_par / gamma * dom_eig
  }
}

# 2D Plot: Filled Contour with Increased Font Sizes
filled.contour(w_vals, p_vals, R0_grid, 
               color.palette = terrain.colors,
               
               # Title settings (cex.main controls size)
               plot.title = {
                 title(main = "Total Population Reproductive Number", 
                       cex.main = 1.8) # 1.8x default size
               },
               
               # Axis and Label settings (cex.lab and cex.axis control size)
               plot.axes = {
                 axis(1, cex.axis = 1.4) # X-axis tick font size
                 axis(2, cex.axis = 1.4) # Y-axis tick font size
                 title(xlab = "Importance Weight (w)", 
                       ylab = "Probability of Network Connection (p)", 
                       cex.lab = 1.5) # Axis titles font size
               },
               
               # Legend tick font size
               key.axes = axis(4, cex.axis = 1.2) 
)

# Define explicit breaking points from 1 to 6
legend_levels <- seq(1.25, 6.25, by = 0.2) # Breaks every 0.5; use seq(1, 6, by = 1) for only whole numbers

# 2D Plot: Filled Contour with Custom Legend Ticks
filled.contour(w_vals, p_vals, R0_grid, 
               color.palette = terrain.colors,
               levels = legend_levels, # Forces the color scale to map between 1 and 6
               
               # Title settings
               plot.title = {
                 title(main = "Total Population Reproductive Number", 
                       cex.main = 1.8) 
               },
               
               # Axis and Label settings
               plot.axes = {
                 axis(1, cex.axis = 1.4) 
                 axis(2, cex.axis = 1.4) 
                 title(xlab = "Importance Weight (w)", 
                       ylab = "Probability of Network Connection (p)", 
                       cex.lab = 1.5) 
               },
               
               # Legend settings: Force ticks to happen at our specified levels
               key.axes = {
                 axis(4, at = seq(1.25, 6.25, by = 1), labels = seq(1.25, 6.25, by = 1), cex.axis = 1.2)
               }
)

# 1. Define a custom palette replacing bright green with steel/teal blue
# Low R0 = Steel Blue -> Muted Teal -> High R0 = Soft Sand Yellow
teal_steel_palette <- colorRampPalette(c("#2B4C7E", "#4A90E2", "#A5CAD2", "#E2F1F4"))
c("#4A7C59", "#688B9A", "#A0B9BF", "#EAD7A1")

teal_orange_palette <- colorRampPalette(c("black","#2B4C7E", "#4A90E2", "#A5CAD2", "#E6AF2E", "#F4EAD4","#A6A6A6"))

subpop_palette <- colorRampPalette(c(
  "#5B9BD5", # 1: Steel Blue
  "#E6AF2E", # 2: Warm Yellow/Orange
  "#5CB89C", # 3: Muted Teal
  "#D597BD", # 4: Muted Pink
  "#A6A6A6"  # 5: Grey
))


# Define your custom breaking points
legend_levels <- seq(1.25, 6.25, by = 0.2) 

# 2D Plot: Filled Contour with Custom Palette and Legend Ticks
filled.contour(w_vals, p_vals, R0_grid, 
               color.palette = teal_orange_palette, # Swapped to your custom muted blues
               levels = legend_levels, 
               
               # Title settings
               plot.title = {
                 title(main = "Total Population Reproductive Number", 
                       cex.main = 1.8) 
               },
               
               # Axis and Label settings
               plot.axes = {
                 axis(1, cex.axis = 1.4) 
                 axis(2, cex.axis = 1.4) 
                 title(xlab = "Importance Weight (w)", 
                       ylab = "Expected Connectance", 
                       cex.lab = 1.5) 
               },
               
               # Legend settings: Ticks placed exactly at your intervals
               key.axes = {
                 axis(4, at = seq(1.25, 6.25, by = 1), 
                      labels = seq(1.25, 6.25, by = 1), 
                      cex.axis = 1.2)
               }
)

