

# IBM code running

movement <- function(inds, xloc = 2, yloc = 3, xmax = 8, ymax = 8){
  total_inds   <- dim(inds)[1]; # Get the number of individuals in inds
  move_dists   <- c(-1, 0, 1);  # Define the possible distances to move
  x_move       <- sample(x = move_dists, size = total_inds, replace = TRUE);
  y_move       <- sample(x = move_dists, size = total_inds, replace = TRUE);
  inds[, xloc] <- inds[, xloc] + x_move;
  inds[, yloc] <- inds[, yloc] + y_move;
  # =========   The reflecting boundary is added below
  for(i in 1:total_inds){               # For each individual i in the array
    if(inds[i, xloc] > xmax){         # If it moved passed the maximum xloc
      inds[i, xloc] <- xmax - 1;    # Then move it back toward the centre
    }
    if(inds[i, xloc] < 1){            # If it moved below 1 on xloc
      inds[i, xloc] <- 2;           # Move it toward the centre (2)
    }
    if(inds[i, yloc] > ymax){         # If it moved passed the maximum yloc
      inds[i, yloc] <- ymax - 1;    # Then move it back toward the centre
    }
    if(inds[i, yloc] < 1){            # If it moved below 1 on yloc
      inds[i, yloc] <- 2;           # Then move it toward the centre (2)
    }
  } 
  # =========  Now all individuals should stay on the landscape
  return(inds);
}

inds           <- array(data = 0, dim = c(1, 3))
colnames(inds) <- c("body_mass", "x_loc", "y_loc")
rownames(inds) <- c("ind_1")
inds[,1]       <- rnorm(n = dim(inds)[1], mean = 23, sd = 3)
inds[,2]       <- sample(x = 1:8, size = dim(inds)[1], replace = TRUE)
inds[,3]       <- sample(x = 1:8, size = dim(inds)[1], replace = TRUE)

plot(
  x = inds[, 2],               # x-coordinate data
  y = inds[, 3],              # y-coordinate data
  pch = 16,                   # Point character (solid circle)
  cex = 2,                    # Point size (enlarged by a factor of 2)
  col = rgb(1, 0.5, 0),
  xlim = c(1, 8),             # X-axis limits (range from 1 to 8)
  ylim = c(1, 8),             # Y-axis limits (range from 1 to 8)
  xlab = "x location",        # X-axis label
  mar = c(5, 5, 1, 1),        # Margins for the plot (bottom, left, top, right)
  ylab = "y location",        # Y-axis label
  cex.lab = 1.5,              # Size of axis labels (enlarged by a factor of 1.5)
  cex.axis = 1.5              # Size of axis values (enlarged by a factor of 1.5)
)

time_steps <- 20

for (i in 1:time_steps){ 
  inds <- movement(inds)
  points(x = inds[, 2],              # x-coordinate data
         y = inds[, 3],              # y-coordinate data
         pch = 16,                   # Point character (solid circle)
         cex = 2,                    # Point size (enlarged by a factor of 2)
         col = rgb(1, 0.5, i/20, alpha = 0.7))
}


 





