############################################
#	Original script from Clark and Lauderdale's 2012 article "The Genealogy of Law"
# https://doi.org/10.1093/pan/mps019
#
#	Code tidied and modified by Matthew Dahl
############################################

## Libs
library(diagram)

#############################################
#	Plotting script
#############################################

plotXi <- function (Xi, file = 'TreePlot.pdf', case_names, case_years, title) {
  
  # Function that calculates the horizonal positioning for each case
  get_horizontal_position <- function (position_matrix, parent_id, child_ids, parent_min_position, parent_max_position) {
    
    # Calculate the appropriate horizontal position for the case given its number of children
    min_position <- (cumsum(offspring_count[child_ids]) - offspring_count[child_ids]) / sum(offspring_count[child_ids]) * (parent_max_position - parent_min_position) + parent_min_position
    max_position <- (cumsum(offspring_count[child_ids])) / sum(offspring_count[child_ids]) * (parent_max_position - parent_min_position) + parent_min_position
    center_position <- (min_position + max_position) / 2
    
    # Set the horizontal position for the case
    position_matrix[child_ids, 2] <- center_position
    
    # Recursively calculate the horizonal position for each child of the case
    for (i in 1:length(child_ids)) {
      next_generation_ids <- which((t(Xi) %*% diag(1,T,T))[child_ids[i],] == 1)
      
      if (length(next_generation_ids) > 0) {
        position_matrix <- get_horizontal_position(
          position_matrix = position_matrix,
          parent_id = i,
          child_ids = next_generation_ids,
          parent_min_position = min_position[i],
          parent_max_position = max_position[i]
        )
      }
    }
    
    return(position_matrix)
  }
  
  # Calculate the total number of offspring per case
  # Offspring value is cumulative, i.e., counts all children of all children etc.
  # Creates a vector where the index is the case ID and the value is the number of that case's offspring
  T <- dim(Xi)[1]
  position_matrix <- matrix(0, T, 2)
  position_matrix[,1] <- case_years
  for (f in Founders) {
    current_generation <- t(Xi) %*% diag(f,T,T)
    generation <- rep(0,T)
    offspring_count <- rep(0,T)
    generation_count <- 0
    
    while (sum(current_generation > 0)) {
      generation_count <- generation_count + 1
      generation[which(current_generation[1,] == 1)] <- generation_count
      offspring_count <- offspring_count + rowSums(current_generation)
      current_generation <- t(Xi) %*% current_generation
    }
    
    offspring_count <- offspring_count + 1
  }
  
  # Calculate the initial horizontal position for the founding case(s)
  founder_count <- length(Founders)
  founder_min_x <- (1:founder_count - 1) / founder_count
  founder_max_x <- (1:founder_count) / founder_count
  founder_center_position <- (founder_min_x + founder_max_x) / 2
  
  # Initiate the subsequent positioning calculations
  for (f in 1:founder_count) {	
    # First, set the founder's position
    position_matrix[Founders[f], 2] <- founder_center_position[f]
    
    # Then, find the first generation of cases descended from the founding precedent
    first_generation_ids <- which((t(Xi) %*% diag(1, T, T))[Founders[f],] == 1)
    first_generation_ids <- intersect(first_generation_ids, Descendents)
    
    # Then, calculate the horizontal position for each case in that generation
    position_matrix <- get_horizontal_position(
      position_matrix = position_matrix,
      parent_id = Founders[f],
      child_ids = first_generation_ids,
      parent_min_position = founder_min_x[f],
      parent_max_position = founder_max_x[f]
    )
  }
  
  # Flip first and second columns in the matrix
  position_matrix[,1:2] <- position_matrix[,2:1] 
  
  # Set file dimensions
  pdf(file = file, width = 13, height = 18)
  
  # Plot points
  plot(
    position_matrix[,1], # x-coordinates, i.e., case horizontal position
    position_matrix[,2], # y-coordinates, i.e., case year
    col = 'grey',
    pch = 8,
    ylim = range(c(case_years - 5), max(case_years + 5)),
    xlim = c(-0.2, 1.2),
    ylab = 'Year',
    xlab = '',
    main = title,
    axes = FALSE
  )
  axis(2)
  
  # Plot lines connecting each case to its children
  for (t in 1:T) {
    for (t2 in 1:t) {
      if (Xi[t, t2] == 1) {
        lines(
          c(position_matrix[t2, 1], position_matrix[t, 1]),
          c(position_matrix[t2, 2], position_matrix[t,2]),
          lwd = 1,
          col = 'grey'
        )
      }
    }
  }
  
  # Plot text labels for each point
  for (t in 1:T) {
    text(position_matrix[t, 1], position_matrix[t, 2], case_names[t], cex = 0.5)
  }	
  
  dev.off()
}