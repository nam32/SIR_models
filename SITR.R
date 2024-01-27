# loading libraries
library(deSolve)
library(reshape2)
library(ggplot2)

# defining sir function:
# β: beta infection rate for untreated individuals
# δ: delta recovery rate for treated infected individuals
# ϵ: epsilon rate at which treated individuals revert to the infected compartment
# γ: gamma recovery rate for untreated infected individuals
# μ: mu infection rate for treated individuals
# θ: theta treatment initiation rate for infected individuals

sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -(beta * I + mu * T) * S 
    dI <- (beta * I + mu * T)* S - (gamma + theta) * I + epsilon * T
    dT <- theta * I - (epsilon + delta) * T
    dR <- gamma * I + delta * T
    return(list(c(dS, dI, dR, dT)))
  })
}

# initial proportions in each compartment
init <- c(S = 0.999999, 
          I = 0.000001, 
          R = 0, 
          T = 0)

# set parameters
parameters <- c(beta = 1.5, 
                delta = 2, 
                epsilon = 0, 
                gamma = 0.1,
                mu = 0.1,
                theta = 0.4 
                )

# set time frame
times <- seq(0, 70, by = 0.01)

# solve using ode
out <- ode(y = init, times = times, func = sir, parms = parameters)

# convert to a data frame
out <- as.data.frame(out)
head(out)
out_long <- melt(out, id.vars = "time", measure.vars = c("S", "I", "T", "R"), 
                 variable.name = "Compartment", value.name = "value")
names(out_long) = c("time","Compartment","value")

##############################
####part 1: plotting graph####
##############################

ggplot(data = out_long, aes(x = time, y = value / max(value))) +
  geom_line(aes(col = Compartment, group = Compartment)
            #, size = 1.2
            ) +
  labs(title = "SITR Model") +
  xlab("Time") +
  ylab("Proportion of Population") +
  scale_color_manual(values = c("#00FF00", "#FF0000", "#56B4E9", "#A020F0")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 1)

###########################################################
####part 2: varying parameters to check for its effects####
###########################################################

plot_list <- list()

parameters <- c(beta = 1.5, 
                delta = 2, 
                epsilon = 0, 
                gamma = 0.1,
                mu = 0.1,
                theta = 0.4 
)

varied_parameter <- "delta"
# "delta"  seq(0.1, 1, by = 0.25)

# varied_parameter <- "mu"
# "mu"  seq(0.1, 1.5, by = 0.1)

# Loop through different values of varied_parameter
for (val in  seq(0.1, 1, by = 0.25)) {
  # Update the delta parameter
  parameters[varied_parameter] <- val
  
  # Solve the differential equations
  out <- ode(y = init, times = seq(0, 70, by = 0.01), func = sir, parms = parameters)
  
  # Convert to a data frame
  out <- as.data.frame(out)
  
  out_long <- melt(out, id.vars = "time", measure.vars = c("S", "I", "T", "R"), 
                   variable.name = "Compartment", value.name = "value")
  
  names(out_long) <- c("time", "Compartment", "value")

  # Plotting graph
  plot_obj <- ggplot(data = out_long, aes(x = time, y = value / max(value))) +
    geom_line(aes(col = Compartment, group = Compartment)#, size = 1.2
              ) + labs(title = paste(varied_parameter, val,"SITR Model")) +
    xlab("Time") +
    ylab("Proportion of Population") +
    scale_color_manual(values = c("#00FF00", "#FF0000", "#56B4E9", "#A020F0")) +
    theme_minimal() +
    theme(#text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    ylim(0, 1)
  
  # Store the plot object in the list
  plot_list[[as.character(val)]] <- plot_obj
}

# Display the plots
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
}

# to delete all the plots:
# dev.off(dev.list()["RStudioGD"])

#############################################################################################################################
####part 3: plotting Recover for different parameter to visualize the effect of changes in parameter to the model dynamic####
#############################################################################################################################

# Function to vary parameters
vary_parameter <- function(param_name, param_values, init, times) {
  results <- data.frame()
  
  for (val in param_values) {
    parameters <- c(beta = 1.5, gamma = 0.1, delta = 2, theta = 0.4, epsilon = 0, mu = 0.1) # Initial parameters for the function
    parameters[param_name] <- val
    
    out <- ode(y = init, times = times, func = sir, parms = parameters)
    out_long <- melt(as.data.frame(out), id.vars = "time", measure.vars = c("S", "I", "T", "R"), 
                     variable.name = "Compartment", value.name = "value")
    out_long[[param_name]] <- val
    results <- rbind(results, out_long)
  }
  
  return(results)
}

# Initial proportions in each compartment
init <- c(S = 0.999999, I = 0.000001, R = 0, T = 0)

# Time frame
times <- seq(0, 70, by = 0.01)

# Varying parameter and its sequence
param_name <- "mu"  # parameter to vary
param_values <- seq(0.1, 2.6, by = 0.5)  # set sequence to vary the parameter by

# param_name <- "delta"  # parameter to vary
# param_values <- seq(0.1, 1, by = 0.25)  # set sequence to vary the parameter by

# Generate results for varying parameter
results <- vary_parameter(param_name, param_values, init, times)

# Plotting
ggplot(results[results$Compartment == "R",], aes(x = time, y = value, color = as.factor(.data[[param_name]]))) +
  geom_line() +
  labs(x = "Time", y = "R", color = paste(param_name, "Value")) +
  theme_minimal() +
  ggtitle(paste("Recovery Over Time for Different", param_name, "Values")) +
  theme(plot.title = element_text(hjust = 0.5))


#############################
####part 4: phase diagram####
#############################

# Plotting Infected (I) against Susceptible (S)
ggplot(data = out, aes(x = S, y = I, color = time)) +
  geom_point() +
  labs(title = "Phase Diagram") +
  xlab("Susceptible (S)") +
  ylab("Infected (I)") +
  scale_color_gradientn(colors = c("blue", "green", "yellow", "red")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0, 1) +
  ylim(0, 1)
