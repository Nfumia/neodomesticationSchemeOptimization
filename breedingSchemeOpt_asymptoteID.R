allSummary <- read.csv("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/allSummary.csv")
allSummary[10:17] <- zoo::na.locf(allSummary[10:17])
wrking <- subset(allSummary,NeStart==50 & Nparents==15 & Nprogeny==20 & h2 == "hLL" & si == 0.75 & sel=="MxAv" & pop=="orphan")

# Sample data (replace with your own dataset)
x <- wrking$Cycle
y <- wrking$meanP2

# Plot the original data for comparison
#plot(x, y, type = "l", main = "Original Data")

# Perform kernel smoothing
#smoothed_data <- ksmooth(x, y, "normal", bandwidth = 3)

# Plot the smoothed data
#plot(smoothed_data$x, smoothed_data$y, type = "l", main = "Smoothed Data")

# Calculate the first derivative
first_derivative <- diff(y) / diff(x)

# Calculate the second derivative
second_derivative <- diff(first_derivative) / diff(x[-1])

### GAIN (INCREASING SIGMOIDAL)
# Fit a sigmoidal curve to the data if GAIN (INCREASING or Upper Asymptote)
fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), start = list(Asym = max(y), xmid = mean(x), scal = 1))
# Extract the parameters of the sigmoidal curve
Asym_estimate <- coef(fit)["Asym"]
xmid_estimate <- coef(fit)["xmid"]
scal_estimate <- coef(fit)["scal"]
# Evaluate the fitted sigmoidal curve at the x-values
predicted_y_values <- predict(fit, newdata = data.frame(x = x))
# Calculate the absolute difference between the predicted y-values and the asymptote
difference <- abs(predicted_y_values - Asym_estimate)
# Find the x-value where the difference becomes negligible (below the threshold)
threshold <- 0.01  # Adjust this threshold as needed
negligible_difference_x <- x[min(which(difference < threshold))]
# Find the corresponding y-value at the negligible difference point
negligible_difference_y <- Asym_estimate
# Print the coordinates of the point with negligible difference
cat("Point with Negligible Difference (x, y):", negligible_difference_x, ",", negligible_difference_y, "\n")
# Print the estimated asymptote threshold
cat("Sigmoid midpoint (x-value):", xmid_estimate, "\n")


### VARIANCE (DECREASING SIGMOIDAL/LOGIT)
# Fit a sigmoidal curve to the data if VARIANCE (DECREASING or Lower Asymptote)
fit <- nls(y ~ SSfpl(x, A, B, xmid, scal), start = list(A = min(y), B = max(y), xmid = mean(x), scal = 1),control = nls.control(minFactor = 0.0001,tol = 1e-4,maxiter = 100))
# Extract the parameters of the sigmoidal curve
AsymU_estimate <- coef(fit)["A"]
AsymL_estimate <- coef(fit)["B"]
xmid_estimate <- coef(fit)["xmid"]
scal_estimate <- coef(fit)["scal"]
# Evaluate the fitted sigmoidal curve at the x-values
predicted_y_values <- predict(fit, newdata = data.frame(x = x))
# Calculate the absolute difference between the predicted y-values and the asymptote
difference <- abs(predicted_y_values - AsymL_estimate)
# Find the x-value where the difference becomes negligible (below the threshold)
threshold <- 0.01  # Adjust this threshold as needed
negligible_difference_x <- x[min(which(difference < threshold))]
# Find the corresponding y-value at the negligible difference point
negligible_difference_y <- AsymL_estimate
# Print the coordinates of the point with negligible difference
cat("Point with Negligible Difference (x, y):", negligible_difference_x, ",", negligible_difference_y, "\n")
# Print the estimated asymptote threshold
cat("Sigmoid midpoint (x-value):", xmid_estimate, "\n")



# Plot the original data and the fitted curve
plot(x, y, type = "l", xlab="Cycle",ylab="Phenotypic Value")
curve(predict(fit, list(x = x)), col = "blue", add = TRUE)

# Find the x-value where the curve reaches the asymptote
asymptote_x <- xmid_estimate

# Print the estimated asymptote threshold
cat("Estimated Asymptote Threshold (x-value):", asymptote_x, "\n")

# Add a vertical line to visualize the asymptote threshold
abline(v = asymptote_x, col = "red")

# Add a horizontal line to visualize the asymptote
abline(h= Asym_estimate,col="green")

##########################################
### Loop Above

# Define parameter combinations
parameter_combinations <- list(
  NeStart = c(25,50,100),  # Replace with your specific values
  Nparents = c(10,15,20,25),  # Replace with your specific values
  Nprogeny = c(1,5,10,20),  # Replace with your specific values
  h2 = c("hHL","hLH","hLL","hMM"),  # Replace with your specific values
  si = c(0.25,0.5,0.75),  # Replace with your specific values
  sel = c("MxAv","PRS","GS"),  # Replace with your specific values
  pop = c("wild","orphan","landrace")  # Replace with your specific values
)

# Create an empty dataframe to store results
results_df <- data.frame(
  NeStart = numeric(),
  Nparents = numeric(),
  Nprogeny = numeric(),
  h2 = character(),
  si = numeric(),
  sel = character(),
  pop = character(),
  ResponseVariable = character(),
  Negligible_x = numeric(),
  Xmid_Estimate = numeric()
)

# Get the names of columns starting with "var" or "mean"
var_mean_columns <- colnames(allSummary[c(11,12,15,16)])

# Loop through different parameter combinations
for (NeStart_val in parameter_combinations$NeStart) {
  for (Nparents_val in parameter_combinations$Nparents) {
    for (Nprogeny_val in parameter_combinations$Nprogeny) {
      for (h2_val in parameter_combinations$h2) {
        for (si_val in parameter_combinations$si) {
          for (sel_val in parameter_combinations$sel) {
            for (pop_val in parameter_combinations$pop) {
              
              # Initialize variables outside the inner loop
              negligible_x <- NA
              xmid_estimate <- NA
              
              # Subset the data based on the current parameter combination
              subset_data <- subset(
                allSummary,
                NeStart == NeStart_val &
                  Nparents == Nparents_val &
                  Nprogeny == Nprogeny_val &
                  h2 == h2_val &
                  si == si_val &
                  sel == sel_val &
                  pop == pop_val
              )
              
              # Loop through columns starting with "var" or "mean"
              for (column_name in var_mean_columns) {
                # Extract the response variable (y)
                x <- subset_data$Cycle
                y <- subset_data[, column_name]
                
                # Calculate the first derivative
                first_derivative <- diff(y) / diff(x)
                
                # Calculate the second derivative
                second_derivative <- diff(first_derivative) / diff(x[-1])
                
                # Code that may generate errors wrapped in tryCatch
                tryCatch({
                if (startsWith(column_name, "var")) {
                  # Apply the VARIANCE (DECREASING SIGMOIDAL/LOGIT) method
                  
                  # Extract the response variable (y)
                  x <- subset_data$Cycle
                  y <- subset_data[, column_name]
                  
                  # Calculate the first derivative
                  first_derivative <- diff(y) / diff(x)
                  
                  # Calculate the second derivative
                  second_derivative <- diff(first_derivative) / diff(x[-1])
                  
                  # Fit a sigmoidal curve to the data if VARIANCE (DECREASING or Lower Asymptote)
                  fit <- nls(y ~ SSfpl(x, A, B, xmid, scal), start = list(A = min(y), B = max(y), xmid = mean(x), scal = 1))
                  # Extract the parameters of the sigmoidal curve
                  AsymU_estimate <- coef(fit)["A"]
                  AsymL_estimate <- coef(fit)["B"]
                  xmid_estimate <- coef(fit)["xmid"]
                  scal_estimate <- coef(fit)["scal"]
                  # Evaluate the fitted sigmoidal curve at the x-values
                  predicted_y_values <- predict(fit, newdata = data.frame(x = x))
                  # Calculate the absolute difference between the predicted y-values and the asymptote
                  difference <- abs(predicted_y_values - AsymL_estimate)
                  # Find the x-value where the difference becomes negligible (below the threshold)
                  threshold <- 0.1  # Adjust this threshold as needed
                  negligible_difference_x <- ifelse(any(difference < threshold), x[min(which(difference < threshold))], "UNKNOWN")
                  # Find the corresponding y-value at the negligible difference point
                  negligible_difference_y <- AsymL_estimate
                  # Print the coordinates of the point with negligible difference
                  cat("Point with Negligible Difference (x, y):", negligible_difference_x, ",", negligible_difference_y, "\n")
                  # Print the estimated asymptote threshold
                  cat("Sigmoid midpoint (x-value):", xmid_estimate, "\n")
                  
                } else if (startsWith(column_name, "mean")) {
                  # Apply the GAIN (INCREASING SIGMOIDAL) method
                  
                  # Extract the response variable (y)
                  y <- subset_data[, column_name]
                  
                  # Calculate the first derivative
                  first_derivative <- diff(y) / diff(x)
                  
                  # Calculate the second derivative
                  second_derivative <- diff(first_derivative) / diff(x[-1])
                  
                  # Fit a sigmoidal curve to the data if GAIN (INCREASING or Upper Asymptote)
                  fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), start = list(Asym = max(y), xmid = mean(x), scal = 1))
                  # Extract the parameters of the sigmoidal curve
                  Asym_estimate <- coef(fit)["Asym"]
                  xmid_estimate <- coef(fit)["xmid"]
                  scal_estimate <- coef(fit)["scal"]
                  # Evaluate the fitted sigmoidal curve at the x-values
                  predicted_y_values <- predict(fit, newdata = data.frame(x = x))
                  # Calculate the absolute difference between the predicted y-values and the asymptote
                  difference <- abs(predicted_y_values - Asym_estimate)
                  # Find the x-value where the difference becomes negligible (below the threshold)
                  threshold <- 0.1  # Adjust this threshold as needed
                  negligible_difference_x <- ifelse(any(difference < threshold), x[min(which(difference < threshold))], "UNKNOWN")
                  # Find the corresponding y-value at the negligible difference point
                  negligible_difference_y <- Asym_estimate
                  # Print the coordinates of the point with negligible difference
                  cat("Point with Negligible Difference (x, y):", negligible_difference_x, ",", negligible_difference_y, "\n")
                  # Print the estimated asymptote threshold
                  cat("Sigmoid midpoint (x-value):", xmid_estimate, "\n")
                }
                }, error = function(e) {
                  # Handle errors here, for example:
                  cat("Error occurred:", conditionMessage(e), "\n")
                })
                # Create a row of results and add it to the dataframe
                result_row <- data.frame(
                  NeStart = NeStart_val,
                  Nparents = Nparents_val,
                  Nprogeny = Nprogeny_val,
                  h2 = h2_val,
                  si = si_val,
                  sel = sel_val,
                  pop = pop_val,
                  ResponseVariable = column_name,
                  Negligible_x = negligible_difference_x,
                  Xmid_Estimate = xmid_estimate,
                  Scal_Estimate = scal_estimate
                )
                
                # Append the result_row to results_df
                results_df <- rbind(results_df, result_row)
              }
            }
          }
        }
      }
    }
  }
}

asymptoteInfo <- results_df
write.csv(asymptoteInfo,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/asymptoteInfo.csv")
asymptoteInfo <- read.csv("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/asymptoteInfo.csv")
#asymptoteInfo <- read.csv("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/asymptoteInfo.csv")
asymptoteInfo[asymptoteInfo == "UNKNOWN"] <- NA
known_data <- asymptoteInfo[complete.cases(asymptoteInfo), ] 
unknown_data <- asymptoteInfo[!complete.cases(asymptoteInfo), ] 

library(Hmisc)
with(asymptoteInfo,rcorr(Negligible_x,Xmid_Estimate))


# Linear model to predict Negligible_x
model <- lm(Negligible_x ~ NeStart + Nparents + Nprogeny + h2 + si + sel + pop + ResponseVariable + Xmid_Estimate, data = known_data)
unknown_predictions <- predict(model, newdata = unknown_data)

# Make the predictions into a dataframe
unknown_predictions <- data.frame(unknown_predictions)
colnames(unknown_predictions)[1] <- "Negligible_x"

# Add the predicted values for Negligible_x to the original dataframe
asymptoteInfo[rownames(unknown_predictions), "Negligible_x"] <- unknown_predictions$Negligible_x
asymptoteInfo$scheme <- paste(asymptoteInfo$pop, asymptoteInfo$sel,asymptoteInfo$h2,asymptoteInfo$si, sep = ".")
asymptoteInfo$Negligible_x <- as.numeric(asymptoteInfo$Negligible_x)

summary_asymptote <- aggregate(Negligible_x ~ ResponseVariable + scheme, data = asymptoteInfo, FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
summary_asymptote <- do.call(data.frame,summary_asymptote)
colnames(summary_asymptote)[3] <- "mean"
colnames(summary_asymptote)[4] <- "se"

# Split 'scheme' column into the 4 identifying criteria
split_strings <- strsplit(summary_asymptote$scheme, "\\.")
# Create a new dataframe with the split values in separate columns
new_df <- data.frame(do.call(rbind, split_strings))
# Combine the two columns with a period between them
new_df$X6 <- paste(new_df$X4, new_df$X5, sep = ".")
new_df <- new_df[,-c(4,5)]
colnames(new_df) <- c("pop", "sel", "h2", "si")
summary_asymptote <- cbind(new_df,summary_asymptote)

write.csv(summary_asymptote,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Dissertation/github/domsim/out/Scheme/summary_asymptote.csv")

meanP1 <- subset(summary_asymptote,ResponseVariable=="meanP1")
meanP2 <- subset(summary_asymptote,ResponseVariable=="meanP2")
var1 <- subset(summary_asymptote,ResponseVariable=="var1")
var2 <- subset(summary_asymptote,ResponseVariable=="var2")


# Plot a facet wrap bar graph around scheme type within each gain and variance
library(ggplot2)
ggplot(var1, aes(x=scheme,y=mean,fill=interaction(pop,sel))) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x="Scheme",y="Cycle Asymptote of Phenotypic Variance",title="Simple Oligogenic Trait (z1)") +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=4),
        strip.text.x=element_text(angle=0,hjust=0.5),
        axis.text.y=element_text(size=5)) +
  scale_y_continuous(limits=c(0,80),breaks=c(0,20,40,60,80))+
  labs(fill="Population and Selection") +
  scale_fill_manual(labels=c("Landrace GS","Orphan GS","Wild GS","Landrace MxAv","Orphan MxAv",
                             "Wild MxAv","Landrace PRS","Orphan PRS","Wild PRS"),
                    values=c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
                             "#AA4499","#44AA99","#999933","#882255")) +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) 

