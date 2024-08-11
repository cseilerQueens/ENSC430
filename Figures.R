# Lotka-Volterra Pearson Correlation
library(deSolve)

interest <- function(t, state, parameters)
{
  with(as.list(c(state)), {
    alpha <- 0.2
    ri <- 1
    m <- 0.2
    gamma <- 0.5
    K <- 100
    
    dY <- ri * Y * (1 - Y / K) - alpha * R * Y
    dR <- alpha * gamma * R * Y - m * R
    
    list(c(dY, dR))
  })
}
state <- c(Y = 1, R = 1)
times <- seq(0, 200, 1)
out <- as.data.frame(ode(state, times, interest))

# Perform the correlation test
correlation_test <- cor.test(out$Y, out$R)

# Get the correlation coefficient
correlation_coefficient <- round(correlation_test$estimate, 2)

# Get the p-value
p_value <- round(correlation_test$p.value, 2)

legend_text <- paste("R = ", correlation_coefficient, ", p-value = ", p_value, sep = "")

png(
  filename = "Lotka-Volterra-correlation.png",
  width = 10,
  height = 10,
  units = "cm",
  res = 500
)
par(pty = "s", mar = c(4, 4, 1, 1))
plot(out$time,
     out$Y,
     type = "l",
     xlab = "Time Step",
     ylab = "Species Density")
lines(out$time, out$R, col = "red")
legend("bottomright",
       c("Predator", "Prey"),
       text.col = c("red", "black"))
legend("topright", legend_text, bty = "n")
dev.off()


# Autoregression
# Forecast the next 100 steps
forecast_length <- 100
pred <- stats::predict(ar_model, n.ahead = forecast_length)

# Combine actual and forecasted values
combined_series <- c(time_series, pred$pred)

# Plot the actual time series and the forecast

png(
  filename = "forecast.png",
  width = 15,
  height = 10,
  units = "cm",
  res = 500
)
par(mar = c(4, 4, 1, 1))
plot(
  combined_series,
  type = "l",
  col = "red",
  lwd = 2,
  xlab = "Time",
  ylab = "Value"
)

lines(1:n, time_series, col = "black", lwd = 2)  # Actual time series in black

# Add a legend
legend(
  "bottomleft",
  legend = c("Actual", "Forecast"),
  col = c("black", "red"),
  lwd = 2
)

dev.off()

# Granger
rm(list = ls())
library(vars)
# Create a sinusoidal time series
t <- seq(0, 20, 0.1) # time axis
n <- length(t) # number of time steps

x <- sin(t)
y <- x[1:(n-9)]
x <- x[10:n]

# Add noise to y
noise <- runif(n=length(y), min = -0.5, max = 0.5)
y <- y + noise

x <- ts(x)
y <- ts(y)

# (2) Create a vector autoregression model
tsDat <- ts.union(x, y)

png(
  filename = "granger.png",
  width = 15,
  height = 10,
  units = "cm",
  res = 500
)
par(mar = c(4, 4, 1, 1))
plot(tsDat, main = "")
dev.off()


