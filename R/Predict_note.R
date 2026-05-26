rm(list=ls())
cat("\f")
library(hydroState)
library(DEoptim)
library(truncnorm)
library(devtools)


data(streamflow_annual_221201)

streamflow_annual_221201_until_2010 <- subset(
  streamflow_annual_221201,
  year <= 2010
)

streamflow_annual_221201_after_2010 <- subset(
  streamflow_annual_221201,
  year > 2010
)


default.model.annual <- hydroState::build(
  input.data = streamflow_annual_221201_until_2010
)

default.model.annual <- fit.hydroState(
  default.model.annual,
  pop.size.perParameter = 10,
  max.generations = 500
)

LLL <- hydroState:::predict(
  default.model.annual,
  t = nrow(streamflow_annual_221201_after_2010)
)

head(LLL)


pdf(
  "C:/Users/parsa/Downloads/Telegram Desktop/Predicted_PDF_with_observed_flow3.pdf",
  width = 10,
  height = 14
)

par(mfrow = c(6, 2), mar = c(4, 4, 3, 1))

for (j in 2:ncol(LLL)) {

  # Prediction step number
  h <- j - 1

  # Corresponding prediction year
  year_number <- streamflow_annual_221201_after_2010$year[h]

  plot(
    LLL$range,
    LLL[[j]],
    type = "l",
    lwd = 2,
    xlab = "Flow",
    ylab = "Density",
    main = paste("Predicted PDF -", year_number)
  )

  obs_flow <- streamflow_annual_221201_after_2010$flow[
    streamflow_annual_221201_after_2010$year == year_number
  ]

  if (length(obs_flow) == 1 && !is.na(obs_flow)) {

    abline(
      v = obs_flow,
      col = "red",
      lwd = 2,
      lty = 2
    )

    legend(
      "topright",
      legend = "Observed flow",
      col = "red",
      lty = 2,
      lwd = 2,
      bty = "n"
    )
  }
}

dev.off()
