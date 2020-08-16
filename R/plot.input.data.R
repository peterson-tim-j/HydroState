##' @export
plot.input.data <- function(input.data, plot.yearRange=numeric()) {

  # Check inputs
  if (!is.data.frame(precipitation))
    stop('The input data must be a data frame.')

  colnames <- names(input.data)
  if (!any(colnames=='flow'))
    stop('The input data must have a column named "flow".')

  if (!any(colnames=='precipitation'))
    stop('The input data must have a column named "precipitation".')

  if (!any(colnames=='year'))
    stop('The input data must have a column named "year".')

  if (any(!is.finite(input.data$precipitation)))
    warning('The precipitation column should not contain gaps. This may compromise QModel objects with auto-correlation.')

  if (all(!is.finite(input.data$precipitation)))
    stop('The flow column contains only non-finite data.')

  if (max(abs(diff(input.data$year)))!=1)
    warning('The year column should not contain gaps. This may compromise QModel objects with auto-correlation.')


  # Get input grapics settings
  op <- par(no.readonly = T);

  # Change graphics settings
  par(mar = rep(2, 4))
  par(oma = rep(3, 4))
  par(mfrow=c(3,1))


  # Handle monthly and yearly time steps
  if (any(names(input.data)=="month")) {
    obsDates.asISO = ISOdate(input.data$year,data$month,1)
    obsDates = cbind(input.data$year,input.data$month)
    plot.units = 'month'
  } else {
    obsDates.asISO = ISOdate(input.data$year,12,1)
    obsDates = input.data$year
    plot.units = 'yr'
  }


  # Calc axis limits.
  ylim.flow <- c(0, ceiling(max( c(max(data$flow,na.rm=T),max(flow.viterbi.est,na.rm=T)))))

  # Setup year range for plotitng
  if (length(plot.yearRange)==2 && all(is.numeric(plot.yearRange)) && all(plot.yearRange>0) && all(plot.yearRange<=as.numeric(format(Sys.Date(), "%Y")))) {
    xlim = c(ISOdate(plot.yearRange[1],1,1), ISOdate(plot.yearRange[2],12,31))
  } else {
    xlim = c(min(obsDates.asISO)-7, max(obsDates.asISO)+7)
  }

  # Plot flow
  plot(obsDates.asISO, input.data$flow, type='l',col='blue', ylim=ylim.flow, xlim=xlim)

  mtext("Flow",side=2,line=3)
  mtext(paste("[mm /",plot.units,"]"),side=2,line=2, cex=0.75)
  mtext(paste("Date [",plot.units,']'),side=1,line=2)

  # Plot precip.
  plot(obsDates.asISO, input.data$precipitation, type='l',col='blue', xlim=xlim)

  mtext("Precipitation",side=2,line=3)
  mtext(paste("[mm /",plot.units,"]"),side=2,line=2, cex=0.75)
  mtext(paste("Date [",plot.units,']'),side=1,line=2)

  # Plot flow divided by precip.
  plot(obsDates.asISO, input.data$flow/input.data$precipitation, type='l',col='blue', xlim=xlim)

  mtext("Flow/Precipitation",side=2,line=3)
  mtext(paste("[-]"),side=2,line=2, cex=0.75)
  mtext(paste("Date [",plot.units,']'),side=1,line=2)

  # Plot precipitation vs flow
  plot(input.data$precipitation, input.data$flow, type='p',col='blue',
       ylim=ylim.flow, xlim=c(floor(min(input.data$precipitation)),ceiling(max(input.data$precipitation))),pch=21)

  # Add axis labels
  mtext("Flow",side=2,line=3)
  mtext(paste("[mm /",plot.units,"]"),side=2,line=2, cex=0.75)
  mtext(paste("Precipitation [mm /",plot.units,']'),side=1,line=1)

  # Put back original plotting options.
  par(op)

}
