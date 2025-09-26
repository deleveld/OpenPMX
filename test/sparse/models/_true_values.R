	trueparams <- list(
		c(tabname="THETA1", 	trueval=1.38629436112,	label="CL",				logv=TRUE,	units="(l/hr)"),
		c(tabname="THETA2", 	trueval=4.24849524205,	label="V",				logv=TRUE,	units="(l)"),
		c(tabname="THETA3", 	trueval=0, 				label="KA",				logv=TRUE,	units="(1/hr)"),
		c(tabname="OMEGA.1.1.",	trueval=0.09, 			label="var(CL)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.2.2.",	trueval=0.09, 			label="var(V)",			logv=FALSE,	units=""),
		c(tabname="OMEGA.3.3.",	trueval=0.09, 			label="var(KA)",		logv=FALSE,	units=""),
		c(tabname="SIGMA.1.1.",	trueval=0.04, 			label="var(err)",		logv=FALSE,	units=""))

	setup_par <- function()
	{
		par(mfrow=c(3, 4), mar=c(2,2,2,0), oma=c(0,0,1,1), cex=0.9)
	}
