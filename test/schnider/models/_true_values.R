	trueparams <- list(
		c(tabname="THETA1", 	trueval=5.85,		label="V1",				logv=FALSE,		units="(l)"),
		c(tabname="THETA2", 	trueval=26,			label="V2",				logv=FALSE,		units="(l)"),
		c(tabname="THETA3", 	trueval=313, 		label="V3",				logv=FALSE,		units="(l)"),
		c(tabname="THETA4", 	trueval=1.92, 		label="CL",				logv=FALSE,		units="(l/min)"),
		c(tabname="THETA5", 	trueval=1.34, 		label="Q2",				logv=FALSE,		units="(l/min)"),
		c(tabname="THETA6", 	trueval=0.863, 		label="Q3",				logv=FALSE,		units="(l/min)"),
		c(tabname="OMEGA.1.1.", trueval=1.04E-01, 	label="var(V1)",		logv=FALSE,		units=""),
		c(tabname="OMEGA.2.2.", trueval=9.71E-02, 	label="var(V2)",		logv=FALSE,		units=""),
		c(tabname="OMEGA.4.4.", trueval=1.94E-02, 	label="var(CL)",		logv=FALSE,		units=""),
		c(tabname="OMEGA.5.5.", trueval=1.06E-01, 	label="var(Q2)",		logv=FALSE,		units=""),
		c(tabname="SIGMA.1.1.", trueval=5.15E-02, 	label="var(err)",		logv=FALSE,		units=""))

	setup_par <- function()
	{
		par(mfrow=c(3, 4), mar=c(2,2,2,0), oma=c(0,0,1,1), cex=0.9)
	}
