	trueparams <- list(
		c(tabname="THETA1", 	trueval=8.04E+00,	label="V",				logv=FALSE,	units="(l)"),
		c(tabname="THETA2", 	trueval=1.35E-01,	label="CL",				logv=FALSE,	units="(l/hr)"),
		c(tabname="THETA3", 	trueval=5.09E-01, 	label="KA",				logv=FALSE,	units="(1/hr)"),
		c(tabname="THETA5", 	trueval=9.65E+01, 	label="E0",				logv=FALSE,	units=""),
		c(tabname="THETA6", 	trueval=1.19E+00, 	label="Emax",			logv=FALSE,	units=""),
		c(tabname="THETA7", 	trueval=2.19E+00, 	label="C50",			logv=FALSE,	units=""),
		c(tabname="THETA8", 	trueval=-3.17E+00, 	label="log(kout)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.1.1.",	trueval=2.61E-02, 	label="var(V)",			logv=FALSE,	units=""),
		c(tabname="OMEGA.2.2.",	trueval=6.94E-02, 	label="var(CL)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.3.3.",	trueval=4.34E-01, 	label="var(Ka)",		logv=FALSE,	units=""),
		c(tabname="SIGMA.1.1.",	trueval=8.91E-03,	label="var(prop err PK)",	logv=FALSE,	units=""),
		c(tabname="SIGMA.2.2.",	trueval=1.26E+00, 	label="var(add err PK)",	logv=FALSE,	units=""),
		c(tabname="SIGMA.3.3.",	trueval=2.25E+01, 	label="var(add err PD)",	logv=FALSE,	units=""))

	setup_par <- function()
	{
		par(mfrow=c(4, 4), mar=c(2,2,2,0), oma=c(0,0,1,1), cex=0.85)
	}
