	# the true omega matrix 
	omegaltvals <- c(4.57E-01,
			2.80E-01,  5.16E-01,
			8.84E-02,  1.15E-01,  1.35E-01,
			4.68E-02,  1.31E-02,  7.97E-03,  9.46E-02,
			2.00E-01,  3.43E-01,  5.89E-02,  2.89E-02,  4.51E-01,
			8.26E-02, -5.14E-02, -1.91E-02,  2.26E-02,  1.61E-01,  2.45E-01)
	n <- as.integer(round((sqrt(8*length(omegaltvals)+1)-1)/2,0))
	omega <- matrix(nrow=n, ncol=n)
	omega[upper.tri(omega, diag=TRUE)] = omegaltvals
	omega[lower.tri(omega)] = t(omega)[lower.tri(omega)]
	corr <- cov2cor(omega)
	corr <- corr * 100.

	trueparams <- list(
		c(tabname="THETA1", 	trueval=1.01E+01,	label="V1",				logv=FALSE,	units="(l)"),
		c(tabname="THETA2", 	trueval=2.63E+01,	label="V2",				logv=FALSE,	units="(l)"),
		c(tabname="THETA3", 	trueval=2.05E+02, 	label="V3",				logv=FALSE,	units="(l)"),
		c(tabname="THETA4", 	trueval=7.03E-01, 	label="CL",				logv=FALSE,	units="(l/min)"),
		c(tabname="THETA5", 	trueval=2.36E+00, 	label="Q2",				logv=FALSE,	units="(l/min)"),
		c(tabname="THETA6", 	trueval=1.49E+00, 	label="Q3",				logv=FALSE,	units="(l/min)"),
		c(tabname="THETA7", 	trueval=1.23E+00, 	label="V scale",		logv=FALSE,	units=""),
		c(tabname="THETA8", 	trueval=3.10E-01, 	label="CL scale",		logv=FALSE,	units=""),
		c(tabname="OMEGA.1.1.", trueval=4.57E-01, 	label="var(V1)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.2.2.", trueval=5.16E-01, 	label="var(V2)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.3.3.", trueval=1.35E-01, 	label="var(V3)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.4.4.", trueval=9.46E-02, 	label="var(CL)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.5.5.", trueval=4.51E-01, 	label="var(Q2)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.6.6.", trueval=2.45E-01, 	label="var(Q3)",		logv=FALSE,	units=""),
		c(tabname="SIGMA.1.1.", trueval=2.62E-02, 	label="var(err)",		logv=FALSE,	units=""),
		c(tabname="corr21", 	trueval=corr[2,1],	label="corr(V2,V1)",	logv=FALSE,	units="(%)"),
		c(tabname="corr31", 	trueval=corr[3,1], 	label="corr(V3,V1)",	logv=FALSE,	units="(%)"),
		c(tabname="corr32", 	trueval=corr[3,2], 	label="corr(V3,V2)",	logv=FALSE,	units="(%)"),
		c(tabname="corr41", 	trueval=corr[4,1], 	label="corr(CL,V1)",	logv=FALSE,	units="(%)"),
		c(tabname="corr42", 	trueval=corr[4,2], 	label="corr(CL,V2)",	logv=FALSE,	units="(%)"),
		c(tabname="corr43", 	trueval=corr[4,3], 	label="corr(CL,V3)",	logv=FALSE,	units="(%)"),
		c(tabname="corr51", 	trueval=corr[5,1], 	label="corr(Q2,V1)",	logv=FALSE,	units="(%)"),
		c(tabname="corr52", 	trueval=corr[5,2], 	label="corr(Q2,V2)",	logv=FALSE,	units="(%)"),
		c(tabname="corr53", 	trueval=corr[5,3], 	label="corr(Q2,V3)",	logv=FALSE,	units="(%)"),
		c(tabname="corr54", 	trueval=corr[5,4], 	label="corr(Q2,CL)",	logv=FALSE,	units="(%)"),
		c(tabname="corr61",		trueval=corr[6,1], 	label="corr(Q3,V1)",	logv=FALSE,	units="(%)"),
		c(tabname="corr62", 	trueval=corr[6,2], 	label="corr(Q3,V2)",	logv=FALSE,	units="(%)"),
		c(tabname="corr63", 	trueval=corr[6,3], 	label="corr(Q3,V3)",	logv=FALSE,	units="(%)"),
		c(tabname="corr64",		trueval=corr[6,4], 	label="corr(Q3,CL)",	logv=FALSE,	units="(%)"),
		c(tabname="corr65", 	trueval=corr[6,5], 	label="corr(Q3,Q2)",	logv=FALSE,	units="(%)"))

	readin <- function(fname)
	{
		d <- read.table(file=fname, header=TRUE)
		for (i in 1:nrow(d)) {
			runres <- d[i, ]

			d[i, "corr21"] <- runres$OMEGA.2.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.2.2.) * 100
			
			d[i, "corr31"] <- runres$OMEGA.3.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.3.3.) * 100
			d[i, "corr32"] <- runres$OMEGA.3.2. / sqrt(runres$OMEGA.2.2. * runres$OMEGA.3.3.) * 100

			d[i, "corr41"] <- runres$OMEGA.4.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.4.4.) * 100
			d[i, "corr42"] <- runres$OMEGA.4.2. / sqrt(runres$OMEGA.2.2. * runres$OMEGA.4.4.) * 100
			d[i, "corr43"] <- runres$OMEGA.4.3. / sqrt(runres$OMEGA.3.3. * runres$OMEGA.4.4.) * 100

			d[i, "corr51"] <- runres$OMEGA.5.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.5.5.) * 100
			d[i, "corr52"] <- runres$OMEGA.5.2. / sqrt(runres$OMEGA.2.2. * runres$OMEGA.5.5.) * 100
			d[i, "corr53"] <- runres$OMEGA.5.3. / sqrt(runres$OMEGA.3.3. * runres$OMEGA.5.5.) * 100
			d[i, "corr54"] <- runres$OMEGA.5.4. / sqrt(runres$OMEGA.4.4. * runres$OMEGA.5.5.) * 100

			d[i, "corr61"] <- runres$OMEGA.6.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.6.6.) * 100
			d[i, "corr62"] <- runres$OMEGA.6.2. / sqrt(runres$OMEGA.2.2. * runres$OMEGA.6.6.) * 100
			d[i, "corr63"] <- runres$OMEGA.6.3. / sqrt(runres$OMEGA.3.3. * runres$OMEGA.6.6.) * 100
			d[i, "corr64"] <- runres$OMEGA.6.4. / sqrt(runres$OMEGA.4.4. * runres$OMEGA.6.6.) * 100
			d[i, "corr65"] <- runres$OMEGA.6.5. / sqrt(runres$OMEGA.5.5. * runres$OMEGA.6.6.) * 100
		}
		d
	}

	setup_par <- function()
	{
		par(mfrow=c(5, 6), mar=c(2,2,2,0), oma=c(0,0,1,1), cex=0.7)
	}
