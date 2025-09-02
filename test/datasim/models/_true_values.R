	poptrue <- matrix(ncol=3, nrow=3)
	poptrue[1, 1] <- 0.218;
	poptrue[2, 2] <- 0.655;
	poptrue[3, 3] <- 0.023;
	poptrue[1, 2] <- 0.906 * sqrt(poptrue[1, 1] * poptrue[2, 2]);
	poptrue[2, 1] <- 0.906 * sqrt(poptrue[2, 2] * poptrue[1, 1]);
	poptrue[1, 3] <- 0.327 * sqrt(poptrue[1, 1] * poptrue[3, 3]);
	poptrue[3, 1] <- 0.327 * sqrt(poptrue[3, 3] * poptrue[1, 1]);
	poptrue[2, 3] <- 0.013 * sqrt(poptrue[2, 2] * poptrue[3, 3]);
	poptrue[3, 2] <- 0.013 * sqrt(poptrue[3, 3] * poptrue[2, 2]);

	tab <- data.frame()
	tab <- read.csv("datasim_paper_results.csv", header=TRUE, stringsAsFactors=FALSE)
	tab <- tab[, !(names(tab) %in% "X")]
	# consider remvoeing SAS FO not as a serious method
	tab <- tab[tab$name != "SAS FO", ]
	tab <- tab[tab$name != "nlme**", ]
	colnames(tab)[which(names(tab) == "SIGMA")] <- "SIGMA.1.1."

	readin <- function(fname)
	{
		d <- read.table(fname, header=TRUE, stringsAsFactors=FALSE)

		d$THETA1 <- exp(d$THETA1)
		d$THETA2 <- exp(d$THETA2)
		d$THETA3 <- exp(d$THETA3)

		for (i in 1:nrow(d)) {
			runres <- d[i, ]

			d[i, "corr21"] <- runres$OMEGA.2.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.2.2.) * 100
			d[i, "corr31"] <- runres$OMEGA.3.1. / sqrt(runres$OMEGA.1.1. * runres$OMEGA.3.3.) * 100
			d[i, "corr32"] <- runres$OMEGA.3.2. / sqrt(runres$OMEGA.2.2. * runres$OMEGA.3.3.) * 100

			popest <- matrix(ncol=3, nrow=3)
			popest[1, 1] <- runres$OMEGA.1.1.
			popest[2, 2] <- runres$OMEGA.2.2.
			popest[3, 3] <- runres$OMEGA.3.3.
			popest[1, 2] <- runres$OMEGA.2.1.
			popest[2, 1] <- runres$OMEGA.2.1.
			popest[1, 3] <- runres$OMEGA.3.1.
			popest[3, 1] <- runres$OMEGA.3.1.
			popest[2, 3] <- runres$OMEGA.3.2.
			popest[3, 2] <- runres$OMEGA.3.2.
			d[i,"det"] <- det(popest %*% solve(poptrue))^(1/3)
		}

		nfail <- sum(is.na(d[["det"]]))
		if (nfail > 0) {
			cat(sprintf("#####################################################\n"))
			cat(sprintf("determinant fails of %i of cases!!!\n", nfail))
			cat(sprintf("#####################################################\n"))
		}
		d
	}

	trueparams <- list(
		c(tabname="THETA1", 	trueval=27.2,			label="V",				logv=FALSE,	units=""),
		c(tabname="THETA2",		trueval=0.232,			label="KE",				logv=FALSE,	units="(1/time)"),
		c(tabname="THETA3",		trueval=0.304,			label="KA-KE",			logv=FALSE,	units="(1/time)"),
		c(tabname="OMEGA.1.1.", trueval=0.218,			label="var(V)",			logv=FALSE,	units=""),
		c(tabname="OMEGA.2.2.", trueval=0.655,			label="var(KE)",		logv=FALSE,	units=""),
		c(tabname="OMEGA.3.3.", trueval=0.023,			label="var(KA-KE)",		logv=FALSE,	units=""),
		c(tabname="SIGMA.1.1.",	trueval=0.064, 			label="var(err)",		logv=FALSE,	units=""),
		c(tabname="corr21", 	trueval=0.906 * 100,	label="corr(V,KE)",		logv=FALSE,	units="(%)"),
		c(tabname="corr31", 	trueval=0.327 * 100,	label="corr(V,KA-KE)",	logv=FALSE,	units="(%)"),
		c(tabname="corr32", 	trueval=0.013 * 100,	label="corr(Ke,KA-KE)",	logv=FALSE,	units="(%)"),
		c(tabname="det", 		trueval=1, 				label="determinant",	logv=FALSE,	units=""))

	setup_par <- function()
	{
		par(mfrow=c(3, 4), mar=c(2,2,2,0), oma=c(0,0,1,1), cex=0.9)
	}
