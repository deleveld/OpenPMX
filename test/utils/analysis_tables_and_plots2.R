	performance_loss <- function(g, gname, baseg, basename)
	{
		perftable <- data.frame()

		for (i in trueparams) {
			paramname <- i["tabname"]
			trueval <- as.numeric(i["trueval"])
			displayname <-i["label"]
			logv <-i["logv"]

			if (logv == TRUE)
				trueval <- exp(trueval)

			v1 <- g[[paramname]]
			if (logv == TRUE)
				v1 <- exp(v1)
			v2 <- baseg[[paramname]]
			if (logv == TRUE)
				v2 <- exp(v2)

			# filter out NA
			valid <- !is.na(v1) & !is.na(v2)
			v1 <- v1[valid]
			v2 <- v2[valid]

			pest1 <- abs((v1 - trueval) / trueval) * 100
			pest2 <- abs((v2 - trueval) / trueval) * 100
###			err1 <- v1 - trueval
###			err2 <- v2 - trueval
###			pest1 <- (sqrt(err1*err1) / trueval) * 100
###			pest2 <- (sqrt(err2*err2) / trueval) * 100
			xlim <- range(c(pest1, pest2))

			v <- pest1 - pest2
			w <- wilcox.test(v, alternative="two.sided")
			pval <- w[["p.value"]]

			plot(0, type="n", yaxt="n", xlab=NA, ylab=NA, xlim=(xlim), ylim=c(0,1.4), main=NA)
			mtext("Density", side=2, line=0.5, cex=par()$cex)

			m <- median(v)
			better <- ""
			if (m < 0 && pval < 0.05) {
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgreen")
				better <- gname
			}
			if (m > 0 && pval < 0.05) {
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey")
				better <- basename
			}

			perftable[displayname, "delta"] <- abs(m)
			perftable[displayname, "pval"] <- pval
			perftable[displayname, "better"] <- better
			perftable <- perftable[order(perftable[["pval"]]), ]

			d1 <- density((pest1), na.rm=TRUE)
			d2 <- density((pest2), na.rm=TRUE)
			lines(x=(d2$x), y=d2$y/max(d2$y), col="darkgrey", lw=3)
			lines(x=(d1$x), y=d1$y/max(d1$y), col="black", lw=1)

			grid()
			box()

			legend("topright",
				legend=c(basename,
						 gname,
						 sprintf("loss=%.1f%%", m),
						 sprintf("p=%.3f", pval)),
				lw=c(3, 1, NA, NA),
				col=c("darkgrey", "black"),
				bty="n", cex=par()$cex*1.2)

			t <- bquote(.(displayname)~"|"~Delta~"|"~"(%)")
#			t <- sprintf("%s Δ(%%)", displayname)
			title(t, line=0.5)
		}
		perftable
	}
	if (exists("gronmem") && exists("nonmem")) {
		setup_par()
		t <- "OpenPMX vs nonmem"
		d <- performance_loss(gronmem, "OpenPMX", nonmem, "NONMEM")
		title(paste0("Distribution of RSE ", t), outer=TRUE)
		print(t)
		print(d)
	}
	if (exists("OpenPMX") && exists("nonmemlaplace")) {
		setup_par()
		t <- "OpenPMX vs laplace"
		d <- performance_loss(gronmem, "OpenPMX", nonmemlaplace, "LAPLACE")
		title(paste0("Distribution of RSE ", t), outer=TRUE)
		print(t)
		print(d)
	}
	if (exists("gronmemicov") && exists("nonmem")) {
		setup_par()
		t <- "OpenPMX_icov vs nonmem"
		d <- performance_loss(gronmemicov, "OpenPMX_icov", nonmem, "NONMEM")
		title(paste0("Distribution of RSE ", t), outer=TRUE)
		print(t)
		print(d)
	}
	if (exists("gronmemicov") && exists("gronmem")) {
		setup_par()
		t <- "OpenPMX_icov vs OpenPMX"
		d <- performance_loss(gronmemicov, "OpenPMX_icov", gronmem, "gronmem")
		title(paste0("Distribution of RSE ", t), outer=TRUE)
		print(t)
		print(d)
	}
	if (exists("gronmemtest") && exists("gronmem")) {
		setup_par()
		t <- "Distribution of RSE gronmemtest vs OpenPMX"
		d <- performance_loss(gronmemtest, "test", gronmem, "OpenPMX")
		title(paste0("Distribution of RSE ", t), outer=TRUE)
		print(t)
		print(d)
	}
	if (exists("gronmem") && exists("validate")) {
		setup_par()
		t <- "OpenPMX validate"
		d <- performance_loss(gronmem, "OpenPMX", validate, "Validate")
		title(paste0("Distribution of RSE ", t), outer=TRUE)
		print(t)
		print(d)
	}

########################################################################

	objfn_compare <- function(g, gname, baseg, basename)
	{
		par(mfrow=c(2, 3), mar=c(4,4,1,0), oma=c(0,0,3,1), cex=1)

		objfn1 <- baseg[["OBJ"]]
		objfn2 <- g[["OBJ"]]

		lim <- range(c(objfn1, objfn2))
		plot(0, type="n", yaxt="n", xlab=NA, ylab=NA, xlim=(lim), ylim=c(0,1.4), main=NA)
		mtext("Density", side=2, line=0.5, cex=par()$cex)
		mtext("Objective function value", side=1, line=2, cex=par()$cex)
		d1 <- density(objfn1, na.rm=TRUE)
		d2 <- density(objfn2, na.rm=TRUE)
		lines(x=(d2$x), y=d2$y/max(d2$y), col="darkgrey", lw=3)
		lines(x=(d1$x), y=d1$y/max(d1$y), col="black", lw=1)

		x <- (objfn1 + objfn2) / 2
		y <- (objfn1 - objfn2)
		xlim <- range(x)
		ylim <- range(c(y,-y))
		plot(x=x, y=y, type="n", xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, main=NA)
		mtext("Difference", side=2, line=2, cex=par()$cex)
		mtext("Average", side=1, line=2, cex=par()$cex)
		points(x=x, y=y, col="black")
		abline(h=0, col="black")

		x <- objfn2
		y <- objfn1
		xlim <- range(c(x,y))
		ylim <- range(c(x,y))
		plot(x=x, y=y, type="n", xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, main=NA)
#		mtext("Difference", side=2, line=2, cex=par()$cex)
#		mtext("Average", side=1, line=2, cex=par()$cex)
		points(x=x, y=y, col="black")
		abline(a=0, b=1, col="black")

		title(sprintf("Objective function %s vs %s", gname, basename), outer=TRUE)
	}
	if (exists("gronmem") && exists("nonmem")) {
		objfn_compare(gronmem, "OpenPMX", nonmem, "NONMEM")
	}
	if (exists("gronmem") && exists("nonmemlaplace")) {
		objfn_compare(gronmem, "OpenPMX", nonmemlaplace, "LAPLACE")
	}
	if (exists("gronmemicov") && exists("nonmem")) {
		objfn_compare(gronmemicov, "OpenPMX_icov", nonmem, "NONMEM")
	}
	if (exists("gronmemicov") && exists("nonmem")) {
		objfn_compare(gronmemicov, "OpenPMX_icov", nonmem, "NONMEM")
	}
	if (exists("gronmemicov") && exists("gronmem")) {
		objfn_compare(gronmemicov, "OpenPMX_icov", gronmem, "OpenPMX")
	}
	if (exists("gronmemtest") && exists("gronmem")) {
		objfn_compare(gronmemtest, "OpenPMX_icov", gronmem, "OpenPMX")
	}
	if (exists("gronmem") && exists("validate")) {
		objfn_compare(gronmem, "OpenPMX", validate, "Validate")
	}
