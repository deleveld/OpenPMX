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

			plot(0, type="n", yaxt="n", xlab=NA, ylab=NA, xlim=c(0,max(xlim)), ylim=c(0,1.4), main=NA, xaxs="i", yaxs="i")
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
			d1x <- d1$x
			d1y <- d1$y/max(d1$y)
			d2x <- d2$x
			d2y <- d2$y/max(d2$y)
			lines(x=d2x, y=d2y, col="darkgrey", lw=3)
			lines(x=d1x, y=d1y, col="black", lw=1)

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
###		title(paste0("Distribution of RSE ", t), outer=TRUE)
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
		objfn1 <- baseg[["OBJ"]]
		objfn2 <- g[["OBJ"]]

		lim <- range(c(objfn1, objfn2))
		plot(0, type="n", yaxt="n", xlab=NA, ylab=NA, xlim=(lim), ylim=c(0,1.4), main=NA)
		mtext("Density", side=2, line=0.5, cex=par()$cex)
		mtext("Objective function value", side=1, line=2, cex=par()$cex)
		grid()
		d1 <- density(objfn1, na.rm=TRUE)
		d2 <- density(objfn2, na.rm=TRUE)
		lines(x=(d2$x), y=d2$y/max(d2$y), col="darkgrey", lw=3)
		lines(x=(d1$x), y=d1$y/max(d1$y), col="black", lw=1)
		legend("topright",
			legend=c(basename, gname),
			col=c("darkgrey", "black"),
			lw=c(3,1),
			bty="n", cex=par()$cex*0.85)

		x <- (objfn1 + objfn2) / 2
		y <- (objfn1 - objfn2)
		xlim <- range(x)
		ylim <- range(c(y,-y))
		plot(x=x, y=y, type="n", xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, main=NA)
		mtext("Objective function difference", side=2, line=2, cex=par()$cex)
		mtext("Objective function average", side=1, line=2, cex=par()$cex)
		grid()
		points(x=x, y=y, col="black")
		abline(h=0, col="black")

		frame()

		x <- objfn2
		y <- objfn1
		xlim <- range(c(x,y))
		ylim <- range(c(x,y))
		plot(x=x, y=y, type="n", xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, main=NA)
		mtext(sprintf("Objective function %s", gname), side=2, line=2, cex=par()$cex)
		mtext(sprintf("Objective function %s", basename), side=1, line=2, cex=par()$cex)
		grid()
		points(x=x, y=y, col="black")
		abline(a=0, b=1, col="black")
	}

	offx <- 0
	if (exists("median_speed_offset_x"))
		offx <- median_speed_offset_x
	offy <- 0
	if (exists("median_speed_offset_y"))
		offy <- median_speed_offset_y
	
	runtime_compare <- function(g, gname, baseg, basename)
	{
		runtime1 <- g[["V2"]]/1000
		runtime2 <- baseg[["V2"]]/1000

		med1 <- median(runtime1)
		med2 <- median(runtime2)

		lim <- range(c(runtime1, runtime2))
		plot(0, type="n", yaxt="n", xlab=NA, ylab=NA, xlim=(lim), ylim=c(0,1.4), main=NA)
		mtext("Density", side=2, line=0.5, cex=par()$cex)
		mtext("Runtime (s)", side=1, line=2, cex=par()$cex)
		grid()
		d1 <- density(runtime1, na.rm=TRUE)
		d2 <- density(runtime2, na.rm=TRUE)
		lines(x=(d2$x), y=d2$y/max(d2$y), col="darkgrey", lw=3)
		lines(x=(d1$x), y=d1$y/max(d1$y), col="black", lw=1)
		abline(v=med1, lt=2, col="black")
		abline(v=med2, lt=2, col="darkgrey")
		text(x=med1+offx, y=1.05+offy, labels=sprintf(" %.1f s", med1), adj=c(0,0))
		text(x=med2, y=1.05, labels=sprintf(" %.1f s", med2), adj=c(0,0), col="darkgrey")
		legend("topright",
			legend=c(basename, gname),
			col=c("darkgrey", "black"),
			lw=c(3,1),
			bty="n", cex=par()$cex*0.85)
	}

	mfrow <- c(2,3)
	if (exists("gronmem") && exists("nonmem")) {
		par(mfrow=mfrow, mar=c(4,4,1,0), oma=c(0,0,3,1), cex=1)
		objfn_compare(gronmem, "OpenPMX", nonmem, "NONMEM")
		runtime_compare(gronmem_ms, "OpenPMX", nonmem_ms, "NONMEM")
	}
	if (exists("validate") && exists("gronmem")) {
		par(mfrow=mfrow, mar=c(4,4,1,0), oma=c(0,0,3,1), cex=1)
		objfn_compare(validate, "Validate", gronmem, "OpenPMX")
		runtime_compare(validate_ms, "Validate", gronmem_ms, "OpenPMX")
	}

