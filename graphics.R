#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

filename <- args[1]
datafile <- args[2]

# the file names we are looking for
extfile <- paste0(filename, ".ext")
phifile <- paste0(filename, ".phi")
covfile <- paste0(filename, ".covariates")
yhatfile <- paste0(filename, ".yhat")
pdffile <- paste0(filename, ".pdf")

# read in files
data <- NA
if (file.exists(datafile)) {
	cat(sprintf("read data %s\n", datafile))
	data <- read.table(datafile, header=TRUE, stringsAsFactors=FALSE)
	if (length(names(data)) == 1) {
		data <- read.csv(datafile, header=TRUE, stringsAsFactors=FALSE)
	}
	names(data) <- toupper(names(data))
}
ext <- NA
if (file.exists(extfile)) {
	cat(sprintf("read ext %s\n", extfile))
	ext <- read.table(extfile, header=TRUE, stringsAsFactors=FALSE)
	names(ext) <- toupper(names(ext))
}
phi <- NA
if (file.exists(phifile)) {
	cat(sprintf("read phi %s\n", phifile))
	phi <- read.table(phifile, header=TRUE, stringsAsFactors=FALSE)
	names(phi) <- toupper(names(phi))
}
cov <- NA
if (file.exists(covfile)) {
	cat(sprintf("read covariates %s\n", covfile))
	cov <- read.table(covfile, header=TRUE, stringsAsFactors=FALSE)
	names(cov) <- toupper(names(cov))
	newid <- c(TRUE, diff(cov[["ID"]]) != 0)
	cov <- cov[newid, ]
	
# if no covariates file, try to construct from data
} else if (is.data.frame(data)) {
	cat(sprintf("make covariates from data\n"))
	cov <- data
	names(cov) <- toupper(names(cov))
	
	# covariates are constant for each individual
	allnames <- names(cov)
	validnames <- c()
	for (j in allnames) {
		covlike <- TRUE
		allids <- unique(cov[["ID"]])
		for (i in allids) {
			v <- cov[cov[["ID"]] == i, j]
			if (length(unique(v)) != 1) {
				covlike <- FALSE
			}
		}
		if (covlike) 
			validnames <- c(validnames, j)
	}

	# covariates have more that one value
	valid2names <- c()
	for (j in validnames) {
		v <- cov[[j]]
		if (length(unique(v)) != 1) {
			valid2names <- c(valid2names, j)
		}
	}
	cov <- cov[ , names(cov) %in% valid2names] 
	newid <- c(TRUE, diff(cov[["ID"]]) != 0)
	cov <- cov[newid, ]
}

yhatdata <- NA
if (file.exists(yhatfile)) {
	cat(sprintf("read yhat %s\n", yhatfile))
	yhatdata <- read.table(yhatfile, header=TRUE, stringsAsFactors=FALSE)
	names(yhatdata) <- toupper(names(yhatdata))
	if (is.data.frame(data)) {
		cat(sprintf("bind data and yhat\n"))
		yhatdata <- cbind(data, yhatdata)
	}
}
pdf(pdffile)

# parameter vs. iteration
if (is.data.frame(ext)) {
	cat(sprintf("plot ext\n"))
	par(mfrow=c(2,2), cex=0.9)
	for (iname in names(ext)) {
		
		v <- ext[["ITERATION"]]
		valid <- ext[["ITERATION"]] > 0

		x <- ext[valid, "ITERATION"]
		y <- ext[valid, iname]

		if (iname != "ITERATION" && !all(y == 0)) {
			plot(x=x, y=y, type="n", xlab="ITERATION", ylab=iname)
			lines(x=x, y=y, col="red")
			grid()
		}
	}
}

# covariate vs. eta
if (is.data.frame(cov) && is.data.frame(phi)) {
	cat(sprintf("plot cov vs. phi\n"))

	# shared eta range
	ylim <- 0
	for (ename in names(phi)) {
		y <- phi[[ename]]

		iseta <- FALSE
		if (length(grep("ETA", ename)) != 0)
			iseta <- TRUE

		if (iseta && !all(y == 0)) 
			ylim <- range(c(y, -y, ylim))
	}

	# eta plots
	plot_continuous_covariate <- function(x, y, ylim, cname, ename)
	{
		plot(x=x, y=y, type="n", ylim=ylim, xlab=cname, ylab=ename)
		points(x=x, y=y, col="grey")
		abline(h=0)
		grid()

		l1 <- lm(y ~ x)
		coef <- summary(l1)$coefficients
		pval_slope <- 1
		if (NROW(coef) >= 2 && NCOL(coef) >= 4)
			pval_slope <- coef[2,4] # the p-value for the slope
		legend <- ""
		p_value_threshold <- 0.05
		if (!is.nan(pval_slope) && pval_slope < p_value_threshold) {
			legend <- sprintf("p=%.4g", pval_slope)
		}
		if (legend != "") {
			abline(l1, col="red")
			legend("topright",
				legend=legend,
				bty="n")
		}
	}

	for (cname in names(cov)) {
		par(mfrow=c(2,2), cex=0.9)

		for (ename in names(phi)) {
			x <- cov[[cname]]
			y <- phi[[ename]]

			iseta <- FALSE
			if (length(grep("ETA", ename)) != 0)
				iseta <- TRUE

			if (iseta && !all(y == 0)) {
				plot_continuous_covariate(x, y, ylim, cname, ename)
			}
		}
	}
}

# yhat vs. yhat var for each dvtype
if (is.data.frame(yhatdata)) {
	cat(sprintf("plot yhatdata\n"))
	pnames <- names(yhatdata)
	
	# do we have a dvid
	if ("DVID" %in% pnames)
		preddvid <- yhatdata[["DVID"]]
	else if ("DVTY" %in% pnames)
		preddvid <- yhatdata[["DVTY"]]
	else 
		preddvid <- rep(1, nrow(yhatdata))
	dvid_types <- unique(preddvid)

	pointcol = "grey"

	if ("DV" %in% pnames &&
		"YHAT" %in% pnames) {

		timevalid <- yhatdata[["YHATVAR"]] != 0
		timelim <- range(yhatdata[timevalid, "TIME"])

		for (dname in dvid_types) {
			par(mfrow=c(2,2), cex=0.9)

			cat(sprintf("plot yhatdata dvid %i\n", dname))
			valid <- yhatdata[["YHATVAR"]] != 0 & preddvid == dname
			if (!all(valid == FALSE)) {
				id <- yhatdata[valid, "ID"]
				dv <- as.numeric(yhatdata[valid, "DV"])
				yhat <- as.numeric(yhatdata[valid, "YHAT"])

				lim <- range(c(dv, yhat))

				if ("PRED" %in% pnames) {
					ppred <- yhatdata[valid, "PRED"]

					x <- ppred
					y <- dv
					lim <- range(c(lim, ppred))
					plot(x=x, y=y, type="n", xlab="PRED", ylab="DV", xlim=lim, ylim=lim)
					points(x=x, y=y, col=pointcol)
					abline(coef = c(0,1))
					grid()
					
				} else
					frame()

				x <- yhat
				y <- dv
				lim <- range(c(x, y))
				plot(x=x, y=y, type="n", xlab="YHAT", ylab="DV", xlim=lim, ylim=lim)
				points(x=x, y=y, col=pointcol)
				abline(coef = c(0,1))
				grid()

				if ("TIME" %in% pnames) {
					time <- yhatdata[valid, "TIME"]

					plot_lines <- function(i, id, time, y)
					{
						thisid <- id == i
						ix <- time[thisid]
						iy <- y[thisid]
						# handle mutiple sessions (time step < 0) within an individual
						# drawn separate line if time jumps backwards
						sess <- cumsum(c(FALSE, diff(ix) < 0))
						for (s in unique(sess)) {
							six <- ix[sess == s]
							siy <- iy[sess == s]
							lines(x=six, y=siy, col="red")
						}
					}

					if ("PRED" %in% pnames) {
						ppred <- yhatdata[valid, "PRED"]
						
						x <- time
						y1 <- dv
						y2 <- ppred
						plot(x=x, y=y1, type="n", xlab="TIME", ylab="DV, PRED", xlim=timelim, ylim=lim)
						points(x=x, y=y1, col=pointcol)
						for (i in unique(id)) 
							plot_lines(i, id, x, y2)
						grid()
						
					} else
						frame()

					x <- time
					y1 <- dv
					y2 <- yhat
					plot(x=x, y=y1, type="n", xlab="TIME", ylab="DV, YHAT", xlim=timelim, ylim=lim)
					points(x=x, y=y1, col=pointcol)
					for (i in unique(id)) 
						plot_lines(i, id, x, y2)
					grid()
				}
				title(sprintf("data type %i", dname), outer=TRUE, line=-1)
			}
		}
	}
}



dummy <- dev.off()
