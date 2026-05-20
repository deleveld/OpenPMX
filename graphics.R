#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

filename <- args[1]
datafile <- args[2]

# the file names we are looking for
configfile <- paste0(filename, ".graphics")
extfile <- paste0(filename, ".ext")
phifile <- paste0(filename, ".phi")
covfile <- paste0(filename, ".covariates")
yhatfile <- paste0(filename, ".yhat")
pdffile <- paste0(filename, ".pdf")

# read in files

cat(sprintf("writing \"%s\"\n", pdffile))

config <- NA
if (file.exists(configfile)) {
	cat(sprintf("read config \"%s\"\n", configfile))
	config <- read.table(configfile, header=TRUE, na.strings = "")
}

data <- NA
# datafile not given, extract from control file
if (!file.exists(datafile)) {
	file_lines <- readLines(filename, warn = FALSE)
	file_content <- paste(file_lines, collapse = "\n")
	data_match <- regmatches(file_content, regexec("\\$DATA\\(([^)]+)\\)", file_content))
	datafile <- data_match[[1]][2]
	if (!is.na(datafile)) {
		datafile <- gsub("['\"]", "", datafile)
		datafile <- trimws(datafile)
	}
}
if (file.exists(datafile)) {
	cat(sprintf("read data \"%s\"\n", datafile))
	data <- read.table(datafile, header=TRUE, stringsAsFactors=FALSE)
	if (length(names(data)) == 1) {
		data <- read.csv(datafile, header=TRUE, stringsAsFactors=FALSE)
	}
	names(data) <- toupper(names(data))
}
ext <- NA
if (file.exists(extfile)) {
	cat(sprintf("read ext \"%s\"\n", extfile))
	ext <- read.table(extfile, header=TRUE, stringsAsFactors=FALSE)
	names(ext) <- toupper(names(ext))
}
phi <- NA
if (file.exists(phifile)) {
	cat(sprintf("read phi \"%s\"\n", phifile))
	phi <- read.table(phifile, header=TRUE, stringsAsFactors=FALSE)
	names(phi) <- toupper(names(phi))
}
cov <- NA
if (file.exists(covfile)) {
	cat(sprintf("read covariates \"%s\"\n", covfile))
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
	cat(sprintf("read yhat \"%s\"\n", yhatfile))
	yhatdata <- read.table(yhatfile, header=TRUE, stringsAsFactors=FALSE)
	names(yhatdata) <- toupper(names(yhatdata))
	if (is.data.frame(data)) {
		cat(sprintf("bind data and yhat\n"))
		yhatdata <- cbind(data, yhatdata)
	}
}

profiles <- NA
profilelist <- list.files(pattern=paste0(filename, ".profile.*.iter"))
for (profile in profilelist) {
	data <- read.table(profile, header=TRUE, na.strings = "")
	if (is.data.frame(profiles))
		profiles <- rbind(profiles, data)
	else 
		profiles <- data
}

pdf(pdffile)

get_config_setting <- function(name_val, sub_val, sub_type, default)
{
	if (!is.data.frame(config))
		return (default)
	
	d <- config[config[["NAME"]] == name_val & 
				config[["VALUE"]] == sub_val &
				config[["TYPE"]] == sub_type, ]
	if (nrow(d) == 1) 
		return (d$SETTING)
	return (default)
}

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
			abline(l1, col="red", lw=3)
			legend("topright",
				legend=legend,
				bty="n")
		}
		grid()
	}

	plot_categorical_covariate <- function(x, y, ylim, cname, ename)
	{
		dall <- density(y)
		plot(dall,
			 xlim = ylim,
			 ylim = range(c(dall$y, dall$y*1.5)),
			 main = NA, 
			 xlab = ename, 
			 col = NA)
		abline(v=0)

		grp_col <- 2
		legend_text <- NULL
		legend_col <- NULL
		legend_lw <- NULL
		
		for (grp in unique(x)) {
			sel <- grp == x 
			dgrp <- density(y[sel])

			g1 <- y[sel]
			g2 <- y[!sel]
			w <- wilcox.test(g1, g2, alternative="two.sided")

			lw <- 1
			p_val <-w$p.value 
			if (p_val < 0.05) {
				lw <- 3
				grp_name <- get_config_setting(cname, grp, "name", grp)

				legend_text <- c(legend_text, sprintf("%s (p=%.3g)", grp_name, p_val))
				legend_col <- c(legend_col, grp_col)
				legend_lw <- c(legend_lw, lw)
			}

			lines(dgrp, col=grp_col, lw=lw)
			grp_col <- grp_col + 1
		}
		grid()

		if (!is.null(legend_text)) {
			legend("topright",
					legend=legend_text,
					col=legend_col,
					lw=legend_lw,
					bty="n", cex=0.8)
		}
	}

	for (cname in names(cov)) {
		n_plots <- 0
		par(mfrow=c(2,2), cex=0.9)

		cat(sprintf("plot cov %s\n", cname))

		for (ename in names(phi)) {
			x <- cov[[cname]]
			y <- phi[[ename]]

			iseta <- FALSE
			if (length(grep("ETA", ename)) != 0)
				iseta <- TRUE

			if (iseta && !all(y == 0)) {

				val <- get_config_setting(cname, ".", "covariate", "missing")
				if (val == "missing") {
					if (length(unique(x)) <= 5) {
						val = "categorical"
					}
				}
				if (val == "categorical") {
					plot_categorical_covariate(x, y, ylim, cname, ename)
				} else {
					plot_continuous_covariate(x, y, ylim, cname, ename)
				}

				if (n_plots == 0) {
					t <- get_config_setting(cname, ".", "name", NA);
					if (!is.na(t))
						title(t, outer=TRUE, line=-1)
					else
						title(cname, outer=TRUE, line=-1)
				}	
				n_plots <- n_plots + 1
				if (n_plots >= prod(par("mfrow")))
					n_plots <- 0
			}
		}		
	}
}

# yhat vs. yhat var for each dvtype
if (is.data.frame(yhatdata)) {
	cat(sprintf("plot yhatdata\n"))
	pnames <- names(yhatdata)
	
	# do we have a DVTY
	if ("DVTY" %in% pnames)
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
			n_plots <- 0
			par(mfrow=c(2,2), cex=0.9)

			logy <- ""
			logxy <- ""
			log_scale <- ifelse(get_config_setting("DVTY", dname, "scale", "") == "log", TRUE, FALSE)
			if (log_scale) {
				logy <- "y"
				logxy <- "xy"
			}
			
			valid <- yhatdata[["YHATVAR"]] != 0 & preddvid == dname
			if (!all(valid == FALSE)) {
				cat(sprintf("plot yhatdata dvid %i\n", dname))

				id <- yhatdata[valid, "ID"]
				dv <- as.numeric(yhatdata[valid, "DV"])
				yhat <- as.numeric(yhatdata[valid, "YHAT"])
				if (log_scale) {
					dv <- exp(dv)
					yhat <- exp(yhat)
				}
				lim <- range(c(dv, yhat))

				if ("PRED" %in% pnames) {
					ppred <- yhatdata[valid, "PRED"]
					if (log_scale) {
						ppred <- exp(ppred)
					}

					x <- ppred
					y <- dv
					lim <- range(c(lim, ppred))
					plot(x=x, y=y, type="n", xlab="PRED", ylab="DV", xlim=lim, ylim=lim, log=logxy)
					points(x=x, y=y, col=pointcol)
					abline(coef = c(0,1))
					grid()
					
				} else
					frame()

				x <- yhat
				y <- dv
				lim <- range(c(x, y))
				plot(x=x, y=y, type="n", xlab="YHAT", ylab="DV", xlim=lim, ylim=lim, log=logxy)
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
						if (log_scale) {
							ppred <- exp(ppred)
						}
						
						x <- time
						y1 <- dv
						y2 <- ppred
						plot(x=x, y=y1, type="n", xlab="TIME", ylab="DV, PRED", xlim=timelim, ylim=lim, log=logy)
						points(x=x, y=y1, col=pointcol)
						for (i in unique(id)) 
							plot_lines(i, id, x, y2)
						grid()
						
					} else
						frame()

					x <- time
					y1 <- dv
					y2 <- yhat
					plot(x=x, y=y1, type="n", xlab="TIME", ylab="DV, YHAT", xlim=timelim, ylim=lim, log=logy)
					points(x=x, y=y1, col=pointcol)
					for (i in unique(id)) 
						plot_lines(i, id, x, y2)
					grid()
				}
				if (n_plots == 0) {
					t <- get_config_setting("DVTY", dname, "name", NA);
					if (!is.na(t))
						title(t, outer=TRUE, line=-1)
					else
						title(sprintf("data type %i", dname), outer=TRUE, line=-1)
				}
				n_plots <- n_plots + 1
				if (n_plots >= prod(par("mfrow")))
					n_plots <- 0
			}
		}
	}
}

# yhat vs. yhat var for each dvtype
if (is.data.frame(profiles)) {

	par(mfrow=c(2,2), cex=0.9)
	
	p <- profiles
	p$profilegroup <- paste0(p[["type"]], p[["index"]])
	groups <- unique(p$profilegroup)

	for (i in groups) {
		cat(sprintf("plot profile %s\n", i))
		
		idata <- p[p$profilegroup == i, ]
		idata <- idata[idata$objfn != 0., ]
		final <- idata[idata$neval == 0, ]

		ox <- order(idata$value)
		idata <- idata[ox, ]
		
		minobjfn <- min(idata$objfn)
		x0 <- idata[idata$objfn - minobjfn < 20., ]
		xlim <- range(x0$value)

		type <- unique(idata[["type"]])
		idx <- unique(idata[["index"]])
		plotname <- paste0(type, "(", idx, ")")

		x <- idata[["value"]]
		y <- idata[["objfn"]]
		plot(x=x, y=y, type="n", xlab=plotname, ylab="Objfn", xlim=xlim, ylim=c(min(y), min(y)+20))

		xl <- seq(from=min(x), to=max(x), length.out=100)
#		sfunc <- splinefun(x=x, y=y, 
#						   method="natural",
#						   ties=mean)
#		lines(x=xl, y=sfunc(xl), col = "blue")

		sobj <- loess(y ~ x, span=0.7)
		pred <- predict(sobj, data.frame(x=xl), se=TRUE)
		upper <- pred$fit + 1.96 * pred$se.fit
		lower <- pred$fit - 1.96 * pred$se.fit
		polygon(c(xl, rev(xl)), c(upper, rev(lower)), col=rgb(1,0,0,0.15), border=NA)
		
		lines(xl, pred$fit, col = "blue")

		points(x=x, y=y)
		
		points(x=final[["value"]], y=final[["objfn"]], col="red", pch=16)
	}
}
dummy <- dev.off()
cat("done\n")

