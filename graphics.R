#!/usr/bin/env Rscript
# Diagnostic plotting for a population PK/PD run (NONMEM-style output files).
# Usage: graphics.R <controlfile/runname> [datafile]
# Design: every input file is optional. Each is read into a data.frame if it
# exists (NA otherwise), and each plot section is guarded by is.data.frame()
# checks, so the script produces whatever plots the available files allow.
# All plots go into a single multi-page A4-landscape PDF.
args = commandArgs(trailingOnly=TRUE)

filename <- args[1]   # run name / control file; all output files share this stem
datafile <- args[2]   # optional: dataset (may be missing or NA)

# the file names we are looking for: all are <stem>.<suffix>
configfile <- paste0(filename, ".graphics")       # plot settings (names, scales, covariate types)
extfile <- paste0(filename, ".ext")               # parameter estimates per iteration
phifile <- paste0(filename, ".phi")               # individual (empirical Bayes) ETAs
covfile <- paste0(filename, ".covariates")        # covariates per record
yhatfile <- paste0(filename, ".yhat")             # individual predictions
npdefile <- paste0(filename, ".npde")             # npde
profilefile <- paste0(filename, ".profile.iter")  # likelihood profile iterations
covlimitfile <- paste0(filename, ".cov.limits")   # covariance-step confidence limits
vpcoutfile <- paste0(filename, ".vpc.out")        # VPC percentile summaries
vpcpointsfile <- paste0(filename, ".vpc.points")  # VPC individual points
pdffile <- paste0(filename, ".pdf")               # output

# read in files
cat(sprintf("output \"%s\"\n", pdffile))

# Optional config: table with columns NAME, VALUE, TYPE, SETTING
# (empty strings read as NA so missing settings are real NAs)
config <- NA
if (file.exists(configfile)) {
	cat(sprintf("read config \"%s\"\n", configfile))
	config <- read.table(configfile, header=TRUE, na.strings = "")
}

data <- NA
# datafile not given, extract from control file
# (file.exists(NA) is FALSE, so a missing argument also lands here)
if (!file.exists(datafile)) {
	# slurp the control file into one string so the regex can span lines
	file_lines <- readLines(filename, warn = FALSE)
	file_content <- paste(file_lines, collapse = "\n")
	# look for a $DATA(<name>) record and capture what is inside the parentheses
	data_match <- regmatches(file_content, regexec("\\$DATA\\(([^)]+)\\)", file_content))
	datafile <- data_match[[1]][2]   # [1] is the full match, [2] the capture group (NA if no match)
	if (!is.na(datafile)) {
		datafile <- gsub("['\"]", "", datafile)  # strip quotes around the filename
		datafile <- trimws(datafile)
	}
}
if (file.exists(datafile)) {
	cat(sprintf("read data \"%s\"\n", datafile))
	data <- read.table(datafile, header=TRUE, stringsAsFactors=FALSE)
	# format sniffing: if whitespace splitting gave one column, the file is comma separated
	if (length(names(data)) == 1) {
		data <- read.csv(datafile, header=TRUE, stringsAsFactors=FALSE)
	}
	# upper-case all column names so later lookups (ID, TIME, DV...) are case-insensitive
	names(data) <- toupper(names(data))
}
# Parameter estimates by iteration (NONMEM .ext layout)
ext <- NA
if (file.exists(extfile)) {
	cat(sprintf("read ext \"%s\"\n", extfile))
	ext <- read.table(extfile, header=TRUE, stringsAsFactors=FALSE)
	names(ext) <- toupper(names(ext))
}
# Individual ETAs (NONMEM .phi layout), one row per individual
phi <- NA
if (file.exists(phifile)) {
	cat(sprintf("read phi \"%s\"\n", phifile))
	phi <- read.table(phifile, header=TRUE, stringsAsFactors=FALSE)
	names(phi) <- toupper(names(phi))
}
# Covariates: reduced to ONE row per individual so rows align with phi
cov <- NA
if (file.exists(covfile)) {
	cat(sprintf("read covariates \"%s\"\n", covfile))
	cov <- read.table(covfile, header=TRUE, stringsAsFactors=FALSE)
	names(cov) <- toupper(names(cov))
	# keep the first row of each run of identical IDs (assumes data sorted by ID)
	newid <- c(TRUE, diff(cov[["ID"]]) != 0)
	cov <- cov[newid, ]
	
# if no covariates file, try to construct from data
} else if (is.data.frame(data)) {
	cat(sprintf("make covariates from data\n"))
	cov <- data
	names(cov) <- toupper(names(cov))
	
	# covariates are constant for each individual
	# i.e. keep only columns that never change within an ID
	# (this drops TIME, DV, AMT, etc. automatically)
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
	# i.e. drop columns that are constant across the whole data set (nothing to plot)
	valid2names <- c()
	for (j in validnames) {
		v <- cov[[j]]
		if (length(unique(v)) != 1) {
			valid2names <- c(valid2names, j)
		}
	}
	cov <- cov[ , names(cov) %in% valid2names] 
	# collapse to one row per individual
	newid <- c(TRUE, diff(cov[["ID"]]) != 0)
	cov <- cov[newid, ]
}

# Individual predictions; if the data is available, columns are bound onto it
# so the plots can use covariates/DV/TIME alongside YHAT (assumes same row order)
yhatdata <- NA
if (file.exists(yhatfile)) {
	cat(sprintf("read yhat \"%s\"\n", yhatfile))
	yhatdata <- read.table(yhatfile, header=TRUE, stringsAsFactors=FALSE)
	names(yhatdata) <- toupper(names(yhatdata))

#	if (is.data.frame(data)) {
#		cat(sprintf("bind data and yhat\n"))
#		yhatdata <- cbind(data, yhatdata)
#	}
	# add some coloums from the data
	if (is.data.frame(data)) {
		yhatdata[["ID"]] <- data[["ID"]]
		yhatdata[["TIME"]] <- data[["TIME"]]
		dnames <- names(data)
		if ("DVTY" %in% dnames) {
			yhatdata[["DVTY"]] <- data[["DVTY"]]
		} else if ("DVID" %in% dnames) {
			yhatdata[["DVID"]] <- data[["DVID"]]
		}
	}
}

npdedata <- NA
if (file.exists(npdefile)) {
	cat(sprintf("read npde \"%s\"\n", npdefile))
	npdedata <- read.table(npdefile, header=TRUE, stringsAsFactors=FALSE)
	names(npdedata) <- toupper(names(npdedata))

#	if (is.data.frame(data)) {
#		cat(sprintf("bind data and npde\n"))
#		npdedata <- cbind(data, npdedata)
#	}
	# add some coloums from the data
	if (is.data.frame(data)) {
		npdedata[["ID"]] <- data[["ID"]]
		npdedata[["TIME"]] <- data[["TIME"]]
		dnames <- names(data)
		if ("DVTY" %in% dnames) {
			npdedata[["DVTY"]] <- data[["DVTY"]]
		} else if ("DVID" %in% dnames) {
			npdedata[["DVID"]] <- data[["DVID"]]
		}
	}
	if (is.data.frame(yhatdata)) {
		npdedata[["YHAT"]] <- yhatdata[["YHAT"]]
	}
}

# Likelihood-profile iterations (one block per parameter)
profiles <- NA
if (file.exists(profilefile)) {
	cat(sprintf("read profile \"%s\"\n", profilefile))
	profiles <- read.table(profilefile, header=TRUE, stringsAsFactors=FALSE)
	names(profiles) <- toupper(names(profiles))
}

# Confidence limits from the covariance step, overlaid on the profiles
covlimits <- NA
if (file.exists(covlimitfile)) {
	cat(sprintf("read cov limits \"%s\"\n", covlimitfile))
	covlimits <- read.table(covlimitfile, header=TRUE, stringsAsFactors=FALSE)
	names(covlimits) <- toupper(names(covlimits))
}

# VPC: percentile summaries (vpcout) and individual points (vpcpoints)
vpcout <- NA
if (file.exists(vpcoutfile)) {
	cat(sprintf("read vpc out \"%s\"\n", vpcoutfile))
	vpcout <- read.table(vpcoutfile, header=TRUE, stringsAsFactors=FALSE)
	names(vpcout) <- toupper(names(vpcout))
}
vpcpoints <- NA
if (file.exists(vpcpointsfile)) {
	cat(sprintf("read vpc points \"%s\"\n", vpcpointsfile))
	vpcpoints <- read.table(vpcpointsfile, header=TRUE, stringsAsFactors=FALSE)
	names(vpcpoints) <- toupper(names(vpcpoints))
}

# do some checking of files
chk <- function(label, x, need=NULL) {
	if (!is.data.frame(x)) {
		cat(sprintf("%-11s MISSING\n", label));
		return(invisible())
	}
	miss <- setdiff(need, names(x))
	cat(sprintf("%-11s %d rows, %d cols%s\n", label, NROW(x), NCOL(x),
		if (length(miss)) paste0("  -- missing: ", paste(miss, collapse=",")) else ""))
}
chk("data",       data,      c("ID", "TIME"))
chk("ext",        ext,       "ITERATION")
chk("phi",        phi,       "ID")
chk("covariates", cov,       "ID")
chk("yhat",       yhatdata,  c("ID", "TIME", "YHATDV", "YHAT", "YHATVAR"))
chk("npde",       npdedata,  c("NPDE"))
chk("profiles",   profiles,  c("TYPE", "INDEX", "VALUE", "OBJFN", "NEVAL", "TARGET"))
chk("vpc.out",    vpcout,    c("NAME", "BINMIN", "BINMAX"))

pdf(pdffile, width=11.69, height=8.27)  # A4 landscape

# Look up a user setting in the config table by (NAME, VALUE, TYPE).
# Returns SETTING if exactly one row matches, else the default.
# Lets the user rename things, mark DVTYs as log-scale, or force covariate type
# without editing this script. Note: reads the global `config`.
get_config_setting <- function(name_val, sub_val, sub_type, default)
{
	if (!is.data.frame(config))
		return (default)
	
	# which() drops NA comparisons (rows with empty NAME/VALUE/TYPE), which would
	# otherwise index as all-NA rows and be returned as a "match"
	d <- config[which(config[["NAME"]] == name_val & 
				config[["VALUE"]] == sub_val &
				config[["TYPE"]] == sub_type), ]
	if (nrow(d) == 1 && !is.na(d$SETTING)) 
		return (d$SETTING)
	return (default)
}

# Scale (transformation) shared by all figures.
# `setting` is the SETTING of a config row with TYPE "scale":
#   "identity"      -> data used as is (default, also for unknown values)
#   "log" or "exp"  -> data are on the log scale: itr = exp maps them back to the
#                      natural scale, tr = log maps natural values to the log scale,
#                      and the axes are drawn logarithmic
# Returns tr/itr plus the `log=` strings for plot() so callers need no own logic.
make_scale <- function(setting)
{
	if (setting %in% c("log", "exp"))
		return (list(tr=log, itr=exp, islog=TRUE, xlog="x", ylog="y", xylog="xy"))
	if (!setting %in% c("identity", ""))
		cat(sprintf("... unknown scale \"%s\", using identity\n", setting))
	list(tr=identity, itr=identity, islog=FALSE, xlog="", ylog="", xylog="")
}

# Scale from the config: row NAME=name_val, VALUE=sub_val, TYPE=scale
get_scale <- function(name_val, sub_val, default="identity")
	make_scale(get_config_setting(name_val, sub_val, "scale", default))

# parameter vs. iteration
# One convergence trace per parameter, 6 per page (2x3).
paramater_versus_iteration <- TRUE
cat("check parameter vs. iteration\n")
if (!is.data.frame(ext) || NROW(ext) == 0) {
	cat("... no ext file\n")
	paramater_versus_iteration <- FALSE
}
if (paramater_versus_iteration) {
	cat("plot parameter vs. iteration\n")
	par(mfrow=c(2,3), cex=1)
	for (iname in names(ext)) {
		
		# iteration <= 0 rows are NONMEM's special rows (final estimates, SE, etc.), not the trace
		valid <- ext[["ITERATION"]] > 0

		x <- ext[valid, "ITERATION"]
		y <- ext[valid, iname]

		# skip the x column itself and parameters that are fixed at zero (e.g. unused OMEGA elements)
		if (iname != "ITERATION" && !all(y == 0)) {
			plot(x=x, y=y, type="n", xlab="ITERATION", ylab=iname)
			lines(x=x, y=y, col="red")
			grid()
		}
	}
}

# covariate vs. eta
# For every covariate, one panel per ETA to look for unexplained covariate relationships.
cat("check cov vs. eta\n")
covariate_versus_eta_plot <- TRUE
if (!is.data.frame(cov) || NROW(cov) == 0) {
	cat("... no cov data\n")
	covariate_versus_eta_plot <- FALSE
}
if (!is.data.frame(phi) || NROW(phi) == 0) {
	cat("... no phi file\n")
	covariate_versus_eta_plot <- FALSE
}
if (covariate_versus_eta_plot) {
	cat("plot cov vs. eta\n")

	# shared eta range
	# Symmetric range around 0 across ALL etas, so panels are visually comparable
	# (range(c(y, -y, ylim)) forces symmetry)
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
	# Continuous covariate: scatter of ETA vs covariate with a linear fit that is
	# only drawn (and the p-value only shown) when the slope is significant at 5%.
	plot_continuous_covariate <- function(x, y, ylim, cname, ename)
	{
		plot(x=x, y=y, type="n", ylim=ylim, xlab=cname, ylab=ename)
		points(x=x, y=y, col="grey")
		abline(h=0)   # reference line at eta = 0

		l1 <- lm(y ~ x)
		coef <- summary(l1)$coefficients
		pval_slope <- 1
		# guard: with degenerate x the slope row may be missing
		if (NROW(coef) >= 2 && NCOL(coef) >= 4)
			pval_slope <- coef[2,4] # the p-value for the slope
		legend <- ""
		p_value_threshold <- 0.05
		if (!is.nan(pval_slope) && pval_slope < p_value_threshold) {
			legend <- sprintf("p=%.4g", pval_slope)
		}
		# only annotate significant relationships, keeping the plots uncluttered
		if (legend != "") {
			abline(l1, col="red", lw=3)
			legend("topright",
				legend=legend,
				bty="n")
		}
		grid()
	}

	# Categorical covariate: one ETA density per group. A group is highlighted
	# (thick line + legend entry with p-value) when a Wilcoxon rank-sum test says
	# its ETAs differ from all other groups pooled (one-vs-rest).
	plot_categorical_covariate <- function(x, y, ylim, cname, ename)
	{
		dall <- density(y)
		# empty frame (col=NA) just to get axes; extra headroom in y leaves room for the legend
		plot(dall,
			 xlim = ylim,
			 ylim = range(c(dall$y, dall$y*1.5)),
			 main = NA, 
			 xlab = ename, 
			 col = NA)
		abline(v=0)

		grp_col <- 2   # start at col 2 (red); 1 is black
		legend_text <- NULL
		legend_col <- NULL
		legend_lw <- NULL
		
		for (grp in unique(x)) {
			sel <- grp == x 
			dgrp <- density(y[sel])

			# one-vs-rest nonparametric test (no normality assumption on ETAs)
			g1 <- y[sel]
			g2 <- y[!sel]
			w <- wilcox.test(g1, g2, alternative="two.sided")

			lw <- 1
			p_val <-w$p.value 
			if (p_val < 0.05) {
				lw <- 3
				# optional human-readable group name from the config file
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
					bty="n", cex=1)
		}
	}

	for (cname in names(cov)) {
		n_plots <- 0   # panels drawn on the current page; used to title only the first panel of a page
		par(mfrow=c(2,3), cex=1)

		cat(sprintf("plot cov %s\n", cname))

		for (ename in names(phi)) {
			x <- cov[[cname]]
			y <- phi[[ename]]
			
			iseta <- FALSE
			if (length(grep("ETA", ename)) != 0)
				iseta <- TRUE

			if (iseta && !all(y == 0)) {

				# covariate type: user setting wins; otherwise auto-detect
				# (<= 5 distinct values => categorical)
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

				# outer page title (covariate name) written on the first panel of each page
				if (n_plots == 0) {
					t <- get_config_setting(cname, ".", "name", NA);
					if (!is.na(t))
						title(t, outer=TRUE, line=-1)
					else
						title(cname, outer=TRUE, line=-1)
				}	
				n_plots <- n_plots + 1
				# wrap the counter when a page (2x3 = 6 panels) is full
				if (n_plots >= prod(par("mfrow")))
					n_plots <- 0
			}
		}		
	}

	# eta vs. eta plot
	phi2 <- phi[, grepl("ETA", names(phi)), drop=FALSE]
	pairs(phi2)
}

# yhat vs. yhat var for each dvtype
# Goodness of fit: DV vs PRED, DV vs YHAT, and DV/PRED/YHAT vs TIME, per data type (DVTY).
cat("check yhat\n")
yhat_diagnostic_plots <- TRUE
if (!is.data.frame(yhatdata) || NROW(yhatdata) == 0) {
	cat("... no yhat data\n")
	yhat_diagnostic_plots <- FALSE
}
if (yhat_diagnostic_plots) {
	cat(sprintf("plot yhat diagnostics\n"))
	pnames <- names(yhatdata)

	# do we have a DVTY
	# observation type column: DVTY, else DVID, else a single type
	if ("DVTY" %in% pnames)
		preddvid <- yhatdata[["DVTY"]]
	else if ("DVID" %in% pnames)
		preddvid <- yhatdata[["DVID"]]
	else 
		preddvid <- rep(1, nrow(yhatdata))
	dvid_types <- unique(na.omit(preddvid))

	pointcol = "grey"

	if ("YHATDV" %in% pnames &&
		"YHAT" %in% pnames) {

		# YHATVAR == 0 marks rows that are not real observations (dose rows etc.)
		timevalid <- yhatdata[["YHATVAR"]] != 0
		timelim <- range(yhatdata[timevalid, "TIME"])  # common time axis for all types
		
		for (dname in dvid_types) {
			n_plots <- 0
			par(mfrow=c(2,3), cex=1)

			# scale of this data type (itr maps the data to the natural scale)
			sc <- get_scale("DVTY", dname)
			itr <- sc$itr
			logy <- sc$ylog
			logxy <- sc$xylog

			valid <- yhatdata[["YHATVAR"]] != 0 & preddvid == dname
			if (!all(valid == FALSE)) {
				cat(sprintf("plot yhatdata dvid %i\n", dname))

				id <- yhatdata[valid, "ID"]
				dv <- itr(as.numeric(yhatdata[valid, "YHATDV"]))
				yhat <- itr(as.numeric(yhatdata[valid, "YHAT"]))
				lim <- range(c(dv, yhat), na.rm=TRUE)
	
				# panel 1: DV vs PRED (population prediction), identity line
				if ("PRED" %in% pnames) {
					ppred <- itr(yhatdata[valid, "PRED"])

					x <- ppred
					y <- dv
					lim <- range(c(lim, ppred))   # same limits on both axes so the identity line is at 45 degrees
					plot(x=x, y=y, type="n", xlab="PRED", ylab="DV", xlim=lim, ylim=lim, log=logxy)
					points(x=x, y=y, col=pointcol)
					abline(coef = c(0,1))
					grid()
					
				} else
					frame()   # keep the 2x3 layout aligned by leaving the slot blank

				# panel 2: DV vs YHAT (individual prediction)
				x <- yhat
				y <- dv
				plot(x=x, y=y, type="n", xlab="YHAT", ylab="DV", xlim=lim, ylim=lim)
				points(x=x, y=y, col=pointcol)
				abline(coef = c(0,1))
				grid()

				# panel 3: CWRES vs YHAT (individual prediction)
				if ("CWRES" %in% pnames) {
					cwres <- itr(yhatdata[valid, "CWRES"])

					x <- yhat
					y <- cwres
					ylim <- range(c(cwres, -cwres), na.rm=TRUE)

					plot(x=x, y=y, type="n", xlab="YHAT", ylab="CWRES", xlim=lim, ylim=ylim)
					points(x=x, y=y, col=pointcol)
					abline(coef = c(0,0))
					grid()
					
				} else 
					frame() 

				if ("TIME" %in% pnames) {
					time <- yhatdata[valid, "TIME"]

					# Draw one individual's predicted curve. If time goes backwards
					# the individual has several sessions, so start a new line
					# rather than connecting across the jump.
					plot_lines <- function(i, id, time, y)
					{
						thisid <- id == i
						ix <- time[thisid]
						iy <- y[thisid]
						# handle mutiple sessions (time step < 0) within an individual
						# drawn separate line if time jumps backwards
						# cumsum of a "went backwards" flag gives a session counter
						sess <- cumsum(c(FALSE, diff(ix) < 0))
						for (s in unique(sess)) {
							six <- ix[sess == s]
							siy <- iy[sess == s]
							lines(x=six, y=siy, col="red")
						}
					}

					# panel 4: DV and PRED vs time
					if ("PRED" %in% pnames) {
						ppred <- itr(yhatdata[valid, "PRED"])
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

					# panel 5: DV and YHAT vs time
					x <- time
					y1 <- dv
					y2 <- yhat
					plot(x=x, y=y1, type="n", xlab="TIME", ylab="DV, YHAT", xlim=timelim, ylim=lim, log=logy)
					points(x=x, y=y1, col=pointcol)
					for (i in unique(id)) 
						plot_lines(i, id, x, y2)
					grid()
					
					# panel 5: CWRES vs time
					if ("CWRES" %in% pnames) {
						cwres <- itr(yhatdata[valid, "CWRES"])

						x <- time
						y <- cwres
						ylim <- range(c(cwres, -cwres), na.rm=TRUE)
						plot(x=x, y=y, type="n", xlab="TIME", ylab="CWRES", xlim=timelim, ylim=ylim)
						points(x=x, y=y, col=pointcol)
						abline(coef = c(0,0))
						grid()
						
					} else 
						frame() 
					
				}
				# page title for this data type: configured name or generic label
				# (note: %i assumes numeric DVTY codes)
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

cat("check npde\n")
npde_diagnostic_plots <- TRUE
if (!is.data.frame(npdedata) || NROW(npdedata) == 0) {
	cat("... no npde data\n")
	npde_diagnostic_plots <- FALSE
}
if (npde_diagnostic_plots) {
	cat(sprintf("plot npde diagnostics\n"))
	pnames <- names(npdedata)

	# do we have a DVTY
	# observation type column: DVTY, else DVID, else a single type
	if ("DVTY" %in% pnames)
		preddvid <- npdedata[["DVTY"]]
	else if ("DVID" %in% pnames)
		preddvid <- npdedata[["DVID"]]
	else 
		preddvid <- rep(1, nrow(npdedata))
	dvid_types <- unique(na.omit(preddvid))

	pointcol = "grey"

	if ("NPDE" %in% pnames) {

		# NPDE == NA marks rows that are not real observations (dose rows etc.)
		timevalid <- npdedata[["NPDE"]] != 0
		timelim <- range(npdedata[timevalid, "TIME"], na.rm=TRUE)  # common time axis for all types
		
		for (dname in dvid_types) {
			n_plots <- 0
			par(mfrow=c(2,3), cex=1)

			# scale of this data type (itr maps the data to the natural scale)
			sc <- get_scale("DVTY", dname)
			itr <- sc$itr
			logy <- sc$ylog
			logxy <- sc$xylog

			valid <- !is.na(npdedata[["NPDE"]]) & preddvid == dname
			if (!all(valid == FALSE)) {
				cat(sprintf("plot npde dvid %i\n", dname))

				npde <- as.numeric(npdedata[valid, "NPDE"])
				lim <- range(c(npde, -npde), na.rm=TRUE)

				# panel 1: NPDE density vs. standard normal N(0, 1)
				npde_clean <- npde[is.finite(npde)]
				if (length(npde_clean) >= 2) {
					npde_mean <- mean(npde_clean)
					npde_sd <- sd(npde_clean)
					d_npde <- density(npde_clean)

					d_xlim <- range(c(d_npde$x, -4, 4), na.rm=TRUE)
					x_grid <- seq(d_xlim[1], d_xlim[2], length.out=200)
					y_norm <- dnorm(x_grid, mean=0, sd=1)
					d_ylim <- c(0, max(c(d_npde$y, y_norm), na.rm=TRUE) * 1.35)

					col <- adjustcolor("blue", 0.25)
					plot(x=d_npde$x, y=d_npde$y, type="n",
						 xlab="NPDE", ylab="Density", xlim=d_xlim, ylim=d_ylim)
					polygon(c(x_grid[1], x_grid, x_grid[length(x_grid)]),
							c(0, y_norm, 0),
							col=col, border=NA)
					abline(v=0, lty=3, col="grey50")
					lines(x=d_npde$x, y=d_npde$y, col="black", lwd=2)
					grid()

					col <- adjustcolor("blue", 0.25)
					legend("topright",
						   legend=c(
							   sprintf("NPDE (%.3g, %.3g)", npde_mean, npde_sd),
							   "N(0, 1)"),
						   border=c(NA, NA),
						   col=c("black", col),
						   lty=c(1, 1),
						   lwd=c(2, 8),
						   bty="n", cex=0.95)
				} else {
					frame()
				}

				# panel 2: Normal Q-Q plot of NPDE vs. N(0, 1)
				if (length(npde_clean) >= 2) {
					npde_sorted <- sort(npde_clean)
					q_theo <- qnorm(ppoints(length(npde_sorted)))
					qq_lim <- range(c(q_theo, npde_sorted), na.rm=TRUE)

					plot(x=q_theo, y=npde_sorted, type="n",
						 xlab="Theoretical Quantiles N(0, 1)", ylab="Sample Quantiles (NPDE)",
						 xlim=qq_lim, ylim=qq_lim)
					points(x=q_theo, y=npde_sorted, col=pointcol)
					abline(0, 1, col="black")
					grid()
				} else {
					frame()
				}

				frame()

				x <- as.numeric(npdedata[valid, "TIME"])
				y <- npde
				plot(x=x, y=y, type="n", xlab="TIME", ylab="NPDE", xlim=timelim, ylim=lim)
				points(x=x, y=y, col=pointcol)
				abline(coef = c(0,0))
				grid()

				if ("YHAT" %in% pnames) {
					x <- itr(npdedata[valid, "YHAT"])
					y <- npde
					plot(x=x, y=y, type="n", xlab="YHAT", ylab="NPDE", ylim=lim)
					points(x=x, y=y, col=pointcol)
					abline(coef = c(0,0))
					grid()
				} else
					frame()

				# page title for this data type: configured name or generic label
				# (note: %i assumes numeric DVTY codes)
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

# likelihood profiles: objective function vs. one parameter value, per parameter,
# with the profile-likelihood and covariance-step confidence intervals shown together
cat("check likelihood profiles\n")
likelihood_profile_plots <- TRUE
if (!is.data.frame(profiles) || NROW(profiles) == 0) {
	cat("... no likelihood profiles\n")
	likelihood_profile_plots <- FALSE
}
if (likelihood_profile_plots) {
	cat("plot likelihood profiles\n")
	par(mfrow=c(2,3), cex=1)

	p <- profiles
	p <- p[!duplicated(p), ]   # the file can repeat rows
	# one group per (parameter TYPE, INDEX), e.g. "theta3", "omega1"
	p$profilegroup <- paste0(p[["TYPE"]], p[["INDEX"]])
	groups <- unique(p$profilegroup)

	# remove infinities
	p <- p[is.finite(p$VALUE), ]

	# load covlimits 
	covlim <- covlimits
	if (is.data.frame(covlim)) 
		covlim$profilegroup <- paste0(covlim[["TYPE"]], covlim[["INDEX"]])

	for (i in groups) {
		idata <- p[p$profilegroup == i, ]

		# omega fixed at zero (unused row/col): nothing to profile
		fv <- idata[idata[["NEVAL"]] == 0 & idata[["OBJFN"]] != 0, "VALUE"]
		if (length(fv) == 0) fv <- idata[["VALUE"]]
		if (all(idata[["TYPE"]] == "omega") && all(fv == 0)) {
			cat(sprintf("skip profile %s (fixed at zero)\n", i))
			next
		}

		cat(sprintf("plot profile %s\n", i))
		# marker rows (OBJFN==0, TARGET!=0) hold the exact profile limits from the root finder;
		# they are not real objective values so they are split off from the curve data
		bounds <- idata[idata$OBJFN == 0. & idata$TARGET != 0., "VALUE"]
		idata <- idata[idata$OBJFN != 0., ]
		final <- idata[idata$NEVAL == 0, ]   # NEVAL==0: the converged estimate itself

		minobjfn <- min(idata$OBJFN)
		# height of the profile CI box: target OFV rise above the minimum
		pheight <- max(final$TARGET - minobjfn)

		# covariance-step (Wald-type) limits for comparison, if available
		covdata <- NA
		covlower <- NA
		covupper <- NA
		cheight <- 0
		if (is.data.frame(covlim)) {
			covdata <- covlim[covlim$profilegroup == i, ]
			covlower <- covdata$LOWER
			covupper <- covdata$UPPER
			# no cov-step limits for this parameter: treat as unavailable
			if (length(covlower) != 1 || length(covupper) != 1) {
				covlower <- NA
				covupper <- NA
			}
			cheight <- max(covlim$DOBJFN)
		}
		dheight <- max(pheight, cheight)   # y scale is set by whichever box is taller

		ox <- order(idata$VALUE)
		idata <- idata[ox, ]

		type <- unique(idata[["TYPE"]])
		idx <- unique(idata[["INDEX"]])
		plotname <- paste0(type, "(", idx, ")")

		# variances (omega/sigma) are positive and skewed, so default to a log x axis for them;
		sc <- make_scale(if (type %in% c("omega", "sigma")) "log" else "identity")
		tr <- sc$tr
		itr <- sc$itr
		xlog <- sc$xlog

		# x range: span the limits plus 20% margin on each side, padded on the transformed scale
		range_lower <- tr(min(bounds, covlower, na.rm=TRUE))
		range_upper <- tr(max(bounds, covupper, na.rm=TRUE))
		range <- range_upper - range_lower
		xlim <- itr(c(range_lower - 0.2 * range, range_upper + 0.2*range))

		x <- idata[["VALUE"]]
		y <- idata[["OBJFN"]]
		# y range shows up to 3x the CI height above the minimum
		plot(x=x, y=y, type="n", xlab=plotname, ylab="Objective function", xlim=xlim, ylim=c(minobjfn, minobjfn+dheight*3), log=xlog)

		# exact limits found by the root finder (OBJFN==0 marker rows)
		# shaded box = profile-likelihood confidence interval
		lower <- tr(min(bounds))
		upper <- tr(max(bounds))
		rect(xleft=itr(lower), xright=itr(upper), ybottom=minobjfn-20, ytop=minobjfn+pheight, col="lightsteelblue", border="lightsteelblue")
		
		# blue outline = covariance-step interval, for comparison with the profile CI
		covlegend <- ""
		if (!is.na(covlower) && !is.na(covupper)) {
			rect(xleft=covlower, xright=covupper, ybottom=minobjfn-20, ytop=minobjfn+cheight, col=NA, border="blue")
			covlegend <- sprintf("cov [%.3g, %.3g]", covlower, covupper)
		}

		# make a smoothed line
		# monotone Hermite interpolation (no overshoot between points), fitted on the
		# transformed x axis so it looks smooth on a log scale; evaluated at 100 points
		sfunc <- splinefun(x=tr(x), y=y, method="monoH.FC", ties=mean)
		xl <- itr(seq(tr(xlim[1]), tr(xlim[2]), length.out=100))
		yl <- sfunc(tr(xl))
		lines(x=xl, y=yl, col="red")

		if (covlegend != "") {
			covmid <- (tr(covlower) + tr(covupper)) / 2
			covhalf <- (tr(covupper) - tr(covlower)) / 2
			yp <- minobjfn + cheight * ((tr(xl) - covmid) / covhalf)^2
			lines(x=xl, y=yp, col="blue", lty=3, lwd=2)
		}

		points(x=x, y=y, pch=21, col="red", bg="white", cex=0.8)   # evaluated points
		value <- final[["VALUE"]]
		points(x=value, y=final[["OBJFN"]], col="black", pch=16)   # the final estimate

		legend("top", 
			legend=c(sprintf("estimate %.3g", value),
					 sprintf("profile [%.3g, %.3g]", itr(lower), itr(upper)),
					 covlegend),
			bty="n", cex=1)
		
		grid()
		box()
	}
}

# VPC plots: one page per VPC name (first column), strata stacked on the page
# Prediction-corrected VPC: observed percentiles (black lines) are drawn over
# simulated confidence bands (red = median, blue = other percentiles).
cat("check vpc plots\n")
vpc_plots <- TRUE
if (!is.data.frame(vpcout) || NROW(vpcout) == 0) {
	cat("... no vpc out data\n")
	vpc_plots <- FALSE
}
if (!is.data.frame(vpcpoints) || NROW(vpcpoints) == 0) {
	cat("... no vpc points data\n")
	vpc_plots <- FALSE
}
if (vpc_plots) {
	cat("plot vpc\n")

	# percentiles present in the file, from the OBS_P<x> columns
	# (so the script adapts to whichever percentiles were computed)
	vpck <- sub("^OBS_", "", grep("^OBS_", names(vpcout), value=TRUE))
	vpcnames <- unique(vpcout$NAME)
	for (nm in vpcnames) {
		cat(sprintf("plot vpc \"%s\"\n", nm))
		v <- vpcout[vpcout$NAME == nm, ]
		p <- vpcpoints[vpcpoints$NAME == nm, ]
		# scale: from the DVTY/DVID column of the VPC file if present, else NAME is taken
		# as the DVTY code; a config row NAME=VPC, VALUE=<vpc name>, TYPE=scale overrides
		dcol <- intersect(c("DVTY", "DVID"), names(v))
		dkey <- if (length(dcol)) v[[dcol[1]]][1] else nm
		sc <- make_scale(get_config_setting("VPC", nm, "scale",
				get_config_setting("DVTY", dkey, "scale", "identity")))
		itr <- sc$itr
		logy <- sc$ylog
		# bin interpolation: config row NAME=DVTY, VALUE=<dvty>, TYPE=vpc.smooth, SETTING=linear
		# (anything else, or no row: one block per bin)
		if (is.na(dkey))
			vpcsmooth <- get_config_setting("DVTY", "-", "vpc.smooth", "none")
		vpcsmooth <- get_config_setting("DVTY", dkey, "vpc.smooth", vpcsmooth)
		strata <- sort(unique(v$STRAT))
		nstrata <- length(strata)
		par(mfrow=c(nstrata, 1), oma=c(0, 0, 2, 0))   # strata stacked vertically, room for outer title
		for (s in strata) {
			vs <- v[v$STRAT == s, ]
			ps <- p[p$STRAT == s, ]
			# simulated lower/upper band limits for every percentile (used for the y range)
			# (data are on the transformed scale; itr maps them to the natural scale)
			pcdv <- itr(ps$PCDV)
			lo <- itr(as.matrix(vs[, paste0("SIM_", vpck, "_LO"), drop=FALSE]))
			hi <- itr(as.matrix(vs[, paste0("SIM_", vpck, "_HI"), drop=FALSE]))
			plot(ps$X, pcdv, type="n", main=NA,
				 xlab="Time", ylab="Prediction-corrected DV", log=logy,
				 ylim=range(pcdv, lo, hi, na.rm=TRUE))
			if (nstrata > 1)
				title(paste("STRAT =", s))
			points(ps$X, pcdv, pch=20, col="grey70")   # observed points, prediction-corrected
			# one rectangle per bin (BINMIN..BINMAX in time, simulated LO..HI in y) for every percentile;
			# observed percentile connected through the bin centres
			vs <- vs[order(vs$BIN), ]
			for (k in vpck) {
				klo <- itr(vs[[paste0("SIM_", k, "_LO")]])
				khi <- itr(vs[[paste0("SIM_", k, "_HI")]])
				col <- adjustcolor(if (k == "P50") "red" else "blue", 0.25)
				if (vpcsmooth == "linear" && nrow(vs) >= 3 && !anyNA(c(klo, khi))) {
					# one polygon with a vertex at every bin centre, edges linearly interpolated;
					# the first half of the first bin and the second half of the last bin are flat
					n <- nrow(vs)
					ctr <- (vs$BINMIN + vs$BINMAX)/2
					px <- c(vs$BINMIN[1], ctr, vs$BINMAX[n])
					pidx <- c(1, 1:n, n)
					polygon(c(px, rev(px)), c(khi[pidx], rev(klo[pidx])), border=NA, col=col)
				} else {
					# bin edges (histogram limits) from the VPC file: neighbouring bins butt up
					rect(vs$BINMIN, klo, vs$BINMAX, khi, border=NA, col=col)
				}
				lines((vs$BINMIN + vs$BINMAX)/2, itr(vs[[paste0("OBS_", k)]]), lwd=2)   # observed percentile through bin centres
			}
		}
		title(nm, outer=TRUE)
	}
	par(mfrow=c(1,1), oma=c(0, 0, 0, 0))   # reset layout
}

# individual plots: one page per individual, all DVTY on the page
# Black dots = observations, red line = individual prediction (YHAT),
# blue dashed = population prediction (PRED).
cat("check individual plots\n")
individual_plots <- TRUE
if (!is.data.frame(yhatdata) || NROW(yhatdata) == 0) {
	cat("... no individual data\n")
	individual_plots <- FALSE
}
if (!all(c("ID", "TIME", "YHATDV", "YHAT", "YHATVAR") %in% names(yhatdata))) {
	cat("... no individual data for ID, TIME, YDV, YHAY, YHATVAR\n")
	individual_plots <- FALSE
}
if (individual_plots) {
	cat("plot individuals\n")
	pnames <- names(yhatdata)

	if ("DVTY" %in% pnames) {
		idvid <- yhatdata[["DVTY"]]
	} else if ("DVID" %in% pnames) {
		idvid <- yhatdata[["DVID"]]
	} else {
		idvid <- rep(1, nrow(yhatdata))
	}
	ihaspred <- "PRED" %in% pnames
	ivalid <- yhatdata[["YHATVAR"]] != 0     # observation rows only
	idvtypes <- sort(unique(idvid[ivalid]))

	# same scaling as the population plots, one scale per DVTY (keyed by DVTY as character)
	iscale <- list()
	for (dname in idvtypes)
		iscale[[as.character(dname)]] <- get_scale("DVTY", dname)
	itrans <- function(v, dname) iscale[[as.character(dname)]]$itr(as.numeric(v))

	# axis limits shared by all individuals (time) and by all individuals per DVTY (y)
	# so pages are directly comparable between individuals
	itimelim <- range(yhatdata[ivalid, "TIME"], na.rm=TRUE)
	iylim <- list()
	for (dname in idvtypes) {
		sel <- ivalid & idvid == dname
		v <- c(itrans(yhatdata[sel, "YHATDV"], dname),
			   itrans(yhatdata[sel, "YHAT"], dname))
		if (ihaspred)
			v <- c(v, itrans(yhatdata[sel, "PRED"], dname))
		iylim[[as.character(dname)]] <- range(v, na.rm=TRUE)   # list keyed by DVTY (as character)
	}

	# one panel per DVTY on each page
	# near-square grid: ncol = ceil(sqrt(n)), nrow = ceil(n / ncol)
	npanel <- length(idvtypes)
	ncol_ <- ceiling(sqrt(npanel))
	nrow_ <- ceiling(npanel / ncol_)

	# draw a profile line, new line when time jumps backwards (multiple sessions)
	iplot_lines <- function(x, y, col, lty=1)
	{
		sess <- cumsum(c(FALSE, diff(x) < 0))
		for (s in unique(sess))
			lines(x=x[sess == s], y=y[sess == s], col=col, lty=lty)
	}

	# individual descriptors for the page title: label = candidate column names (first found is used)
	# Each label maps to a list of possible column names; the first one present in the data is kept
	# ([1] gives NA if none are present, and those labels are then dropped)
	idesc <- list(STY=c("STY", "STUDY"),
				  AGE=c("AGE"),
				  WT=c("WGT", "WT", "WEIGHT", "BW"),
				  HT=c("HGT", "HT", "HEIGHT"),
				  SEX=c("SEX", "GENDER"))
	idesc <- lapply(idesc, function(n) n[n %in% pnames][1])
	idesc <- idesc[!is.na(idesc)]
	# Build a title such as "ID 12, AGE 34, WT 70.5" from the individual's first row
	ititle <- function(i)
	{
		row <- yhatdata[yhatdata[["ID"]] == i, , drop=FALSE][1, ]
		s <- sprintf("ID %s", i)
		for (lab in names(idesc)) {
			v <- row[[idesc[[lab]]]]
			# numbers with 4 significant digits, text (e.g. SEX) as is
			s <- paste0(s, sprintf(", %s %s", lab, if (is.numeric(v)) sprintf("%.4g", v) else v))
		}
		s
	}

	# one page per individual
	for (i in unique(yhatdata[["ID"]])) {
		par(mfrow=c(nrow_, ncol_), cex=1, oma=c(0, 0, 2, 0))

		for (dname in idvtypes) {
			logy <- iscale[[as.character(dname)]]$ylog
			ylim <- iylim[[as.character(dname)]]

			t <- get_config_setting("DVTY", dname, "name", NA)
			if (is.na(t))
				t <- sprintf("data type %s", dname)

			# this individual's observation rows for this DVTY
			sel <- ivalid & idvid == dname & yhatdata[["ID"]] == i
			# empty frame with shared limits first, so axes are identical on every page
			plot(x=NA, y=NA, type="n", main=t, xlab="TIME", ylab="DV, YHAT, PRED",
				 xlim=itimelim, ylim=ylim, log=logy)
			grid()

			if (any(sel, na.rm=TRUE)) {
				x <- yhatdata[sel, "TIME"]
				dv <- itrans(yhatdata[sel, "YHATDV"], dname)
				yhat <- itrans(yhatdata[sel, "YHAT"], dname)
				if (ihaspred)
					iplot_lines(x, itrans(yhatdata[sel, "PRED"], dname), col="blue", lty=2)
				iplot_lines(x, yhat, col="red")
				points(x=x, y=dv, col="black", pch=16)
			} else
				# keep the panel but say so when this individual has no data for this DVTY
				# (geometric mean of ylim is the visual centre of a log axis)
				text(x=mean(itimelim), y=if (logy == "y") sqrt(prod(ylim)) else mean(ylim),
					 labels="no data")
		}
		title(ititle(i), outer=TRUE)
	}
	par(mfrow=c(1,1), oma=c(0, 0, 0, 0))
}

# close the PDF device
dummy <- dev.off()   
cat("done\n")
