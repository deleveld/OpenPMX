########################################################################
	# add methods to table with empty places
	i <- 1
	allmethods <- list()
	if (exists("nonmem")) {
		allmethods[[i]] <-list(data=nonmem,			name="NONMEM",		col="grey",			lw=3)
		i <- i + 1
	}
	if (exists("nonmemlaplace")) {
		allmethods[[i]] <- list(data=nonmemlaplace,	name="LAPLACE",		col="lightblue",	lw=3)
		i <- i + 1
	}
	if (exists("nonmemsaem")) {
		allmethods[[i]] <- list(data=nonmemsaem,	name="NONMEM-SAEM",		col="lightblue",	lw=3)
		i <- i + 1
	}
	if (exists("nonmemitsb")) {
		allmethods[[i]] <- list(data=nonmemitsb,	name="NONMEM-ITSB",		col="gold",	lw=1)
		i <- i + 1
	}
	if (exists("validate")) {
		allmethods[[i]] <- list(data=validate, 	name="Validate",		col="lightgreen",	lw=5)
		i <- i + 1
	}
	if (exists("gronmem")) {
		allmethods[[i]] <- list(data=gronmem, 		name="OpenPMX",		col="black",		lw=1)
		i <- i + 1
	}
	if (exists("gronmemicov")) {
		allmethods[[i]] <- list(data=gronmemicov, 	name="Icov",	col="blue",			lw=1)
		i <- i + 1
	}
	if (exists("gronmemtest")) {
		allmethods[[i]] <- list(data=gronmemtest, 	name="TEST",		col="lightgreen",	lw=3)
		i <- i + 1
	}

########################################################################
	if (!exists("tab"))
		tab <- data.frame()
		
	nmethods <- length(allmethods)	
	for (i in 1:nmethods) {
		d <- allmethods[[i]]
		name <- d[["name"]]

		newr <- tab[1, ]
		newr[["name"]] <- name
		newr[["type"]] <- "BIAS"
		for (i in trueparams) 
			newr[[i["tabname"]]] <- 0
		tab <- rbind(tab, newr)

		newr[["type"]] <- "RMSE"
		tab <- rbind(tab, newr)
	}

	bias <- data.frame()
	rmse <- data.frame()
	pe <- data.frame()
	ape <- data.frame()
	mdifabspe <- data.frame()
	correlation <- data.frame()

	setup_par()
	for (i in trueparams) {
		paramname <- i["tabname"]
		trueval <- as.numeric(i["trueval"])
		displayname <-i["label"]
		logv <-i["logv"]
		units <-i["units"]

		if (logv == TRUE)
			trueval <- exp(trueval)

		cat("--------------------------------------------------------\n")
		print(paramname)
		
		nonmemdata <- data.frame()
		
		allv <- c()
		for (i in 1:nmethods) {
			d <- allmethods[[i]]
			data <- d[["data"]]
			v <- data[[paramname]]
			if (logv == TRUE)
				v <- exp(v)
			allv <- c(allv, v)
			name <- d[["name"]]
			if (name == "NONMEM")
				nonmemdata <- d[["data"]]
		}
		
		calc_bias <- function(est, trueval)
		{
			est <- (est - trueval) / trueval * 100
			ret <- mean(est)
			ret
		}
		calc_pe <- function(est, trueval)
		{
			est <- (est - trueval) / trueval * 100
			ret <- median(est)
			ret
		}
		calc_ape <- function(est, trueval)
		{
			est <- abs(est - trueval) / trueval * 100
			ret <- median(est)
			ret
		}
		calc_rmse <- function(est, trueval)
		{
#			err <- (est - trueval) / trueval
#			ret <- sqrt(mean(err*err)) * 100
#			ret

### corrected equation but the same results as above
			err <- (est - trueval)
			ret2 <- sqrt(mean(err*err) / (trueval*trueval)) * 100
			ret2
		}
		calc_medrerr <- function(est1, est2, trueval)
		{
			est1 <- abs(est1 - trueval) / trueval * 100
			est2 <- abs(est2 - trueval) / trueval * 100

			ret <- quantile(est1 - est2, 0.5)
			ret
		}
		calc_perf <- function(est, trueval)
		{
			ret <- abs(est - trueval) / trueval * 100
			ret
		}
		for (i in 1:nmethods) {
			d <- allmethods[[i]]
			data <- d[["data"]]
			v <- data[[paramname]]
			if (logv == TRUE)
				v <- exp(v)

			# remove invalid cases
			valid <- !is.na(v)
			v_ <- v[valid]
		
			bias_ <- calc_bias(v_, trueval)
			rmse_ <- calc_rmse(v_, trueval)
			pe_ <- calc_pe(v_, trueval)
			ape_ <- calc_ape(v_, trueval)

			name <- d[["name"]]
			bias[name, displayname] <- bias_
			rmse[name, displayname] <- rmse_
			pe[name, displayname] <- pe_
			ape[name, displayname] <- ape_

			tabr <- which(tab[["name"]] == name & tab[["type"]] == "BIAS")
			stopifnot(length(tabr) == 1)
			tab[tabr, paramname] <- bias_
			tab[tabr, paramname] <- round(bias_, digits=1)

			tabr <- which(tab[["name"]] == name & tab[["type"]] == "RMSE")
			stopifnot(length(tabr) == 1)
			tab[tabr, paramname] <- rmse_
			tab[tabr, paramname] <- round(rmse_, digits=1)

			# remove invalid cases from NONMEM too
			nm1 <- nonmemdata[[paramname]]
			if (logv == TRUE)
				nm1 <- exp(nm1)
			valid <- !is.na(v) & !is.na(nm1)
			v_ <- v[valid]
			nm1_ <- nm1[valid]

			mdifabspe[name, displayname] <- calc_medrerr(v_, nm1_, trueval)

			correlation[name, displayname] <- cor(v_, nm1_)
		}

		tallv <- allv
		xlim <- range(c(trueval, allv), na.rm=TRUE)
		dupper <- max(xlim)-trueval
		dlower <- trueval-min(xlim)
		xlim <- range(c(xlim), trueval+dlower, trueval-dupper)
		ttt <- bquote(.(displayname) ~ .(units))

		log <- ""
		
		xaxs <- "r"
		if (any(grepl("var", ttt) == TRUE)) {
			xaxs <- "i"
			rdelta <- xlim[2] - trueval
			if (trueval - rdelta < 0)
				xlim[1] <- 0.
		}
		if (any(grepl("corr", ttt) == TRUE)) {
			xaxs <- "i"
			if (max(allv) > 90)
				xlim[2] <- 100
			if (min(allv) < -90)
				xlim[1] <- -100
		}
		if (any(grepl("deter", ttt) == TRUE)) {
			xaxs <- "i"
			if (min(xlim) < 0)
				xlim[1] <- 0
		}
		
		d0 <- density(tallv, na.rm=TRUE)
		plot(d0, type="n", yaxt="n", xlab=NA, ylab=NA, xlim=(xlim), ylim=c(0,1.4), log=log, main=NA, xaxs=xaxs, yaxs="i")
		mtext("Density", side=2, line=0.5, cex=par()$cex)

		for (i in 1:nmethods) {
			d <- allmethods[[i]]
			data <- d[["data"]]
			col <- d[["col"]]
			lw <- d[["lw"]]
			v <- data[[paramname]]
			if (logv == TRUE)
				v <- exp(v)	
			d1 <- density(v, na.rm=TRUE)
			lines(x=(d1$x), y=d1$y/max(d1$y), col=col, lw=lw)
		}
		lines(x=(c(trueval,trueval)), y=c(0,1), col="black", lw=2)

		grid()
		box()
		
#		ttt <- sprintf("%s %s", displayname, units)
		title(ttt, line=0.5, cex=par()$cex)

		allname <- NA
		allcol <- NA
		alllw <- NA
		for (i in 1:nmethods) {
			d <- allmethods[[i]]
			name <- d[["name"]]
			namerow <- sprintf("%s (bias %.1f%% rmse %.1f%%)", name, bias[name, displayname], rmse[name, displayname]) 
			col <- d[["col"]]
			lw <- d[["lw"]]
			allname <- c(allname, namerow)
			allcol <- c(allcol, col)
			alllw <- c(alllw, lw)
		}
		allname <- allname[!is.na(allname)]
		allcol <- allcol[!is.na(allcol)]
		alllw <- alllw[!is.na(alllw)]

		legend("topleft",
			legend=allname,
			col=allcol,
			lw=alllw,
			bty="n", cex=par()$cex*0.85)
	}
	title("Distribution of parameter estimates", outer=TRUE)
	
	print("bias")
	print(bias)

	print("rmse")
	print(rmse)
	print("ape")
	print(ape)
	print("median diff abs pe")
	print(mdifabspe)
	print("correlation")
	print(correlation)

	# rank performance
	tabbias <- tab[tab[["type"]] == "BIAS", ]
	tabrmse <- tab[tab[["type"]] == "RMSE", ]
	row.names(tabbias) <- tabbias[["name"]]
	row.names(tabrmse) <- tabrmse[["name"]]
	tabbias <- tabbias[ , !(colnames(tabbias) %in% c("type", "name"))]
	tabrmse <- tabrmse[ , !(colnames(tabrmse) %in% c("type", "name"))]

	calculate_performance_rank <- function()
	{
		tabbias_rank <- tabbias
		tabrmse_rank <- tabrmse
		for (i in trueparams) {
			paramname <- i["tabname"]

			v <- tabbias[[paramname]]
			tabbias_rank[[paramname]] <- rank(abs(v))
			
			v <- tabrmse[[paramname]]
			tabrmse_rank[[paramname]] <- rank(abs(v))
		}
		tabbias_rank[["avgrank"]] <- rowMeans(tabbias_rank)
		tabrmse_rank[["avgrank"]] <- rowMeans(tabrmse_rank)

		tabbias[["avgrank"]] <<- tabbias_rank[["avgrank"]]
		tabrmse[["avgrank"]] <<- tabrmse_rank[["avgrank"]]
	}
	calculate_performance_rank()

	calculate_performance_score <- function(tabv)
	{
		absv <- abs(tabv)
		minv <- apply(absv, 2, min)
		maxv <- apply(absv, 2, max)
		delta <- maxv - minv
		t1 <- sweep(absv, 2, minv, FUN="-")
		t2 <- sweep(t1, 2, delta, FUN="/") * 100
		t2 <- rowMeans(t2)
		t2
	}
	tabbias[["avgperf"]] <- calculate_performance_score(tabbias)
	tabrmse[["avgperf"]] <- calculate_performance_score(tabrmse)
	
#	tabbias <- tabbias[order(tabbias[["avgrank"]]), ]
#	tabrmse <- tabrmse[order(tabrmse[["avgrank"]]), ]

	tabbias <- tabbias[order(tabbias[["avgperf"]]), ]
	tabrmse <- tabrmse[order(tabrmse[["avgperf"]]), ]

	print(tabbias)
	print(tabrmse)
