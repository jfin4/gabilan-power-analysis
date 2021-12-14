# set options
options(browser = "C:/Users/jinman/firefox/firefox.exe")
library(nlme)
set.seed(1)

# clear memory
rm(list =ls())

# define functions
get_sd <- function(x) {
    sqrt(var(x))
}

# load data
data <- read.csv("data.csv")

# transform data
data$SampleDate <- as.Date(data$SampleDate, format = "%m/%d/%Y")
data$Year <- as.numeric(format(data$SampleDate, "%Y"))
data$Month <- as.numeric(format(data$SampleDate, "%m"))
data$Result <- log10(data$Result + 0.01)

# subset data
# subset by year
# data <- data[data$Year > 2016, ]
# subset by month
start <- 10
stop <- 4
data <- data[data$Month >= start | data$Month <= stop, ]
data <- data[c("ParentProject", "StationCode", "SampleDate", "Result", "Year")]
# Removed 309TEM because no data
sites <- c("309GAB", "309NAD", "309ALG", "309ALD", "309JON", "309TEH", 
           "309TDW", "309OLD", "309ASB", "309MER", "309ESP", "309RTA")
data <- data[data$StationCode %in% sites, ]

# run model
mod <- lme(Result ~ Year, data = data, random = ~ 1 | StationCode)

# check residuals
# plot(mod)

# create parameter tables
ns_wide <- tapply(data$Result, list(data$StationCode, data$Year), length)
ns_wide <- as.data.frame(ns_wide)
ns_wide <- cbind(StationCode = row.names(ns_wide), ns_wide)

means_wide <- tapply(data$Result, list(data$StationCode, data$Year), mean)
means_wide <- as.data.frame(means_wide)
means_wide <- cbind(StationCode = row.names(means_wide), means_wide)

sds_wide <- tapply(data$Result, list(data$StationCode, data$Year), get_sd)
sds_wide <- as.data.frame(sds_wide)
sds_wide <- cbind(StationCode = row.names(sds_wide), sds_wide)

#>>>>>>> Tuning knobs >>>>>>>
# tuning knob for baseline
baseline_years <- 5
means_wide_sub <- means_wide[, (ncol(means_wide) - baseline_years):ncol(means_wide)]
baseline_means <- apply(means_wide_sub, 1, mean, na.rm = TRUE)
baseline <- mean(baseline_means)

# tuning knobs for final goals
end_goal <- 11 # NTU
total_years <- 20
log_end_goal <- log10(end_goal + 0.01)
rate <- (log_end_goal / baseline)^(1 / total_years) - 1

# tuning knob for interim goals 
no_years <- 5

# tuning knob for sample size
n_per_month <- 1
n <- ((12 - start + 1) + stop) * n_per_month

# tuning knob for number of iterations
N <- 10000

#<<<<<<< Tuning knobs <<<<<<<<
pvals <- rep(NA, N)
for (i in 1:N) {
    # calculate interim goal
    log_goal <- baseline * (1 + rate)^no_years
    
    # generate means
    gen_means_wide <- sapply(baseline_means, seq, log_goal, length.out = no_years)
    gen_means_wide <- as.data.frame(t(gen_means_wide))
    names(gen_means_wide) <- paste0("y", 1:no_years)
    gen_means_wide <- cbind(StationCode = means_wide[, 1], gen_means_wide)
    
    # generate sds 
    sds_wide_sub <- sds_wide[, (ncol(sds_wide) - baseline_years):ncol(sds_wide)]
    baseline_sds <- apply(sds_wide_sub, 1, mean, na.rm = TRUE)
    gen_sds_wide <- sapply(baseline_sds, rep, no_years)
    gen_sds_wide <- as.data.frame(t(gen_sds_wide))
    names(gen_sds_wide) <- paste0("y", 1:no_years)
    gen_sds_wide <- cbind(StationCode = sds_wide[, 1], gen_sds_wide)
    # attenuate sd over time proportional to mean (trickiest part)
    normalization <- apply(gen_means_wide[-1], 2, "/", gen_means_wide$y1)
    gen_sds_wide[-1] <- gen_sds_wide[-1] * normalization
    
    # reshape generated parameters to long
    gen_means <- reshape(gen_means_wide,
                direction = "long",
                v.names = "mean",
                times = names(gen_means_wide)[2:ncol(gen_means_wide)],
                timevar = "Year",
                varying = list(2:ncol(gen_means_wide)))
    gen_means$id <- NULL
    gen_sds <- reshape(gen_sds_wide,
                direction = "long",
                v.names = "sd",
                times = names(gen_sds_wide)[2:ncol(gen_sds_wide)],
                timevar = "Year",
                varying = list(2:ncol(gen_sds_wide)))
    gen_sds$id <- NULL
    
    # generate data
    stub <- merge(gen_means, gen_sds, by = c("StationCode", "Year"))
    stub$Year <- as.numeric(substr(stub$Year, 2, 4))

    if (n == 1) {
        gen_data_wide <- mapply(rnorm, n, stub$mean, stub$sd)
        gen_data_wide <- cbind(stub, "1" = gen_data_wide)
    } else {
        gen_data_wide <- t(mapply(rnorm, n, stub$mean, stub$sd))
        gen_data_wide <- cbind(stub, gen_data_wide)
    }
    
    # reshape generated data to long
    gen_data <- reshape(gen_data_wide, 
                        direction = "long",
                        v.names = "Result",
                        varying = list(5:ncol(gen_data_wide)))
    
    # run model on generated data
    gen_mod <- lme(Result ~ Year, data = gen_data, random = ~ 1 | StationCode)
    pval <- summary(gen_mod)$tTable["Year", "p-value"]
    pval <- round(pval, 4)
    
    pvals[i] <- pval
}
# calculate power
P <- sum(pvals <= 0.05) / N
P * 100

