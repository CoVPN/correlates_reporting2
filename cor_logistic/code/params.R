.mfrow <- c(1, 1)
numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows", 1, future::availableCores()))

