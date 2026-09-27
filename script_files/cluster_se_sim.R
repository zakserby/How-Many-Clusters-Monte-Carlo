#===============================================================================
# AEM 6850
# How many clusters do you need for clustered standard errors?
#===============================================================================

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 1). Preliminary -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Clean up workspace and load or install necessary packages if necessary
rm(list=ls())
want <- c("multiwayvcov")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only=TRUE))
rm(want, need)

# Working directories
dir <- list()
dir$root <- dirname(getwd())
dir$output_figure <- paste(dir$root,"/output_figure",sep="")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 2). Main code -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Generate new data
set.seed(123)
nlist <- c(100,200,500, 10^3, 10^4) # sample sizes
glist <- c(2,4,5,10,20,25,50,100)   # clusters
r <- 1000 # number of simulations

d <- lapply(nlist, function(n){

  d <- lapply(glist, function(g){

    # Generate data (Note: BC's X's fixed we simulate them outside the MC simulation)
    x1 <- rnorm(n) # generate the X that is correlated within clusters
    x2 <- rnorm(n) + rep(rnorm(g), each=n/g) # generate the X that is correlated within clusters
    data <- data.frame(x1=x1,x2=x2)
    data$clu <- rep(1:g, each=n/g) # cluster id
    beta0 <- 1
    beta1 <- 1
    beta2 <- 1

    # Run MC simulation
    out3 <- lapply(1:r, function(sim){

      # Generate data
      e1 <- rnorm(n) # unique variance for each obs
      e2 <- rep(rnorm(g), each=n/g) # variance common within each cluster
      e  <- e1 + e2
      data$y <- beta0 + beta1*data$x1 + beta2*data$x2 + e

      # Run OLS
      reg <- lm(y~x1+x2, data)

      # Get coefficients
      b <- coef(reg)

      # Get SEs with different approaches
      se.naive  <- sqrt(diag(vcov(reg)))
      se.clu    <- sqrt(diag(cluster.vcov(reg, data$clu, df_correction=F)))

      # Export
      o <- c(b, se.naive, se.clu)
      names(o) <- paste(rep(c("b.ols","se.naive","se.clu"), each=3), names(b))
      o
    })

    out3 <- do.call("rbind", out3)

    # Calculate the ratio to see whether your clustered SE works well for SE(beta2_hat)
    ratio <- mean(out3[,"se.clu x2"]) / sd(out3[,"b.ols x2"])

    data.frame(n=n, g=g, ratio=ratio)
  })

  do.call("rbind", d)
})

d <- do.call("rbind", d)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 3). Plot -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Save as png
png(paste0(dir$output_figure, "/How_Many_Clusters.png"), width=1400, height=850, res = 150)

# Settings
xticks <- 1:length(glist)
ylim <- c(0.0, 1.1)
yticks <- seq(0.0, 1.0, by = 0.2)
ytick_labels <- format(yticks, nsmall = 1) # To format 0.0 and 1.0

# Colors for different lines/# of observations
line_colors <- c("#c7e9b4", "#7fcdbb", "#44b6c4", "#2c7fb8", "#253494")

# Plot figure
plot(0, type="n", axes=FALSE, xlim=range(xticks), ylim=ylim,
     main = "How many clusters do you need?", font.main = 2, cex.main = 1.5,
     xlab="Number of clusters", ylab="Std Err estimate / Std Dev of sampling distribution")

axis(1, at=xticks, labels=glist)
axis(2, at=yticks, labels=ytick_labels, las=2)
abline(h=1, lty=3)

# Loop over sample sizes to plot lines and points
for(i in seq_along(nlist)){
  subset <- d[d$n == nlist[i], ] # subset for current sample size
  lines(xticks, subset$ratio, col=line_colors[i], lwd=2)
  points(xticks, subset$ratio, col=line_colors[i], pch=19)
}

# Legend
legend("bottomright", title = "Number of observations", legend=nlist,
       col=line_colors, lwd=2, pch=19, inset=c(0.025, 0.05))
box()

dev.off()


