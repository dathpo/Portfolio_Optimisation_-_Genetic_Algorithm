
library(quantmod)
library(GA)
library(ggplot2)
library(scales)
library(reshape2)
library(PortfolioAnalytics)
library(ROI)
library(ROI.plugin.quadprog)

options(repr.plot.width=6, repr.plot.height=4)

assets <- c("FCT.MI", "BP.L", "BA.L",
            "BATS.L", "GOOGL", "AAPL",
            "GSK", "XLNX", "CCL.L",
            "HL.L", "LLOY.L", "MKS.L",
            "PSON.L", "RR.L", "CNA.L", "SBRY.L")
assets <- sort(assets)
assets

stocks = lapply(assets, function(sym) {
  monthlyReturn(na.omit(getSymbols(sym,
                                   src="yahoo",
                                   from="2015-10-01",
                                   to="2017-10-01",
                                   auto.assign=FALSE)))
})

df <- as.data.frame(do.call(merge.xts, stocks))
df <- cbind(month = rownames(df), df)
names(df)[2:17] <- assets
df

fh.df <- df[1:9]
sh.cols <- c(1, 10:17)
sh.df <- df[sh.cols]

melted.fh.df <- melt(fh.df, id="month")
melted.sh.df <- melt(sh.df, id="month")


ggplot(data=melted.fh.df, aes(x=as.Date(melted.fh.df$month), y=melted.fh.df$value,
                           group=melted.fh.df$variable, color=melted.fh.df$variable)) +
geom_line() + geom_point() +
xlab("Date") + ylab("Return") + scale_y_continuous(labels=percent) +
scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
scale_color_discrete("Stock") + ggtitle("Stock Performance Oct 2015 - Oct 2017 (pt. 1)") +
theme(axis.text.x = element_text(angle = 45))
ggsave("stockperf.pdf", width=8, height=4)


ggplot(data=melted.sh.df, aes(x=as.Date(melted.sh.df$month), y=melted.sh.df$value,
                           group=melted.sh.df$variable, color=melted.sh.df$variable)) +
geom_line() + geom_point() +
xlab("Date") + ylab("Return") + scale_y_continuous(labels=percent) +
scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
scale_color_discrete("Stock") + ggtitle("Stock Performance Oct 2015 - Oct 2017 (pt. 2)") +
theme(axis.text.x = element_text(angle = 45))
ggsave("stockperf2.pdf", width=8, height=4)


mean.returns <- c(sapply(df[2:17], function(x) mean(x)))
sort(mean.returns, decreasing = TRUE)

mean.mean.returns <- sum(mean.returns)/length(mean.returns)
mean.mean.returns*100

cov.matrix <- cov(df[2:17])

obj.weights <- seq(from = 0, to = 1, by = 0.01)
#obj.weights <- obj.weights[-1]
hof <- vector("list", length(obj.weights))
hof.fitness <- vector("list", length(obj.weights))
lower.bounds <- c(rep(0, length(assets)))
upper.bounds <- c(rep(1, length(assets)))

normalise <- function(chromosome){
    total <- sum(chromosome)
    return(chromosome / total)
}

getMeanReturn <- function(mean.returns) {
    return(sum(mean.returns))
}

getTotalReturn <- function(mean.returns, weights) {
    return(sum(mean.returns * weights))
}

getTotalRisk <- function(cov.matrix, weights) {
    m <- cov.matrix
    m <- mapply("*", as.data.frame(m), weights)
    m <- t(m)
    m <- mapply("*", as.data.frame(m), weights)
    return(sum(m))
}

getFitness <- function(chromosome) {
    weights = normalise(chromosome)
    f = getTotalReturn(mean.returns, weights)
    g = getTotalRisk(cov.matrix, weights)
    fitness = w * f + (1 - w) * (-g)
    return(fitness)
}

for(i in seq_along(obj.weights)){
    w <- obj.weights[i]
    print(paste0("Run ", i, ", Ret. w: ", w))
    GA <- GA::ga(type='real-valued',
                 fitness = getFitness,
                 maxiter = 2000,
                 popSize = 100,
                 pcrossover = 0.8,
                 pmutation = 0.05,
                 lower = lower.bounds,
                 upper = upper.bounds,
                 seed = 123)
    hof[[i]] <- normalise(GA@solution)
    hof.fitness[[i]] <- GA@fitnessValue
}

summary(GA)
plot(GA)

hof

hof.fitness
sum(unlist(hof[1]))

returns <- vector('list', length(obj.weights))
risks <- vector('list', length(obj.weights))
sharpe_ratios <- vector('list', length(obj.weights))
points <- vector('list', length(obj.weights))

for(i in 1:length(obj.weights)) {
  sol.weights <- (unlist(hof[i]))   
  risk <- getTotalRisk(cov.matrix, sol.weights)
  exp.return <- getTotalReturn(mean.returns, sol.weights)
  returns[i] <- exp.return
  risks[i] <- risk
  sharpe_ratios[i] <- exp.return / sqrt(risk)
  points[[i]] <- c(risk, exp.return)
}

even.weighted.solution <- normalise(rep(1, length(assets)))
e.w.sol.return <- getTotalReturn(mean.returns, even.weighted.solution)
e.w.sol.risk <- getTotalRisk(cov.matrix, even.weighted.solution)

rand.solutions <- 1000
rand.obj.weights <- vector('list', length(assets))
rand.returns <- vector('list', rand.solutions)
rand.risks <- vector('list', rand.solutions)
rand.ratios <- vector('list', rand.solutions)

for(j in 1:rand.solutions) {
  rand.obj.weights[[j]] <- normalise(runif(length(assets)))
  rand.risk <- getTotalRisk(cov.matrix, rand.obj.weights[[j]])
  rand.exp.return <- getTotalReturn(mean.returns, rand.obj.weights[[j]])
  rand.returns[j] <- rand.exp.return
  rand.risks[j] <- rand.risk
  rand.ratios[j] <- rand.exp.return / sqrt(rand.risk)
}

p <- portfolio.spec(assets)
p <- add.constraint(portfolio = p, type = 'weight_sum', 
                    min_sum = 1, max_sum = 1)
p <- add.constraint(portfolio = p, type = 'box', min = 0, max = 1)
p <- add.objective(portfolio = p, type = 'risk', name = 'StdDev')
p <- add.objective(portfolio=p, type="return", name="mean")
pa.df <- df
pa.df["month"] <- NULL
opt <- optimize.portfolio(pa.df, portfolio=p, optimize_method='quadprog',
                          search_size=10000)
print(opt)

pa.weights <- extractWeights(opt)
pa.weights

pa.rand.weights <- extractWeights(opt)
pa.rand.weights

#print(unlist(sharpe_ratios))
#print(max(sapply(sharpe_ratios, max)))
#print(points[2])
#print(unlist(risks))
#print(min(sapply(risks, min)))
print(unlist(returns))
print(max(sapply(returns, max)))

pa.risk <- getTotalRisk(cov.matrix, pa.weights)
pa.exp.return <- getTotalReturn(mean.returns, pa.weights)

pa.rand.risk <- getTotalRisk(cov.matrix, pa.rand.weights)
pa.rand.exp.return <- getTotalReturn(mean.returns, pa.rand.weights)

rand.df <- as.data.frame(cbind(as.numeric(rand.risks), as.numeric(rand.returns)))
ga.df <- as.data.frame(cbind(as.numeric(risks), as.numeric(returns)))

plotReturns <- function(hof) {
    ggplot() +
      geom_line(data=ga.df, aes(x=ga.df$V1, y=ga.df$V2, linetype="Efficient Frontier"), size=0.3) +
      geom_point(data=rand.df, aes(x=rand.df$V1, y=rand.df$V2, fill=unlist(rand.ratios)), color=unlist(rand.ratios), shape=21) +
      geom_point(data=ga.df, aes(x=ga.df$V1, y=ga.df$V2, fill=unlist(sharpe_ratios)), shape=23, size=2) +
      geom_point(aes(x=e.w.sol.risk, y=e.w.sol.return, colour="Weighted Evenly"), size=3) +
      geom_point(aes(x=ga.df$V1[8], y=ga.df$V2[8], colour="Highest Sharpe Ratio"), size=3) +
      geom_point(aes(x=ga.df$V1[1], y=ga.df$V2[1], colour="Minimum Risk"), size=3) +
      geom_point(aes(x=ga.df$V1[99], y=ga.df$V2[99], colour="Maximum Return"), size=3) +
      geom_point(aes(x=pa.rand.risk, y=pa.rand.exp.return, colour="PortfolioAnalytics Random"), size=3) +
      geom_point(aes(x=pa.risk, y=pa.exp.return, colour="PortfolioAnalytics ROI"), size=3) +
      ggtitle("Portfolio Solutions") +
      scale_linetype("Line Type") + scale_color_discrete("Portfolios") +
      scale_fill_gradient("Sharpe Ratio", low="red", high="blue") +
      scale_y_continuous(name="Return", labels=percent) +
      scale_x_continuous(name="Risk", labels=percent)
}


plotReturns(hof) 
ggsave("GAcurv2.pdf", width=8, height=4.3)    
print(hof[8])
print(paste0(ga.df$V1[8], ", ", ga.df$V2[8]))

plotShareAlph <- function(title, solution) {
    ggplot(solution, aes(x=sort(solution$asset, decreasing = TRUE), y=solution$pf.share, fill=solution$asset)) + 
    geom_bar(stat = "identity", show.legend=FALSE) + 
    scale_y_continuous(labels=percent) + scale_x_discrete(labels=sort(solution$asset, decreasing = TRUE)) +
    ggtitle(title) + labs(x = "Stock", y = "Portfolio Share") + coord_flip()
}

plotShare <- function(title, solution) {
    ggplot(solution, aes(x=reorder(solution$asset, solution$pf.share), y=solution$pf.share, fill=solution$asset)) +
    geom_bar(stat = "identity", show.legend=FALSE) +
    scale_y_continuous(labels=percent) +
    labs(x = "Stock", y = "Portfolio Share") +
    ggtitle(title) + coord_flip()
}

create.sol.df <- function(solution) {
    sol.df <- data.frame(
        asset <- assets,
        pf.share <- solution
    )
    return(sol.df)
}

best.sharpe.sol <- create.sol.df(as.numeric(unlist(hof[8])))
plotShareAlph("Portfolio with Highest Sharpe Ratio", best.sharpe.sol)
ggsave("max_sharpe_sh.pdf", width=4, height=3)

max.ret.sol <- create.sol.df(as.numeric(unlist(hof[99])))
plotShareAlph("Portfolio with Maximum Return", max.ret.sol)
ggsave("max_ret_sh.pdf", width=4, height=3)

min.risk.sol <- create.sol.df(as.numeric(unlist(hof[1])))
plotShareAlph("Portfolio with Minimum Risk", min.risk.sol)
ggsave("min_risk_sh.pdf", width=4, height=3)

pa.sol <- create.sol.df(pa.weights)
plotShareAlph("PortfolioAnalytics ROI", pa.sol)
ggsave("pa_roi_sh.pdf", width=4, height=3)

pa.rand.sol <- create.sol.df(pa.rand.weights)
plotShareAlph("PortfolioAnalytics Random", pa.rand.sol)
ggsave("pa_rand_sh.pdf", width=4, height=3)

future.stocks = lapply(assets, function(sym) {
  monthlyReturn(na.omit(getSymbols(sym,
                                   src="yahoo",
                                   from="2017-10-01",
                                   to="2018-10-01",
                                   auto.assign=FALSE)))
})

future.df <- as.data.frame(do.call(merge.xts, future.stocks))
future.df <- cbind(month = rownames(future.df), future.df)
names(future.df)[2:17] <- assets
future.df

future.df[9,8] <- future.df[8,8]
future.df[14,8] <- -0.11346154
future.df <- na.omit(future.df)
future.df

future.mean.returns <- c(sapply(future.df[2:17], function(x) mean(x)))
sort(future.mean.returns, decreasing = TRUE)

mean.fut.mean.ret <- sum(future.mean.returns)/length(future.mean.returns)
mean.fut.mean.ret*100

getTotalReturn(future.mean.returns, as.numeric(unlist(hof[8])))*100 # Highest Sharpe Ratio PF

getTotalReturn(future.mean.returns, even.weighted.solution)*100 # Even weighted PF

getTotalReturn(future.mean.returns, as.numeric(unlist(hof[1])))*100 # Min Risk PF

getTotalReturn(future.mean.returns, as.numeric(unlist(hof[99])))*100 # Max Return PF

getTotalReturn(future.mean.returns, pa.weights)*100 # PA ROI

getTotalReturn(future.mean.returns, pa.rand.weights)*100 # PA Random

fh.fut.df <- future.df[1:9]
sh.fut.df <- future.df[sh.cols]

melted.fh.fut.df <- melt(fh.fut.df, id="month")
melted.sh.fut.df <- melt(sh.fut.df, id="month")


ggplot(data=melted.fh.fut.df, aes(x=as.Date(melted.fh.fut.df$month), y=melted.fh.fut.df$value,
                           group=melted.fh.fut.df$variable, color=melted.fh.fut.df$variable)) +
geom_line() + geom_point() +
xlab("Date") + ylab("Return (%)") + scale_y_continuous(labels=percent) +
scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
scale_color_discrete("Stock") + ggtitle("Stock Performance Oct 2017 - Oct 2018 (pt. 1)") +
theme(axis.text.x = element_text(angle = 45))
ggsave("stockperf2019.pdf", width=8, height=4)


ggplot(data=melted.sh.fut.df, aes(x=as.Date(melted.sh.fut.df$month), y=melted.sh.fut.df$value,
                           group=melted.sh.fut.df$variable, color=melted.sh.fut.df$variable)) +
geom_line() + geom_point() +
xlab("Date") + ylab("Return (%)") + scale_y_continuous(labels=percent) +
scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
scale_color_discrete("Stock") + ggtitle("Stock Performance Oct 2017 - Oct 2018 (pt. 2)") +
theme(axis.text.x = element_text(angle = 45))
ggsave("stockperf2_2019.pdf", width=8, height=4)