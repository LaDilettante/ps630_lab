#######
#Replication File for "Voter Buying: Shaping the Electorate Through Clientelism"
#######

#Exact version of R and packages used is listed at the bottom of this file

#To install packages, use the "install.packages" command. Example: "install.packages('reshape2')"
#Load packages
library(reshape2); library(ggplot2); library(Hmisc); library(scales); library(car); library(sandwich); library(AER); library(stargazer); library(rdd); library(reporttools); library(tidyr); library(dplyr)

#Set Directory to save tables and plots with all results
#Tables will be saved as Latex (tex) files
#Plots will be saved as PDFs
setwd('.')

#Set File Path to Replication Data
load("replication_data.RData")

###
#Functions
####


#Create function that calculates RD estimates and standard errors
rd.estim <- function(y, data, bw.dm, bw.ll, forcevar.name = "elecpopratio", cutpoint = .8){
  force.var <- data[, forcevar.name]
  data$y <- y
  data$force.var <- force.var
  data$cutpoint <- cutpoint
  data$above <- ifelse(data$force.var > cutpoint, 1, 0)
  data <- data[is.na(data$force.var) == FALSE, ]
  
  rfdm <- (lm(y ~ above , data = data, subset = (force.var > cutpoint - bw.dm) & (force.var < cutpoint + bw.dm)))
  rfdm.coef <- coef(rfdm)[2]
  rfdm.se <- sqrt(vcovHC(rfdm,type="HC1")[2,2])
  dm.baseline <- coef(rfdm)[1]
  
  ivdm <- ivreg(y ~ treat | above , data = data[(data$force.var > cutpoint - bw.dm) & (data$force.var < cutpoint + bw.dm), ])
  ivdm.coef <- coef(ivdm)[2]
  ivdm.se <- sqrt(vcovHC(ivdm, type = "HC1")[2,2])
  
  rfll <- lm(y ~ above * I(force.var - cutpoint), data = data, subset = (force.var > cutpoint - bw.ll) & (force.var < cutpoint + bw.ll))
  rfll.coef <- coef(rfll)[2]
  rfll.se <- sqrt(vcovHC(rfll,type="HC1")[2,2])
  ll.baseline <- coef(rfll)[1]
  
  ivll <- ivreg(y ~ treat * I(force.var - cutpoint) | above * I(force.var - cutpoint) , data = data[(data$force.var > cutpoint - bw.ll) & (data$force.var < cutpoint + bw.ll), ])
  ivll.coef <- coef(ivll)[2]
  ivll.se <- sqrt(vcovHC(ivll, type = "HC1")[2,2])
  
  force.var <- na.omit(force.var)
  
  results <- list(rfdm.coef = rfdm.coef, rfdm.se = rfdm.se, ivdm.coef = ivdm.coef, ivdm.se = ivdm.se, dm.n = sum((force.var > cutpoint - bw.dm) & (force.var < cutpoint + bw.dm)), dm.baseline = dm.baseline,
                  rfll.coef = rfll.coef, rfll.se = rfll.se, ivll.coef = ivll.coef, ivll.se = ivll.se, ll.n = sum((force.var > cutpoint - bw.ll) & (force.var < cutpoint + bw.ll)), ll.baseline = ll.baseline)
  results
}

#Create function that estimates heterogenous treatment effects with a binary covariate

dummy.interaction.results <- function(data, dv, inter){
  data$inter <- data[,inter]
  data$dm.weights <- 1 - abs((data$fv) / bw.dm)
  dm.model0 <- ivreg(as.formula(paste(dv, "~ treat | above", sep = "")), data = data, subset = inter == 0 & (fv < bw.dm & fv > -bw.dm))
  dm.model1 <- ivreg(as.formula(paste(dv, "~ treat | above", sep = "")), data = data, subset = inter == 1 & (fv < bw.dm & fv > -bw.dm))
  dm.model.inter <- ivreg(as.formula(paste(dv, "~ treat * ", inter, "| above * ", inter, sep = "")), data = data, subset = (fv < bw.dm & fv > -bw.dm))
  ll.model0 <- ivreg(as.formula(paste(dv, "~ treat * fv | above * fv", sep = "")), data = data, subset = inter == 0 & (fv < bw.ll & fv > -bw.ll))
  ll.model1 <- ivreg(as.formula(paste(dv, "~ treat * fv | above * fv", sep = "")), data = data, subset = inter == 1 & (fv < bw.ll & fv > -bw.ll))
  ll.model.inter <- ivreg(as.formula(paste(dv, "~ treat * fv *", inter, "| above * fv *", inter, sep = "")), data = data, subset = (fv < bw.ll & fv > -bw.ll))
  list(dm.coef0 = coef(dm.model0)[2], dm.se0  = sqrt(vcovHC(dm.model0, type = "HC1")[2,2]), dm.coef1 = coef(dm.model1)[2], dm.se1  = sqrt(vcovHC(dm.model1, type = "HC2")[2,2]),
       ll.coef0 = coef(ll.model0)[2], ll.se0  = sqrt(vcovHC(ll.model0, type = "HC1")[2,2]), ll.coef1 = coef(ll.model1)[2], ll.se1  = sqrt(vcovHC(ll.model1, type = "HC2")[2,2]),
       dm.inter.coef = coef(dm.model.inter)[4], dm.se.inter = sqrt(vcovHC(dm.model.inter, type = "HC1")[4,4]),
       ll.inter.coef = coef(ll.model.inter)[6], ll.se.inter = sqrt(vcovHC(ll.model.inter, type = "HC1")[6,6])
  )
}

#Create function for combining plots
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}


########
### Figure 3
########

#Reshape transfers data from wide to long format
netTransfers <- melt(data[ ,c('codigo', 'region', 'dist.statecap', 'pop2000', 'votingage.pop07', 'nb.estpop2008', grep('^transfers.inflow', names(data), value = TRUE), grep('^transfers.outflow', names(data), value = TRUE))], id.vars = c('codigo', 'region', 'dist.statecap', 'votingage.pop07', 'pop2000', 'nb.estpop2008'))
netTransfers <- cbind(netTransfers[, c(-7)], colsplit(netTransfers$variable, "_", names = c('variable', 'year')))
netTransfers <- dcast(netTransfers, codigo + region + votingage.pop07 + pop2000 + dist.statecap + nb.estpop2008 +  year ~ variable)
#Normalize transfers data by voting age population
netTransfers$inflowpct <- 100 * netTransfers$transfers.inflow / netTransfers$votingage.pop07
netTransfers$outflowpct <- 100 * netTransfers$transfers.outflow / netTransfers$votingage.pop07
netTransfers$close <- ifelse(netTransfers$dist.statecap < 50, 'Close', 'Far')
netTransfers$state_cap <- as.factor(ifelse(netTransfers$dist.statecap == 0 , "State Capital", "Non State-Capital"))
netTransfers$nettransfers <- netTransfers$transfers.inflow - netTransfers$transfers.outflow
netTransfers$nettransferspct <- 100 * netTransfers$nettransfers / netTransfers$votingage.pop07
netTransfers$importexportratio <- 100 * netTransfers$transfers.inflow / netTransfers$transfers.outflow
netTransfers_plot_data <- melt(netTransfers, id.vars =  c('state_cap', 'close', 'year', 'region'), measure.vars = c('importexportratio', 'nettransferspct', 'inflowpct'))

#Plot transfers data
pdf(file = "figure_3.pdf", height = 7, width = 11.5)
ggplot(netTransfers_plot_data[netTransfers_plot_data$close == "Close" & netTransfers_plot_data$variable == 'nettransferspct', ], aes(x = year, y= value, color = state_cap, lty = state_cap)) +
  stat_summary(fun.y = 'mean', geom = 'line', size = 2) +
  theme_bw() +
  stat_summary(fun.y = 'mean', geom = 'point', size = 2) +
  scale_x_continuous(breaks = seq(min(netTransfers_plot_data$year, na.rm = TRUE) - 1, max(netTransfers_plot_data$year, na.rm = TRUE), 2)) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  coord_cartesian(ylim = c(-1,3.5)) + xlab("Year") +
  ylab("Net Transfers (%)") + labs(color = "", linetype = "") +
  scale_colour_brewer(palette = "Set1")
dev.off()

#######
#Figure 4a
#######

#Plot First Stage
data$bin.equal2 <- NA
data$bin.equal2[data$fv >= 0] <- as.numeric(as.character(cut2(data$elecpopratio[data$fv >= 0], m = 50, levels.mean = TRUE))) * 100
data$bin.equal2[data$fv <= 0] <- as.numeric(as.character(cut2(data$elecpopratio[data$fv <= 0], m = 50, levels.mean = TRUE)))  * 100
pdf(file = "figure_4_left_panel.pdf",  height = 5, width = 7)
ggplot(data[data$elecpopratio <= 1.1 & data$elecpopratio > .30,], aes(y = treat, x = elecpopratio.pct))  +
  stat_smooth(method = loess, aes(fill = above), se = FALSE, size = 1.5, color = "black", span = .5, degree = 1) +
  scale_x_continuous("Electorate as % of Population") +  scale_y_continuous("Percent of Municipalities Audited", labels = percent, limits = c(0,1)) +
  stat_summary(aes(x = bin.equal2, color = above), alpha = .8, fun.y = mean, geom = "point", size = 3) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(-.01, 1))  +
  geom_vline(xintercep = 80, lty = 2)
dev.off()


#######
#Figure 4b
#######

#Plot relationship between forcing variable and population 
pdf(file = "figure_4_right_panel.pdf", height = 5, width = 7)
#Scatter Plot of Population and Elec/Pop Ratio
ggplot(data, aes(y = pop2008, x = elecpopratio.pct)) +
  scale_x_continuous("Electorate as a % of Population", limits = c(30, 110), breaks = c(25, 50, 75, 100, 125)) +
  scale_y_continuous("Population (2008)", trans = "log", breaks = c(1000, 5000, 25000, 100000, 300000, 1000000, 3000000, 9000000), labels = comma_format(digits = 10)) +
  geom_point(shape=16, size = .85)  +
  geom_vline(xintercep = 80, lty = 2) +
  theme_bw()
dev.off()

#######
#Figure 5
#######

#Set bandwidths
bw.dm <- .015
bw.ll <- .04

#Choose variables to check balance on
covar <- c("incumb.partywon.04", "turnout.04.perpop", "pct.70aboveelect.06", "pov2000", "theil2000", "incumb.reelected.in.2004", "rural.pop.07", "dist.statecap", "numrevisions.0205", "transfers.inflow0506", "transfers.outflow0506", "nb.transfers.outflow0506", "nb.transfers.inflow0506", "branconulo.perpop04", "nordeste", "sudeste", "num.neighbors", "pop2008", "electorate2006", "vote.margin.04.perpop", "gov_pref_sameparty", "pres_pref_sameparty", "pct.immigrants00", "pct.emigrants00", "num_vereadores04", 'electoratechange0708.perpop', 'muni_resources2006_percap', 'pub.employee.percap.2006', 'average_section_size', 'radio')
bal.data <- data[ , c("above", "elecpopratio", "treat", covar)]
bal.data$above <- as.numeric(as.character(bal.data$above))
#Supply median for missing values
bal.data <- data.frame(apply(bal.data, 2, function(x){x[is.na(x)] <- median(x, na.rm = TRUE); x}))
covar.names <- c("Incumbent Party Won (04)", "Percent Turnout (04)", "% of Electorate Above Age 70 (06)", "% Poverty (00)", "Theil Index (00)", "Incumbent Reelected (04)", "Rural Population (07)", "Distance to State Capital", "# of Previous Audits", "Transfers In (05-06)", "Transfers Out (05-06)", "Neighbor Munis Transfers Out (05-06)", "Neighbor Munis Transfers In (05-06)", "% Votes Blank and Invalid (04)", "Northeast Region", "Southeast Region", "# of Neighboring Municipalities", "Population (08)", "Electorate (06)", "Win Margin (04)", "Mayor, Governor Same Party", "Mayor, President Same Party", "Immigration Rate (00)", "Emigration Rate (00)", "Councilors in Municipality", 'Change in Electorate (02-07)', 'Log Municipal Budget Per Capita (06)', 'Public Employees Per Capita (06)', 'Average Precinct Size', 'Radio in Municipality')

#Calculate balance statistics
rd.balresults <- apply(bal.data[,covar], 2, function(x) rd.estim(y = x, data = bal.data, bw.dm = bw.dm, bw.ll = bw.ll))
rd.baldm.p <- pt(sapply(rd.balresults, function(x) abs(x$rfdm.coef/x$rfdm.se)), df = sum(data$elecpopratio > .8 - bw.dm & data$elecpopratio <.8 + bw.dm, na.rm = TRUE) ,lower.tail = FALSE)
rd.baldm.stddiff <- sapply(rd.balresults, function(x) x$rfdm.coef) / apply(bal.data[, covar], 2, sd)
rd.balll.p <- pt(sapply(rd.balresults, function(x) abs(x$rfll.coef/x$rfll.se)), df = sum(data$elecpopratio > .8 - bw.ll & data$elecpopratio <.8 + bw.ll, na.rm = TRUE) ,lower.tail = FALSE)
rd.balll.stddiff <- sapply(rd.balresults, function(x) x$rfll.coef) / apply(bal.data[, covar], 2, sd)
fs.bal.p <- apply(bal.data[ , covar], 2, function(x) t.test(x[bal.data$above == 1], x[bal.data$above == 0])$p.value)
fs.bal.stddiff <- apply(bal.data[, covar], 2, function(x) (mean(x[bal.data$above == 1]) - mean(x[bal.data$above == 0]))/sd(x))
bal.table <- data.frame(Covariate = rep(covar.names, 3), 
                        pvalue = c(rd.baldm.p, rd.balll.p, fs.bal.p), 
                        stddiff = c(rd.baldm.stddiff, rd.balll.stddiff, fs.bal.stddiff), 
                        Specification = c(rep("Difference-in-Means", length(covar.names)), rep("Local Linear", length(covar.names)), rep("Full Sample Difference-in-Means", length(covar.names))))
bal.table <- melt(bal.table, id.vars = c("Covariate", "Specification"))
names(bal.table)[3] <- "Statistic"
bal.table$Statistic <- ifelse(as.character(bal.table$Statistic) == "pvalue", "t-test p-value", "Standardized Difference")
bal.table$Covariate <- as.factor(as.character(bal.table$Covariate))
bal.table$Covariate <- factor(bal.table$Covariate, levels = rev(levels(bal.table$Covariate)))

bal.table$Statistic <- factor(bal.table$Statistic)
bal.table$vline <- ifelse(bal.table$Statistic == "t-test p-value", .05, 0)

pdf(file = "figure_5.pdf", height = 5, width = 9)
ggplot(bal.table, aes(y = Covariate, x = (value), color = Specification, shape = Specification)) +
  geom_point(alpha = .7, size = 2) + scale_colour_brewer(palette = "Set1", guide = guide_legend(nrow = 1)) +
  theme_bw() +
  scale_x_continuous("Statistic") +
  facet_grid(. ~ Statistic, scales = "free_x") +
  geom_vline(aes(xintercept = vline), linetype = 2, size = .2) +
  theme(legend.position="bottom")
dev.off()


#######
#Figure 6
#######

#Plot discontinuity at threshold

rdplot.electorate <- ggplot(data[data$elecpopratio < .8401 & data$elecpopratio > .759,] , aes(x=elecpopratio, y =  (electorate.perpop08 - electorate.perpop07)/100)) +
  geom_point(aes(color = above), alpha = .7, size = 1.5) +
  geom_vline(xintercept = .80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = 1, degree = 1, alpha = .8)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(.76, .84, .01), labels = percent_format()) +
  scale_y_continuous("Change in Registration (% of Population)", limits= c(-.20, .20), label = percent) +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1")  +
  theme(axis.text.x = element_text(size=rel(1.2))) +
  theme(axis.text.y = element_text(size=rel(1.2)) , axis.title.x = element_text(vjust = -.5), axis.title.y = element_text(vjust = .4))


num.obsperbin <- 20
data$bin.equal <- NA
data$bin.equal[data$fv >= 0] <- as.numeric(as.character(cut2(data$elecpopratio[data$fv >= 0], m = num.obsperbin, levels.mean = TRUE))) * 100
data$bin.equal[data$fv <= 0] <- as.numeric(as.character(cut2(data$elecpopratio[data$fv <= 0], m = num.obsperbin, levels.mean = TRUE)))  * 100

rdplot.incumbwins <- ggplot(data[data$elecpopratio < .8401 & data$elecpopratio > .759 & data$incumb.reelected.in.2004 == 0,], aes(y = incumbent.won08, x = elecpopratio.pct / 100)) +
  geom_vline(xintercept = .80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = 1, degree = 1) +
  scale_x_continuous("Electorate as % of Population", breaks = seq(.76, .84, .01), label = percent, limits = c(.76,.84)) +
  scale_y_continuous("Incumbent Reelection Rate", labels = percent) +
  stat_summary(aes(x = bin.equal/100, color = bin.equal/100 > .80), fun.y = mean, geom = "point", size = 3) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(.2, .8), xlim = c(.759, .8401)) +
  theme(axis.text.x = element_text(size=rel(1.2))) +
  theme(axis.text.y = element_text(size=rel(1.2)) , axis.title.x = element_text(vjust = -.5), axis.title.y = element_text(vjust = .4))


pdf(file = "figure_6a.pdf", height = 5, width = 9)
rdplot.electorate
dev.off()

pdf(file = "figure_6b.pdf", height = 5, width = 9)
rdplot.incumbwins
dev.off()

#######
#Table 1
#######

#Estimate main results of paper

main.results <- list()
main.results$electorate.results <- rd.estim(y = data$electorate.perpop08 - data$electorate.perpop07, data = data, bw.dm = bw.dm, bw.ll = bw.ll)
main.results$turnout.results <- rd.estim(y = data$turnout.08.perpop - data$turnout.04.perpop, data = data, bw.dm = bw.dm, bw.ll = bw.ll)
main.results$incumbwins.results <- rd.estim(y = data$incumbent.won08[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = bw.dm, bw.ll = bw.ll)
main.results$incumbpartywins.results <- rd.estim(y = data$incumb.partywon.08, data = data, bw.dm = bw.dm, bw.ll = bw.ll)
main.results$runreelect.results <- rd.estim(y = data$run.reelec[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = bw.dm, bw.ll = bw.ll)
main.results$incumbvoteshare.results <- rd.estim(y = data$incumbent.vote.2008.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1] - data$winner.vote.2004.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1], data = data[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1, ], bw.dm = bw.dm, bw.ll = bw.ll)


#Create table with estimates
main_results_table <- cbind('$\\hat \\tau_A$' = c(sapply(main.results, function(x) x$rfll.coef),
                                                  sapply(main.results, function(x) x$rfdm.coef)),
                            'SE$_{\\tau_A}$' = c(sapply(main.results, function(x) x$rfll.se),
                                                 sapply(main.results, function(x) x$rfdm.se)),
                            '$\\hat \\tau_R$' = c(sapply(main.results, function(x) x$ivll.coef),
                                                  sapply(main.results, function(x) x$ivdm.coef)),
                            'SE$_{\\tau_R}$' = c(sapply(main.results, function(x) x$ivll.se),
                                                 sapply(main.results, function(x) x$ivdm.se)),
                            'Baseline' = c(sapply(main.results, function(x) x$ll.baseline),
                                           sapply(main.results, function(x) x$dm.baseline)),
                            '$n$' = c(sapply(main.results, function(x) x$ll.n),
                                      sapply(main.results, function(x) x$dm.n)))
main_results_caption <- "Effect of Voter Audits on Brazil's 2008 Municipal Elections (RDD Estimates)"
main_results_bottomnote <- "$\\hat \\tau_A$ is the local average effect of a municipality having an electorate more than 80\\% of its population (``reduced form'') and SE$_{\\tau_A}$ is its associated standard error. $\\hat \\tau_R$ is the estimated local average effect of the audit and SE$_{\\tau_R}$ is its  standard error. Standard errors are ``robust''.  Baseline estimates value of the dependent variable among controls at $E_i = 80$ For  \\textit{Incumbent Reelected}, \\textit{Incumbent Runs for Reelection}, and \\textit{Change in Incumbent Vote Share} variables, sample includes only municipalities with incumbents eligible for re-election."

#Save latex table with main results
main_results_tex <- latex(main_results_table, cdec = c(2, 2, 2, 2, 2, 0), col.just = rep('c', 6),
                          rowname = c('Change in Registration (\\%)', 'Change in Turnout (\\%)', 'Incumbent Reelected', 'Incumbent Party Reelected', 'Incumbent Runs for Reelection', 'Change in Incumbent Vote Share (\\%)'),
                          rowlabel = 'Outcome', ctable = TRUE, booktabs = FALSE, dcolumn = TRUE,
                          n.rgroup = c(6, 6), rgroup = c('Local Linear Specification', 'Difference-in-Means Specification'),
                          caption = main_results_caption, insert.bottom = main_results_bottomnote,
                          file ="table_1.tex")


#######
#Figure 7
#######

#Estimate heterogenous treatment effects

hetero.incumbwins <- list()
hetero.incumbwins$transfersin <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.transfersin.dummy")
hetero.incumbwins$nbtransfers <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.nbtransfers.dummy")
hetero.incumbwins$nboutin <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0 , ], dv = "incumbent.won08", inter = "inter.nbout.inflow.dummy")

hetero.registration <- list()
hetero.registration$transfersin <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.transfersin.dummy")
hetero.registration$nbtransfers <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.nbtransfers.dummy")
hetero.registration$nboutin <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0 , ], dv = "electoratechange0708.perpop", inter = "inter.nbout.inflow.dummy")

#Create dataset to use for plotting heterogenous treatment efefcts

hetero.incumb.plotdata <- data.frame(Estimate = unlist(sapply(hetero.incumbwins[c(1:3)], function(x) x[grep("coef", names(x))])),
                                     SE = unlist(sapply(hetero.incumbwins[c(1:3)], function(x) x[grep("se", names(x))]))
                                     ,
                                     Strata = rep(c("Below\nMedian",  "Above\nMedian",  "Below\nMedian", "Above\nMedian", "Difference", "Difference"), 3)
                                     ,
                                     Specification = rep(c(rep("Difference-in-Means", 2), rep("Local Linear", 2), "Difference-in-Means", "Local Linear"), 3),
                                     Covariate = c(rep("(d) Transfers Into Municipality", 6), rep("(e) Transfers Out of Neighbors", 6), rep("(f) Neighbors' Outflows & Transfers In", 6)))
hetero.incumb.plotdata$Strata <- factor(hetero.incumb.plotdata$Strata, levels = c("Below\nMedian", "Above\nMedian", "Difference"))
hetero.incumb.plotdata$Covariate <- factor(hetero.incumb.plotdata$Covariate, levels = c("(d) Transfers Into Municipality", "(e) Transfers Out of Neighbors", "(f) Neighbors' Outflows & Transfers In"))


#Plot heterogenous treatment effects

hetero.incumb.plot <- ggplot(hetero.incumb.plotdata[hetero.incumb.plotdata$Covariate != "Population", ], aes(x = Strata, y = Estimate, color = Specification, shape = Specification)) +
  #  geom_pointrange(aes(ymin = Estimate - 1.65 * SE, ymax = Estimate + 1.65 * SE), size = .8, position=position_dodge(width=0.4)) +
  facet_wrap( ~ Covariate, ncol  = 3, scales = "free_x") + scale_colour_brewer(palette = "Set1", guide = guide_legend(nrow = 1)) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position="bottom") +
  scale_y_continuous("Effect on Incumbent Reelection\n(proportion)") +
  geom_pointrange(aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 * SE), size = .8, position=position_dodge(width=0.4)) +
  xlab("") + 
  guides(color = FALSE, shape = FALSE)

#Create dataset to use for plotting heterogenous treatment efefcts


hetero.registration.plotdata <- data.frame(Estimate = unlist(sapply(hetero.registration[c(1:3)], function(x) x[grep("coef", names(x))])),
                                           SE = unlist(sapply(hetero.registration[c(1:3)], function(x) x[grep("se", names(x))]))
                                           ,
                                           Strata = rep(c("Below\nMedian",  "Above\nMedian",  "Below\nMedian", "Above\nMedian", "Difference", "Difference"), 3)
                                           ,
                                           Specification = rep(c(rep("Difference-in-Means", 2), rep("Local Linear", 2), "Difference-in-Means", "Local Linear"), 3),
                                           Covariate = c(rep("(a) Transfers Into Municipality", 6), rep("(b) Transfers Out of Neighbors", 6), rep("(c) Neighbors' Outflows & Transfers In", 6))
) 
hetero.registration.plotdata$Strata <- factor(hetero.registration.plotdata$Strata, levels = c("Below\nMedian", "Above\nMedian", "Difference"))
hetero.registration.plotdata$Covariate <- factor(hetero.registration.plotdata$Covariate, levels = c("(a) Transfers Into Municipality", "(b) Transfers Out of Neighbors", "(c) Neighbors' Outflows & Transfers In"))

#Plot heterogenous treatment effects

hetero.reg.plot <- ggplot(hetero.registration.plotdata[hetero.registration.plotdata$Covariate != "Population", ], aes(x = Strata, y = Estimate, color = Specification, shape = Specification)) +
  #  geom_pointrange(aes(ymin = Estimate - 1.65 * SE, ymax = Estimate + 1.65 * SE), size = .8, position=position_dodge(width=0.4)) +
  facet_wrap( ~ Covariate, ncol  = 3, scales = "free_x") + scale_colour_brewer(palette = "Set1", guide = guide_legend(nrow = 1)) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position="bottom") + scale_y_continuous("Effect on Registration\n(%)") +
  geom_pointrange(aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 * SE), size = .8, position=position_dodge(width=0.4)) +
  xlab("") + 
  guides(color = FALSE, shape = FALSE)

pdf(file = "figure_7_top_panel.pdf", height = 4, width = 8.2)
hetero.reg.plot
dev.off()

pdf(file = "figure_7_bottom_panel.pdf", height = 4, width = 8.2)
hetero.incumb.plot
dev.off()

#######
#Figure 8
#######

#Create dataset to use for plotting heterogenous treatment efefcts

hetero.incumbwins$emigration <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.emigration.dummy")
hetero.incumbwins$elderly <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.elderly.dummy")
hetero.incumbwins$branconulo <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.branconulo.dummy")
hetero.incumbwins$faltosos <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.faltosos.dummy")
hetero.incumbwins$transfers.out9707 <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.transfersout9707.dummy")
hetero.incumbwins$transfers.out0207 <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "inter.transfersout0207.dummy")
hetero.incumbwins$radio <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "incumbent.won08", inter = "radio")

hetero.registration$emigration <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.emigration.dummy")
hetero.registration$elderly <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.elderly.dummy")
hetero.registration$branconulo <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.branconulo.dummy")
hetero.registration$transfers.out9707 <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.transfersout9707.dummy")
hetero.registration$transfers.out0207 <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "inter.transfersout0207.dummy")
hetero.registration$radio <- dummy.interaction.results(data = data[data$incumb.reelected.in.2004 == 0, ], dv = "electoratechange0708.perpop", inter = "radio")


hetero.alternative.plotdata <- data.frame(Estimate = unlist(sapply(hetero.incumbwins[c('emigration', 'elderly', 'branconulo', 'radio')], function(x) x[grep("coef", names(x))])),
                                          SE = unlist(sapply(hetero.incumbwins[c('emigration', 'elderly', 'branconulo', 'radio')], function(x) x[grep("se", names(x))])),
                                          Strata = rep(c("Below\nMedian",  "Above\nMedian",  "Below\nMedian", "Above\nMedian", "Difference", "Difference"), 4),
                                          Specification = rep(c(rep("Difference-in-Means", 2), rep("Local Linear", 2), "Difference-in-Means", "Local Linear"), 4),
                                          Covariate = c(rep("(a) Emigration", 6), rep("(b) Elderly Citizens", 6), rep("(c) Blank and Invalid Votes", 6),  rep('(d) Radio Station in Municipality', 6)))

hetero.alternative.plotdata$Strata <- factor(hetero.alternative.plotdata$Strata, levels = c("Below\nMedian", "Above\nMedian", "Difference"))
hetero.alternative.plotdata$Covariate <- factor(as.character(hetero.alternative.plotdata$Covariate), levels = c('(a) Emigration', "(b) Elderly Citizens", "(c) Blank and Invalid Votes", '(d) Radio Station in Municipality'))

#Plot heterogenous treatment effects to test alternative mechanisms

alt.incumb.plot <- ggplot(hetero.alternative.plotdata[hetero.alternative.plotdata$Covariate %in% c('(a) Emigration', "(b) Elderly Citizens", "(c) Blank and Invalid Votes", "(d) Radio Station in Municipality"), ],
                          aes(x = Strata, y = Estimate, color = Specification, shape = Specification))+
  facet_wrap( ~ Covariate, ncol  = 4, scales = "free_x") +
  scale_colour_brewer(palette = "Set1", guide = guide_legend(nrow = 1)) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position="bottom") +
  scale_y_continuous("Audit Effect on Incumbent Reelection") +
  geom_pointrange(aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 * SE), size = .8, position=position_dodge(width=0.4)) +
  xlab("") +
  guides(color = FALSE, shape = FALSE)

hetero.reg.alternative.plotdata <- data.frame(Estimate = unlist(sapply(hetero.registration[c('emigration', 'elderly', 'branconulo',  'radio')], function(x) x[grep("coef", names(x))])),
                                              SE = unlist(sapply(hetero.registration[c('emigration', 'elderly', 'branconulo',  'radio')], function(x) x[grep("se", names(x))])),
                                              Strata = rep(c("Below\nMedian",  "Above\nMedian",  "Below\nMedian", "Above\nMedian", "Difference", "Difference"), 4),
                                              Specification = rep(c(rep("Difference-in-Means", 2), rep("Local Linear", 2), "Difference-in-Means", "Local Linear"), 4),
                                              Covariate = c(rep("(a) Emigration", 6), rep("(b) Elderly Citizens", 6), rep("(c) Blank and Invalid Votes", 6),  rep('(d) Radio Station in Municipality', 6)))

hetero.reg.alternative.plotdata$Strata <- factor(hetero.reg.alternative.plotdata$Strata, levels = c("Below\nMedian", "Above\nMedian", "Difference"))
hetero.reg.alternative.plotdata$Covariate <- factor(as.character(hetero.reg.alternative.plotdata$Covariate), levels = c('(a) Emigration', '(a) Transfers Out (02-07)',"(b) Elderly Citizens", "(c) Blank and Invalid Votes", '(b) Transfers Out (97-07)', '(c) Missing Voters', '(d) Radio Station in Municipality'))

alt.reg.plot <- ggplot(hetero.reg.alternative.plotdata[hetero.reg.alternative.plotdata$Covariate %in% c('(a) Emigration', "(b) Elderly Citizens", "(c) Blank and Invalid Votes", "(d) Radio Station in Municipality"), ],
                       aes(x = Strata, y = Estimate, color = Specification, shape = Specification))+
  facet_wrap( ~ Covariate, ncol  = 4, scales = "free_x") +
  scale_colour_brewer(palette = "Set1", guide = guide_legend(nrow = 1)) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position="bottom") +
  scale_y_continuous("Audit Effect on Registration", breaks=c(0, -5, -10, -15, -20, -25), labels = c("0", "-5.0", "-10.0", "-15.0", "-20", "-25")) +
  geom_pointrange(aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 * SE), size = .8, position=position_dodge(width=0.4)) +
  xlab("") +
  guides(color = FALSE, shape = FALSE)


pdf(file = "figure_8_top_panel.pdf", height = 4, width = 8.2)
alt.reg.plot
dev.off()

pdf(file = "figure_8_bottom_panel.pdf", height = 4, width = 8.2)
alt.incumb.plot
dev.off()


#######
#Table 2
#######

#Estimate Bahia Results

reg_mod1 <- lm((pct_revision_cancellations) ~ (nbtransfer_dummy)  + rural + pct_incumb_party_filiados + pct_totalafilliated   + log(reg_voters_2004) + mun - 1, data = secao_incumb)
reg_mod1$se <- sqrt(diag(vcovHC(reg_mod1, type = "HC1")))
reg_mod2 <- lm((pct_revision_cancellations) ~ (log_nbtransfer_pctreg) + rural + pct_incumb_party_filiados + pct_totalafilliated   + log(reg_voters_2004) + mun - 1, data = secao_incumb)
reg_mod2$se <- sqrt(diag(vcovHC(reg_mod2, type = "HC1")))

incumb_mod1 <- lm(incumbchange_pct_0408 ~  pct_revision_cancellations   + rural + pct_incumb_party_filiados + pct_totalafilliated + log(reg_voters_2004) + mun - 1, data = secao_incumb)
incumb_mod1$se <- sqrt(diag(vcovHC(incumb_mod1, type = "HC1")))
incumb_mod2 <- lm(incumbchange_pct_0408 ~  pct_revision_cancellations + pct_revision_cancellations * nbtransfer_dummy  + rural + pct_incumb_party_filiados + pct_totalafilliated + log(reg_voters_2004) + mun - 1, data = secao_incumb)
incumb_mod2$se <- sqrt(diag(vcovHC(incumb_mod2, type = "HC1")))
incumb_mod3 <- lm(incumbchange_pct_0408 ~  pct_revision_cancellations + pct_revision_cancellations * log_nbtransfer_pctreg  + rural + pct_incumb_party_filiados + pct_totalafilliated + log(reg_voters_2004) + mun -1, data = secao_incumb)
incumb_mod3$se <- sqrt(diag(vcovHC(incumb_mod2, type = "HC1")))

#Create Table with Bahia Results

stargazer(reg_mod1, reg_mod2, incumb_mod1, incumb_mod2, incumb_mod3,
          omit = c("mun", "rural|filiados|total|reg_voters"),
          type = "latex",
          omit.labels = c("Muni Fixed Effects", "Controls"),
          covariate.labels= c("Transfers from Neighbors (Above Median)", "Transfers from Neighbors (Logged)", "\\% Audit Removals", "\\% Audit Removals x Transfers from Neighbors (Above Median)", "\\% Audit Removals x Transfers from Neighbors (Logged)"),
          dep.var.labels = c("Audit Removals", "Change in Incumbent Vote Share (04-08)"),
          keep.stat = "n",
          digits = 2, 
          notes = c("Robust standard errors in parentheses."),
          notes.align = "r",
          title = c("Precinct Level Results for the State of Bahia"),
          label = "tab:precinct_results",
          font.size = "footnotesize",
          out = 'table_2.tex')

#################
#Online Appendix
#################


#######
#Figure A1
#######

#Plot Results of McCrary Density Test

pdf(file = "figure_A1.pdf", height = 5, width = 9)
par(mfrow=c(1,2))
mccrary_pvalue <- DCdensity(runvar = data$elecpopratio[data$elecpopratio > .4  & data$elecpopratio < 1.2 ], .8)
abline(v = .8, lwd = 2)
title(main = "Full Sample", ylab = "Frequency Count", xlab = "Electorate to Population Ratio")
mccrary_pvalue <- DCdensity(runvar = data$elecpopratio[data$elecpopratio > .4  & data$elecpopratio < 1.2 &  data$incumb.reelected.in.2004 == 1], .8)
abline(v = .8, lwd = 2)
title(main = "Municipalities with Mayors\nEligible for Re-election",  ylab = "Frequency Count", xlab = "Electorate to Population Ratio")
dev.off()


#######
#Table A1
#######

#Estimate First Stage

first_stage <- list()
first_stage$dm <- lm(treat ~ above, data = data, subset = abs(fv) < bw.dm)
first_stage$dm$se <- sqrt(diag(vcovHC(first_stage$dm, type = "HC1")))
first_stage$ll <- lm(treat ~ above * fv, data = data, subset = abs(fv) < bw.ll)
first_stage$ll$se <- sqrt(diag(vcovHC(first_stage$ll, type = "HC1")))

first_stage_table <- stargazer(first_stage$dm,
                               first_stage$ll,
                               se = list(first_stage$dm$se,
                                         first_stage$ll$se),
                               dep.var.labels = c("Revision"),
                               covariate.labels = c("Above .8", "Baseline"),
                               omit = "fv",
                               omit.stat = c("rsq", "adj.rsq", "f", "ser"),
                               digits = 2,
                               float = FALSE)
first_stage_table <- c(first_stage_table[1:18],
                       "Specification & Difference-in-Means & Local Linear \\\\",
                       first_stage_table[19:24])
cat(paste(first_stage_table, collapse = "\n"), file = 'table_A1.tex')

#######
#Figure A2
#######

# Show Robustness to alternative bandwidths (local linear)

robust.bw <- seq(.005, .06, by = .005)
robust.incumbwins <- lapply(robust.bw, function(x) rd.estim(y = data$incumbent.won08[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = x, bw.ll = x))
robust.electorate <- lapply(robust.bw, function(x) rd.estim(y =  data$electorate.perpop08 - data$electorate.perpop07, data = data, bw.dm = x, bw.ll = x))
robust.turnout <- lapply(robust.bw, function(x) rd.estim(y =  data$turnout.08.perpop - data$turnout.04.perpop, data = data, bw.dm = x, bw.ll = x))
robust.incumbparty <- lapply(robust.bw, function(x) rd.estim(y =  data$incumb.partywon.08, data = data, bw.dm = x, bw.ll = x))
robust.runreelection <- lapply(robust.bw, function(x) rd.estim(y =  data$run.reelec[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = x, bw.ll = x))
robust.incumbvoteshare <- lapply(robust.bw, function(x) rd.estim(y =  data$incumbent.vote.2008.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1] - data$winner.vote.2004.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1], data = data[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1, ], bw.dm = x, bw.ll = x))


robust.plot.data <- rbind(data.frame(bw = robust.bw, est = sapply(robust.incumbwins, function(x) x$ivll.coef), se = sapply(robust.incumbwins, function(x) x$ivll.se), Variable = "Incumbent Reelected"),
                          data.frame(bw = robust.bw, est = sapply(robust.electorate, function(x) x$ivll.coef), se = sapply(robust.electorate, function(x) x$ivll.se), Variable = "Change in Electorate (%)"),
                          data.frame(bw = robust.bw, est = sapply(robust.turnout, function(x) x$ivll.coef), se = sapply(robust.turnout, function(x) x$ivll.se), Variable = "Change in Turnout (%)"),
                          data.frame(bw = robust.bw, est = sapply(robust.incumbparty, function(x) x$ivll.coef), se = sapply(robust.incumbparty, function(x) x$ivll.se), Variable = "Incumbent Party Reelected"),
                          data.frame(bw = robust.bw, est = sapply(robust.runreelection, function(x) x$ivll.coef), se = sapply(robust.runreelection, function(x) x$ivll.se), Variable = "Incumbent Runs for Reelection"),
                          data.frame(bw = robust.bw, est = sapply(robust.incumbvoteshare, function(x) x$ivll.coef), se = sapply(robust.incumbvoteshare, function(x) x$ivll.se), Variable = "Change in Incumbent Vote Share (%)")
)
robust.plot.data$Variable <- factor(robust.plot.data$Variable, levels = c("Incumbent Reelected", "Incumbent Party Reelected", "Incumbent Runs for Reelection", "Change in Incumbent Vote Share (%)", "Change in Electorate (%)", "Change in Turnout (%)"))
robust.plot.data$Preferred <- ifelse(robust.plot.data$bw == bw.ll, 1, 0)

pdf(file = "figure_A2.pdf", height = 10, width = 8.5)
ggplot(robust.plot.data, aes(x = bw * 100, y = est)) +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se,  color = Preferred)) +
  facet_wrap(~ Variable, ncol = 1, scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  scale_x_continuous("Bandwidth", breaks = robust.bw * 100) +
  scale_y_continuous("Estimate") +
  theme(legend.position="none")
dev.off()

#######
#Figure A3
#######

# Show Robustness to alternative bandwidths (difference in means)

# Collect results with alternative bandwidths
robust.bw <- seq(.0025, .0275, by = .0025)
robust.incumbwins <- lapply(robust.bw, function(x) rd.estim(y = data$incumbent.won08[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = x, bw.ll = x))
robust.electorate <- lapply(robust.bw, function(x) rd.estim(y =  data$electorate.perpop08 - data$electorate.perpop07, data = data, bw.dm = x, bw.ll = x))
robust.turnout <- lapply(robust.bw, function(x) rd.estim(y =  data$turnout.08.perpop - data$turnout.04.perpop, data = data, bw.dm = x, bw.ll = x))
robust.incumbparty <- lapply(robust.bw, function(x) rd.estim(y =  data$incumb.partywon.08, data = data, bw.dm = x, bw.ll = x))
robust.runreelection <- lapply(robust.bw, function(x) rd.estim(y =  data$run.reelec[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = x, bw.ll = x))
robust.incumbvoteshare <- lapply(robust.bw, function(x) rd.estim(y =  data$incumbent.vote.2008.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1] - data$winner.vote.2004.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1], data = data[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1, ], bw.dm = x, bw.ll = x))

#Plot alternative bandwidth results
robust.plot.data <- rbind(data.frame(bw = robust.bw, est = sapply(robust.incumbwins, function(x) x$ivdm.coef), se = sapply(robust.incumbwins, function(x) x$ivdm.se), Variable = "Incumbent Reelected"),
                          data.frame(bw = robust.bw, est = sapply(robust.electorate, function(x) x$ivdm.coef), se = sapply(robust.electorate, function(x) x$ivdm.se), Variable = "Change in Electorate (%)"),
                          data.frame(bw = robust.bw, est = sapply(robust.turnout, function(x) x$ivdm.coef), se = sapply(robust.turnout, function(x) x$ivdm.se), Variable = "Change in Turnout (%)"),
                          data.frame(bw = robust.bw, est = sapply(robust.incumbparty, function(x) x$ivdm.coef), se = sapply(robust.incumbparty, function(x) x$ivdm.se), Variable = "Incumbent Party Reelected"),
                          data.frame(bw = robust.bw, est = sapply(robust.runreelection, function(x) x$ivdm.coef), se = sapply(robust.runreelection, function(x) x$ivdm.se), Variable = "Incumbent Runs for Reelection"),
                          data.frame(bw = robust.bw, est = sapply(robust.incumbvoteshare, function(x) x$ivdm.coef), se = sapply(robust.incumbvoteshare, function(x) x$ivdm.se), Variable = "Incumbent Vote Share (%)")
)
robust.plot.data$Variable <- factor(robust.plot.data$Variable, levels = c("Incumbent Reelected", "Incumbent Party Reelected", "Incumbent Runs for Reelection", "Incumbent Vote Share (%)", "Change in Electorate (%)", "Change in Turnout (%)"))
robust.plot.data$Preferred <- ifelse(as.character(robust.plot.data$bw) == bw.dm, 1, 0)
robust.plot.data$Preferred[robust.plot.data$Preferred == bw.dm]

pdf(file = "figure_A3.pdf", height = 10, width = 8.5)
ggplot(robust.plot.data, aes(x = bw * 100, y = est)) +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se,  color = Preferred)) +
  facet_wrap(~ Variable, ncol = 1, scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  scale_x_continuous("Bandwidth", breaks = robust.bw * 100) +
  scale_y_continuous("Estimate") +
  theme(legend.position="none")
dev.off()

#######
#Table A2
#######

#Descriptive Statistics
descriptive.var <- c("electorate.perpop08", "turnout.08.perpop", "incumb.partywon.08", "run.reelec","incumbent.vote.2008.perpop",  "branconulo.perpop08", "incumb.partywon.04", "turnout.04.perpop", "pct_age70aboveelect.06", "pov2000", "theil2000", "incumb.reelected.in.2004", "rural.pop.07", "dist.statecap", "numrevisions.0205", "transfers.inflow0506", "transfers.outflow0506", "nb.transfers.outflow0506", "nb.transfers.inflow0506", "branconulo.perpop04", "nordeste", "sudeste", "num.neighbors", "pop2008", "electorate.perpop06", "vote.margin.04.perpop")

descriptive.data <- data[, (descriptive.var)]
descriptive.data$transfers.inflow0506 <- descriptive.data$transfers.inflow0506 * 100

names(descriptive.data) <- c("Electorate as \\% of Pop", "Turnout (08) as \\% of Pop", "Incumbent Party Won (08)", "Ran for Re-election (08)", "Incumbent Vote (08) as \\% of Pop",  "Blank and Invalid (08) as \\% of  Pop", "Incumbent Party Won (04)", "Turnout (04) as \\% of Pop", "Electorate Above 70 (06) as \\% of Pop", "Pct Poverty (00)", "Theil Index (00)", "Incumbent Re-elected (04)", "Rural Population (07)", "Distance to State Capital", "Num of Previous Audits", "Transfers In (05-06)", "Transfers Out (05-06)", "Neighbor Munis Transfers Out (05-06)", "Neighbor Munis Transfers In (05-06)", "Blank and Invalid (04) as \\% of  Pop", "Northeast Region", "Southeast Region", "Num of Neighboring Municipalities", "Population (08)", "Electorate (06) as \\% of  Pop", "Win Margin (04) as \\% of Pop")
data$Revision <- ifelse(data$treat == 1, "Revision", "No Revision")

tableContinuous(descriptive.data[, ], stats = list("n",  "\\textbf{Mean}" = function(x){return(mean(x))}, "\\textbf{Median}" = function(x){return(median(x))}, "\\textbf{SD}" = function(x){return(sd(x))}),
                    prec = 2, cap = "Descriptive Statistics", longtable = FALSE, file = "table_A2.tex")

#######
#Figure A4
#######

#Plot histogram of forcing variable

data$Status <- as.factor(ifelse(data$treat == 1, "Audit", "No Audit"))
data$elecpopratio.pct <- data$elecpopratio * 100

pdf(file = "figure_A4.pdf", height = 5, width = 7)
ggplot(data[data$elecpopratio >.3 & data$elecpopratio<1.3,], aes(x=elecpopratio.pct)) +
  geom_bar(binwidth = 1, aes(fill = Status), color = "black", size = .1) +
  geom_vline(xintercept = 80, col = "black") +
  scale_fill_brewer(palette = "Set3")  +
  theme_bw() +
  scale_x_continuous("Electorate as % of Population")
dev.off()


#######
#Table A3
#######

#Estimate Alternative Specifications (5th order polynomial and local linear with covariates)

alt_spec <- list()
alt_spec$poly_electorate <- lm(electorate.perpop08 - electorate.perpop07 ~ poly(fv, degree = 5, raw = TRUE) * above, data = data, subset =  abs(fv) < .3)
alt_spec$poly_electorate$se <- sqrt(diag(vcovHC(alt_spec$poly_electorate, type = "HC1")))
alt_spec$cov_electorate <- (lm(electorate.perpop08 - electorate.perpop07 ~ poly(fv, degree = 1, raw = TRUE) * above + uf + turnout.04.perpop + vote.margin.04.perpop + pov2000 + incumb.partywon.04 + theil2000 + rural.pop.07, data = data, subset = abs(fv) < bw.ll, weights = 1 - abs(fv/bw.ll)))
alt_spec$cov_electorate$se <- sqrt(diag(vcovHC(alt_spec$cov_electorate, type = "HC1")))

alt_spec$poly_turnout <- lm(turnout.08.perpop - turnout.04.perpop ~ poly(fv, degree = 5, raw = TRUE) * above, data = data, subset =  abs(fv) < .3)
alt_spec$poly_turnout$se <- sqrt(diag(vcovHC(alt_spec$poly_turnout, type = "HC1")))
alt_spec$cov_turnout <- (lm(turnout.08.perpop - turnout.04.perpop ~ poly(fv, degree = 1, raw = TRUE) * above + uf + turnout.04.perpop + vote.margin.04.perpop + pov2000 + incumb.partywon.04 + theil2000 + rural.pop.07, data = data, subset = abs(fv) < bw.ll, weights = 1 - abs(fv/bw.ll)))
alt_spec$cov_turnout$se <- sqrt(diag(vcovHC(alt_spec$cov_turnout, type = "HC1")))

alt_spec$poly_incumb <- (lm(incumbent.won08 ~ poly(fv, degree = 5, raw = TRUE) * above, data = data, subset = incumb.reelected.in.2004 == 0 & abs(fv) < .3))
alt_spec$poly_incumb$se <- sqrt(diag(vcovHC(alt_spec$poly_incumb, type = "HC1")))
alt_spec$cov_incumb <-(lm(incumbent.won08 ~ poly(fv, degree = 1, raw = TRUE) * above + uf + turnout.04.perpop + vote.margin.04.perpop + pov2000 + incumb.partywon.04 + theil2000 + rural.pop.07, data = data, subset = incumb.reelected.in.2004 == 0 & abs(fv) < bw.ll, weights = 1 - abs(fv/bw.ll)))
alt_spec$cov_incumb$se <- sqrt(diag(vcovHC(alt_spec$cov_incumb, type = "HC1")))

alt_spec$poly_partyincumb <- (lm(incumb.partywon.08 ~ poly(fv, degree = 5, raw = TRUE) * above, data = data, subset = abs(fv) < .3))
alt_spec$poly_partyincumb$se <- sqrt(diag(vcovHC(alt_spec$poly_partyincumb, type = "HC1")))
alt_spec$cov_partyincumb <-(lm(incumb.partywon.08 ~ poly(fv, degree = 1, raw = TRUE) * above + uf + turnout.04.perpop + vote.margin.04.perpop + pov2000 + incumb.partywon.04 + theil2000 + rural.pop.07, data = data, subset = abs(fv) < bw.ll, weights = 1 - abs(fv/bw.ll)))
alt_spec$cov_partyincumb$se <- sqrt(diag(vcovHC(alt_spec$cov_partyincumb, type = "HC1")))

alt_spec$poly_run <- (lm(run.reelec ~ poly(fv, degree = 5, raw = TRUE) * above, data = data, subset = incumb.reelected.in.2004 == 0 & abs(fv) < .3))
alt_spec$poly_run$se <- sqrt(diag(vcovHC(alt_spec$poly_run, type = "HC1")))
alt_spec$cov_run <-(lm(run.reelec ~ poly(fv, degree = 1, raw = TRUE) * above + uf + turnout.04.perpop + vote.margin.04.perpop + pov2000 + incumb.partywon.04 + theil2000 + rural.pop.07, data = data, subset = incumb.reelected.in.2004 == 0 & abs(fv) < bw.ll, weights = 1 - abs(fv/bw.ll)))
alt_spec$cov_run$se <- sqrt(diag(vcovHC(alt_spec$cov_run, type = "HC1")))

alt_spec$poly_voteshare <- (lm(incumbent.vote.2008.perpop - winner.vote.2004.perpop ~ poly(fv, degree = 5, raw = TRUE) * above, data = data, subset = incumb.reelected.in.2004 == 0 & run.reelec == 1 & abs(fv) < .3))
alt_spec$poly_voteshare$se <- sqrt(diag(vcovHC(alt_spec$poly_voteshare, type = "HC1")))
alt_spec$cov_voteshare <-(lm(incumbent.vote.2008.perpop - winner.vote.2004.perpop ~ poly(fv, degree = 1, raw = TRUE) * above + uf + turnout.04.perpop + vote.margin.04.perpop + pov2000 + incumb.partywon.04 + theil2000 + rural.pop.07, data = data, subset = incumb.reelected.in.2004 == 0 & run.reelec == 1 & abs(fv) < bw.ll, weights = 1 - abs(fv/bw.ll)))
alt_spec$cov_voteshare$se <- sqrt(diag(vcovHC(alt_spec$cov_voteshare, type = "HC1")))

alt_spec_table_cov <- stargazer(
  alt_spec$cov_incumb,
  alt_spec$cov_partyincumb,
  alt_spec$cov_run,
  alt_spec$cov_voteshare,
  alt_spec$cov_electorate,
  alt_spec$cov_turnout,
  se = list(
    alt_spec$cov_incumb$se,
    alt_spec$cov_partyincumb$se,
    alt_spec$cov_run$se,
    alt_spec$cov_voteshare$se,
    alt_spec$cov_electorate$se,
    alt_spec$cov_turnout$se
  ),
  dep.var.labels = c("Incumb. Reelected", "Incumb. Party Reelected", "Incumb. Re-runs", "Incumb. Vote Share", "Electorate Change", "Turnout Change"),
  covariate.labels = c("$\\hat \\tau_A$", "Baseline"),
  omit = c("uf", "poly", "turnout", "vote", "transfers", "incumb", "pib", "theil", "rural","Constant"),
  omit.stat = c("rsq", "adj.rsq", "f", "ser"),
  digits = 2,
  float = FALSE,
  out = 'table_A3.tex')


#######
#Table A4
#######

#Table with estimates using a polynomial specification

alt_spec_table_poly <- stargazer(
  alt_spec$poly_incumb,
  alt_spec$poly_partyincumb,
  alt_spec$poly_run,
  alt_spec$poly_voteshare,
  alt_spec$poly_electorate,
  alt_spec$poly_turnout,
  se = list(
    alt_spec$poly_incumb$se,
    alt_spec$poly_partyincumb$se,
    alt_spec$poly_run$se,
    alt_spec$poly_voteshare$se,
    alt_spec$poly_electorate$se,
    alt_spec$poly_turnout$se
  ),
  dep.var.labels = c("Incumb. Reelected", "Incumb. Party Reelected", "Incumb. Re-runs", "Incumb. Vote Share", "Electorate Change", "Turnout Change"),
  covariate.labels = c("$\\hat \\tau_A$", "Baseline"),
  omit = c("uf", "poly", "turnout", "vote", "transfers", "incumb", "pib"),
  omit.stat = c("rsq", "adj.rsq", "f", "ser"),
  digits = 2,
  float = FALSE,
  out = 'table_A4.tex')

#######
#Figure A5
#######

#Estimate effects at placebo bandwidths

placebo.incumbwins <- list()
placebo.cut <- seq(.7, .9, by = .02)
placebo.incumbwins <- lapply(placebo.cut, function(x) rd.estim(y = data$incumbent.won08[data$incumb.reelected.in.2004 == 0] - data$incumb.partywon.04[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = bw.dm, bw.ll = bw.dm, cutpoint = x))
placebo.electorate <- lapply(placebo.cut, function(x) rd.estim(y =  data$electoratechange0608.perpop, data = data, bw.dm = bw.dm, bw.ll = bw.dm, cutpoint = x))
placebo.turnout <- lapply(placebo.cut, function(x) rd.estim(y =  data$turnout.08.perpop - data$turnout.04.perpop, data = data, bw.dm = bw.dm, bw.ll = bw.dm,  cutpoint = x))
placebo.incumbparty <- lapply(placebo.cut, function(x) rd.estim(y =  data$incumb.partywon.08 - data$incumb.partywon.04, data = data, bw.dm = bw.dm, bw.ll = bw.dm,  cutpoint = x))
placebo.runreelection <- lapply(placebo.cut, function(x) rd.estim(y =  data$run.reelec[data$incumb.reelected.in.2004 == 0], data = data[data$incumb.reelected.in.2004 == 0, ], bw.dm = bw.dm, bw.ll = bw.dm,  cutpoint = x))
placebo.incumbvoteshare <- lapply(placebo.cut, function(x) rd.estim(y =  data$incumbent.vote.2008.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1] - data$winner.vote.2004.perpop[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1], data = data[data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1, ], bw.dm = bw.dm, bw.ll = bw.dm, cutpoint = x))


placebo.plot.data <- rbind(data.frame(bw = placebo.cut, est = sapply(placebo.incumbwins, function(x) x$rfll.coef), se = sapply(placebo.incumbwins, function(x) x$rfll.se), Variable = "Incumbent Wins"),
                           data.frame(bw = placebo.cut, est = sapply(placebo.electorate, function(x) x$rfll.coef), se = sapply(placebo.electorate, function(x) x$rfll.se), Variable = "Electorate Change"),
                           data.frame(bw = placebo.cut, est = sapply(placebo.turnout, function(x) x$rfll.coef), se = sapply(placebo.turnout, function(x) x$rfll.se), Variable = "Turnout Change"),
                           data.frame(bw = placebo.cut, est = sapply(placebo.incumbparty, function(x) x$rfll.coef), se = sapply(placebo.incumbparty, function(x) x$rfll.se), Variable = "Incumbent Party Wins"),
                           data.frame(bw = placebo.cut, est = sapply(placebo.runreelection, function(x) x$rfll.coef), se = sapply(placebo.runreelection, function(x) x$rfll.se), Variable = "Incumbent Runs Again"),
                           data.frame(bw = placebo.cut, est = sapply(placebo.incumbvoteshare, function(x) x$rfll.coef), se = sapply(placebo.incumbvoteshare, function(x) x$rfll.se), Variable = "Incumbent Vote Share Change")
)

placebo.plot.data$Variable <- factor(placebo.plot.data$Variable, levels = c("Electorate Change", "Turnout Change", "Incumbent Wins",  "Incumbent Party Wins", "Incumbent Runs Again", "Incumbent Vote Share Change"))

pdf(file = "figure_A5.pdf", height = 10, width = 8.5)
ggplot(placebo.plot.data, aes(x = bw * 100, y = est)) +
  geom_pointrange(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  facet_wrap(~ Variable, ncol = 1, scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  scale_x_continuous("Placebo Thresholds", breaks = placebo.cut * 100) +
  scale_y_continuous("Estimate")
dev.off()


#######
#Figure A6
#######

#Additional heterogenous treatment effects analysis

hetero.alternative.appendix.plotdata <- data.frame(Estimate = unlist(sapply(hetero.incumbwins[c('transfers.out0207', 'transfers.out9707', 'faltosos')], function(x) x[grep("coef", names(x))])),
                                          SE = unlist(sapply(hetero.incumbwins[c('transfers.out0207', 'transfers.out9707', 'faltosos')], function(x) x[grep("se", names(x))])),
                                          Strata = rep(c("Below\nMedian",  "Above\nMedian",  "Below\nMedian", "Above\nMedian", "Difference", "Difference"), 3),
                                          Specification = rep(c(rep("Difference-in-Means", 2), rep("Local Linear", 2), "Difference-in-Means", "Local Linear"), 3),
                                          Covariate = c(rep('(a) Transfers Out (02-07)', 6), rep('(b) Transfers Out (97-07)', 6), rep('(c) Missing Voters', 6)))

hetero.alternative.appendix.plotdata$Strata <- factor(hetero.alternative.appendix.plotdata$Strata, levels = c("Below\nMedian", "Above\nMedian", "Difference"))
hetero.alternative.appendix.plotdata$Covariate <- factor(as.character(hetero.alternative.appendix.plotdata$Covariate), levels = c('(a) Transfers Out (02-07)', '(b) Transfers Out (97-07)', '(c) Missing Voters'))


pdf(file = "figure_A6.pdf", height = 5, width = 8)
ggplot(hetero.alternative.appendix.plotdata[hetero.alternative.appendix.plotdata$Covariate %in% c('(a) Transfers Out (02-07)', '(b) Transfers Out (97-07)', '(c) Missing Voters'), ],
       aes(x = Strata, y = Estimate, color = Specification, shape = Specification))+
  facet_wrap( ~ Covariate, ncol  = 3, scales = "free_x") + scale_colour_brewer(palette = "Set1", guide = guide_legend(nrow = 1)) + 
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(legend.position="bottom") +
  scale_y_continuous("Audit Effect on Incumbent Reelection") +
  geom_pointrange(aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 * SE), size = .8, position=position_dodge(width=0.4)) +
  xlab("")
dev.off()


#######
#Figure A8
#######

#Show trajectory of registration across years
reg_traj <- data[, c('codigo', 'elecpopratio', 'electorate2002', 'electorate2003', 'electorate2004', "electorate2005",  'electorate2006', 'electorate2007', 'electorate2008', 'electorate2009', 'electorate2010', 'electorate2011', 'electorate2012', 'pop2007')]
names(reg_traj) <- c('codigo', 'elecpopratio', 'electorate_2002', 'electorate_2003', 'electorate_2004', "electorate_2005",  'electorate_2006', 'electorate_2007', 'electorate_2008', 'electorate_2009', 'electorate_2010', 'electorate_2011', 'electorate2012', 'pop2007')
reg_traj <- gather(reg_traj, year, electorate, -codigo, -pop2007, -elecpopratio)
reg_traj$year <- gsub("electorate_", reg_traj$year, replacement = "")
reg_traj$electorate_pct <- 100 * reg_traj$electorate / reg_traj$pop2007
reg_traj <- reg_traj[order(reg_traj$codigo, reg_traj$year), ]
reg_traj <- group_by(reg_traj, codigo) %>%
  mutate(lag_electorate_pct = lag(electorate_pct), electorate_change = electorate_pct - lag_electorate_pct)
reg_traj <- reg_traj[reg_traj$elecpopratio < .8401 & reg_traj$elecpopratio > .759, ]
reg_traj$above <- as.factor(ifelse(reg_traj$elecpopratio > .8, 1, 0))
reg_traj <- reg_traj[reg_traj$year > 2002, ]

pdf(file = "figure_A8_top_panel.pdf", height = 7, width = 11.5)
ggplot(reg_traj[reg_traj$year < 2008, ], aes(x = elecpopratio * 100, y = electorate_pct)) +
  geom_point(aes(color = above), alpha = .7, size = .5) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = .7, degree = 1, alpha = 1) +
  facet_wrap(~year, nrow = 1) +
  scale_x_continuous("Electorate as % of Population", breaks = seq(76, 84, 1), labels = comma_format()) +
  scale_y_continuous("Registration (% of Population)", limits = c(50, 100)) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Before Audits")
dev.off()


pdf(file = "figure_A8_bottom_panel.pdf", height = 7, width = 11.5)
ggplot(reg_traj[reg_traj$year > 2007, ], aes(x = elecpopratio * 100, y = electorate_pct)) +
  geom_point(aes(color = above), alpha = .7, size = .5) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = .7, degree = 1, alpha = 1) +
  facet_wrap(~year, nrow = 1) +
  scale_x_continuous("Electorate as % of Population", breaks = seq(76, 84, 1), labels = comma_format()) +
  scale_y_continuous("Registration (% of Population)", limits= c(50, 100)) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("After Audits")
dev.off()

#######
#Figure A9
#######

#Show distribution of candidate characteristics across broad range of the forcing variable

exval_data <- data[data$elecpopratio < .92 & data$elecpopratio > .67,]
exval_data$bin <- (as.character(cut(exval_data$elecpopratio, breaks = seq(.67,.92,.05))))

exval_winmargin_plot <- ggplot(exval_data, aes( x = elecpopratio * 100, y = winmargin.04)) +
  geom_point(aes(color = above), alpha = .6, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, se = TRUE, size = 1.5, color = "black", span = .3)  +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Vote Margin (2004)") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle('Vote Margin (2004)')


exval_donation_plot <-ggplot(exval_data, aes( x = elecpopratio * 100, y = donation_diff)) +
  geom_point(aes(color = above), alpha = .6, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, se = TRUE, size = 1.5, color = "black", span = .3)  +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous(limits = c(-50000, 50000), "Difference in Total Donations (BRL)") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle('Difference in Winner / Runner Up Donations (2008)')

exval_incumbage_plot <- ggplot(exval_data, aes( x = elecpopratio * 100, y = incumb_age)) +
  geom_point(aes(color = above), alpha = .6, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, se = TRUE, size = 1.5, color = "black", span = .3)  +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Incumbent Age") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle('Incumbent Age')

exval_chall_plot <- ggplot(exval_data, aes( x = elecpopratio * 100, y = chall_age)) +
  geom_point(aes(color = above), alpha = .6, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, se = TRUE, size = 1.5, color = "black", span = .3)  +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Challenger Age") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle('Challenger Age')


pdf(file = "figure_A9.pdf", height = 10, width = 8)
arrange(exval_winmargin_plot,
        exval_donation_plot,
        exval_incumbage_plot,
        exval_chall_plot,
        ncol = 1)
dev.off()


#######
#Figure A10
#######

#Show discontinuity plots for selected covariates


smoothness_data <- data[data$elecpopratio < .90 & data$elecpopratio > .70 & data$incumb.reelected.in.2004 == 0,]
smoothness_data$bin <-  100 * as.numeric(as.character(cut(smoothness_data$elecpopratio, breaks = seq(.7,.9,.015), labels = seq(.705, .8901,.015))))

smoothness_vars <- c('vote.margin.04.perpop', 'incumb_age', 'num_vereadores04', 'muni_resources2006_percap', 'pub.employee.percap.2006', 'radio')

incumbage_plot_data <- group_by(smoothness_data[, c('bin', smoothness_vars)], bin) %>% summarise(bin_mean = mean(incumb_age, na.rm = TRUE), bin_se = sd(incumb_age, na.rm = TRUE) / sqrt(n()))

smooth_plot_incumbage <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  incumb_age)) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Incumbent Age") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Incumbent Age")

smooth_plot_challage <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  chall_age)) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Challenger Age") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Challenger Age")


smooth_plot_votemargin04 <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699 & smoothness_data$winmargin.04 < 100,] , aes(x=elecpopratio.pct, y =  winmargin.04)) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Vote Margin (2004)") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle('Vote Margin (2004)')


smooth_plot_numver <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  electorate2006 / (num_vereadores04))) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Electorate to Number\nof City Councilors Ratio", limits = c(0,30000)) +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Electorate to\nNumber of City Councilors Ratio")


smooth_plot_budget <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  log(muni_resources2006_percap))) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Log Budget Per Capita (2006)") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Log Municipal Budget Per Capita (2006)")


smooth_plot_pubempl <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  pub.employee.percap.2006)) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("# of Public Employees Per Capita", limits = c(0,.3)) +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Number of Public Employees Per Capita")


smooth_plot_radio <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,], aes(y = radio, x = elecpopratio.pct / 100)) +
  geom_vline(xintercept = .80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = 1, degree = 1, alpha = 1) +
  scale_x_continuous("Electorate as % of Population", breaks = seq(.70, .90, .02), label = percent) +
  scale_y_continuous("Radio in Municipality", labels = percent) +
  stat_summary(aes(x = bin.equal/100, color = bin.equal/100 > .80), fun.y = mean, geom = "point", size = 1, alpha = .7) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Radio in Municipality")

smooth_plot_seccsize <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  average_section_size)) +
  geom_point(aes(color = above), alpha = .5, size = 1) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Average Precinct Size") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Average Precinct Size")

smooth_plot_turnout <- ggplot(smoothness_data[smoothness_data$elecpopratio < .901 & smoothness_data$elecpopratio > .699,] , aes(x=elecpopratio.pct, y =  turnout.04.perpop)) +
  geom_point(aes(color = above), alpha = .5, size = .8) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, alpha = 1)  +
  scale_x_continuous("Electorate as % of Population", breaks = seq(70, 90, 2), labels = comma_format()) +
  scale_y_continuous("Turnout (2004)") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Turnout (2004)")


pdf(file = "figure_A10.pdf", height = 10, width = 8)
arrange(smooth_plot_votemargin04,
        smooth_plot_budget,
        smooth_plot_pubempl,
        smooth_plot_incumbage,
        smooth_plot_challage,
        smooth_plot_radio,
        smooth_plot_numver,
        smooth_plot_seccsize,
        ncol = 2)
dev.off()


#######
#Figure A11
#######

#Incumbent reelection discontinuity plot showing full range of the data

pdf(file = "figure_A11.pdf", height = 7, width = 11.5)
ggplot((data[data$elecpopratio < 1.1 & data$elecpopratio > .5 & data$incumb.reelected.in.2004 == 0  ,]), aes(y = incumbent.won08, x = elecpopratio.pct)) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = .4, degree = 1) +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Incumbent Reelected in 2008", limit = c(0,1)) +
  stat_summary(aes(x = bin.equal, color = bin.equal > 80), fun.y = mean, geom = "point", size = 3) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1")
dev.off()

#######
#Figure A12
#######

#Incubment run for reelection discontinuity plot showing full range of the data

pdf(file = "figure_A12.pdf", height = 7, width = 11.5)
ggplot(data[data$elecpopratio < 1.1 & data$elecpopratio > .5 & data$incumb.reelected.in.2004 == 0,], aes(y = run.reelec, x = elecpopratio.pct )) + geom_vline(xintercept = 80)  +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Incumbent Reelected in 2008", limit = c(0,1)) +
  stat_summary(aes(x = bin.equal, color = bin.equal > 80), fun.y = mean, geom = "point", size = 3) +
  theme_bw() + theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set1") +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = .4, degree = 1)
dev.off()

#######
#Figure A13
#######

#Incubment vote share discontinuity plot showing full range of the data

pdf(file = "figure_A13.pdf", height = 7, width = 11.5)
ggplot((data[data$elecpopratio < 1.1 & data$elecpopratio > .5  & data$incumb.reelected.in.2004 == 0 & data$run.reelec == 1,]), aes(y = incumbent.vote.2008.perpop - winner.vote.2004.perpop, x = elecpopratio.pct))  +
  geom_point(aes(color = above), alpha = .7, size = 1.5) +
  geom_vline(xintercept = 80) +
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", degree = 1, span = .6)  +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Change in Vote Share (2008-2004)") +
  theme_bw() + theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(-40, 80))
dev.off()

#######
#Figure A14
#######

#Electorate discontinuity plot showing full range of the data


pdf(file = "figure_A14.pdf", height = 7, width = 11.5)
ggplot((data[data$elecpopratio < 1.1 & data$elecpopratio > .5 & data$incumb.reelected.in.2004 == 0  ,]), aes(y = electorate.perpop08 - electorate.perpop06, x = elecpopratio.pct)) +
  geom_point(aes(color = above), alpha = .7, size = 1.5) +
  geom_vline(xintercept = 80) + stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = .4, degree = 1) +
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) + 
  scale_y_continuous("Change in Electorate (% of Population)") +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_colour_brewer(palette = "Set1")
dev.off()

#######
#Figure A15
#######

#Turnout discontinuity plot showing full range of the data

pdf(file = "figure_A15.pdf", height = 7, width = 11.5)
ggplot((data[data$elecpopratio < 1.1 & data$elecpopratio > .5 & data$incumb.reelected.in.2004 == 0  ,]), aes(y = turnout.08.perpop - turnout.04.perpop, x = elecpopratio.pct)) +
  geom_point(aes(color = above), alpha = .7, size = 1.5) +
  geom_vline(xintercept = 80) + 
  stat_smooth(method = loess, aes(fill = above), se = TRUE, size = 1.5, color = "black", span = .4, degree = 1) + 
  scale_x_continuous("Electorate as % of Population", labels = comma_format()) +
  scale_y_continuous("Change in Turnout (% of Population)") + 
  theme_bw() + theme(legend.position = "none")  + 
  scale_colour_brewer(palette = "Set1")
dev.off()


#############################
#Version of R and Packages used are listed below
#############################

# R version 3.2.0 (2015-04-16)
# Platform: x86_64-apple-darwin14.3.0 (64-bit)
# Running under: OS X 10.10.3 (Yosemite)
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_0.4.1       tidyr_0.2.0       reporttools_1.1.1 xtable_1.7-4      rdd_0.56          stargazer_5.1     AER_1.2-3        
# [8] lmtest_0.9-33     zoo_1.7-12        sandwich_2.3-3    car_2.0-25        scales_0.2.4      Hmisc_3.15-0      Formula_1.2-1    
# [15] survival_2.38-1   lattice_0.20-31   ggplot2_1.0.1     reshape2_1.4.1   
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.11.5         RColorBrewer_1.1-2  nloptr_1.0.4        plyr_1.8.2          tools_3.2.0         rpart_4.1-9        
# [7] digest_0.6.8        lme4_1.1-7          gtable_0.1.2        nlme_3.1-120        mgcv_1.8-6          Matrix_1.2-0       
# [13] DBI_0.3.1           parallel_3.2.0      SparseM_1.6         proto_0.3-10        stringr_0.6.2       cluster_2.0.1      
# [19] nnet_7.3-9          foreign_0.8-63      latticeExtra_0.6-26 minqa_1.2.4         magrittr_1.5        MASS_7.3-40        
# [25] splines_3.2.0       assertthat_0.1      pbkrtest_0.4-2      colorspace_1.2-6    quantreg_5.11       acepack_1.3-3.3    
# [31] munsell_0.4.2 
