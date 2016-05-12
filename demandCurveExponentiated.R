# 
#    Copyright 2016 Shawn Gilroy
#
#    This file is part of Small N Stats.
#
#    Small N Stats is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 2.
#
#    Small N Stats is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Small N Stats.  If not, see <http://www.gnu.org/licenses/gpl-2.0.html>.
#
# summary
# R script for individual/batched demand curve analyses.
#
# dependencies
# @ggplot2 = individual/aggregated visual of logged demand curves (Copyright 2016 - Hadley Wickham - GPLv2)
# link @ https://cran.r-project.org/web/packages/ggplot2/index.html
# license @ https://cran.r-project.org/web/licenses/GPL-2
#
# @nlstools = confint2 for bootstrapped confidence intervals (Copyright 2015 - Baty and Delignette-Muller - GPLv2+)
# link @ https://cran.r-project.org/web/packages/nlstools/index.html
# license @ https://cran.r-project.org/web/licenses/GPL-2
#
# @nlmrt (R package) = Nash's customized optimization of L-M residual reduction (Copyright 2016 - John C. Nash - GPLv2)
# link @ https://cran.r-project.org/web/packages/nlmrt/index.html
# license @ https://cran.r-project.org/web/licenses/GPL-2
#
# params = Demand Function fittings
# @p <- participant id  {e.g.,  c(1,1,1,  2, ... )    }
# @y <- consumption     {e.g.,  c(10,7,3, 9, ... )    }
# @x <- pricing         {e.g.,  c(0.5,1,  5, ... )    }
# @k <- k value         {e.g.,  c(4,4,4,  4, ... )    }
#

if (!require(ggplot2)) { install.packages('ggplot2', repos = 'http://cran.us.r-project.org') }

if (!require(nlstools)) { install.packages('nlstools', repos = 'http://cran.us.r-project.org') }

if (!require(nlmrt)) { install.packages('nlmrt', repos = 'http://cran.us.r-project.org') }

SourceFrame <- data.frame(
  p=pLoad,
  y=yLoad,
  x=xLoad,
  k=kLoad)

nSimulated <- max(SourceFrame$p)

fitFrame <- data.frame(
  p=seq(1,nSimulated,1),
  q0=rep(NA,nSimulated),
  q0err=rep(NA,nSimulated),
  alpha=rep(NA,nSimulated),
  alphaerr=rep(NA,nSimulated),
  k=rep(NA,nSimulated),
  r2=rep(NA,nSimulated),
  absSS=rep(NA,nSimulated),
  sdResid=rep(NA,nSimulated),
  q0low=rep(NA,nSimulated),
  q0high=rep(NA,nSimulated),
  alow=rep(NA,nSimulated),
  ahigh=rep(NA,nSimulated))

for (i in 1:nSimulated)
{
  fit <- wrapnls(data=SourceFrame[SourceFrame$p==i,], y ~ q0 * 10^(k * (exp(-alpha*q0*x)-1)), start=c(q0=3, alpha=0.000000001), control = list(maxiter = 1000))
  
  fitFrame[fitFrame$p==i,]$q0 <- as.numeric(coef(fit)["q0"])
  fitFrame[fitFrame$p==i,]$alpha <- as.numeric(coef(fit)["alpha"])
  fitFrame[fitFrame$p==i,]$k <- min(SourceFrame$k)
  fitFrame[fitFrame$p==i,]$q0err <- summary(fit)[[10]][1,2]
  fitFrame[fitFrame$p==i,]$alphaerr <- summary(fit)[[10]][2,2]
  fitFrame[fitFrame$p==i,]$r2 <- 1.0 -(deviance(fit)/sum((SourceFrame[SourceFrame$p==i,]$y-mean(SourceFrame[SourceFrame$p==i,]$y))^2))
  fitFrame[fitFrame$p==i,]$absSS <- deviance(fit)  
  fitFrame[fitFrame$p==i,]$sdResid <- sqrt(deviance(fit)/df.residual(fit))
  fitFrame[fitFrame$p==i,]$q0low <- confint2(fit)[1]
  fitFrame[fitFrame$p==i,]$q0high <- confint2(fit)[3]
  fitFrame[fitFrame$p==i,]$alow <- confint2(fit)[2]
  fitFrame[fitFrame$p==i,]$ahigh <- confint2(fit)[4]
}

xDraw <- seq(min(SourceFrame$x), max(SourceFrame$x), 0.01)

p.rep <- seq(1,max(SourceFrame$p),1)

graphFrame<-data.frame(
  Individual=rep(seq(min(p.rep), max(p.rep),1),each=length(xDraw)),
  DemandSeries=rep(seq(1:length(xDraw)-1),length(p.rep)),
  YSeries=rep(seq(1:length(xDraw)-1),length(p.rep)),
  XSeries=rep(seq(1:length(xDraw)-1),length(p.rep))
)

for (j in 1:max(SourceFrame$p))
{
  for (i in 1:length(xDraw))
  {
    qTemp <- fitFrame[fitFrame$p==j,]$q0
    aTemp <- fitFrame[fitFrame$p==j,]$alpha
    kTemp <- fitFrame[fitFrame$p==j,]$k
    
    graphFrame[ graphFrame$Individual==j & graphFrame$DemandSeries==as.numeric(i),]$YSeries <- log10(qTemp) + kTemp * (exp(-aTemp*qTemp*xDraw[i])-1)
    graphFrame[ graphFrame$Individual==j & graphFrame$DemandSeries==as.numeric(i),]$XSeries <- xDraw[i]
  }
}

custom_axis <- function(l) { 
  l <- paste("10^", l, sep = "")  
  parse(text=l) 
} 

pointFrame <- data.frame(X=SourceFrame$x, Y=SourceFrame$y, Individual=SourceFrame$p)

logChart <- ggplot() +
  geom_line(data=graphFrame, aes(x=XSeries, y=YSeries, group=Individual, colour = factor(Individual))) + 
  geom_point(data=pointFrame, aes(x=pointFrame$X, y=log10(pointFrame$Y), shape=factor(Individual))) +
  expand_limits(y=0) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("Fitted Demand Curves\n") +
  ylab("log(Consumption)") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  scale_y_continuous(labels=custom_axis) +
  annotation_logticks(sides = "b") +
  xlab("log(Price)") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(legend.direction = "vertical") + 
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  guides(col = guide_legend(ncol = 3))
