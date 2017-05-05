## Adam Chaffee
## Analysis of Flint lead data using different kriging methods
## March 2017

##############################
## Packages used
##############################

## geoR, gstat, and sp are geospatial analysis packages
library(geoR)
library(gstat)
library(sp)

## Maps package can provide maps. 
library(maps)

## companion to applied regression, for model diagnostics
library(car)

## mosaic is helpful for visualizations. Also loads the populat ggplot2 package
library(mosaic)

########################################################################
## Load the data
## Data notes: 
## (1) Public data. URL: http://www.michigan.gov/flintwater/0,6092,7-345-76292_76294_76297---,00.html
##     Clicked residential testing Jan 1 - March 31, 2017
## (2) data was provided in Excel format which lended better to cleaning in Excel
##      (a) Looked at only two bottle tests. Some of the single bottle tests were
##          missing data. Both types of tests are supposed to produce same result
##      (b) Found a geocoder online to convert addresses to lat/lon. Removed about
##          50 locations which the geocoder could not convert
########################################################################
setwd("C:/Users/achaf/OneDrive/Documents/UCLA/Winter 17/Stat C273/Project")
flint = read.csv("Flint data final.csv", header = T)

## Needed to fix the latitude column name, it loads up weird for me
names(flint) = c("Latitude", "Longitude", "Pb", "Cu")

################################
## Initial plot of data - bubble plot
################################
x = flint$Longitude
y = flint$Latitude
xmin = min(flint$Latitude)
xmax = max(flint$Latitude)
ymin = min(flint$Longitude)
ymax = max(flint$Longitude)

## Establish colors and point sizes for different levels
## 15 PPB is the EPA mandatory action level
lead_colors = c("green", "goldenrod", "red")
lead_pch = c(1, 1.5, 2)
lead_levels = cut(flint$Pb, c(-1, 0.5, 15, 10000))

## Plot of the data
plot(x,y, col = "white", main = "Flint Lead Data",
     xlab = "Latitude", ylab = "Longitude", xlim = c(-83.81, -83.6), "n")
points(x, y, col = lead_colors[as.numeric(lead_levels)], 
       pch = 19, cex = lead_pch[as.numeric(lead_levels)])
legend("topleft", c("0 PPB", "<15 PPB", ">15 PPB"), pch = 19, 
       col = lead_colors, cex = 1.2)


###########################################
## Now assess potential for co-kriging
## Scatter plot for Pb, Cu and regression run
###########################################
plot(flint$Cu ~flint$Pb,  col = "blue")
flint_lm = lm(flint$Cu ~ flint$Pb)
abline(flint_lm, col = "red", lwd = 2)

## Weak linear association. R-squared at about 0.25
summary(flint_lm)

## Diagnostics. Concern over normality of residuals
## Hard to visualize in a histogram due to the large spread and few data points
## Large spike at -37 because we had lots of 0 values and the regression equation has
## intercept at 37
qqPlot(flint_lm)
histogram(flint_lm$residuals, breaks = c(-600,-100:15, 20, 25, 30, 2000),
          xlim = c(-50,50))


###########################################
## Variogram fitting
###########################################
names(flint) = c("y", "x", "Pb", "Cu")

# Create gstat object
g1 = gstat(id = "Pb", formula = Pb~1, locations = ~x+y, data = flint)

# Append copper for cross variogram
g1 = gstat(g1, id = "Cu", formula = Cu~1, locations = ~x+y, data = flint)

#Plot the variograms and cross-variograms:
plot(variogram(g1))

#Fit a model variogram to the target variogram:
g = gstat(id="Pb", formula = Pb~1, locations = ~x+y, data = flint) 
v.fit = fit.variogram(variogram(g), vgm(80, "Exp", 0.03, 80)) 

## Basically pure nugget effect
plot(variogram(g),v.fit)
v.fit$range
vm = variogram(g1) 

# Fit a model to all the variograms
vm.fit = fit.lmc(vm, g1, model=v.fit) 
vm.fit

#Plot the fitted variograms to all the sample variograms:
plot(variogram(g1),vm.fit)

#######################################
## Analysis of ordinary, simple, and universal kriging methods
#######################################

## Perform co-kriging predictions:
ck1 = predict(vm.fit, grd)
ck1

## Perform leave one out CV and compute sum squared error (PRESS):
cv_ck = gstat.cv(vm.fit)
PRESS_ck = sum(cv_ck$residual^2)

## CV using ordinary kriging. Apply similar steps as above
gg = gstat(id = "Pb", formula = Pb~1, locations = ~x+y, data = flint)

plot(variogram(gg))

v.fit = fit.variogram(variogram(gg), vgm(80, "Exp", .03, 80))

plot(variogram(gg), v.fit)
mean(flint$Pb)

## CV for ordinary, simple, and universal kriging
## Simple kriging assumes mean is known. Mean is set to sample average
cv_pr_ok = krige.cv(Pb~1, data = flint, locations = ~x+y, model = v.fit)

cv_pr_sk = krige.cv(Pb~1, data = flint, locations = ~x+y, model = v.fit, 
                    beta = mean(flint$Pb))

cv_pr_uk = krige.cv(Pb~x+y, data = flint, locations = ~x+y, 
                   model = v.fit, nfold=nrow(flint)) 

## Calculate more sum square errors
PRESS_ok = sum(cv_pr_ok$residual^2)
PRESS_uk = sum(cv_pr_uk$residual^2)
PRESS_sk = sum(cv_pr_sk$residual^2)

## Look at the cross validated sum square error. Best model has lowest PRESS
PRESS_ok
PRESS_ck
PRESS_uk
PRESS_sk

## Move forward with prediction for co-kriging

#######################################
## Generation of prediction raster/contour map for co-kriging
######################################

## Credit to Dr. Nicolas Christou for sample code

## Create the grid for predictions:
x.range = round(range(flint[,2]), 2) 
y.range = round(range(flint[,1]), 2) 
grd = expand.grid(x=seq(from=x.range[1], to=x.range[2], by=.0025), 
                   y=seq(from=y.range[1], to=y.range[2], by=.0025)) 
y = seq(from = y.range[1], to = y.range[2], by = 0.0025)

## Collapse the predicted values into a matrix:
qqq = matrix(ck1$Pb.pred,
              length(seq(from=x.range[1], to=x.range[2], by=.0025)),
              length(seq(from=y.range[1], to=y.range[2], by=.0025)))

## Create image without extrapolating
X = seq(from=x.range[1], to=x.range[2], by=.0025)
Y = seq(from=y.range[1], to=y.range[2], by=.0025)
dz.mask = function(grid, pts, x.dom, y.dom, window, mitre=2) 
{
  N = length(pts[ ,1]) ; mask <- array( NA, dim(grid) )
  for(j in 1:N) 
  {
    dx = abs(x.dom - pts$x[j])
    dy = abs(y.dom - pts$y[j])
    d.Mx = tcrossprod( dx , rep(1, length(dy)) )^mitre +
      tcrossprod( rep(1, length(dx)), dy )^mitre
    mask[ d.Mx < window^mitre ] = FALSE
  }
  return(mask+grid)
}
qqq.masked <- dz.mask(qqq, flint, X, Y, 0.01)
image(X,Y,qqq.masked, xlab="West to East", ylab="South to North", main="Raster map - Pb concentration")
points(flint$x, flint$y, pch = 1, cex = .7)

contour(seq(from=x.range[1], to=x.range[2], by=0.0025), 
        seq(from=y.range[1],to=y.range[2], by=0.0025), qqq.masked, add=TRUE, col="black", 
        levels=seq(0, 10, by=1), labcex=1)

## Additional contour for >15PPB in blue to identify most serious areas of concern
contour(seq(from=x.range[1], to=x.range[2], by=0.0025), 
        seq(from=y.range[1],to=y.range[2], by=0.0025), qqq.masked, add=TRUE, col="blue", 
        levels=15, labcex=1, lwd = 2)

#################################################
## Indicator Kriging analysis
#################################################
image.orig = image

## Set up indicator data
c1 = ifelse(flint$Pb > 15, 1 , 0)
flint$c1 = c1

indic = gstat(id="c1", formula=c1~1, locations=~x+y, data=flint)

plot(variogram(indic))

## From the plot spherical variogram is likely the best
v.fit = fit.variogram(variogram(indic), vgm(0.2,"Sph", .02, 0.03))
plot(variogram(indic),v.fit, main = "Variogram - Indicator Data")

## Cross validate to find best universal, simple, and indicator kriging results
cv_pr_uk_ind = krige.cv(c1~x+y, data = flint, locations = ~x+y, 
                    model=v.fit, nfold=nrow(flint)) 
cv_pr_sk_ind = krige.cv(c1~1, data = flint, locations = ~x+y, 
                        model = v.fit, beta = mean(c1))
cv_pr_ok_ind = krige.cv(c1~1, data = flint, locations = ~x+y, model = v.fit)

## Compare sum squared residuals to find best fit
sum(cv_pr_ok_ind$residual^2)
sum(cv_pr_sk_ind$residual^2)
sum(cv_pr_uk_ind$residual^2)

## Simple kriging gave the best model. Create a simple kriging object
pr_sk_ind = krige(id="c1",c1~1, locations=~x+y, model=v.fit, 
                   data=flint, newdata=grd, beta = mean(c1))

## Set values < 0 equal to zero, and values > 1 equal to 1:
for(i in 1:length(pr_sk_ind$c1.pred) ){if(pr_sk_ind$c1.pred[i] <0) {pr_sk_ind$c1.pred[i]=0}}

for(i in 1:length(pr_sk_ind$c1.pred) ){if(pr_sk_ind$c1.pred[i] >1) {pr_sk_ind$c1.pred[i]=1}}


##################################################
## Raster/contour map for indicator kriging
##################################################

#Collapse the predicted values into a matrix: 
qqq_ind = matrix(pr_sk_ind$c1.pred, 
              length(seq(from=x.range[1], to=x.range[2], by=0.0025)), 
              length(seq(from=y.range[1], to=y.range[2], by=0.0025))) 

qqq.masked <- dz.mask(qqq_ind, flint, X, Y, 0.01)

image.orig(X,Y,qqq.masked, main = "Raster Map using Indicator Kriging",
           xlab = "West to East", ylab = "South to North")
points(flint$x, flint$y, pch = 1, cex = .7)

## Add contour lines:
contour(seq(from=x.range[1], to=x.range[2], by=0.0025), 
        + seq(from=y.range[1],to=y.range[2], by=0.0025), qqq.masked, 
        levels=seq(0, 1, by=0.05), add=TRUE, col="black", labcex=1) 

