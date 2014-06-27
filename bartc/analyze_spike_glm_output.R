# load up useful vars
source('helpers.R')

# load analysis output
load(file=paste(ddir, 'spkfitdata', sep='/'))

######## code to plot heatmap of regression coefficients ##########
betas <- ldply(fitobjs, .fun = function(x) {data.frame(t(x$beta))}) 
effects <- exp(betas) * 100 - 100
effects <- cbind(data.frame(unit=1:dim(effects)[1]), effects)
df <- melt(effects, id.vars=c('unit'))

plt <- plot_spike_coefficient_grid(df)
print(plt)
