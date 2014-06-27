# load up useful vars
source('helpers.R')

# load analysis output
load(file=paste(ddir, 'spkfitdata', sep='/'))

######## code to plot heatmap of regression coefficients ##########
df <- ldply(fitobjs, .fun = function(x) {data.frame(t(x$beta))}) 

plt <- plot_coefficient_grid(df)
print(plt)
