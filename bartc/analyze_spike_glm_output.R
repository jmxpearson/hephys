# load up useful vars
source('helpers.R')

# load analysis output
load(file=paste(ddir, 'spkfitdata', sep='/'))

######## code to plot heatmap of regression coefficients ##########
df <- ldply(fitobjs, .fun = function(x) {data.frame(t(x$beta))}) 

# useful summary stats
band_stats <- ddply(df, ~band, summarize, mean = mean(value), 
    mean_abs = mean(abs(value)), std = sd(value), std_abs = sd(abs(value)))

plt <- plot_coefficient_grid(df)
print(plt)
