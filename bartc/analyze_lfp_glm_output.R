# load up useful vars
source('helpers.R')

# load analysis output
load(file=paste(ddir, 'lfpfitdata', sep='/'))

# now plot an ROC curve for a given subject
ind <- 10  # numbered consecutively, starting at 1

# get performance dataframe
perf <- get_performance(ind)

# plot to file
pdf(file='~/Dropbox/hephys/media/figs/roc.pdf', paper='USr', width=11, height=8.5)
plotroc(perf)
dev.off()

######## code to plot heatmap of regression coefficients ##########
df <- extract_coeffs(fitobjs[[ind]])

# useful summary stats
band_stats <- ddply(df, ~band, summarize, mean = mean(value), 
    mean_abs = mean(abs(value)), std = sd(value), std_abs = sd(abs(value)))

plt <- plot_coefficient_grid(df)
print(plt)