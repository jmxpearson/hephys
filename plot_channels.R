# miscellaneous lfp channel plots
adir <- '~/Dropbox/hephys/analysis'  # analysis directory
ddir <- '~/data/bartc'  # data directory
cfile <- 'lfp_channel_file.csv'
chanlist <- read.csv(paste(adir, cfile, sep='/'), header=F)
chanlist <- chanlist[, 1:2]
chanlist <- chanlist[!duplicated(chanlist),]
numunits <- dim(chanlist)[1]

library(ggplot2)
library(reshape)

ind = 10
fname = paste(chanlist[ind, 1], chanlist[ind, 2], 'chanmeans', 'csv', sep='.')
# load analysis output
dat <- read.csv(file=paste(ddir, fname, sep='/'), header = F)
names(dat) = c('time', 1:(dim(dat)[2] - 1))
dat <- melt(dat, id.vars=c('time'))

pdf(file='~/Dropbox/hephys/media/figs/chanplot.pdf', paper='USr', width=11, height=8.5)
plt <- ggplot(dat, aes(time, value))
plt + geom_line(aes(color=variable)) + scale_colour_hue(h=c(90, 180), guide=F) + labs(x='\nTime From Stopping (s)', y='Normalized Channel Power\n') +theme_bw() + theme(plot.title = element_text(lineheight=1, size=36), axis.text=element_text(color='black', size=12),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20)
  )
dev.off()