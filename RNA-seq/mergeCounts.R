# get work directory
wd.curr <- getwd()
# assume run in <project_dir>/src
# set work directory as <project_dir>
setwd(dirname(wd.curr))
# L3_rep1-3, WP_rep1-3
id <- paste0(rep(c('L3_rep','WP_rep'), each=3), 1:3)
paths <- paste0('results/featurecounts/',id,'_counts.txt')

counts_df <- Reduce(cbind,lapply(paths, function(i) {
  read.table(i,sep = '\t', row.names = 1, skip = 2,
             colClasses = c('character',rep('NULL',5),'integer'))
}))
colnames(counts_df) <- id
write.table(counts_df, file='results/Counts.csv', sep=',', col.names=NA)
