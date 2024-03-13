rm(list = ls(all.names = TRUE))

library(ggplot2)

summary_info <- read.table('~/Drive/ncsu/archaea_chapter/data/cml_out_parsed/summary_tab.tsv', header = TRUE,
                           stringsAsFactors = FALSE, sep = '\t')

# Remove one gene which has a really high evolutionary rate:
summary_info <- summary_info[-which(summary_info$locus == 'fam2396'), ]

# Number of polymorphisms vs omega.
num.poly_vs_omega <- ggplot(data = summary_info, aes(x = summary_info$num_polymorphic_sites, y = log10(omega))) +
  geom_point(shape = 1) +
  theme_bw() +
  xlab('Number of polymorphic sites') +
  ylab(expression(paste("log"[10], "(", omega, ")"))) +
  ggtitle('') +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(filename = '~/Drive/ncsu/archaea_chapter/num.poly_vs_omega.pdf',
       plot = num.poly_vs_omega,
       width = 4, height = 4, units = 'in',
       dpi = 600)
