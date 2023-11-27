library(reshape2)
library(philentropy)
file_path <- '/mnt/snRNA_seq/snrna/my_fig/altas'
test.clean <- load('/mnt/snRNA_seq/snrna/my_fig/rda/test.clean.rda')
seu <- test.clean
seu.df <- test.clean@meta.data

probMatrix <- data.frame(table(seu.df[, c('orig.ident', 'subtype')]))
probMatrix <- probMatrix %>% group_by(subtype) %>% mutate( prob =  Freq / sum(Freq)) %>%
  dcast(orig.ident~subtype, value.var = 'prob') 
row.names(probMatrix) <- probMatrix[,1]
probMatrix <- t(probMatrix[,-1])

baseline <- data.frame(table(seu.df[, 'orig.ident'] ))
baseline <- baseline %>% mutate(prob = Freq / sum(Freq))
row.names(baseline) <- baseline[, 1]
baseline <- t(baseline[,-c(1,2)])

spercific <- c()
for(i in 1:length(unique(seu.df$subtype))){
  spercific <- append(spercific, JSD(rbind(probMatrix[i,], baseline)))
}

sper.df <- data.frame(celltype = factor(row.names(probMatrix), levels = row.names(probMatrix)), 'spercific' = spercific)

f1_g <- ggplot(sper.df, aes(x= celltype, y = spercific)) + geom_bar(stat = 'identity') +  
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        # axis.text.y = element_blank(),
        # axis.text.x = element_blank(),
        # #axis.text = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.line.x = element_blank(),
        # axis.line.y = element_blank(),
        # plot.title = element_blank(),
        axis.text = element_text(size = 10)
  ) + coord_flip()
ggsave('test.clean.spercific_barplot.pdf', f1_g, width = 8, height = 12, units = 'cm', path = file_path)


f1_f <- ggplot(seu.df, aes(x = subtype, fill = orig.ident)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c('#8D574C', '#3D6CB4', '#2FA148', '#ff7f0e', '#D52A29'))+
  scale_x_discrete(position = 'top') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 0, vjust = 0, size = 9, lineheight = 0),
        #axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_blank(),
        legend.position = 'bottom',
        legend.margin = margin(0,0,0,0, unit = 'cm'),
        legend.box.margin = margin(0,0,0,0, unit = 'cm')
        
  ) +
  guides(fill = guide_legend(title = 'Age     ', 
                             title.theme = element_text(size = 9), 
                             direction = 'horizontal', nrow = 1, 
                             keywidth = unit(4, 'mm'),
                             keyheight = unit(4, 'mm'),
                             
                             label.theme = element_text(size = 9))) +
  NoLegend() +
  coord_flip()

ggsave('test.clean.barplot.NL.pdf', f1_f, width = 8, height = 12, units = 'cm', path = file_path)
