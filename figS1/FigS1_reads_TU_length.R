# this part is to find the degree of reads number contributes to annotated TU length
# Rui Shao
# 2020 Feb

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
system("mkdir figs")
system("mkdir data")

# ----------------------------------------------------------------------------------
files = list.files("../data/TU_anno/other/WT26", ".gtf", full.names = T)
file.names = c("10%","20%", "40%", "60%", "80%", "100%")

ReadsNum = read.table('../data/TU_anno/other/WT26/result.out')
ReadsNum = ReadsNum$V1
ReadsNum = c(ReadsNum, ReadsNum[5] / 8*10)

dat$ReadsNum = ReadsNum / 1e+6 # to million reads

dat = data.frame()
for( i in seq_along(files)){
  temp.gr = importRanges(files[i])
  temp.table = aggregate(width(temp.gr), 
                         list(temp.gr$location),
                         sum)
  colnames(temp.table) = c("TU_location", "Total_width")
  temp.table$ReadsNum = ReadsNum[i] / 1e+6 # to million reads
  dat = rbind(dat, temp.table)
}
# get percentage of TU width comparing to the full size bam file
dat$Width_percent = 0
for( i in unique(dat$TU_location))
{
  dat[dat$TU_location == i, "Width_percent"] = 
    dat[dat$TU_location == i, "Total_width"] / 
    dat[dat$TU_location == i & dat$ReadsNum == max(ReadsNum/1e+6), "Total_width"] * 100
}

write.table(dat, file = "data/FigS1_reads_TU_length.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

dat %>% filter(TU_location %in%
                 c("protein_coding", "intergenic", "antisense", "uaRNA", "conRNA")) %>%
  ggplot(aes(x = ReadsNum, y = Width_percent, color = TU_location)) + 
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3),
              se = FALSE, lwd=0.4) +
  scale_color_viridis_d(option = "E", direction = -1) + 
  xlab('Million reads') +
  ylab('Annotated length (%)') +
  ggtitle("Total reads splits") +
  theme_minimal() +
  theme_setting -> g

ggsave(g, filename = "FigS1.Bam_splits_Annotated_percentage.pdf", path = "figs",
       device = "pdf", width = 5, height = 4 )

