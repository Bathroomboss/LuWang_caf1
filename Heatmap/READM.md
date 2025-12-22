# a whole pipeline generate ChIPseq profile and heatmap ---
# step 1 mapping to TE
sh TEChip.pbs

# step 2 compare bw (optional)
bigwigCompare -p 24 -b1 ${IP}_TE.bw -b2 ${Input}_TE.bw --skipZeroOverZero --operation log2 -o ${sample}_WT_vs_Input_log2.bw

# step 3 generate matrix for plot 
computeMatrix scale-regions -R TE.bed -S ${sample}_WT_vs_Input_log2.bw --missingDataAsZero --skipZeros -o ${sample}_WT_vs_Input_log2_atTE.gz 
zcat ${sample}_WT_vs_Input_log2_atTE.gz  |awk -v OFS=',' '{$1=$1}1' > ${sample}_WT_vs_Input_log2_atTE.csv

# step 4 plot profile / heatmap
# in R 
data <- read.csv('${sample}_WT_vs_Input_log2_atTE.csv',header = F)
# plot profile
plot(smooth.spline(colMeans(data[,7:106]),spar=0.1),type="l",lwd=2,ylab="CHIP intensity",
     main="IP_VS_Input_atTE",xaxt='n',xlab="",cex.axis=1.3,cex.lab=1.3,col="#e31a1c")

# plot heatmap
library(pheatmap)
my_colors <- c(colorRampPalette(c('blue', 'white'))(50),  
               colorRampPalette(c('white', 'red'))(50))  
breaks_custom <- c(seq(-3, 0, length.out = 51),  
                   seq(0, 2, length.out = 51)[-1])  
hm_data <- data[, 7:106]  
pheatmap(hm_data, 
         color = my_colors,
         border_color = "NA", 
         breaks = breaks_custom,
         scale = "none", 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         show_rownames = FALSE, show_colnames = FALSE, 
         main = "")
