#ne moe ne nado
#(from Volans)



library(ggplot2)

### Monogyne maf 0.05 in monogyne group
## r2 > 0.05
Cn_ant_vcf <- read.table("/mono_maf0.05_r20.05_chr.ld", header = F, stringsAsFactors = F)
head(Cn_ant_vcf)

LDsub2 <- subset(Cn_ant_vcf,Cn_ant_vcf$V9>=0.81)

ggplot(Cn_ant_vcf, aes(x=V2/1000,y=V6/1000, color=V9)) + geom_point() + 
  geom_point(data=LDsub2, aes(x=V6/1000,y=V2/1000,color=V9 ))+
  scale_color_gradient(low = "#f7f8e8", high = "#D2603A", na.value = 'white', name = expression(paste("R"^2))) +
  facet_wrap(vars(V1),scales= "free") +
  theme_bw() + labs(title =expression(paste("Monogyne maf 0.05 in group, R"^2,">=0.05 v.s. R"^2,">=0.81")), x = "SNP position 1 (Kbp)", y = "SNP position 2 (Kbp)")

# Monogyne only chromosome 2, r2 >=0.05 v.s. r2>=0.64
LDchr2_sub <- subset(Cn_ant_vcf,Cn_ant_vcf$V1==2)
LDchr2_sub1 <- subset(LDchr2_sub,LDchr2_sub$V9>=0.64)

ggplot(LDchr2_sub, aes(x=V2/1000,y=V6/1000, color=V9)) + geom_point() + 
  geom_point(data=LDchr2_sub1, aes(x=V6/1000,y=V2/1000,color=V9 ))+
  scale_color_gradient(low = "#f7f8e8", high = "#D2603A", na.value = 'white', name = expression(paste("R"^2))) +
  theme_bw() + labs(title =expression(paste("Monogyne maf 0.05 in group, R"^2,">=0.05 v.s. R"^2,">=0.64, Chr2")), x = "SNP position 1 (Kbp)", y = "SNP position 2 (Kbp)")

#=========================================
# Polygyne maf 0.05 in polygyne group

Cn_ant_vcf <- read.table("/poly_maf0.05_r20.05_chr.ld", header = F, stringsAsFactors = F)
head(Cn_ant_vcf)

LDsub2 <- subset(Cn_ant_vcf,Cn_ant_vcf$V9>=0.81)

ggplot(Cn_ant_vcf, aes(x=V2/1000,y=V6/1000, color=V9)) + geom_point() + 
  geom_point(data=LDsub2, aes(x=V6/1000,y=V2/1000,color=V9 ))+
  scale_color_gradient(low = "#f7f8e8", high = "#D2603A", na.value = 'white', name = expression(paste("R"^2))) +
  facet_wrap(vars(V1),scales= "free") +
  theme_bw() + labs(title =expression(paste("Polygyne maf 0.05 in group, R"^2,">=0.05 v.s. R"^2,">=0.81")), x = "SNP position 1 (Kbp)", y = "SNP position 2 (Kbp)")

#Plygyne only chromosome 2,  
LDchr2_sub <- subset(Cn_ant_vcf,Cn_ant_vcf$V1==2)
LDchr2_sub1 <- subset(LDchr2_sub,LDchr2_sub$V9>=0.64)

ggplot(LDchr2_sub, aes(x=V2/1000,y=V6/1000, color=V9)) + geom_point() + 
  geom_point(data=LDchr2_sub1, aes(x=V6/1000,y=V2/1000,color=V9 ))+
  scale_color_gradient(low = "#f7f8e8", high = "#D2603A", na.value = 'white', name = expression(paste("R"^2))) +
  theme_bw() + labs(title =expression(paste("Polygyne maf 0.05 in group, R"^2,">=0.05 v.s. R"^2,">=0.64, Chr2")), x = "SNP position 1 (Kbp)", y = "SNP position 2 (Kbp)")

