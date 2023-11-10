# Predict age 
require(lumi)
require(wateRmelon)
require(glmnet) 
library(readxl)
library(ENmix)
library(dplyr)
setwd("Daniel/Biomed_AI_Jun2022")



# File locations for reading in each dataset
lbc21 = dir("../../Biomedical_AI_2022/LBC21_Norm_RDS_Files/", pattern="Methyl21", full.names=T)
lbc21 = lbc21[-grep(".Noob|.Raw", lbc21)] # Remove earlier version of Noob normalisation and Raw

lbc36 = dir("../../Biomedical_AI_2022/LBC36_Norm_RDS_Files/", pattern="Methyl36", full.names=T)
lbc36 = lbc36[-grep(".Raw", lbc36)]

lbccomb = dir("../../Biomedical_AI_2022/LBCComb_Norm_RDS_Files/", pattern="MethylComb", full.names=T)
lbccomb = lbccomb[-grep(".Raw", lbccomb)]

methods_21 = methods_36 = methods_comb = r2_21 = r2_36 = r2_comb = output_final_21 = output_final_36 = output_final_comb = list()

targets = read.csv("../../LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")

library(foreign)


# Predict in LBC36
dat <- readRDS("../../Biomedical_AI_2022/LBC36_Norm_RDS_Files/Methyl36.dasen.rds")

for(i in lbc36){
  print(i)
  norm = gsub(".*//|Methyl36.|.rds", "", i)

  methods_36[[norm]] = readRDS(i)
  print(paste0(i, " ", nrow(methods_36[[norm]])))

  if(norm=="BMIQ"|norm=="fuks"|norm=="Funnorm"){
    betas = t(methods_36[[norm]])
  } else if(norm=="tost") {
    betas = t(methods_36[[norm]]$beta)
    } else {
  betas = t(getBeta(methods_36[[norm]])) # Get betas and transpose
}


  if(any(is.infinite(betas)) | any(is.na(betas))){ # change infinites to NA
  for(j in 1:ncol(betas)){
    betas[which(is.infinite(betas[,j])),j] = NA
  #   betas[which(is.na(betas[,j])),j] = median(betas[,j], na.rm=T)
    }
  }


## predict age
betas = betas[colnames(dat),]
stopifnot(isTRUE(all.equal(colnames(dat), rownames(betas))))
dat$methylAge <- methyAge(beta=t(betas), fastImputation=F,  normalize=F)

## Format Output
output <- data.frame(dat$ID , dat$Basename , dat$methylAge)
output_merge = merge(targets , output , by.x = 'Basename' , by.y = 'dat.Basename')
output_final_36[[norm]] = output_merge[,c('Basename','ID', 'age' ,'sex' ,'mAge_Hovath')]
colnames(output_final_36[[norm]])[5] <- 'Predicted_Age'
}
saveRDS(output_final_36, file="lbc36_predicted_values_Age_postreview.rds")


# correlations between predicted and actual age

# Predict in LBC21
dat <- readRDS("../../Biomedical_AI_2022/LBC21_Norm_RDS_Files/Methyl21.dasen.rds")
for(i in lbc21){

  norm = gsub(".*//|Methyl21.|.rds", "", i)

  methods_21[[norm]] = readRDS(i)
    print(paste0(i, " ", nrow(methods_21[[norm]])))

if(norm=="BMIQ"| norm=="fuks"| norm=="Funnorm"){
    betas = t(methods_21[[norm]])
  } else if(norm=="tost") {
    betas = t(methods_21[[norm]]$beta)
    } else {
  betas = t(getBeta(methods_21[[norm]])) # Get betas and transpose
}  
# Get betas and transpose
  if(any(is.infinite(betas)) | any(is.na(betas))){ # change infinites to NA
  for(j in 1:ncol(betas)){
    betas[which(is.infinite(betas[,j])),j] = NA
   #  betas[which(is.na(betas[,j])),j] = median(betas[,j], na.rm=T)
    }
    }
  


## predict age
betas = betas[colnames(dat),]
stopifnot(isTRUE(all.equal(colnames(dat), rownames(betas))))
dat$methylAge <- methyAge(beta=t(betas), fastImputation=F,  normalize=F)

## Format Output
output <- data.frame(dat$ID , dat$Basename , dat$methylAge)
output_merge = merge(targets , output , by.x = 'Basename' , by.y = 'dat.Basename')
output_final_21[[norm]] = output_merge[,c('Basename','ID', 'age' ,'sex' ,'mAge_Hovath')]
colnames(output_final_21[[norm]])[5] <- 'Predicted_Age'

# m1 = summary(lm(log(bmi) ~ age + as.factor(sex) + Predicted_Age, data=output_final_21[[norm]]))$r.squared * 100
# m0 = summary(lm(log(bmi) ~ age + as.factor(sex), data=output_final_21[[norm]]))$r.squared * 100
# r2_21[[norm]] = m1-m0
}

saveRDS(output_final_21, file="lbc21_predicted_values_Age_postreview.rds")
# saveRDS(r2_21, file="lbc21_r2_postreview.rds")

# Predict in LBCComb
dat <- readRDS("../../Biomedical_AI_2022/LBCComb_Norm_RDS_Files/MethylComb.dasen.rds")
for(i in lbccomb){
  norm = gsub(".*//|MethylComb.|.rds", "", i)
  methods_comb[[norm]] = readRDS(i)
    print(paste0(i, " ", nrow(methods_comb[[norm]])))
if(norm=="BMIQ"|norm=="fuks"|norm=="Funnorm"){
    betas = t(methods_comb[[norm]])
  } else if(norm=="tost") {
    betas = t(methods_comb[[norm]]$beta)
    } else {
  betas = t(getBeta(methods_comb[[norm]])) # Get betas and transpose
} 
  if(any(is.infinite(betas)) | any(is.na(betas))){ # change infinites to NA
  for(j in 1:ncol(betas)){
    betas[which(is.infinite(betas[,j])),j] = NA
  #   betas[which(is.na(betas[,j])),j] = median(betas[,j], na.rm=T)
    }
  }


## predict age
betas = betas[colnames(dat),]
stopifnot(isTRUE(all.equal(colnames(dat), rownames(betas))))
dat$methylAge <- methyAge(beta=t(betas), fastImputation=F,  normalize=F)

## Format Output
output <- data.frame(dat$ID , dat$Basename , dat$methylAge)
output_merge = merge(targets , output , by.x = 'Basename' , by.y = 'dat.Basename')
output_final_comb[[norm]] = output_merge[,c('Basename','ID', 'age' ,'sex' ,'mAge_Hovath')]
colnames(output_final_comb[[norm]])[5] <- 'Predicted_Age'

# m1 = summary(lm(log(bmi) ~ age + as.factor(sex) + Predicted_Age, data=output_final_comb[[norm]]))$r.squared * 100
# m0 = summary(lm(log(bmi) ~ age + as.factor(sex), data=output_final_comb[[norm]]))$r.squared * 100
# r2_comb[[norm]] = m1-m0
}
saveRDS(output_final_comb, file="lbccomb_predicted_values_Age_postreview.rds")
# saveRDS(r2_comb, file="lbccomb_r2_postreview.rds")


cor21 = cor36 = cor_comb = list()
for(i in names(output_final_36)){
  cor21[[i]] = cor.test(output_final_21[[i]][,'Predicted_Age'], output_final_21[[i]][,'age'])
  cor36[[i]] = cor.test(output_final_36[[i]][,'Predicted_Age'], output_final_36[[i]][,'age'])
  cor_comb[[i]] = cor.test(output_final_comb[[i]][,'Predicted_Age'], output_final_comb[[i]][,'age'])
}

cor21 = lapply(cor21, function(x){c(x$estimate, lci=x$conf.int[1], uci=x$conf.int[2])})
cor36 = lapply(cor36, function(x){c(x$estimate, lci=x$conf.int[1], uci=x$conf.int[2])})
cor_comb = lapply(cor_comb, function(x){c(x$estimate, lci=x$conf.int[1], uci=x$conf.int[2])})

cor21 = do.call("rbind", cor21) %>% as.data.frame
cor36 = do.call("rbind", cor36) %>% as.data.frame
cor_comb = do.call("rbind", cor_comb) %>% as.data.frame

cor21$Cohort = 'LBC1921'
cor36$Cohort = 'LBC1936'
cor_comb$Cohort = 'Combined'

cor21$Method = rownames(cor21)
cor36$Method = rownames(cor36)
cor_comb$Method = rownames(cor_comb)


fig2 = rbind(cor21, cor36, cor_comb)
fig2$Method = gsub("swan", "SWAN", fig2$Method)
fig2$Method = gsub("noob", "Noob", fig2$Method)
fig2$Method = gsub("tost", "Tost", fig2$Method)
fig2$Method = gsub("fuks", "PBC", fig2$Method)
fig2$Cohort = factor(fig2$Cohort, levels=c("LBC1921","LBC1936","Combined"))
saveRDS(fig2, file='figure2_plotdata_age.rds')

p1 = ggplot(fig2, aes(x=factor(Method), y=cor, group=Cohort)) + 
     geom_point(size=2, aes(colour=Cohort), position = position_dodge(0.5)) + 
    scale_colour_manual(values=pal[1:3]) +
    geom_errorbar(aes(ymin = lci, ymax = uci, colour=Cohort), width = 0.2, position = position_dodge(0.5), color = "black")  +
    theme_bw() + ylim(0,1) + ylab('Correlation between predicted and actual age') + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab('Normalisation Method') +  
 theme(plot.title = element_text(hjust = 0.5) ,axis.title=element_text(size=14)) +
 ggtitle('Correlation between predicted age and actual age')

ggsave("Figure2_postreview_age.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')
print(p1)
dev.off()

# Get correlations between BMI and Episcore
cor_21 = cor_36 = cor_comb = list()

for(i in names(output_final_21)){
  cor_21[[i]] = cor(output_final_21[[i]]$age, output_final_21[[i]]$Predicted_Age, use="pairwise.complete.obs")
}

for(i in names(output_final_36)){
  cor_36[[i]] = cor(output_final_36[[i]]$age, output_final_36[[i]]$Predicted_Age, use="pairwise.complete.obs")
}

for(i in names(output_final_comb)){
  cor_comb[[i]] = cor(output_final_comb[[i]]$age, output_final_comb[[i]]$Predicted_Age, use="pairwise.complete.obs")
}

table_s1a = data.frame(Cohort = rep("LBC21", 16), 'Normalisation Method' = names(cor_21), Correlation=unlist(cor_21))
table_s1a$rank = rank(1/table_s1a$Correlation)

table_s1b = data.frame(Cohort = rep("LBC36", 16), 'Normalisation Method' = names(cor_36), Correlation=unlist(cor_36))
table_s1b$rank = rank(1/table_s1b$Correlation)

table_s1c = data.frame(Cohort = rep("Combined", 16), 'Normalisation Method' = names(cor_comb), Correlation=unlist(cor_comb))
table_s1c$rank = rank(1/table_s1c$Correlation)

table_s1 = rbind(table_s1a, table_s1b, table_s1c)
table_s1$Normalisation.Method = gsub("swan", "SWAN", table_s1$Normalisation.Method)
table_s1$Normalisation.Method = gsub("noob", "Noob", table_s1$Normalisation.Method)
table_s1$Normalisation.Method = gsub("fuks", "PBC", table_s1$Normalisation.Method)
table_s1$Normalisation.Method = gsub("tost", "Tost", table_s1$Normalisation.Method)

names(table_s1)[grep("rank", names(table_s1))] = 'Rank (Within Cohort)'
write.csv(table_s1, file="2023-03-03_Table_S1_age_postreview.csv", quote=F, row.names=F)

new_s1 = data.frame("Normalisation Method" = unique(table_s1$Normalisation.Method),
                    LBC1921 = NA,
                    LBC1936 = NA,
                    Combined = NA)

new_s1$LBC1921 = signif(table_s1[which(table_s1$Cohort=="LBC21"), ][match(new_s1$Normalisation.Method, table_s1[which(table_s1$Cohort=="LBC21"), "Normalisation.Method"]), "Correlation"], 3)
new_s1$LBC1936 = signif(table_s1[which(table_s1$Cohort=="LBC36"), ][match(new_s1$Normalisation.Method, table_s1[which(table_s1$Cohort=="LBC36"), "Normalisation.Method"]), "Correlation"], 3)
new_s1$Combined = signif(table_s1[which(table_s1$Cohort=="Combined"), ][match(new_s1$Normalisation.Method, table_s1[which(table_s1$Cohort=="Combined"), "Normalisation.Method"]), "Correlation"], 3)
# write.csv(new_s1, file="2023-03-03_Table_S1_postreview.csv", quote=F, row.names=F)

# Mean/SD correlations for main text
for(i in c("LBC21", "LBC36", "Combined")) {
  avg = mean(table_s1[which(table_s1$Cohort==i), "Correlation"])
  sdev =  sd(table_s1[which(table_s1$Cohort==i), "Correlation"])
  print(paste0(signif(avg, 3), " (SD = ", 
               signif(sdev, 3), ")"))
}
# [1] "0.121 (SD = 0.0236)"
# [1] "0.0537 (SD = 0.0301)"
# [1] "0.123 (SD = 0.0248)"



## Density comparison plots ###
#LBC21
comparison_21 = data.frame(ID = output_final_21[[1]]$ID, 
                           Sex = output_final_21[[1]]$sex,
                           Actual.Age = output_final_21[[1]]$age)

for(i in sort(names(output_final_21))){
  tmp = data.frame(ID = output_final_21[[i]]$ID,
             Pred = output_final_21[[i]]$Predicted_Age)
  names(tmp)[2] = paste0("Predicted ", i)
  comparison_21 = merge(comparison_21, tmp, by="ID")
}

#LBC36
comparison_36 = data.frame(ID = output_final_36[[1]]$ID, 
                           Sex = output_final_36[[1]]$sex,
                           Actual.Age = output_final_36[[1]]$age)

for(i in sort(names(output_final_36))){
  tmp = data.frame(ID = output_final_36[[i]]$ID,
             Pred = output_final_36[[i]]$Predicted_Age)
  names(tmp)[2] = paste0("Predicted ", i)
  comparison_36 = merge(comparison_36, tmp, by="ID")
}


# Combined
comparison_comb = data.frame(ID = output_final_comb[[1]]$ID, 
                           Sex = output_final_comb[[1]]$sex,
                           Actual.Age = output_final_comb[[1]]$age)

for(i in sort(names(output_final_comb))){
  tmp = data.frame(ID = output_final_comb[[i]]$ID,
             Pred = output_final_comb[[i]]$Predicted_Age)
  names(tmp)[2] = paste0("Predicted ", i)
  comparison_comb = merge(comparison_comb, tmp, by="ID")
}





# Figure 3 MAD between episcores separately/combined
library(dplyr)
diff21 = list()
cols = names(comparison_comb)[4:19]

# for(i in cols){
#   id = gsub("Predicted ", "", i)
#     id = gsub("fuks", "PBC", i)

# diff21[[id]] = stats::mad(comparison_21[,i], comparison_comb[match(comparison_21$ID, comparison_comb$ID),i], na.rm=T)
# }
# diff21 = do.call("rbind", diff21) %>% as.data.frame
# names(diff21) = "MAD"
# diff21$Cohort = "LBC1921"
# diff21$Method = rownames(diff21)


# diff36 = list()
# cols = names(comparison_comb)[4:19]

# for(i in cols){
#   id = gsub("Predicted ", "", i)
#       id = gsub("fuks", "PBC", i)

# diff36[[id]] = stats::mad(comparison_36[,i], comparison_comb[match(comparison_36$ID, comparison_comb$ID),i], na.rm=T)
# }
# diff36 = do.call("rbind", diff36) %>% as.data.frame
# names(diff36) = "MAD"
# diff36$Cohort = "LBC1936"
# diff36$Method = rownames(diff36)

# fig3 = rbind(diff21, diff36)
# fig3$Method = gsub("swan", "SWAN", fig3$Method)
# fig3$Method = gsub("noob", "Noob", fig3$Method)
# fig3$Method = gsub("tost", "Tost", fig3$Method)

# fig3$Method = gsub("Predicted ", "", fig3$Method)

library(wesanderson)
pal <- wes_palette("Darjeeling1", 5)
pal2 <- wes_palette("Darjeeling2", 5)

# col_LBC21=pal[1]
# col_LBC36=pal[2]     #tst2
# col_actual=pal2[5]   #tst1
# col_combined=pal[3]  #tst3
# line_width=3
#   palette = c(col_LBC21, col_LBC36)

# p2 =  ggplot(fig3, aes(fill=Cohort, y=MAD, x=reorder(Method, MAD , median))) + theme_bw() + scale_fill_manual(values=palette) + scale_color_manual(values=palette) + theme(plot.title = element_text(hjust = 0.5) ,axis.title=element_text(size=14)) +
#     geom_bar(position="dodge", stat="identity") + xlab('Normalisation Method') + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     ylab("Mean Absolute Difference") + ggtitle('Age')+
#     geom_hline(yintercept = 0, size=0.8) 
# ggsave('Figure3_age_postreview.tiff',plot=p2, width=7, height=7, dpi=300)

for(i in cols){
  id = gsub("Predicted ", "", i)
    id = gsub("fuks", "PBC", i)
diff21[[id]] = abs(comparison_21[,i] - comparison_comb[match(comparison_21$ID, comparison_comb$ID),i]) 
}
diff21 = t(do.call("rbind", diff21)) %>% as.data.frame %>%  melt
diff21$Cohort = 'LBC1921'

diff36 = list()
cols = names(comparison_comb)[4:19]

for(i in cols){
  id = gsub("Predicted ", "", i)
    id = gsub("fuks", "PBC", i)
diff36[[id]] = abs(comparison_36[,i] - comparison_comb[match(comparison_36$ID, comparison_comb$ID),i])
}
diff36 = t(do.call("rbind", diff36)) %>% as.data.frame %>%  melt
diff36$Cohort = 'LBC1936'

fig3 = rbind(diff21, diff36)
fig3$variable = gsub("Predicted ", "", fig3$variable)
fig3$variable = gsub("swan", "SWAN", fig3$variable)
fig3$variable = gsub("noob", "Noob", fig3$variable)
fig3$variable = gsub("tost", "Tost", fig3$variable)

p2 = ggplot(fig3, aes(x=reorder(variable, value , mean, na.rm=T), y=value, fill=Cohort)) + 
geom_boxplot() + 
ylab("Absolute Difference") + ggtitle('Absolute Difference between Cohort and Combined (Age)') + 
theme_bw() + 
xlab('Normalisation Method') + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
theme(plot.title = element_text(hjust = 0.5) ,axis.title=element_text(size=14)) +
scale_fill_manual(values=pal[1:2])

saveRDS(fig3, file='fig3_barplot_age.rds')
ggsave('Figure3_age_postreview.tiff',plot=p2, width=7, height=7, dpi=300)



# Repair ranking heatmap
# hm = read.csv("2022-02-06_heatmap_values.csv")
hm21 = read.table("Biomedical_AI_2022/R_Code/Normalisation_performance_LBC1921.txt", header=T) %>% mutate(Cohort="LBC1921")
hm36 = read.table("Biomedical_AI_2022/R_Code/Normalisation_performance_LBC1936.txt", header=T) %>% mutate(Cohort="LBC1936")
hmcomb = read.table("Biomedical_AI_2022/R_Code/Normalisation_performance_LBCCombined.txt", header=T) %>% mutate(Cohort="LBC Combined")
hm = rbind(hm21, hm36, hmcomb)


hm = melt(hm, measure.vars = c("Rank.of.DMRSE", "Rank.of.GCOSE", "Rank.of.Seabird", "Ranked.Means"))# %>% select(-Mean.of.Ranks)
hm$variable = gsub("Rank.of.", "", hm$variable)
hm$variable = gsub("Ranked.Means", "Overall Rank", hm$variable)
hm$Normalization = gsub("tost", "Tost", hm$Normalization)
hm$Normalization = gsub("fuks", "PBC", hm$Normalization)
hm$variable = factor(hm$variable, levels = c("DMRSE", "GCOSE", "Seabird" ,"Overall Rank"))
hm$Normalization = factor(hm$Normalization, levels=c("Raw", unique(sort(hm$Normalization[-which(hm$Normalization=="Raw")]))))
fig_s1 = ggplot(hm, aes(x = Normalization, y = Cohort, fill = value)) +
  geom_tile() + 
  xlab("Normalisation Method") +
  ylab("Dataset") +
  coord_fixed() +
 geom_text(aes(label=value), size=3, colour="white") + 
 facet_wrap(~variable, nrow=4) + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


ggsave('FigureS1_postreview.tiff',plot=fig_s1, width=7, height=5, dpi=300)



# projection
library(haven)
#### TO BE CHECKED ###
descale <- function(test_bmi, test_age, test_sex){
lm_tst <- lm(log(pheno$bmi) ~ pheno$age + pheno$sex, na.action = "na.exclude")
mean_lm <- mean(resid(lm_tst) , na.rm = TRUE)
std_lm <- sd(resid(lm_tst) , na.rm = TRUE)
descale <- mean_lm + test_bmi*std_lm
pred <- coef(lm_tst)[1] + coef(lm_tst)[2]*test_age +   
        coef(lm_tst)[3]*test_sex   # Check this?
return(exp(descale + pred))
}

for(i in names(output_final_21)){
  output_final_21[[i]]$sex2 = gsub("M", "1", output_final_21[[i]]$sex)
  output_final_21[[i]]$sex2 = gsub("F", "0", output_final_21[[i]]$sex2)
  output_final_21[[i]]$sex2 = as.numeric(output_final_21[[i]]$sex2)
  output_final_21[[i]]$bmi_descaled = descale(output_final_21[[i]]$Predicted_Age, test_age=output_final_21[[i]]$age, test_sex=output_final_21[[i]]$sex2)
  output_final_21[[i]]$Norm.Method = i

}


for(i in names(output_final_36)){
  output_final_36[[i]]$sex2 = gsub("M", "1", output_final_36[[i]]$sex)
  output_final_36[[i]]$sex2 = gsub("F", "0", output_final_36[[i]]$sex2)
  output_final_36[[i]]$sex2 = as.numeric(output_final_36[[i]]$sex2)
  output_final_36[[i]]$bmi_descaled = descale(output_final_36[[i]]$Predicted_Age, test_age=output_final_36[[i]]$age, test_sex=output_final_36[[i]]$sex2)
  output_final_36[[i]]$Norm.Method = i
}

for(i in names(output_final_comb)){
  output_final_comb[[i]]$sex2 = gsub("M", "1", output_final_comb[[i]]$sex)
  output_final_comb[[i]]$sex2 = gsub("F", "0", output_final_comb[[i]]$sex2)
  output_final_comb[[i]]$sex2 = as.numeric(output_final_comb[[i]]$sex2)
  output_final_comb[[i]]$bmi_descaled = descale(output_final_comb[[i]]$Predicted_Age, test_age=output_final_comb[[i]]$age, test_sex=output_final_comb[[i]]$sex2)
  output_final_comb[[i]]$Norm.Method = i
}



# Figure S2 (BMI Plot):
i = 1
pheno = readRDS("phenos.rds")

bmi_21 = data.frame(id=output_final_21[[i]]$ID, bmi=output_final_21[[i]]$bmi, Cohort="LBC1921")
bmi_36 = data.frame(id=output_final_36[[i]]$ID, bmi=output_final_36[[i]]$bmi, Cohort="LBC1936")
bmi_gs = data.frame(id=pheno$id, bmi=pheno$bmi, Cohort="GS")

s20 = rbind(bmi_21, bmi_36, bmi_gs)

p_s20 = ggplot(s20, aes(x=bmi, fill=Cohort)) + 
 xlab(expression('BMI (kg/'*m^2*')')) + 
 ylab("Density") +
geom_density(alpha=0.5) 

ggsave("Figure_S2_postreview.tiff", plot = p_s20, units="in", width=8, height=8, dpi=300, compression = 'lzw')



# Appendix Figure 3 MAD between episcores separately/combined
library(dplyr)
descale <- function(test_bmi, test_age, test_sex){
lm_tst <- lm(log(pheno$bmi) ~ pheno$age + pheno$sex, na.action = "na.exclude")
mean_lm <- mean(resid(lm_tst) , na.rm = TRUE)
std_lm <- sd(resid(lm_tst) , na.rm = TRUE)
descale <- mean_lm + test_bmi*std_lm
pred <- coef(lm_tst)[1] + coef(lm_tst)[2]*test_age +   
        coef(lm_tst)[3]*test_sex   # Check this?
return(exp(descale + pred))
}



# Supplementary figure (performance metrics)
library(dplyr)
fig_s3_a = read.table("../../Biomedical_AI_2022/R_Code/Normalisation_performance_plotdata_LBC1921.txt", header=T) %>% mutate(Group="LBC1921")
fig_s3_b = read.table("../../Biomedical_AI_2022/R_Code/Normalisation_performance_plotdata_LBC1936.txt", header=T) %>% mutate(Group="LBC1936")
fig_s3_c = read.table("../../Biomedical_AI_2022/R_Code/Normalisation_performance_plotdata_LBCCombined.txt", header=T) %>% mutate(Group="Combined")
fig_s3 = rbind(fig_s3_a, fig_s3_b, fig_s3_c)

# DMRSE = x10-3
# gcose = x10-5
# SEABIRD = x10-2

fig_s3$DMRSE = fig_s3$DMRSE * 1000
fig_s3$GCOSE = fig_s3$GCOSE * 100000
fig_s3$SEABIRD = fig_s3$SEABIRD * 100

s3_pd = melt(fig_s3)


s3_pd$facets = factor(s3_pd$variable, labels = c(
    "DMRSE~(x10^{-4})", 
    "GCOSE~(x10^{-5})", 
    "Seabird~(x10^{-2})"))
s3_pd$Method = gsub("swan", "SWAN", s3_pd$Method)
s3_pd$Method = gsub("noob", "Noob", s3_pd$Method)
s3_pd$Method = gsub("fuks", "PBC", s3_pd$Method)
s3_pd$Method = gsub("tost", "Tost", s3_pd$Method)


figure_s3 = ggplot(s3_pd, aes(x=Method, y=value, colour=Group)) + 
geom_point() + 
facet_wrap(.~facets, labeller=label_parsed) + 
ylim(0,7) + 
ylab("Score") +
xlab("Normalisation Method") + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave('Metrics_plot_postreview.tiff',plot=figure_s3, width=8, height=6, dpi=300)
