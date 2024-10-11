library(readxl)
library(tidyverse)
library(pROC)
library(modEvA)
library(caret)
library(ggplot2)
library(ggsci)
library(scales)
library(vegan)
library(ape)
library(Rtsne)
library(xlsx)
library(randomForest)
library(multiROC)
library(datasets)
library(dplyr)
library(reshape2)
library(ggpicrust2)

data_sp = species[species$tax %in% mark$tax,]
rownames(data_sp) = gsub("^.*\\;s__", "", rownames(data_sp))
data_sp = select(data_sp, -tax)
data_sp = as.data.frame(t(data_sp))
rownames(data_sp) = gsub("\\.", "-", rownames(data_sp))
data_sp = merge(data_sp, dplyr::select(metadata,c("SampleID","response")), by.x = "row.names", by.y = "SampleID")
rownames(data_sp) = data_sp$Row.names
data_sp = dplyr::select(data_sp, -Row.names)

#######Divide into training cohort and validation cohort######
train_data_index = createDataPartition(data_sp$response, p = 0.8, list = F)
train_data = data_sp[train_data_index,]
test_data = data_sp[-train_data_index,]

######Constructing a random forest model and determining the optimal model through  five fold cross test#####
set.seed(13)
fit.rf = train(response~., data = train_data, method = "rf",trControl = trainControl(method = "CV", number = 5))

rf.test = predict(fit.rf$finalModel, newdata = test_data, type = "class")
rf.cf = caret::confusionMatrix(as.factor(rf.test1),as.factor(test_data$response))

#######ROC#####
rf.test_ROC1 = predict(fit.rf$finalModel, newdata = test_data, type = "prob")
rf.test_ROC1 = data.frame(rf.test_ROC1)
colnames(rf.test_ROC1) <- paste0(colnames(rf.test_ROC1), "_pred_RF")

dmy = dummyVars(~ response, data = test_data)
true_label = data.frame(predict(dmy, newdata = test_data))
colnames(true_label) = c("PD","PR", "SD")
colnames(true_label) <- paste0(colnames(true_label), "_true")
final_df <- cbind(true_label, rf.test_ROC1)


roc_res <- multi_roc(final_df, force_diag = F)
pr_res <- multi_pr(final_df, force_diag = F)
plot_roc_df <- plot_roc_data(roc_res)
plot_pr_df <- plot_pr_data(pr_res)

p = ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))
#####AUC#####
AUC = c()
for(i in unique(plot_roc_df$Group)){
  tmp = plot_roc_df[plot_roc_df$Group == i,]
  a = mean(tmp$AUC)
  AUC = c(AUC, i, a)
}
AUC
########Importance of Features#####
importance = as.data.frame(importance(fit.rf$finalModel, type = 2))
importance$feature = rownames(importance)
importance = importance %>%mutate(MeanDecreaseGini = round(MeanDecreaseGini, 2))
importance = importance %>%top_n(20, MeanDecreaseGini)
p1 = ggplot()+geom_bar(data = importance,mapping=aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini, fill = "grey" ,color = "grey"),size = 1.2,alpha=0.7,
                       position="dodge", stat="identity",width = 0.85)+ 
  scale_fill_manual(values=c("grey")) +
  scale_color_manual(values=c("grey")) +
  geom_text(
    data = importance,
    aes(x = feature, y =  MeanDecreaseGini,label = format(MeanDecreaseGini, nsmall = 2)),
    position = position_dodge(width = 0.85),
    hjust = -0.15,
    size = 3
    
  ) +
  coord_flip() +
  theme_bw()+ labs(y="Num of MAGs", x="")+
  scale_y_continuous(expand = c(0, 0), limit = c(0, 5))+
  #annotate("rect", xmin = 0, xmax =6.48,  ymin = -0.1, ymax = 3.2, alpha = 0.2,fill="#FAE3AD") +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(strip.background = element_rect(fill=c("#FFF6E1")))+
  theme(axis.text=element_text(colour='black',size=8, face = "bold"))+
  theme(axis.text=element_text(colour='black',size=8, face = "bold"))+
  theme(legend.position = "none")






