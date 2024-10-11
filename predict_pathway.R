kegg_abundance = ko2kegg_abundance(data = ko_abundance)
  
#########Analysis of Abundance Differences#######
daa_results_df <- pathway_daa(kegg_abundance, metadata = metadata,
                              group = "response", daa_method = "DESeq2")
  
daa_results_df = daa_results_df[!is.na(daa_results_df$p_values),]
daa_results_df = daa_results_df[daa_results_df$p_adjust <= 0.05,]
daa_annotated_results = pathway_annotation(pathway = "KO",
                                            daa_results_df = daa_results_df, ko_to_kegg = TRUE)
  
data_ko = kegg_abundance[daa_annotated_results$feature,]
data_ko = merge(data_ko, select(daa_annotated_results,c("feature","pathway_name")), by.x = "row.names", by.y = "feature")
data_ko= data_ko[-c(3,6,9),]
rownames(data_ko) = data_ko$pathway_name
data_ko = select(data_ko, -c("Row.names","pathway_name"))
data_ko = as.data.frame(t(data_ko))
  
  
data_ko = merge(data_ko, dplyr::select(metadata,c("SampleID","response")), by.x = "row.names", by.y = "SampleID")
rownames(data_ko) = data$Row.names
data_ko = select(data_ko, -Row.names)
colnames(data_ko) = gsub(" ", "_", colnames(data_ko))
colnames(data_ko) = gsub("-", "_", colnames(data_ko))

######Constructing a random forest model and determining the optimal model through  five fold cross test#####
train_data_index = createDataPartition(data_ko$response, p = 0.8, list = F)
train_data = data_ko[train_data_index,]
test_data = data_ko[-train_data_index,]

set.seed(13)
fit.rf = train(response~., data = train_data, method = "rf",trControl = trainControl(method = "CV", number = 5))
fit.rf$finalModel
rf.test1 = predict(fit.rf$finalModel, newdata = test_data, type = "class")
rf.cf1 = caret::confusionMatrix(as.factor(rf.test1),as.factor(test_data$response))
rf.cf1

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
p
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
  scale_y_continuous(expand = c(0, 0), limit = c(0, 10))+
  #annotate("rect", xmin = 0, xmax =6.48,  ymin = -0.1, ymax = 3.2, alpha = 0.2,fill="#FAE3AD") +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(strip.background = element_rect(fill=c("#FFF6E1")))+
  theme(axis.text=element_text(colour='black',size=8, face = "bold"))+
  theme(axis.text=element_text(colour='black',size=8, face = "bold"))+
  theme(legend.position = "none")
