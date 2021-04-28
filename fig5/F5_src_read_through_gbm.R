library(caret)
library(gbm)
# library(xgboost)

library(pROC)


set.seed(9527)
test.idx <- sample(nrow(dat_test_gbm), nrow(dat_test_gbm) / 5) # 20% genes as the test set

fitControl <- trainControl(method = "cv", number = 10)

if (F) {
  gbmGrid <-  expand.grid(interaction.depth = c(10, 20, 40), 
                          n.trees = (1:20)*50, 
                          shrinkage = 0.1,
                          n.minobsinnode = 20)
  
  gbmFit <- train(Read.through.length ~ .,
                  data = dat_test_gbm[-test.idx, ], 
                  method = "gbm", 
                  trControl = fitControl, 
                  tuneGrid = gbmGrid,
                  verbose = F)
  gbmFit
  plot(gbmFit) 
}

# use tuned parameters
gbmGrid2 <-  expand.grid(interaction.depth = 40, 
                         n.trees = 20*50, 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)

gbm.Test <- predict(gbmFit, 
                    newdata = dat_test_gbm[test.idx, ],
                    type="raw") # correlation 0.7543092

data.frame(x = dat_test_gbm[test.idx, 1] / 1e3, y = gbm.Test / 1e3) %>%
  ggplot(aes(x = x, y = y, color = get_dens(x, y))) +
  geom_point(cex = 0.5) +
  annotate("text", x = Inf, y = Inf,
           hjust = 1.2, vjust = 1.2, 
           label = paste0("     r = ", 0.75, "\nn = ", length(test.idx))) +
  scale_color_viridis(direction = -1) +
  xlab("Read-through distance (kb)") + ylab("Predicted read-through (kb)") +
  theme_setting +
  theme(legend.position = 'none') 
  
ggsave(filename = "FigS5_gbm_predicted_read_though.png",
       device = "png", path = "../figS5/figs",
       width = 3.5, height = 3.5)

# feature importance
gbmImp <- varImp(gbmFit, scale = FALSE)

dat_gbmImp <- gbmImp$importance 
dat_gbmImp$Feature <- rownames(dat_gbmImp)
dat_gbmImp$Feature_class <- c(rep("Gene", 5),#
                              rep("Speed", 2),#
                              rep("Remodeler", 7), #
                              rep("Accessibility", 3), #
                              rep("Domain", 4), #
                              rep("RNA/DNA", 3), #
                              "Histone acetylation", "Histone variant", "Histone acetylation", 
                              rep("Histone methylation", 5), "Histone acetylation", "Histone methylation",
                              "Histone ubiquitylation", rep("Histone variant", 2), "Histone ubiquitylation",
                              "Histone phosphorylation", rep("Histone acetylation", 2))
dat_gbmImp <- dat_gbmImp[order(dat_gbmImp$Overall, decreasing = T), ]
dat_gbmImp$Feature <- factor(dat_gbmImp$Feature, levels = dat_gbmImp$Feature)

ggplot(dat_gbmImp, aes(x = Feature, y = Overall, fill = Feature_class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Overall importance",
                     labels = function(x) format(x, scientific = TRUE)) +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "FigS9_Read_thru_gbm_feature_importance.png", device = "png",
       width = 11, height = 4, path = "../figS5/figs")

# bin read-through distance to class vector
dat_test_gbm$Read.through.length.class <- cut(dat_test_gbm$Read.through.length,
                                              breaks = c(0, 3000, 5000, 8000, 15000), 
                                              labels = letters[1:4])

fitControl <- trainControl(method = "cv", number = 10, 
                           classProbs = TRUE)

gbmGrid2 <-  expand.grid(interaction.depth = 40, 
                        n.trees = 20*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)

gbmFit2 <- train(Read.through.length.class ~ .,
                 data = dat_test_gbm[-test.idx, -1], 
                 method = "gbm", 
                 trControl = fitControl, 
                 tuneGrid = gbmGrid,
                 metric = "ROC",
                 verbose = F)


gbm.probs <- predict(gbmFit2,
                   newdata = dat_test_gbm[test.idx, ],
                   type="prob")

png("../figS5/figs/FigS5_gbm_prediction_ROC.png", width = 350, height = 350)
plot(roc(dat_test_gbm[test.idx, "Read.through.length.class"],
         gbm.probs[, "a"]), col = 1) # Area under the curve: 0.8743
plot(roc(dat_test_gbm[test.idx, "Read.through.length.class"],
         gbm.probs[, "b"]), col = 2, add = T) # Area under the curve: 0.7401
plot(roc(dat_test_gbm[test.idx, "Read.through.length.class"],
         gbm.probs[, "c"]), col = 3, add = T) # Area under the curve: 0.846
plot(roc(dat_test_gbm[test.idx, "Read.through.length.class"],
         gbm.probs[, "d"]), col = 4, add = T) # Area under the curve: 0.8137
legend(0.5, 0.3, legend = c("0-3 kb AUC 0.8743", "3-5 kb AUC 0.7401", 
                            "5-8 kb AUC 0.846", "8-15 kb AUC 0.8137" ),
       col = 1:4, lty=1, cex=0.8)
dev.off()
