library('pls')
library("tidyr")

data <- read.csv("data/clean_data.csv")

df_sub <- data[data$Tissue.Type == 'cortex' & data$PND == 30,]

df_sub$logNormalizedArea <- log10(df_sub$Normalized.Peak.Area..Peak.Area.of.Lipid.Species....Peak.Area.of.Internal.Standard...Tissue.Weight..)

df_subset <- df_sub[c('Sample.Number','Treatment','Individual.Lipid.Species','logNormalizedArea')]

# reshaping the data
#df_subset_reshaped <- reshape(df_subset, idvar = "Individual.Lipid.Species",timevar = "logNormalizedArea", direction = "wide")

df_subset_reshaped <- spread(df_subset, key = "Individual.Lipid.Species", value = 'logNormalizedArea')

sapply(df_subset_reshaped, anyNA)

df_subset_reshaped_clean = df_subset_reshaped[, colSums(is.na(df_subset_reshaped)) == 0]

df_subset_reshaped_clean$Treatment <- as.factor(df_subset_reshaped_clean$Treatment)

xl.cal = df_subset_reshaped_clean$Treatment
yl.cal = df_subset_reshaped_clean[,3:ncol(df_subset_reshaped_clean)]

#model1 <- plsr(xl.cal ~ yl.cal, data=df_subset_reshaped_clean, scale=TRUE)
model2 = plsda(yl.cal, xl.cal, ncomp = 2, cv = 1, info = 'Lypydomics data example')

summary(model2)
plot(model2)

## 2. Show performance plots for a model
par(mfrow = c(2, 2))
plotSpecificity(model2)
plotSensitivity(model2)
plotMisclassified(model2)
plotMisclassified(model2, nc = 2)
par(mfrow = c(1, 1))

## 3. Show both class and y values predictions
par(mfrow = c(2, 2))
plotPredictions(model2)
plotPredictions(model2, res = "cal", ncomp = 2, nc = 2)
plotPredictions(structure(model2, class = "regmodel"))
plotPredictions(structure(model2, class = "regmodel"), ncomp = 2, ny = 2)
par(mfrow = c(1, 1))

## 4. All plots from ordinary PLS can be used, e.g.:
par(mfrow = c(2, 2))
plotXYScores(model2, ncomp = 2)
plotYVariance(model2, ncomp = 2)
plotXResiduals(model2, ncomp = 2)
plotRegcoeffs(model2, ny = 2)
par(mfrow = c(1, 1))

### Examples for PLS-DA model class

library(mdatools)

## 1. Make a PLS-DA model with full cross-validation and show model overview

# make a calibration set from iris data (3 classes)
# use names of classes as class vector
x.cal = iris[seq(1, nrow(iris), 2), 1:4]
c.cal = iris[seq(1, nrow(iris), 2), 5]

c.cal
x.cal
iris

model = plsda(x.cal, c.cal, ncomp = 3, cv = 1, info = 'IRIS data example')
model = selectCompNum(model, 1)

# show summary and basic model plots
# misclassification will be shown only for first class
summary(model)
plot(model)

# summary and model plots for second class
summary(model, nc = 2)
plot(model, nc = 2)

# summary and model plot for specific class and number of components
summary(model, nc = 3, ncomp = 3)
plot(model, nc = 3, ncomp = 3)

## 2. Show performance plots for a model
par(mfrow = c(2, 2))
plotSpecificity(model)
plotSensitivity(model)
plotMisclassified(model)
plotMisclassified(model, nc = 2)
par(mfrow = c(1, 1))

## 3. Show both class and y values predictions
par(mfrow = c(2, 2))
plotPredictions(model)
plotPredictions(model, res = "cal", ncomp = 2, nc = 2)
plotPredictions(structure(model, class = "regmodel"))
plotPredictions(structure(model, class = "regmodel"), ncomp = 2, ny = 2)
par(mfrow = c(1, 1))

## 4. All plots from ordinary PLS can be used, e.g.:
par(mfrow = c(2, 2))
plotXYScores(model)
plotYVariance(model)
plotXResiduals(model)
plotRegcoeffs(model, ny = 2)
par(mfrow = c(1, 1))