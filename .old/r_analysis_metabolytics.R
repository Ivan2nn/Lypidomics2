library('openxlsx')
library(mdatools)
library("tidyr")
library('pls')

setwd('~/programming/python/matteoUni/Lypidomics_01_2023/')

data <- read.xlsx('data/2022_12_Laezza_Bernabucci_Mouse_Brain_COMBINED_METABOLOMICS.xlsx')

length(unique(data$Sample.Submission.Date))

df_sub <- data[data$Tissue.Type == 'cortex' & data$PND == 30,]

df_sub$logNormalizedArea <- log10(df_sub$`Normalized.Peak.Area.(Peak.Area.of.Metabolite./.(Peak.Area.of.Internal.Standard.*.Tisue.weight))`)

# Rimuoviamo i campioni che hanno Normalized area  = 0.000000

df_sub_clean = df_sub[df_sub$`Normalized.Peak.Area.(Peak.Area.of.Metabolite./.(Peak.Area.of.Internal.Standard.*.Tisue.weight))` > 0.000000,]

# Prendiamo un subset del dataframe,solo le colonne che ci possono servire

#df_subset_clean <- df_sub_clean[c('Sample.Name', 'Treatment', 'Metabolite', 'logNormalizedArea')]
df_subset_clean <- df_sub_clean[c('Sample.Name', 'Treatment', 'Metabolite', 'Normalized.Peak.Area.(Peak.Area.of.Metabolite./.(Peak.Area.of.Internal.Standard.*.Tisue.weight))')]

# Facciamo un reshapeing dei dati. I metaboliti sono ripetuti, a noi interessa ogni metabolite (come fosse una feature).

#df_subset_reshaped <- spread(df_subset_clean, key = "Metabolite", value = 'logNormalizedArea')
df_subset_reshaped <- spread(df_subset_clean, key = "Metabolite", value = 'Normalized.Peak.Area.(Peak.Area.of.Metabolite./.(Peak.Area.of.Internal.Standard.*.Tisue.weight))')
# Tohliamo le colonne, cioè i metaboliti con dei NAN nelle colonne

df_subset_reshaped_clean = df_subset_reshaped[, colSums(is.na(df_subset_reshaped)) == 0]

# La colonna Treatment non deve contnedere Stringhe ma in R devono esserci i tipi FACTOR, cos' da poter costurire il mdoello

df_subset_reshaped_clean$Treatment <- as.numeric(as.factor(df_subset_reshaped_clean$Treatment))

# Costruiamo la matrice dei valori dei metaboliti e il vettore delle classi da passare al modello

yl.cal = df_subset_reshaped_clean$Treatment
xl.cal = df_subset_reshaped_clean[,3:ncol(df_subset_reshaped_clean)]

model2 = plsda(xl.cal, yl.cal, ncomp = 2, cv = 1, info = 'Lypydomics data example')
model <- plsr(Treatment~. - Sample.Name , data=df_subset_reshaped_clean, scale=TRUE, validation="CV")

summary(model)

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

plotScores(model2, c(1,2))

par(mfrow = c(2, 2))
plotXYScores(model2, ncomp = 2)
plotYVariance(model2, ncomp = 2)
plotXResiduals(model2, ncomp = 2)
plotRegcoeffs(model2, ncomp = 2, ny = 2)
par(mfrow = c(1, 1))
plot.plsda(model2)

summary(model2)

p <- ggplot(dfscores, aes(x=`Comp 1`,y=`Comp 2`, color = Ycalib)) + geom_point()
