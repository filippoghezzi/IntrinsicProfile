library(lmerTest)
library(car)
library(ggplot2)
library(FSA)

data = read.table("C://Users//Butt Lab//Documents//GitHub//IntrinsicProfile//IntrinsicProfiles_Results.csv", sep = ",", header = TRUE)

data$CellID = as.factor(data$CellID)
data$MouseDevelopment = factor(data$MouseDevelopment, levels = c("P5-P8", "P9-P13","P14-P18"))
data$MouseID = as.factor(data$MouseID)
data$BrainArea = as.factor(data$BrainArea)
data$CellTarget = as.factor(data$CellTarget)
data$CellType = as.factor(data$CellType)
data$CellLayer = as.factor(data$CellLayer)
summary(data)

###################### L4 Pyr analysis #######################################################
df = subset(data, BrainArea=="V1" & CellType=="Pyr" & CellLayer=="4")
summary(df)


### Rin ###
boxplot(Rin ~ MouseDevelopment, data = df)
model = aov(Rin ~ MouseDevelopment, data = df)

shapiro.test(df$Rin)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Rin ~ MouseDevelopment, data=df)

kruskal.test(Rin ~ MouseDevelopment, data=df)
dunnTest(Rin ~ MouseDevelopment, data=df, method="bh") 


### RMP ###
boxplot(RMP ~ MouseDevelopment, data = df)
model = aov(RMP ~ MouseDevelopment, data = df)

shapiro.test(df$RMP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(RMP ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### Tau ###
boxplot(Tau ~ MouseDevelopment, data = df)
model = aov(Tau ~ MouseDevelopment, data = df)

shapiro.test(df$Tau)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Tau ~ MouseDevelopment, data=df)

summary(model)


### Rheobase ###
boxplot(Rheobase ~ MouseDevelopment, data = df)
model = aov(Rheobase ~ MouseDevelopment, data = df)

shapiro.test(df$Rheobase)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Rheobase ~ MouseDevelopment, data=df)

kruskal.test(Rheobase ~ MouseDevelopment, data=df)


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

kruskal.test(MaxFiringFrequency ~ MouseDevelopment, data=df)
dunnTest(MaxFiringFrequency ~ MouseDevelopment, data=df, method="bh") 


### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

kruskal.test(AdaptationIndex ~ MouseDevelopment, data=df)


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

kruskal.test(VoltageSag ~ MouseDevelopment, data=df)


### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

### AP Half-Width ###
boxplot(AP_HalfWidth ~ MouseDevelopment, data = df)
model = aov(AP_HalfWidth ~ MouseDevelopment, data = df)

shapiro.test(df$AP_HalfWidth)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_HalfWidth ~ MouseDevelopment, data=df)

kruskal.test(AP_HalfWidth ~ MouseDevelopment, data=df)
dunnTest(AP_HalfWidth ~ MouseDevelopment, data=df, method="bh") 


### AP Height ###
boxplot(AP_Height ~ MouseDevelopment, data = df)
model = aov(AP_Height ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Height)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Height ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### AP AHP ###
boxplot(AP_AHP ~ MouseDevelopment, data = df)
model = aov(AP_AHP ~ MouseDevelopment, data = df)

shapiro.test(df$AP_AHP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_AHP ~ MouseDevelopment, data=df)

summary(model)
###################################################################################################