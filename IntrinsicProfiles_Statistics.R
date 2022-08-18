library(lmerTest)
library(car)
library(ggplot2)
library(FSA)

data = read.table("C://Users//Butt Lab//Documents//GitHub//IntrinsicProfile//IntrinsicProfiles_Results_withMaturity.csv", sep = ",", header = TRUE)

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

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

###################################################################################################
###################### L2/3 Pyr analysis #######################################################
df = subset(data, BrainArea=="V1" & CellType=="Pyr" & CellLayer=="23")
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
dunnTest(Rheobase ~ MouseDevelopment, data=df, method="bh") 


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

kruskal.test(AdaptationIndex ~ MouseDevelopment, data=df)
dunnTest(AdaptationIndex ~ MouseDevelopment, data=df, method="bh") 


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

kruskal.test(VoltageSag ~ MouseDevelopment, data=df)
dunnTest(VoltageSag ~ MouseDevelopment, data=df, method="bh") 


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

kruskal.test(AP_AHP ~ MouseDevelopment, data=df)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

kruskal.test(MaturityIndex ~ MouseDevelopment, data=df)
dunnTest(MaturityIndex ~ MouseDevelopment, data=df, method="bh") 

###################################################################################################
###################### L5 Pyr analysis #######################################################
df = subset(data, BrainArea=="V1" & CellType=="Pyr" & CellLayer=="5")
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

summary(model)
TukeyHSD(model)


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

summary(model)
TukeyHSD(model)


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

kruskal.test(AP_Height ~ MouseDevelopment, data=df)
dunnTest(AP_Height ~ MouseDevelopment, data=df, method="bh") 


### AP AHP ###
boxplot(AP_AHP ~ MouseDevelopment, data = df)
model = aov(AP_AHP ~ MouseDevelopment, data = df)

shapiro.test(df$AP_AHP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_AHP ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model) 

###################################################################################################
###################### L6 Pyr analysis #######################################################
df = subset(data, BrainArea=="V1" & CellType=="Pyr" & CellLayer=="6")
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

kruskal.test(RMP ~ MouseDevelopment, data=df)
dunnTest(RMP ~ MouseDevelopment, data=df, method="bh") 


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


### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

kruskal.test(AdaptationIndex ~ MouseDevelopment, data=df)
dunnTest(AdaptationIndex ~ MouseDevelopment, data=df, method="bh") 


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

summary(model)


### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

kruskal.test(AP_Threshold ~ MouseDevelopment, data=df)

### AP Half-Width ###
boxplot(AP_HalfWidth ~ MouseDevelopment, data = df)
model = aov(AP_HalfWidth ~ MouseDevelopment, data = df)

shapiro.test(df$AP_HalfWidth)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_HalfWidth ~ MouseDevelopment, data=df)

kruskal.test(AP_HalfWidth ~ MouseDevelopment, data=df)


### AP Height ###
boxplot(AP_Height ~ MouseDevelopment, data = df)
model = aov(AP_Height ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Height)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Height ~ MouseDevelopment, data=df)

kruskal.test(AP_Height ~ MouseDevelopment, data=df)


### AP AHP ###
boxplot(AP_AHP ~ MouseDevelopment, data = df)
model = aov(AP_AHP ~ MouseDevelopment, data = df)

shapiro.test(df$AP_AHP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_AHP ~ MouseDevelopment, data=df)

kruskal.test(AP_AHP ~ MouseDevelopment, data=df)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

kruskal.test(MaturityIndex ~ MouseDevelopment, data=df)
dunnTest(MaturityIndex ~ MouseDevelopment, data=df, method="bh") 

###################################################################################################
###################### L4 FS analysis #######################################################


df = subset(data, BrainArea=="V1" & (CellType=="FS" | CellType=="PV") & CellLayer=="4")
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
dunnTest(Rheobase ~ MouseDevelopment, data=df, method="bh") 


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

summary(model)


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

kruskal.test(VoltageSag ~ MouseDevelopment, data=df)
dunnTest(VoltageSag ~ MouseDevelopment, data=df, method="bh") 


### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

kruskal.test(AP_Threshold ~ MouseDevelopment, data=df)

### AP Half-Width ###
boxplot(AP_HalfWidth ~ MouseDevelopment, data = df)
model = aov(AP_HalfWidth ~ MouseDevelopment, data = df)

shapiro.test(df$AP_HalfWidth)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_HalfWidth ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

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

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

kruskal.test(MaturityIndex ~ MouseDevelopment, data=df)
dunnTest(MaturityIndex ~ MouseDevelopment, data=df, method="bh") 

###################################################################################################
###################### V1 SST analysis #######################################################

df = subset(data, BrainArea=="V1" & CellType=="SST" )
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

### Tau ###
boxplot(Tau ~ MouseDevelopment, data = df)
model = aov(Tau ~ MouseDevelopment, data = df)

shapiro.test(df$Tau)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Tau ~ MouseDevelopment, data=df)

kruskal.test(Tau ~ MouseDevelopment, data=df)


### Rheobase ###
boxplot(Rheobase ~ MouseDevelopment, data = df)
model = aov(Rheobase ~ MouseDevelopment, data = df)

shapiro.test(df$Rheobase)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Rheobase ~ MouseDevelopment, data=df)

kruskal.test(Rheobase ~ MouseDevelopment, data=df)
dunnTest(Rheobase ~ MouseDevelopment, data=df, method="bh") 


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

summary(model)


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

summary(model)
 


### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

kruskal.test(AP_Threshold ~ MouseDevelopment, data=df)

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
TukeyHSD(model)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


###################################################################################################
data = read.table("C://Users//Butt Lab//Documents//GitHub//IntrinsicProfile//IntrinsicProfiles_Results_S1PyrOnly.csv", sep = ",", header = TRUE)

data$CellID = as.factor(data$CellID)
data$MouseDevelopment = factor(data$MouseDevelopment, levels = c("P5-P8", "P9-P13","P14-P18"))
data$MouseID = as.factor(data$MouseID)
data$BrainArea = as.factor(data$BrainArea)
data$CellTarget = as.factor(data$CellTarget)
data$CellType = as.factor(data$CellType)
data$CellLayer = as.factor(data$CellLayer)
summary(data)

###################### L4 Pyr analysis #######################################################
df = subset(data, BrainArea=="S1BF" & CellType=="Pyr" & CellLayer=="4")
summary(df)


### Rin ###
boxplot(Rin ~ MouseDevelopment, data = df)
model = aov(Rin ~ MouseDevelopment, data = df)

shapiro.test(df$Rin)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Rin ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


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
dunnTest(Rheobase ~ MouseDevelopment, data=df, method="bh") 


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

summary(model)

### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

summary(model)


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

summary(model)


### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

kruskal.test(AP_Threshold ~ MouseDevelopment, data=df)
dunnTest(AP_Threshold ~ MouseDevelopment, data=df, method="bh") 

### AP Half-Width ###
boxplot(AP_HalfWidth ~ MouseDevelopment, data = df)
model = aov(AP_HalfWidth ~ MouseDevelopment, data = df)

shapiro.test(df$AP_HalfWidth)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_HalfWidth ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### AP Height ###
boxplot(AP_Height ~ MouseDevelopment, data = df)
model = aov(AP_Height ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Height)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Height ~ MouseDevelopment, data=df)

summary(model)

### AP AHP ###
boxplot(AP_AHP ~ MouseDevelopment, data = df)
model = aov(AP_AHP ~ MouseDevelopment, data = df)

shapiro.test(df$AP_AHP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_AHP ~ MouseDevelopment, data=df)

summary(model)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

kruskal.test(MaturityIndex ~ MouseDevelopment, data=df)
dunnTest(MaturityIndex ~ MouseDevelopment, data=df, method="bh") 

###################### L4 FS analysis #######################################################
df = subset(data, BrainArea=="S1BF" & (CellType=="FS" |  CellType=="PV") )
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

summary(model)


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

summary(model)



### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

summary(model)


### AP Half-Width ###
boxplot(AP_HalfWidth ~ MouseDevelopment, data = df)
model = aov(AP_HalfWidth ~ MouseDevelopment, data = df)

shapiro.test(df$AP_HalfWidth)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_HalfWidth ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### AP Height ###
boxplot(AP_Height ~ MouseDevelopment, data = df)
model = aov(AP_Height ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Height)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Height ~ MouseDevelopment, data=df)

summary(model)

### AP AHP ###
boxplot(AP_AHP ~ MouseDevelopment, data = df)
model = aov(AP_AHP ~ MouseDevelopment, data = df)

shapiro.test(df$AP_AHP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_AHP ~ MouseDevelopment, data=df)

summary(model)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

kruskal.test(MaturityIndex ~ MouseDevelopment, data=df)
dunnTest(MaturityIndex ~ MouseDevelopment, data=df, method="bh") 

###################### L4 SST analysis #######################################################
df = subset(data, BrainArea=="S1BF" & CellType=="SST")
summary(df)


### Rin ###
boxplot(Rin ~ MouseDevelopment, data = df)
model = aov(Rin ~ MouseDevelopment, data = df)

shapiro.test(df$Rin)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(Rin ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### RMP ###
boxplot(RMP ~ MouseDevelopment, data = df)
model = aov(RMP ~ MouseDevelopment, data = df)

shapiro.test(df$RMP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(RMP ~ MouseDevelopment, data=df)

summary(model)

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

summary(model)


### Max Firing Frequency ###
boxplot(MaxFiringFrequency ~ MouseDevelopment, data = df)
model = aov(MaxFiringFrequency ~ MouseDevelopment, data = df)

shapiro.test(df$MaxFiringFrequency)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaxFiringFrequency ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

### Adaptation Index ###
boxplot(AdaptationIndex ~ MouseDevelopment, data = df)
model = aov(AdaptationIndex ~ MouseDevelopment, data = df)

shapiro.test(df$AdaptationIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AdaptationIndex ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### Voltage Sag ###
boxplot(VoltageSag ~ MouseDevelopment, data = df)
model = aov(VoltageSag ~ MouseDevelopment, data = df)

shapiro.test(df$VoltageSag)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(VoltageSag ~ MouseDevelopment, data=df)

kruskal.test(VoltageSag ~ MouseDevelopment, data=df)
dunnTest(VoltageSag ~ MouseDevelopment, data=df, method="bh") 

### AP Threshold ###
boxplot(AP_Threshold ~ MouseDevelopment, data = df)
model = aov(AP_Threshold ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Threshold)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Threshold ~ MouseDevelopment, data=df)

summary(model)


### AP Half-Width ###
boxplot(AP_HalfWidth ~ MouseDevelopment, data = df)
model = aov(AP_HalfWidth ~ MouseDevelopment, data = df)

shapiro.test(df$AP_HalfWidth)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_HalfWidth ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)


### AP Height ###
boxplot(AP_Height ~ MouseDevelopment, data = df)
model = aov(AP_Height ~ MouseDevelopment, data = df)

shapiro.test(df$AP_Height)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_Height ~ MouseDevelopment, data=df)

summary(model)

### AP AHP ###
boxplot(AP_AHP ~ MouseDevelopment, data = df)
model = aov(AP_AHP ~ MouseDevelopment, data = df)

shapiro.test(df$AP_AHP)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(AP_AHP ~ MouseDevelopment, data=df)

summary(model)
TukeyHSD(model)

### Maturity index###
boxplot(MaturityIndex ~ MouseDevelopment, data = df)
model = aov(MaturityIndex ~ MouseDevelopment, data = df)

shapiro.test(df$MaturityIndex)
qqnorm(model$residuals)
qqline(model$residuals)
leveneTest(MaturityIndex ~ MouseDevelopment, data=df)

kruskal.test(MaturityIndex ~ MouseDevelopment, data=df)
dunnTest(MaturityIndex ~ MouseDevelopment, data=df, method="bh") 
