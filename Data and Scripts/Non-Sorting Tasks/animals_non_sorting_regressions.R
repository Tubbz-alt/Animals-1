# Load packages
library(rstan);
require(dplyr);
require(lme4);
require(ggplot2);
require(tidyr);
require(MASS);
library(ordinal);
library(car);
library(tables);
library(vcd);

# Set working directory
setwd("~/Documents/Animals Paper/Data/Final Analyses")

######## Familiarity Ratings ########
# Ordered Logistic Regression (cumulative link)
# Fixed: group
# Random: subject, item
# (Using CLMM + Ordinal package because polr doesn't allow random effects) 

d.fam <- read.csv('animals_familiarity.csv', header=TRUE, sep=",")
d.fam$rating<- factor(d.fam$rating, ordered=TRUE, levels=c("1","2","3","4"))
m.fam <- clmm(rating ~ group + (1|subject) + (1|item), data=d.fam)

# The proportional odds assumption is violated with our data...
# But using linear model gives basically same result:

xtabs(~ group + rating, d.fam)
d.fam$rating<- as.numeric(d.fam$rating)
m.fam2 <- lmer(rating ~ group + (1|subject) + (1|item), data=d.fam)
Anova(m.fam2)

######## Texture Forced Choice ########

# Logistic Regression
# Fixed: group 
# Random: subject, item 
# Maximal s

d.texture <- read.csv('animals_texture.csv',header=TRUE,sep=",")
m.texture <- glmer(expected ~ group + (1|subject) + (1|animal), data=d.texture, family=binomial(link="logit"))
Anova(m.texture)

# convert to odds ratio
exp(fixef(m.texture))

# CI 
se<-sqrt(diag(vcov(m.texture)))
tab<-cbind(Est=fixef(m.texture),LL=fixef(m.texture)-1.96*se, UL=fixef(m.texture)+1.96*se)
exp(tab)

######## Shape Odd-one-out ########

# Logistic Regression
# Fixed: group 
# Random: subject, item 

d.shape <- read.csv('animals_shape.csv',header=TRUE,sep=",")
d.shape$item <- factor(d.shape$item)
m.shape <- glmer(expected ~ group + (1|subject) + (1|item), data=d.shape, family=binomial(link="logit"))
Anova(m.shape)

# convert to odds ratio
exp(fixef(m.shape))

# CI 
se<-sqrt(diag(vcov(m.shape)))
tab<-cbind(Est=fixef(m.shape),LL=fixef(m.shape)-1.96*se, UL=fixef(m.shape)+1.96*se)
exp(tab)
