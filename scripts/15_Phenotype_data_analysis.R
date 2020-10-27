library(survminer)
library(survival)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(popbio)
library(car)
library(lme4)
library(sjPlot)
library(emmeans)
library(agricolae)








#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Survival Data analysis ##################################################################################
setwd("C:/Users/james/Documents/Grad_school/OA_Project/Cost_experiments/")


# Read in data table
SurvData <- fread("SurvDataTransplant_complete_sexratio.txt")

# Create uniqe sorting column to group each technical and biological replicate

SurvData <- unite(SurvData,#the data frame
                  unique, #the name of the new column
                  c(Treatment, Rep, Beak), #the existing columns you want to combine
                  remove = FALSE) 

# change number of individuals to numeric format
SurvData$nx <- as.numeric(SurvData$nx)


# Modify dataframe to repeat number of rows by the number of remaining individuals
# Other vectors now weighted to reflect the number of remaining individuals for that specific time point
SurvData.kaplan <- SurvData[rep(seq(nrow(SurvData)), 
                                       SurvData$nx),#the column which holds the value you wish to have everything repeated by
                                   ]
# Designate "Survival" for the individuals that survived to adulthood
# help from: https://stackoverflow.com/questions/51642416/changing-number-of-rows-based-on-changing-vector-in-r

SurvData.kaplan <- SurvData.kaplan %>%
  group_by(Treatment, Rep, Beak) %>%
  mutate(
    survival = case_when(
      Food == 0 ~ 1,# no starved treatments survived
      nx == min(nx) ~ 0, #need to use the "==" and VALUE OF 1 EQUALS DEATH HAPPENING!!!!!!
      nx > min(nx) ~ 1 )) %>% 
  group_by(unique) %>%
  mutate(keep = max(nx) - nx) %>% 
  group_by(unique, keep) %>% #need to group by unique and keep to get all the animals that died in between the start and the lowest amount of dead animals
  filter(survival == 1 & nx == min(nx) |
           survival == 0) %>% 
  filter(row_number() %in% 1:unique(keep) |
           survival == 0) %>% 
  select(-keep) %>% 
  ungroup()


# Remove starting point since they all started alive
SurvData.kaplan <- SurvData.kaplan[!SurvData.kaplan$time == 0,] 

# Create a category for unique treatments and food
SurvData.kaplan <- unite(SurvData.kaplan,#the data frame
                         Treat.food, #the name of the new column
                         c(Treatment, Food), #the existing columns you want to combine
                         remove = FALSE) 

## Plot only food limited and food replete treatments

# Remove starved individuals to plot 
SurvData.kaplan.no.starve <- SurvData.kaplan[!SurvData.kaplan$Food == 0,]

# New grouping vectors for line and developmental treatment
SurvData.kaplan.no.starve$line <- if_else(SurvData.kaplan.no.starve$Treatment == 1 | SurvData.kaplan.no.starve$Treatment == 6, 
                                          "AA", "HH")

SurvData.kaplan.no.starve$Treatment2 <- if_else(SurvData.kaplan.no.starve$Treatment == 1 | SurvData.kaplan.no.starve$Treatment == 5, 
                                                "AA", "HH")

# Create a survival object
surv_object <- Surv(time = SurvData.kaplan.no.starve$time, event = SurvData.kaplan.no.starve$survival)

# Create a model to fit the object
surv_object_fit <- survfit(surv_object~Treatment+Food, data = SurvData.kaplan.no.starve)

surv_object

surv_object_fit


# log-rank test results for resulting model
surv_pvalue(surv_object_fit, SurvData.kaplan.no.starve)

# Two-way ANOVA comparing day-specific survivorship by Food abundance and Line
surv.anova <- aov(lx~Food*line, data = SurvData.kaplan.no.starve)
summary(surv.anova)



# Plot the fits
surv.plot <- ggsurvplot(surv_object_fit,
                        data = SurvData.kaplan.no.starve,
                        facet.by = "Food",
                        #colors are the "lines"
                        palette = c("cornflowerblue", #AAAA
                                    "brown2", #HHHH
                                    "cornflowerblue", #AAHH
                                    "brown2" #HHAA
                        ),
                        #linetypes are the "treatments"
                        linetype =c("solid", #AAAA 250
                                    "solid", #AAAA 800
                                    "dashed", #HHHH 250
                                    "dashed", #HHHH 800
                                    "dashed", #AAHH 250
                                    "dashed", #AAHH 800
                                    "solid", #HHAA 250
                                    "solid" #HHAA 800
                        ),
                        conf.int = TRUE,
                        break.time.by = 10
                        
)

surv.plot.large <- ggpar(surv.plot, font.x = 48, font.y = 48, font.tickslab = 48)

surv.plot.large

aspect.ratio <- 1.15
ggsave(filename = "survplot_tall2.pdf", surv.plot.large, height = 8.5*aspect.ratio, width = 8.5, units = "in", device = "pdf")



#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Egg Production and Hatching #############################################################################



cost.dir <- "C:/Users/james/Documents/Grad_school/OA_Project/Cost_experiments/"
Cost.epr <- fread(paste(cost.dir, "Cost_EPR_total.txt", sep = ""))




# Pairwise comparison of EPR results
Cost.epr$Food.f <- as.factor(as.numeric(Cost.epr$Food))
Cost.epr$Treatment.f <- as.factor(as.numeric(Cost.epr$Treatment))
#pairwise.t.test(Cost.epr$EPR, Cost.epr$Food.f:Cost.epr$Treatment.f, p.adj = "none")

# test if Bins are significant
#m <- lm(EPR~Treatment*Bin*Food, Cost.epr)
#m2 <- Anova(m)
#m2

m <- lmer(EPR~factor(Treatment)*factor(Food)+(1|Bin), 
          # REML = FALSE,
          # family = gaussian,
          data = Cost.epr)
tab_model(m)


m2 <- Anova(m)
m2$factors <- rownames(m2)

fwrite(m2, "Cost_EPR_rand_effects_anova.txt", sep = "\t")

Cost.epr$Line <- case_when(Cost.epr$Treatment == 1 ~ "AM",
                                Cost.epr$Treatment == 4 ~ "GH",
                                Cost.epr$Treatment == 5 ~ "GH",
                                Cost.epr$Treatment == 6 ~ "AM")

Cost.epr$Environment <- case_when(Cost.epr$Treatment == 1 ~ "AM",
                                      Cost.epr$Treatment == 4 ~ "GH",
                                      Cost.epr$Treatment == 5 ~ "AM",
                                      Cost.epr$Treatment == 6 ~ "GH")



m3 <- lmer(EPR~factor(Line)*factor(Environment)*factor(Food)+(1|Bin), data = Cost.epr)

tab_model(m3)

m4 <- Anova(m3)
m4
m4$factors <- rownames(m4)

fwrite(m4, "Cost_EPR_3way_rand_effects_anova.txt", sep = "\t")



m.emmeans <- emmeans(m, pairwise ~ Treatment | Food)

EPR.pw <- summary(m.emmeans, adjust = "none")
EPR.contrasts <- EPR.pw$contrasts
fwrite(EPR.contrasts, file = "EPR_pw_random_effects.txt", sep = "\t")


Cost.epr.low.food <- Cost.epr %>% 
  filter(Food == 250)

m5 <- lmer(EPR~Line*Environment+(1|Bin), data = Cost.epr.low.food)

tab_model(m5)
m6 <- Anova(m5)
m6$factors <- rownames(m6)

fwrite(m6, "EPR_low_food_GxE.txt", sep = "\t")


plot_model(m5, 'int')

Cost.epr.high.food <- Cost.epr %>% 
  filter(Food == 800)


m7 <- lmer(EPR~Line*Environment+(1|Bin), data = Cost.epr.high.food)
tab_model(m7)

m8 <- Anova(m7)
m8$factors <- rownames(m8)


fwrite(m8, "EPR_high_food_GxE.txt", sep = "\t")



plot_model(m7, 'int')

## Plot EPR results


Cost.epr.mean <- Cost.epr %>% 
  group_by(Treatment, Food) %>% 
  summarise(mean = mean(EPR, na.rm = TRUE),
            sd = sd(EPR, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)



Cost.epr.mean$Line <- case_when(Cost.epr.mean$Treatment == 1 ~ "AM",
                                Cost.epr.mean$Treatment == 4 ~ "GH",
                                Cost.epr.mean$Treatment == 5 ~ "GH",
                                Cost.epr.mean$Treatment == 6 ~ "AM")

Cost.epr.mean$Treatment2 <- case_when(Cost.epr.mean$Treatment == 1 ~ "AM",
                                      Cost.epr.mean$Treatment == 4 ~ "GH",
                                      Cost.epr.mean$Treatment == 5 ~ "AM",
                                      Cost.epr.mean$Treatment == 6 ~ "GH")


Cost.epr.mean$Food <-as.numeric(Cost.epr.mean$Food)

# Plot only food limited and food replete conditions
Cost.epr.mean.no.starve <- filter(Cost.epr.mean, Food > 1)


Cost.epr.mean.no.starve$Food <- as.factor(Cost.epr.mean.no.starve$Food)


# Create variables for grouping when plotting
Cost.epr.mean.no.starve$Order <- case_when(Cost.epr.mean.no.starve$Treatment == 1 ~ "A",
                                           Cost.epr.mean.no.starve$Treatment == 4 ~ "C",
                                           Cost.epr.mean.no.starve$Treatment == 5 ~ "D",
                                           Cost.epr.mean.no.starve$Treatment == 6 ~ "B")


font.size <- 48

Cost.epr.point.plot <- ggplot(data = Cost.epr.mean.no.starve, aes(x=Food, 
                                                                  y=mean, 
                                                                  color = factor(Treatment2), 
                                                                  shape = factor(Line), 
                                                                  group = paste(Food, Order)))+
  
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), colour = "black", size=2, width = 0.4, position = position_dodge(width = 0.68))+
  
  geom_point(position = position_dodge(width = 0.68), size = 13.5)+
  
  scale_shape_manual(name = "Line",
                     values = c(16, 17))+
  
  scale_color_manual(name = "Treatment",
                     values = c("cornflowerblue", "brown2"))+
  
  scale_x_discrete(limits = c("800", "250"))+
  
  theme_classic()+
  theme(legend.text.align = 0,
        legend.text = element_text(size = font.size),
        legend.title = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(size = font.size, colour = "black"),
        axis.text.y = element_text(size = font.size, colour = "black"),
        legend.position = "bottom")+
  xlab(expression ("Food Concentration ("*mu*"g C/L)"))+
  ylab(expression (atop("Egg Production Rate",(day^-1~female^-1))))

Cost.epr.point.plot

ggsave(filename = "Cost.epr.point.plot3.pdf", Cost.epr.point.plot, height = 8.5*aspect.ratio, width = 8.5, units = "in", device = "pdf")




#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Calculating Malthusian Parameter ########################################################################


SurvData$lx <- as.numeric(SurvData$lx)

## Extract the sex ratio data for scaling fecundity
SurvData$M.Ratio <- as.numeric(SurvData$M.Ratio)
SurvData$F.Ratio <- as.numeric(SurvData$F.Ratio)

SurvDataSex <- SurvData %>%
  filter(F.Ratio >= 0 | M.Ratio >= 0) %>%
  group_by(Treatment, Rep, Beak) %>%
  summarise(F.Ratio = last(F.Ratio))

SurvDataSex.Mean <- SurvDataSex %>%
  group_by(Treatment) %>%
  summarise(F.Ratio = mean(F.Ratio, na.rm = TRUE)) %>%
  as.data.frame(SurvDataSex.Mean)

SurvDataSex.Mean




# Create a vector formatted for day-specific survivorship as it fits along the off-diagonal of a Leslie Matrix
# See Caswell, H. 2001. Matrix Population Models for further details


SurvData1 <- SurvData %>%
  group_by(Treatment) %>%
  mutate(lx_diag = if_else(time == 0, 1, as.numeric(lx/lag(lx, default = first(lx))))) %>% ## always start with 100%
  mutate(lx_diag = if_else(lx_diag <= 1.000000, lx_diag, lag(lx_diag, n=1, default =last(lx_diag)))) %>%
  mutate(days = if_else(time == 0, 1, as.numeric(time-lag(time)))) # create a new column for the number of days spent at the respective survivorships

# Survivorship is reflective of prior day. Therefore, there can be no 0 survivorship in this vector.
# If all animals die, then the vector ends and the matrix is truncated at the appropriate time
SurvData1 <- filter(SurvData1, lx_diag > 0) 


# Check if there is any super survivorship. There can be no survivorship > 1.
# No individuals can be lost and reappear


any(SurvData1$lx_diag > 1)



SurvData1 <- as.data.frame(SurvData1)


# elongate the data frame to make it reflect actual days spent over the experiment. This essentially changes the matrix to a daily matrix.

SurvData1 <- SurvData1[rep(seq(nrow(SurvData1)), SurvData1$days),] 

# remove starved samples
SurvData1 <- SurvData1[!SurvData1$Food == 0,]





#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Development Time ########################################################################################


Dev.time <- filter(SurvData, SurvData$Cdev > 0)
Dev.time <- Dev.time[rep(seq(nrow(Dev.time)), Dev.time$Cdev),]

Dev.time.sum <- Dev.time %>%
  group_by(Treatment, Rep, Food) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


Dev.time.sum <- Dev.time.sum[!Dev.time.sum$Food == 0,]

Dev.time.sum <- unite(Dev.time.sum,#the data frame
                      unique, #the name of the new column
                      c(Treatment, Rep, Food), #the existing columns you want to combine
                      remove = FALSE)
Dev.time.sum.short <- Dev.time.sum[,-c(6:10)]

## Add the dev.time to the survival table
SurvData1 <- inner_join(SurvData1,
                       Dev.time.sum.short,
                       by = c("Treatment", "Rep", "Food"))

SurvData1 <- SurvData1 %>% rename(dev.time = mean) # the new name of the column has to come first

SurvData1 <- unite(SurvData1,#the data frame
                   unique2, #the name of the new column
                   c(Treatment, Food), #the existing columns you want to combine
                   remove = FALSE)

Survival.list <- split(SurvData1, f = SurvData1$unique2) # don't create the list with a list of organizers. Then you create a complete matrix of samples which include some tables with no data



# lists within lists
Survival.list <- lapply(Survival.list, function (x) split(x, f = x$Rep))

View(Survival.list)

##### EPR data #####


Cost.epr$fecundity <- Cost.epr$EPR*Cost.epr$HF

# These sex ratio values are based on previously calculated values in survival data

Cost.epr$F.ratio <- case_when(Cost.epr$Treatment == 4 ~ 0.45,
                            Cost.epr$Treatment == 5 ~ 0.44,
                            Cost.epr$Treatment == 1 ~ 0.54,
                            Cost.epr$Treatment == 6 ~ 0.52)


# Caluclate per capita sex-specific fecundity

Cost.epr$sex.spec.fecundity <- Cost.epr$fecundity*Cost.epr$F.ratio

# remove starved animals
Cost.epr <- Cost.epr[!Cost.epr$Food == 0,]

Cost.epr <- unite(Cost.epr,#the data frame
                      unique, #the name of the new column
                      c(Treatment, Food), #the existing columns you want to combine
                      remove = FALSE)


### NOTE: lists of survivorship and EPR data MUST have same names and same length of items

EPR.list <- split(Cost.epr, f = Cost.epr$unique)

View(EPR.list)
EPR.list$`1_250`

names(Survival.list)
names(EPR.list)


names.surv <- names(Survival.list)
names.surv <- as.list(names.surv)
names.surv

# Create a dummy table to add data to

lambda.results <- data.frame(Variables = "dummy", Replicate = 0, lambda = 0)
#i=1
#j=1



for (i in 1:length(Survival.list)) { #list of names for each gen and treatment
  
  diag.ls <- Survival.list[[i]]
  
  epr.df <- EPR.list[[i]] # select the appropriate data frame to pull epr values from
  
  for (j in 1: length(diag.ls)) {
    
    diag.df <-  diag.ls[[j]] # temporary list used to index at each point
    
    diag.v <- diag.df$lx_diag # extract the vector of interest for the diagonal
    
    leslie.matrix <- diag(diag.v) # create the matrix without the top row added
    
    surv.rep <- mean(diag.df$Rep)
    
    dev.time.value <- as.integer(mean(diag.df$dev.time)) # calculate the development time
    
    zero <- matrix(0, nrow = 1, ncol = (dev.time.value)) # create an empty, one-row matrix to add the epr values to. This will be added to the leslie matrix as the top row.
    
    epr.vector <- epr.df$sex.spec.fecundity # select the vector with sex spec fecundity
    
    epr.count <- as.integer(dim(leslie.matrix)[1]-dev.time.value) # calculate the number of days to incorporate egg production based on development time
    
    if (epr.count < 1) {
      
      epr.count <- 1
      
    } 
    
    
    for(k in 1:length(epr.vector)) {
      
      epr.value <- epr.vector[k] # use the sex specific fecundity values in sequence with the same survival matrix
      
      if (is.na(epr.value) == TRUE) {
        
        epr.value <- 0 ## use an arrow when assigning numbers, not a "=="
        
      } 
      
      
      
      fecundity.row <- t(c(zero, rep(epr.value, epr.count))) # combine the zero row and the epr.value for the matrix and transpose it to make it a row
      
      if (ncol(fecundity.row) > ncol(leslie.matrix)) { # if the dev.time is somehow greater than the survival matrix, then only make the last column reflective of epr
        
        delete <- ncol(fecundity.row)-ncol(leslie.matrix) # find the number of days that the dev time is greater than the survivorship
        
        fecundity.row <- fecundity.row[,-c(1:delete)] # delete those days from the fecundity row to make it the same size as the leslie matrix
        
      }
      
      leslie.matrix1 <- rbind(fecundity.row, leslie.matrix) # add the fecundity row to the matrix
      
      matrix.final <- leslie.matrix1[-nrow(leslie.matrix1),] # eliminate the last row to make it a square
      
      eigen.calcs <- eigen.analysis(matrix.final, zero = FALSE) # calculate the eigen values
      
      lambda.value <- eigen.calcs$lambda1 # extract the dominant eigen values
      
      lambda.row <- data.frame(Variables = names.surv[i], Replicate = surv.rep, lambda = lambda.value) # create a 1x2 data frame to add to the end of the final data frame
      # data frame has to have the same colnames in order to rbind
      
      colnames(lambda.row) <- colnames(lambda.results) # make sure the data frames have the same names
      
      lambda.results <- bind_rows(lambda.results, lambda.row) # append the data frame with new results
      
      
    }
    
    
    
  }
}

######################################################################################################################################################


lambda.results <- separate(lambda.results, "Variables", into = c("Treatment","Replicate", "Food"))


# remove the first row as a last step
lambda.results <- lambda.results[-1,] 


View(lambda.results)
# calculate malthusian parameter
lambda.results$Malthusian <- log10(lambda.results$lambda)


# create continuous food vector for anova and plotting
lambda.results$Food.c <- as.numeric(as.character(lambda.results$Food)) 

lambda.results$Food <- as.factor(as.numeric(lambda.results$Food))
lambda.results$Treatment <- as.numeric(lambda.results$Treatment)

# evaluate line x environment interactions
lambda.results$Line <- case_when(lambda.results$Treatment == 1 ~ "AM",
                                 lambda.results$Treatment == 4 ~ "GH",
                                 lambda.results$Treatment == 5 ~ "GH",
                                 lambda.results$Treatment == 6 ~ "AM")


lambda.results$Environment <- case_when(lambda.results$Treatment == 1 ~ "AM",
                                 lambda.results$Treatment == 4 ~ "GH",
                                 lambda.results$Treatment == 5 ~ "AM",
                                 lambda.results$Treatment == 6 ~ "GH")









##### Lambda statistics #####



#setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/Cost_experiment/")
lambda.results <- fread(paste(cost.dir, "lambda_results_cost_devtime_lines_sep.txt", sep = ""))


l <- lmer(lambda~factor(Treatment)*factor(Food)+(1|Replicate), 
           # REML = FALSE,
           # family = gaussian,
            data = lambda.results)

summary(l)
tab_model(l)


lambda.emm <- emmeans(l, pairwise ~ Treatment | Food)
lambda.tukey.groups <- CLD(lambda.emm)

pairs(lambda.emm)

summary(lambda.emm, adjust = "none")
lambda.pw <- summary(lambda.emm, adjust = "none")
lambda.contrasts <- lambda.pw$contrasts
fwrite(lambda.contrasts, file = "lambda_pw_random_effects.txt", sep = "\t")


fwrite(lambda.tukey.groups, file = "Lambda_cost_devtime_lines_tukey_groups.txt", sep = "\t")


## find grouping letters for lambda values
lambda.results <- unite(lambda.results,
                        unique,
                        c(Treatment, Food),
                        remove = F)

l1 <- lm(lambda~unique, lambda.results)


test <- LSD.test(l1,trt = "unique")
lambda.lsd.groups <- test$groups
lambda.lsd.groups$unique <- rownames(lambda.lsd.groups)

fwrite(lambda.lsd.groups, file = "lambda_pairwise_devtime_groups.txt", sep = "\t")



## two-way ANOVA for Line x Environment interactions without separating by food treatment
l2 <- lmer(lambda~factor(Line)*factor(Environment) +(1|Replicate), data = lambda.results)
summary(l2)
tab_model(l2)
Anova(l2)


## three-way ANOVA for Line x Environment x Food interaction
l3 <- lmer(lambda~factor(Line)*factor(Environment)*factor(Food.c) +(1|Replicate), data = lambda.results)
tab_model(l3)
Anova(l3)
lambda.cost.3way.anova <- Anova(l3)
lambda.cost.3way.anova$factors <- rownames(lambda.cost.3way.anova)
fwrite(lambda.cost.3way.anova, file = "lambda_cost_3way_anova.txt", sep = "\t")



## analyzing different food groups individually for Line x Environment interactions

# Low food
lambda.results.low.food <- lambda.results %>% 
  filter(Food == 250)


l4 <- lmer(lambda~Line*Environment+(1|Replicate), data = lambda.results.low.food)

l5 <- Anova(l4)
l5$factors <- rownames(l5)



plot_model(l4, 'int')

fwrite(l5, "lambda_low_food_GxE.txt", sep = "\t")


# High food
lambda.results.high.food <- lambda.results %>% 
  filter(Food == 800)


l6 <- lmer(lambda~Line*Environment+(1|Replicate), data = lambda.results.high.food)

l7 <- Anova(l6)
l7$factors <- rownames(l7)

plot_model(l6, 'int')


fwrite(l7, "lambda_high_food_GxE.txt", sep = "\t")


###### Lambda Point Plot #####

lambda.list.mean <- lambda.results %>%
  group_by(Treatment, Food) %>%
  summarise(mean = mean(lambda, na.rm = TRUE),
            sd = sd(lambda, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

fwrite(lambda.list.mean, file = "lambda_mean_devtime.txt", sep = "\t")

lambda.list.mean$Line <- case_when(lambda.list.mean$Treatment == 1 ~ "AM",
                                           lambda.list.mean$Treatment == 4 ~ "GH",
                                           lambda.list.mean$Treatment == 5 ~ "GH",
                                           lambda.list.mean$Treatment == 6 ~ "AM")

lambda.list.mean$Treatment2 <- case_when(lambda.list.mean$Treatment == 1 ~ "AM",
                                                lambda.list.mean$Treatment == 4 ~ "GH",
                                                lambda.list.mean$Treatment == 5 ~ "AM",
                                                lambda.list.mean$Treatment == 6 ~ "GH")


lambda.list.mean$Order <- case_when(lambda.list.mean$Treatment == 1 ~ "A",
                                            lambda.list.mean$Treatment == 4 ~ "C",
                                            lambda.list.mean$Treatment == 5 ~ "D",
                                            lambda.list.mean$Treatment == 6 ~ "B")

font.size <- 48

lambda.list.mean$Food<- as.factor(lambda.list.mean$Food)

Lambda.point.plot <- ggplot(data = lambda.list.mean, aes(x=Food, 
                                                                  y=mean, 
                                                                  color = factor(Treatment2), 
                                                                  shape = factor(Line), 
                                                                  group = paste(Food, Order)))+
  
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), colour = "black", size=2, width = 0.4, position = position_dodge(width = 0.68))+
  geom_point(position = position_dodge(width = 0.68), size = 13.5)+
  
  scale_shape_manual(name = "Line",
                     values = c(16, 17))+
  
  scale_color_manual(name = "Treatment",
                     values = c("cornflowerblue", "brown2"))+
  
  scale_x_discrete(limits = c("800", "250"))+
  
  theme_classic()+
  theme(legend.text.align = 0,
        legend.text = element_text(size = font.size),
        legend.title = element_text(size = font.size),
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(size = font.size, colour = "black"),
        axis.text.y = element_text(size = font.size, colour = "black"),
        legend.position = "bottom")+
  xlab(expression ("Food Concentration ("*mu*"g C/L)"))+
  ylab(expression (atop("lambda",(generation^-1))))



Lambda.point.plot
setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/Cost_experiment/")
ggsave(filename = "Malthusian.point.plot.tall_new_analysis_devtime.pdf", Malthusian.point.plot, height = 8.5*aspect.ratio, width = 8.5, units = "in", device = "pdf")

###### Old stats methods #####





malthusian.model <- lm(lambda~Food.c*Treatment, data=lambda.results)


# use a continuous generation model for the anova
m3 <- Anova(malthusian.model)
m3

#setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/Cost_experiment/")
fwrite(m3, file = "lambda_cost_Anova_devtime.txt", sep = "\t")

shapiro.test(malthusian.model$residuals)

hist(malthusian.model$residuals)

leveneTest(malthusian.model)

foo <- plot_model(malthusian.model, 'int')

foo

# Pairwise comparison of malthusian parameters

# create the model to use for pairwise comparisons -- generation must be a factor, not a continuous variable
malthusian.model.2 <- lm(lambda~Food*Treatment, data=lambda.results)


emmip(malthusian.model.2, Treatment ~ Food)

cost.emm_1 <- emmeans(malthusian.model.2, pairwise ~ Treatment | Food)

# Tukey comparison
p <- pairs(cost.emm_1)
p


lambda.results <- unite(lambda.results,#the data frame
                        unique, #the name of the new column
                        c(Treatment, Food), #the existing columns you want to combine
                        remove = FALSE)

# non-adjusted pairwise t-test
p2 <- tidy(pairwise.t.test(lambda.results$lambda, lambda.results$unique, p.adjust.method = "none"))


fwrite(lambda.results, file = "lambda_results_cost_devtime_lines.txt", sep = "\t")
fwrite(p2, file = "Lambda_cost_pairwise_devtime.txt", sep = "\t")

# model for writing the tukey results
malthusian.model.3 <- aov(lambda~unique, data = lambda.results)
p3 <- TukeyHSD(malthusian.model.3, "unique")
p3.tukey <- tidy(p3)
fwrite(p3.tukey, file = "lambda_cost_devtime_Tukey.txt", sep = "\t")


