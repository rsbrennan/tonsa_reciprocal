library(survminer)
library(survival)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)




#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Survival Data analysis ##################################################################################


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
                        conf.int = TRUE
                        
)

surv.plot



#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Egg Production and Hatching #############################################################################


Cost.epr <- read.table(file = "051018_HH_Feeding_cost_epr.txt", sep = "\t", header = TRUE)
Cost.epr$Temp <- c(22)
Cost.epr$pH <- c(7.5)
Cost.epr$Treatment <- c(4)

Cost.epr.trans <- read.table(file = "051018_AAHH_Feeding_cost_epr.txt", sep = "\t", header = TRUE)
Cost.epr.trans$Temp <- c(18)
Cost.epr.trans$pH <- c(8.2)
Cost.epr.trans$Treatment <- c(5)

Cost.epr.trans2 <- read.table(file = "050818_AA_Feeding_cost_epr.txt", sep = "\t", header = TRUE)
Cost.epr.trans2$Temp <- c(18)
Cost.epr.trans2$pH <- c(8.2)
Cost.epr.trans2$Treatment <- c(1)

Cost.epr.trans3 <- read.table(file = "050818_HHAA_Feeding_cost_epr.txt", sep = "\t", header = TRUE)
Cost.epr.trans3$Temp <- c(22)
Cost.epr.trans3$pH <- c(7.5)
Cost.epr.trans3$Treatment <- c(6)


Cost.epr <- rbind(Cost.epr, Cost.epr.trans, Cost.epr.trans2, Cost.epr.trans3)


# Pairwise comparison of EPR results
pairwise.t.test(Cost.epr$EPR, Cost.epr$Food:Cost.epr$Treatment, p.adj = "none")


## Plot EPR results

# Calculate means and 95% confidence intervals
detach("package:dplyr", unload = TRUE)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=TRUE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

Cost.epr.mean <- summarySE(Cost.epr, measurevar = "EPR", groupvars = c("Treatment", "Food"))


library(dplyr)
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


font.size <- 32

Cost.epr.point.plot <- ggplot(data = Cost.epr.mean.no.starve, aes(x=Food, 
                                                                  y=EPR, 
                                                                  color = factor(Treatment2), 
                                                                  shape = factor(Line), 
                                                                  group = paste(Food, Order)))+
  
  geom_errorbar(aes(ymin=EPR-ci, ymax=EPR+ci), colour = "black", size=2, width = 0.4, position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5), size = 13.5)+
  
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



#################################################################################################################################################################################
#################################################################################################################################################################################
####################################################################### Calculating Malthusian Parameter ########################################################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(popbio)




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
  mutate(lx_diag = if_else(lx_diag <= 1.000000, lx_diag, lag(lx_diag, n=1, default =last(lx_diag))))

# Survivorship is reflective of prior day. Therefore, there can be no 0 survivorship in this vector.
# If all animals die, then the vector ends and the matrix is truncated at the appropriate time
SurvData1 <- filter(SurvData1, lx_diag > 0) 


# Check if there is any super survivorship. There can be no survivorship > 1.
# No individuals can be lost and reappear
if (any(SurvData1$lx_diag > 1) == TRUE ){
  SurvData1$lx_diag <- if_else(SurvData1$lx_diag > 1, lag(SurvData1$lx_diag), SurvData1$lx_diag)
}

any(SurvData1$lx_diag > 1)



##### EPR data #####


Cost.epr$fecundity <- Cost.epr$EPR*Cost.epr$HF


Cost.epr$F.ratio <- case_when(Cost.epr$Treatment == 4 ~ 0.45,
                            Cost.epr$Treatment == 5 ~ 0.44,
                            Cost.epr$Treatment == 1 ~ 0.54,
                            Cost.epr$Treatment == 6 ~ 0.52)


# Caluclate per capita sex-specific fecundity

Cost.epr$sex.spec.fecundity <- Cost.epr$fecundity*Cost.epr$F.ratio


# Create a list of lists to use for matrix creations.


### NOTE: lists of survivorship and EPR data MUST have same names and same length of items

EPR.list <- split(Cost.epr, f = list(Cost.epr$Treatment,
                                   Cost.epr$Bin,
                                   Cost.epr$Food))

EPR.list <- lapply(EPR.list, function (x) split(x, f = list(x$Bin,
                                                            x$Food)))
names(EPR.list)
View(EPR.list)



Survival.list <- split(SurvData1, f = list(SurvData1$Treatment,
                                          SurvData1$Rep,
                                          SurvData1$Food))


names(Survival.list)
names.surv <- names(Survival.list)
names.surv <- as.list(names.surv)
names.surv



# Create a dummy table to add data to

lambda.results <- data.frame(Variables = "dummy", lambda = 0)


for (i in 1:length(Survival.list)) { #list of names for each gen and treatment
  
  diag.df <- Survival.list[[i]] # temporary list used to index at each point
  
  diag.v <- diag.df$lx_diag # extract the vector of interest for the diagonal
  
  leslie.matrix <- diag(diag.v) # create the matrix without the top row added
  
  # create a matrix of zeros if the randomly vector chosen is a vector of no survival
  if(dim(leslie.matrix) == 0) {
    leslie.matrix <- matrix(0, 10, 10)
  }
  
  zero <- matrix(0, nrow = 1, ncol = (dim(leslie.matrix)-1)) # make a zero vector to add to the top of the matrix
  
  
  
  
  epr.df <- EPR.list[[i]] # select the appropriate data frame to pull epr values from
  
  epr.vector <- epr.df$sex.spec.fecundity # select the vector with sex spec fecundity
  
  for(j in 1:length(epr.vector)) {
    
    epr.value <- epr.vector[j] # use the sex specific fecundity values in sequence with the same survival matrix
    
    if(is.na(epr.value) == TRUE || nrow(epr.df) == 0) { ## use two || for an "or" statement 
      epr.value <- 10 ## use an arrow when assigning numbers, not a "=="
      
    }
    
    fecundity.row <- t(c(zero, epr.value)) # combine the zero row and the epr.value for the matrix and transpose it to make it a row
    
    
    leslie.matrix1 <- rbind(fecundity.row, leslie.matrix) # add the fecundity row to the matrix
    
    matrix.final <- leslie.matrix1[-nrow(leslie.matrix1),] # eliminate the last row to make it a square
    
    eigen.calcs <- eigen.analysis(matrix.final, zero = FALSE) # calculate the eigen values
    
    lambda.value <- eigen.calcs$lambda1 # extract the dominant eigen values
    
    lambda.row <- data.frame(Variables = names.surv[i], lambda = lambda.value) # create a 1x2 data frame to add to the end of the final data frame
    # data frame has to have the same colnames in order to rbind
    
    colnames(lambda.row) <- colnames(lambda.results) # make sure the data frames have the same names
    
    lambda.results <- bind_rows(lambda.results, lambda.row) # append the data frame with new results
    
    
  }
  
  
  
}


lambda.results <- separate(lambda.results, "Variables", into = c("Treatment", "Rep", "Food"))


# remove the first row as a last step
lambda.results <- lambda.results[-1,] 


# remove all data where lambda < 0 (i.e. no data present)
lambda.results <- filter(lambda.results, lambda > 0)

# calculate malthusian parameter
lambda.results$Malthusian <- log10(lambda.results$lambda)

# Remove results of starved individuals
lambda.results.no.starve <- filter(lambda.results, lambda.results$Food > 0)



lambda.results.no.starve <- unite(lambda.results.no.starve,#the data frame
                                  unique, #the name of the new column
                                  c(Treatment, Food), #the existing columns you want to combine
                                  remove = FALSE)


# Pairwise comparison of malthusian parameters
malthusian.pairwise <- pairwise.t.test(lambda.results.no.starve$Malthusian, lambda.results.no.starve$unique)


# Two-way ANOVA
malthusian.model <- aov(Malthusian~Treatment*Food, data=lambda.results.no.starve)




###### Malthusian boxplot #####

lambda.results.no.starve$Line <- case_when(lambda.results.no.starve$Treatment == 1 ~ "AM",
                                           lambda.results.no.starve$Treatment == 4 ~ "GH",
                                           lambda.results.no.starve$Treatment == 5 ~ "GH",
                                           lambda.results.no.starve$Treatment == 6 ~ "AM")

lambda.results.no.starve$Treatment <- case_when(lambda.results.no.starve$Treatment == 1 ~ "AM",
                                                lambda.results.no.starve$Treatment == 4 ~ "GH",
                                                lambda.results.no.starve$Treatment == 5 ~ "AM",
                                                lambda.results.no.starve$Treatment == 6 ~ "GH")



## Create a specific column to order everything by
lambda.results.no.starve$Order <- case_when(lambda.results.no.starve$Culture == 1 ~ "A",
                                            lambda.results.no.starve$Culture == 4 ~ "C",
                                            lambda.results.no.starve$Culture == 5 ~ "D",
                                            lambda.results.no.starve$Culture == 6 ~ "B")



Malthusian.boxplot <- ggplot(data = lambda.results.no.starve, aes(factor(Food), Malthusian, 
                                                                  
                                                                  group = paste(Food, Order) # use this method to get everything in order based on the group
                                                                  ## help from: https://stackoverflow.com/a/14529172/9335733
                                                                  
))+
  geom_boxplot(lwd=1.15, aes(color=Treatment, fill = Line))




Malthusian.boxplot+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("cornflowerblue", "brown2"))+
  scale_fill_manual(values = alpha(c("cornflowerblue", "brown2"), 0.6))+ # make the fill of the boxes slightly transparent
  scale_shape_manual(values = c(21,21,21,21))+
  theme_classic()+
  labs(y="Malthusian Parameter", x="Food")+
  scale_x_discrete(limits = c("800", "250"))+ # add quotes to make the x axis order in the order you desire. Also use "limits" instead of "breaks"
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.text = element_text(size = font.size), 
        axis.title.x = element_text(size = font.size), 
        axis.text.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.y = element_text(size = font.size),
        legend.position = "bottom",
        legend.title = element_text(size = font.size))





