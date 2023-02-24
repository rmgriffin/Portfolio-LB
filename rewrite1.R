rm(list=ls())

# Setting up replicable environment
install.packages("renv")
library(renv)
renv::init(bare = TRUE)

# ######################################
# This code is based on code provided by Geret DePiper: "portfolio_efficiency_code_JIN.R"
# to implement the decay factor
# and Howard Townsend: "EfficientFrontier2.R"
# ######################################

PKG <- c("tidyverse","readxl","priceR","RColorBrewer","viridis","openxlsx","grDevices","extrafont","boot",
         "kernlab","reshape","matrixcalc")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}
rm(p,PKG)

renv::snapshot() # Save package snapshot

# Options
theme_set(theme_bw() + 
            theme(text = element_text(size = 20),
                  axis.title = element_text(size = 20),
                  legend.position="top",
                  panel.background = element_rect(fill='transparent'), #transparent panel bg
                  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                  # panel.grid.major = element_blank(), #remove major gridlines
                  # panel.grid.minor = element_blank(), #remove minor gridlines
                  legend.background = element_rect(fill='transparent'), #transparent legend bg
                  legend.box.background = element_rect(fill='transparent', color='NA') #transparent legend panel
            ))

options(scipen=999) # Reduce the use of scientific notation

# Load data
#setwd("~/Portfolio project/Lauran Frontier Analysis")
raw<-read_excel("./Data/foss_landings_NMFSRegion_Year_Species_State_UpTo2021Data_StandardisedPrice.xlsx") 
raw<-raw %>% mutate_if(is.character,as.factor)
# Standardized price seems inflation adjusted before import
#setwd("~/Portfolio project/Howard Frontier Analysis/Howard Efficient Frontier2")

# Data analysis
Pub.Man<-raw %>%
  #filter(!(Management.Group=="Other") & Confidentiality=="Public") %>%
  filter(Confidentiality=="Public") %>%
  arrange(Management.Group, Scientific.Name, Year) %>%
  filter(!grepl("SEAWEED", NMFS.Name))

YrSps <- Pub.Man %>%
  select(-c("State", "Source")) %>%
  aggregate(cbind(Standardized.Price, Dollars, Metric.Tons, Pounds) ~ #R Wouldn't want to aggregate Dollars as it is not inflation adjusted
              Year+NMFS.Name+Confidentiality+Collection+Scientific.Name+Tsn+Management.Group+NMFS.Namev2, sum)

YrSps <- YrSps %>%
  group_by(Year) %>%
  drop_na(Metric.Tons) %>%
  mutate(PercentAnnualCatch = Metric.Tons / sum(Metric.Tons)*100) %>% 
  dplyr::rename("IYear"="Year",   #rename to match Geret's code
                "Value"="Standardized.Price",
                "Catch"="Metric.Tons",
                "Taxonkey"="Tsn") %>%
  as.data.frame()

NEFMC <- YrSps%>%
  group_by(NMFS.Namev2) %>%
  filter(Management.Group=="NE Multispecies") %>%
  droplevels()

colSums(is.na(NEFMC))

ggplot(NEFMC, aes(x=NMFS.Namev2, y=IYear)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust=1), 
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

NEFMC$NMFS.Name %>% 
  fct_count() %>%
  print(n = nrow(.))

SCALING_F<-1.0e6 #R Hundreds of millions

NEFMCv2<-NEFMC %>%
  filter(!NMFS.Namev2 %in% c("POUT, OCEAN", "WOLFFISH, ATLANTIC")) %>%
  filter(IYear >= 1981) %>%
  mutate(Value=Value/SCALING_F)

colSums(is.na(NEFMCv2))

#Aggregating to ignore species. Calculating total annual revenue and catch
LB_SUM<-aggregate(cbind(Catch, Value)~IYear, NEFMCv2, sum) %>%
  dplyr::rename(SumCT=Catch,
                SumVAL=Value)

LB_SUM2<-merge(NEFMCv2, LB_SUM, by=c("IYear")) %>%
  mutate(Price=Value/Catch) #R Hundreds of millions $ per metric ton

data<-LB_SUM2
#Use this info (YEAR through N) if you need to go line-by-line to look at how the code is functioning
YEAR=1993
YEARS<-as.integer(unlist(unique(data %>% filter(IYear<=YEAR) %>% distinct(IYear))))
j<-length(YEARS)
N<-j

#Run this function to add the decay and run and plot the frontiers generated
#from quantiles of the distribution of annual revenue up until time t.
#See "2_RiskGapCalc_Dec22.R" function to calculate actual portfolios and risk gaps.
portfolio_var <- function(data,YEAR,YEARS,TAXA,N,j,EBFM) {
  d_2 <- data %>% 
    filter(IYear %in% c(YEARS)) %>%
    aggregate(cbind(Value,Catch)~IYear+Taxonkey, FUN='sum') %>%
    arrange(IYear, Taxonkey, Value, Catch)

  ALL_YEAR <- unique(d_2$IYear)
  ALL_TAXA <- unique(d_2$Taxonkey)
  ALL_COMB <- expand.grid(ALL_YEAR,ALL_TAXA) %>%
    dplyr::rename(IYear=Var1,
           Taxonkey=Var2)

  d_2 <- merge(d_2,ALL_COMB, all=TRUE,by=c('IYear','Taxonkey')) %>%
    dplyr::mutate(Value = replace_na(Value, 0),
                  Catch = replace_na(Catch, 0),
                  Price = replace_na(Value/Catch,0))
  
  P_LAMBDA <- 0.549 #decay factor. If this is 1, each year is given equal weight.
  P_GAMMA <- 1 #sustainability parameter
  S_I <- 0.5 #R ?
  z = length(unique(d_2$Taxonkey))
  
  #Adding the decay factor into the covariance matrix (Eqn. 2)
  d_2 <- d_2 %>% #R Why "YEAR - IYear + 1"? Treats reference year as in the past. 
    mutate(WCatch = P_LAMBDA^(YEAR-IYear+1)*Catch, #check this with Geret because this was d_2$WCatch <- (P_LAMBDA^(YEAR-d_2$IYear+1))*d_2$Catch
           WAvgVal = WCatch*Price, #Value used to estimate omega. Omega= weighted average of catches over time with decay
           Lambda = P_LAMBDA^(YEAR-IYear+1),  #this is the decay value applied to each year.#ggplot (data=d_2, aes(x=desc(IYear), y=Lambda)) + geom_point() #check this with Geret because this was d_2$Lambda <- P_LAMBDA^(YEAR-d_2$IYear+1)
           WPrice = Lambda*Price,
           WVal = Lambda*Value)
  WAVG_CATCH <- aggregate(cbind(Value,WCatch,Lambda,WAvgVal,WPrice,WVal)~Taxonkey,data=d_2,FUN='sum') %>%
    mutate(WCatch = WCatch/Lambda, #R Annual average decay-weighted catch
           Omega = WAvgVal/WPrice, #R Looks good. average biomass. Omega= weighted average of catches over time with decay. Eqn. 3
           WVal = WVal/Lambda) #R Eq. 3
  MAX_CATCH <- aggregate(Catch~Taxonkey,data=d_2,FUN='max') %>%
    mutate(MaxCatch = Catch*P_GAMMA) %>%
    select(Taxonkey,MaxCatch)
  WAVG_CATCH <- merge(MAX_CATCH,WAVG_CATCH,by='Taxonkey') %>%
    mutate(MaxWeight = MaxCatch/Omega, #R Eq. 4. maximum weight = maximum biomass/average biomass)
           Year = YEAR,
           ImpWeight = Value/WVal, #R Summed net revenue divided by weighted average annual revenue. implicit weight from actual return.
           WRatio = ImpWeight/MaxWeight, #ratio of implicit to max weight
           CReturn = pmin(ImpWeight,MaxWeight)*WVal)
  DIFF_RETURNS <-cbind(sum(WAVG_CATCH$CReturn[which(WAVG_CATCH$Year==YEAR)]), sum(WAVG_CATCH$Value[which(WAVG_CATCH$Year==YEAR)]))
  
  #Estimating weighted Variance-Covariance matrix
  VAR_DATA <- merge(d_2,WAVG_CATCH, by='Taxonkey') %>%
    select('Taxonkey','IYear','Lambda.x','Value.x','WVal.y') %>%
    dplyr::rename(Year=IYear, #R Year
           Lambda=Lambda.x, #R Decay weight
           Value=Value.x, #R Annual rev
           WVal=WVal.y) %>% #R Eq. 3
    mutate(Value = Value-WVal) %>% #R Parentheticals eq. 4 
      select(!WVal)
  
  VAR_DATA1 <- subset(VAR_DATA, select=c('Taxonkey','Year','Value')) %>%
    cast(Taxonkey~Year)
  
  VAR_DATA <- VAR_DATA %>%
    mutate(Value = Value*Lambda) %>% #R Numerator eq. 4 (first two terms)
    select(!Lambda) %>%
    cast(Taxonkey~Year)
  
  VAR_names <- VAR_DATA$Taxonkey
  
  VAR_DATA <- matrix(unlist(subset(VAR_DATA, select=-Taxonkey)),nrow=z, ncol=j) #z=number of spp. in portfolio, j=years in portfolio
  VAR_DATA1 <- matrix(unlist(subset(VAR_DATA1, select=-Taxonkey)),nrow=z, ncol=j)
  VAR_COVAR <- VAR_DATA%*%t(VAR_DATA1) #R Eq. 4 numerator
  VAR_COVAR <- VAR_COVAR/WAVG_CATCH$Lambda #R Eq. 4
  pmean <- function(x,y) (x+y)/2 #R Parallel mean function for two vectors
  VAR_COVAR <- pmean(VAR_COVAR,matrix(VAR_COVAR,nrow(VAR_COVAR),byrow=TRUE)) ##????
  
  #Estimate the variance, which will be minimized
  VAR_COVAR_SS <- diag(diag(VAR_COVAR[])) # puts the vector above into a diagonal matrix for use in optimizing based on single species assumptions 
  MEAN_RETURNS <- matrix(WAVG_CATCH$WVal) #R Eq. 3, mean annual decayed returns
  TOT_REVENUE <- aggregate(cbind(Value)~IYear, data=d_2, FUN='sum') #R Tot rev/yr, all species
  #TARGET_RETURNS <- quantile(TOT_REVENUE$Value, seq(.00,1,by=.05)) #R Quantiles of total rev/yr, all species
  TARGET_RETURNS <- seq(1,max(TOT_REVENUE$Value),1) #R Quantiles of total rev/yr, all species
  MAX_WEIGHTS <- matrix(WAVG_CATCH$MaxWeight,nrow=z) #R Eq. 4, max harvest weight per species over time horizon 
  #CONSTRAINT_MATRIX <- diag(z)
  #CONSTRAINT_MATRIX <- -(rbind(CONSTRAINT_MATRIX,t(MEAN_RETURNS)))
  #CONSTRAINTS <- -c(MAX_WEIGHTS,TARGET_RETURNS)
  MIN_WEIGHTS <- matrix(0,nrow=z)
  
  OPTIMAL_WEIGHTS <- NULL
  OPTIMAL_REVENUE <- NULL
  OPTIMAL_LANDINGS <- NULL
  OPTIMAL_VAR <- NULL
  
  #Changed the margins from 10^-4 because of the computational singularity issue
  #with some of the years
  #TARGET<-50 
  for (TARGET in TARGET_RETURNS) {
    if(EBFM == TRUE){
      SOLUTION <- ipop(-MEAN_RETURNS,VAR_COVAR,t(-MEAN_RETURNS),-TARGET,
                       MIN_WEIGHTS, MAX_WEIGHTS, TARGET, margin=10^-2, sig=6,
                       verb=FALSE)
    } else {
      SOLUTION <- ipop(-MEAN_RETURNS,VAR_COVAR_SS,t(-MEAN_RETURNS),-TARGET,
                       MIN_WEIGHTS, MAX_WEIGHTS, TARGET, margin=10^-2, sig=6,
                       verb=FALSE)
    }
    
    
    #Pulling the optimum weights from the solution
    TEMP_WEIGHTS <- matrix(SOLUTION@primal) 
    colnames(TEMP_WEIGHTS) <- paste("Optimal_weights_",TARGET)
    
    #Calculating optimal revenue, by stock complex
    # TEMP_REVENUE <- MEAN_RETURNS*TEMP_WEIGHTS*SCALING_F  
    TEMP_REVENUE <- MEAN_RETURNS*TEMP_WEIGHTS
    colnames(TEMP_REVENUE) <- paste("Optimized_Revenue_",TARGET)
    
    #Calculate optimal landings, by stock complex
    # TEMP_LANDINGS <- WAVG_CATCH$WCatch*TEMP_WEIGHTS*SCALING_F
    TEMP_LANDINGS <- WAVG_CATCH$WCatch*TEMP_WEIGHTS
    colnames(TEMP_LANDINGS) <- paste("Optimized_Landings_",TARGET)
    
    #Calculate the realized variance for each optimization solution
    ###!!!! use full COV matirx here
    TEMP_VAR<-t(TEMP_WEIGHTS)%*%VAR_COVAR%*%TEMP_WEIGHTS
    # t(OPTIMAL_WEIGHTS)%*%VAR_COVAR%*%OPTIMAL_WEIGHTS)
    rownames(TEMP_VAR)<- paste("Optimized_Variance_",TARGET)
    
    # Accumulating the results for each optimization solution
    OPTIMAL_WEIGHTS <- cbind(OPTIMAL_WEIGHTS,TEMP_WEIGHTS)
    OPTIMAL_REVENUE <- cbind(OPTIMAL_REVENUE,TEMP_REVENUE)
    OPTIMAL_LANDINGS <- cbind(OPTIMAL_LANDINGS, TEMP_LANDINGS)
    OPTIMAL_VAR <- rbind(OPTIMAL_VAR,TEMP_VAR)
    
    print(TARGET) # check on progress of quadratic programming optimizer
  }
  
  # OPTIMAL_VAR <- (OPTIMAL_VAR*SCALING_F)^2
  #rownames(OPTIMAL_WEIGHTS) <- colnames(REVENUE[,1:stockNum])
  sum(OPTIMAL_WEIGHTS)
  #rownames(OPTIMAL_LANDINGS) <- colnames(REVENUE[,1:stockNum])
  #rownames(OPTIMAL_REVENUE) <- colnames(REVENUE[,1:stockNum])
  #colnames(OPTIMAL_VAR) <- "Optimized_Variance"
  MEAN_VAR <- cbind(OPTIMAL_VAR,colSums(OPTIMAL_REVENUE))
  colnames(MEAN_VAR) <- c("OptimizedVariance",'OptimizedRevenue')
  MEAN_VAR <- as.data.frame(MEAN_VAR) %>%
    mutate(OptimizedRevenue = OptimizedRevenue,
           OptimizedStDev = sqrt(OptimizedVariance))
  
  #You would need the below if calculating the actual portfolios and risk gaps
  #See "XXXXXXX" code to do that.
  IMPLICIT_WEIGHTS = cbind(WAVG_CATCH$ImpWeight, WAVG_CATCH$Taxonkey)
  colnames(IMPLICIT_WEIGHTS) <- c("ImpWeight",'Taxonkey')

  AP<-matrix(t(WAVG_CATCH$ImpWeight)%*%VAR_COVAR%*%WAVG_CATCH$ImpWeight)
  colnames(AP) <- 'ActualP'
  AP<-as.data.frame(cbind(AP, TARGET)) #you need to sqrt the AP to get SD

  FP<-matrix(t(OPTIMAL_WEIGHTS)%*%VAR_COVAR%*%OPTIMAL_WEIGHTS)
  colnames(FP) <- 'FrontierP'
  FP<-as.data.frame(cbind(FP, TARGET)) #you need to sqrt the FP to get SD
  
  
  ############### Returns a list of data frames ############
  optimums <- list ("MEAN_VAR"=MEAN_VAR,"OPTIMAL_WEIGHTS" =OPTIMAL_WEIGHTS, "COVARIANCE" =VAR_COVAR, "IMPLICIT_WEIGHTS" = IMPLICIT_WEIGHTS, "ACTUAL_PORTOLIO"= AP, "FRONTIER_PORTFOLIO"=FP)
  return(optimums)
} 

test<-portfolio_var(data=LB_SUM2,YEAR=1993,YEARS=YEARS,TAXA=Taxonkey,N=N,j=j, EBFM = TRUE)

## Problems with frontier portfolio, too many results (probably due to different optimal weights for different targets) - start here





#Loop through the years so the frontier is going up until time t. 
YEARS<-unique(LB_SUM2$IYear) #use if you are good looping through the whole time-series
# otherwise specify the timeframe below
# YEARS = 2012:2020
FrontierResults<-NULL
ImpResults<-NULL
Portfolio<-NULL

for (k in YEARS) {
  
  yr = min(YEARS):k
  N = as.numeric(length(yr))
  j = N
  LAST_YEAR = max(yr)
  
  a<-portfolio_var(data=LB_SUM2,YEAR= LAST_YEAR,YEARS=yr,TAXA=Taxonkey,N=N,j=j, EBFM = TRUE)
  b<-portfolio_var(data=LB_SUM2,YEAR= LAST_YEAR,YEARS=yr,TAXA=Taxonkey,N=N,j=j, EBFM = FALSE)
  
  
  EBFM_EF <- a$MEAN_VAR
  Type <- rep("EBFM",nrow(a$MEAN_VAR))
  EBFM_EF <- cbind(Type, EBFM_EF)
  
  SS_EF <- b$MEAN_VAR
  Type <- rep("SS",nrow(b$MEAN_VAR))
  SS_EF <- cbind(Type, SS_EF)
  
  All_EF <- rbind(EBFM_EF, SS_EF)
  # All_EF <- rbind(All_EF, SS_EF)
  All_EF$Iteration<-LAST_YEAR
  Iteration<-LAST_YEAR
  
  FrontierResults <- rbind(All_EF, FrontierResults)
  
  ImpWeights <- a$IMPLICIT_WEIGHTS
  Iteration<-LAST_YEAR
  # Taxonkey<-
  ImpWeights <- cbind(Iteration, ImpWeights)

  ImpResults <- rbind(ImpWeights, ImpResults)
  ImpResults <-as.data.frame(ImpResults)
  
  Act_Portfolio <- a$ACTUAL_PORTOLIO
  Front_Portolio<- a$FRONTIER_PORTFOLIO
  AF<-cbind(Iteration, Act_Portfolio, Front_Portolio)
  
  Portfolio<-rbind(AF, Portfolio)
  
  # ImpWeights <- a$IMPLICIT_WEIGHTS
  # Iteration<-LAST_YEAR
  # # Taxonkey<-
  # ImpWeights <- cbind(Iteration, ImpWeights)
  # 
  # ImpResults <- rbind(ImpWeights, ImpResults)
  # ImpResults <-as.data.frame(ImpResults)
  # 
  # Act_Portfolio <- a$ACTUAL_PORTOLIO
  # Front_Portolio<- a$FRONTIER_PORTFOLIO
  # AF<-cbind(Iteration, Act_Portfolio, Front_Portolio)
  # 
  # Portfolio<-rbind(AF, Portfolio)
  
  print(Iteration) #This prints out the Year at the end of each iteration of the loop
  #so that if you have problems with convergence you can see which year it occurs in
}

money<-LB_SUM2%>%
  group_by(IYear)%>%
  mutate(Dollars2 = Value) %>%
  summarise(val = sum(Dollars2 , na.rm =TRUE),
            sd = sd(Dollars2 , na.rm = TRUE)) %>%
  mutate(Iteration = IYear)

ggplot(data=money, aes(x=IYear, y=val)) + geom_point()

### plot to check, but removing the first year (min(Iteration))
# ggplot() +
#   geom_line(data =FrontierResults %>%
#               filter(!Iteration==min(Iteration)), aes(x= OptimizedStDev, y=OptimizedRevenue, color=Type))+
#   # filter(!Iteration==min(Iteration)), aes(x= OptimizedStDev/1000000, y=OptimizedRevenue/1000000, color=Type))+
#   # geom_point(data =money, aes(x = sd/1000000, y = val/1000000))+
#   # geom_text(data = money, aes(x = sd/1000000, y = val/1000000,label = as.factor(IYear)), hjust=0.25, vjust=2) +
#   labs(y ="Revenue (Hundred Million Dollars)", x = "Risk (SD of Revenue)") +
#   facet_wrap(~Iteration, scales="free")
# # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "NE Multispecies 11 spp with decay.png",
# #        width = 25, height = 12, dpi = 120) 

ggplot() +
  geom_line(data=FrontierResults %>%
              filter(Iteration!=min(Iteration)), aes(x= OptimizedStDev, y=OptimizedRevenue, color=Type))+
  # filter(!Iteration==min(Iteration)), aes(x= OptimizedStDev/1000000, y=OptimizedRevenue/1000000, color=Type))+
  # geom_point(data =money, aes(x = sd/1000000, y = val/1000000))+
  # geom_text(data = money, aes(x = sd/1000000, y = val/1000000,label = as.factor(IYear)), hjust=0.25, vjust=2) +
  labs(y ="Revenue (Million Dollars)", x = "Risk (SD of Revenue)") +
  facet_wrap(~Iteration, scales="free")
