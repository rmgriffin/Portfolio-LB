portfolio_risk <- function(data,YEAR,YEARS,TAXA,N,j,EBFM) {
  d_2 <- data %>% 
    filter(IYear %in% c(YEARS)) %>%
    aggregate(cbind(Value,Catch)~IYear+Taxonkey, FUN='sum') %>%
    arrange(IYear, Taxonkey, Value, Catch)
  
  ALL_YEAR <- unique(d_2$IYear)
  ALL_TAXA <- unique(d_2$Taxonkey)
  ALL_COMB <- expand.grid(ALL_YEAR,ALL_TAXA) %>%
    rename(IYear=Var1,
           Taxonkey=Var2)
  
  d_2 <- merge(d_2,ALL_COMB, all=TRUE,by=c('IYear','Taxonkey')) %>%
    dplyr::mutate(Value = replace_na(Value, 0),
                  Catch = replace_na(Catch, 0),
                  Price = replace_na(Value/Catch,0))
  
  P_LAMBDA <- 0.549 #decay factor. If this is 1, each year is given equal weight.
  P_GAMMA <- 1 #sustainability parameter
  S_I <- 0.5
  z = length(unique(d_2$Taxonkey))
  
  #Adding the decay factor into the covariance matrix (Eqn. 2)
  d_2 <- d_2 %>%
    mutate(WCatch = P_LAMBDA^(YEAR-IYear+1)*Catch, #check this with Geret because this was d_2$WCatch <- (P_LAMBDA^(YEAR-d_2$IYear+1))*d_2$Catch
           WAvgVal = WCatch*Price, #Value used to estimate omega. Omega= weighted average of catches over time with decay
           Lambda = P_LAMBDA^(YEAR-IYear+1),  #this is the decay value applied to each year.#ggplot (data=d_2, aes(x=desc(IYear), y=Lambda)) + geom_point() #check this with Geret because this was d_2$Lambda <- P_LAMBDA^(YEAR-d_2$IYear+1)
           WPrice = Lambda*Price,
           WVal = Lambda*Value)
  WAVG_CATCH <- aggregate(cbind(Value,WCatch,Lambda,WAvgVal,WPrice,WVal)~Taxonkey,data=d_2,FUN='sum') %>%
    mutate(WCatch = WCatch/Lambda,
           Omega = WAvgVal/WPrice, #average biomass. Omega= weighted average of catches over time with decay. Eqn. 3
           WVal = WVal/Lambda)
  MAX_CATCH <- aggregate(Catch~Taxonkey,data=d_2,FUN='max') %>%
    mutate(MaxCatch = Catch*P_GAMMA) %>%
    select(Taxonkey,MaxCatch)
  WAVG_CATCH <- merge(MAX_CATCH,WAVG_CATCH,by='Taxonkey') %>%
    mutate(MaxWeight = MaxCatch/Omega, #maximum weight = maximum biomass/average biomass)
           Year = YEAR,
           ImpWeight = Value/WVal, #implicit weight from actual return.
           WRatio = ImpWeight/MaxWeight, #ratio of implicit to max weight
           CReturn = pmin(ImpWeight,MaxWeight)*WVal)
  DIFF_RETURNS <-cbind(sum(WAVG_CATCH$CReturn[which(WAVG_CATCH$Year==YEAR)]), sum(WAVG_CATCH$Value[which(WAVG_CATCH$Year==YEAR)]))
  
  #Estimating weighted Variance-Covariance matrix
  VAR_DATA <- merge(d_2,WAVG_CATCH, by='Taxonkey') %>%
    select('Taxonkey','IYear','Lambda.x','Value.x','WVal.y') %>%
    rename(Year=IYear,
           Lambda=Lambda.x,
           Value=Value.x,
           WVal=WVal.y) %>%
    mutate(Value = Value-WVal) %>%
    select(!WVal)
  
  VAR_DATA1 <- subset(VAR_DATA, select=c('Taxonkey','Year','Value')) %>%
    cast(Taxonkey~Year)
  
  
  VAR_DATA <- VAR_DATA %>%
    mutate(Value = Value*Lambda) %>%
    select(!Lambda) %>%
    cast(Taxonkey~Year)
  
  VAR_names <- VAR_DATA$Taxonkey
  
  VAR_DATA <- matrix(unlist(subset(VAR_DATA, select=-Taxonkey)),nrow=z, ncol=j) #z=number of spp. in portfolio, j=years in portfolio
  VAR_DATA1 <- matrix(unlist(subset(VAR_DATA1, select=-Taxonkey)),nrow=z, ncol=j)
  VAR_COVAR <- VAR_DATA%*%t(VAR_DATA1)
  VAR_COVAR <- VAR_COVAR/WAVG_CATCH$Lambda
  pmean <- function(x,y) (x+y)/2
  VAR_COVAR <- pmean(VAR_COVAR,matrix(VAR_COVAR,nrow(VAR_COVAR),byrow=TRUE)) ##????
  
  ###########
  
  VAR_COVAR_SS <- diag(diag(VAR_COVAR[])) # puts the vector above into a diagonal matrix for use in optimizing based on single species assumptions 
  
  MEAN_RETURNS <- matrix(WAVG_CATCH$WVal)
  
  TOT_REVENUE <- aggregate(cbind(Value)~IYear, data=d_2, FUN='sum')
  
  TARGET_RETURNS <- sum(d_2$Value[which(d_2$IYear==YEAR)]) #Geret
  # TARGET_RETURNS2 <- d_2 %>%
  #   group_by(IYear)
  #   summarise(sum(Value))
  
  MAX_WEIGHTS <- matrix(WAVG_CATCH$MaxWeight,nrow=z)
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
  # TARGET<-0.5903352 
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
  
  #rownames(OPTIMAL_WEIGHTS) <- colnames(REVENUE[,1:stockNum])
  sum(OPTIMAL_WEIGHTS)
  #rownames(OPTIMAL_LANDINGS) <- colnames(REVENUE[,1:stockNum])
  #rownames(OPTIMAL_REVENUE) <- colnames(REVENUE[,1:stockNum])
  #colnames(OPTIMAL_VAR) <- "Optimized_Variance"
  MEAN_VAR <- cbind(OPTIMAL_VAR,colSums(OPTIMAL_REVENUE))
  colnames(MEAN_VAR) <- c("OptimizedVariance",'OptimizedRevenue')
  MEAN_VAR <- as.data.frame(MEAN_VAR)
  MEAN_VAR$OptimizedRevenue <- MEAN_VAR$OptimizedRevenue
  MEAN_VAR$OptimizedStDev <- sqrt(MEAN_VAR$OptimizedVariance)
  
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
  
  #Loop through the years so the frontier is going up until time t.
  YEARS<-unique(LB_SUM2$IYear)
  # YEARS = 1981:2016
  RiskResults<-NULL
  ImpResults<-NULL
  Portfolio<-NULL
  
  for (k in YEARS) {
    
    yr = min(YEARS):k
    N = as.numeric(length(yr))
    j = N
    LAST_YEAR = max(yr)
    
    c<-portfolio_risk(data=LB_SUM2,YEAR=LAST_YEAR,YEARS=yr,TAXA=Taxonkey,N=N,j=j, EBFM = TRUE)
    d<-portfolio_risk(data=LB_SUM2,YEAR=LAST_YEAR,YEARS=yr,TAXA=Taxonkey,N=N,j=j, EBFM = FALSE)
    
    
    EBFM_EF_Risk <- c$MEAN_VAR
    Type <- rep("EBFM",nrow(c$MEAN_VAR))
    EBFM_EF_Risk <- cbind(Type, EBFM_EF_Risk)
    
    SS_EF_Risk <- d$MEAN_VAR
    Type <- rep("SS",nrow(d$MEAN_VAR))
    SS_EF_Risk <- cbind(Type, SS_EF_Risk)
    
    All_EF_Risk <- rbind(EBFM_EF_Risk, SS_EF_Risk)
    # All_EF <- rbind(All_EF, SS_EF)
    All_EF$Iteration<-LAST_YEAR
    Iteration<-LAST_YEAR
    
    RiskResults <- rbind(All_EF_Risk, RiskResults)
    
    ImpWeights <- c$IMPLICIT_WEIGHTS
    Iteration<-LAST_YEAR
    # Taxonkey<-
    ImpWeights <- cbind(Iteration, ImpWeights)
    
    ImpResults <- rbind(ImpWeights, ImpResults)
    ImpResults <-as.data.frame(ImpResults)
    
    Act_Portfolio <- c$ACTUAL_PORTOLIO
    Front_Portolio<- c$FRONTIER_PORTFOLIO
    AF<-cbind(Iteration, Act_Portfolio, Front_Portolio)
    
    Portfolio<-rbind(AF, Portfolio)
    
    print(Iteration)
    
  }
  
  #TARGET here is the annual revenue/SCALING_F for each year of the timeseries
  #Drop the second TARGET column otherwise you will get an error when you try to plot TARGET
  Portfolio<-Portfolio %>%
    select(c(1:4)) 
  
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
  
  ### plot to check. You get FrontierResults from running the code in "1_FrontierCalc_Quantiles_Decay_Dec22"
  ### Removes the first year (min(Iteration))ggplot() +
  ggplot() +
    geom_line(data =FrontierResults %>%
                # filter(Type=='EBFM') %>% #Use this line if you want to drop the SS blue line.
                filter(!Iteration==min(Iteration)), 
              aes(x= OptimizedStDev, y=OptimizedRevenue, color=Type))+
    geom_point(data =Portfolio %>%
                 filter(!Iteration==min(Iteration)), aes(x = sqrt(ActualP), y = TARGET))+
    # geom_point(data =Portfolio, aes(x = sqrt(FrontierP), y = TARGET*SCALING_F/1000000), col="red")+
    # geom_text(data = Portfolio, aes(x = sqrt(FrontierP), y = TARGET*SCALING_F/1000000),label = as.factor(Portfolio$Iteration), hjust=0.25, vjust=2) +
    labs(y ="Revenue (Hundred Million Dollars)", x = "Risk (SD of Revenue)") +
    facet_wrap(~Iteration, scales="free")
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "NE Multispecies 11 spp with decay and actual portfolio.png",
  # width = 20, height = 15, dpi = 120, bg='transparent')
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting/SSremoved", filename = "NE Multispecies 11 spp with decay and actual portfolio.png",
  #        width = 20, height = 15, dpi = 120, bg='transparent')
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "Top 30 by MetricLandings with decay and actual portfolio.png",
  # width = 20, height = 15, dpi = 400, bg='transparent')
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "Top 30 by Revenue with decay and actual portfolio.png",
  #        width = 20, height = 15, dpi = 400, bg='transparent')
  
  ###############################################
  # Risk gap calculation: #This is the numerator of Eqn. 6 in Jin et al. 2016
  Portfolio <- Portfolio %>%
    mutate(Riskgap=sqrt(ActualP)-sqrt(FrontierP)) 
  
  #Now you need to calculate the denominator of Eqn. 6. Edit for the timeperiod in question.
  LRB<-LB_SUM2 %>%
    # filter(!IYear %in% c(2017:2020)) %>%
    group_by(IYear, Taxonkey) %>%
    summarise(AnnualRev=sum(Value)) %>%
    dplyr::rename(Iteration=IYear)
  # mutate(lb=lb*SCALING_F) %>%
  
  #ImpResults contains the implicit weights (w~,t) for each species, each year, from the function above. 
  #Merge those with the dataframe containing annual revenue by species (μt).
  LRB<-merge(LRB, ImpResults, by=c("Iteration", "Taxonkey"))
  
  #Loop over the years to calculate the denominator for Eqn. 6
  #This is doing matrix algebra using species implicit weights and revenues for each year.
  IT=min(YEARS):max(YEARS)
  output<-NULL
  
  for (m in IT) {
    
    df<-LRB %>%
      filter(Iteration==m)
    
    denom<-t(df$ImpWeight)%*%df$AnnualRev
    colnames(denom) <- c("Denominator")
    # denom$V1 <-colnames("XXXX")
    Iteration<-unique(df$Iteration)
    
    denom2<-cbind(denom, Iteration)
    
    output<-rbind(denom2, output)
  }
  
  df2<-merge(Portfolio, output, by=c("Iteration")) %>%
    # mutate(Riskgap_PerDollar = (Riskgap*1.0e8)/(Denominator*1.0e8))
  mutate(Riskgap_PerDollar = Riskgap/Denominator)
  
  ggplot(data=df2, aes(x=Iteration, y=Riskgap_PerDollar)) +
    geom_line() +
    geom_point() +
    xlab("Year") +
    ylab("Riskgap") +
    scale_x_continuous(breaks = seq(min(df2$Iteration), max(df2$Iteration), by = 5))
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "NE Multispecies 11 spp Riskgap.png",
  # width = 12, height = 10, dpi = 120, bg='transparent')
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "Top 30 by MetricLandings Riskgap.png",
  # width = 12, height = 10, dpi = 400, bg='transparent')
  # ggsave(path = "Figures for Dec 22 Steering Committee meeting", filename = "Top 30 by Revenue Riskgap.png",
  #        width = 12, height = 10, dpi = 400, bg='transparent')
  
  
  