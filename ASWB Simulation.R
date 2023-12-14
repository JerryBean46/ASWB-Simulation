## Set Working Directory
setwd("E:/ASWB Exam Project")

library(mirt)
library(psych)
library(lessR)

## Get DIF items function
get.dif.items <- function(f.data,p.val=.05,parms){
  r.warnings = ""
  keep.vars <- c("X2", "df", "p") # just keep these variables
  f.data <- f.data[keep.vars]
  f.data$p = round(f.data$p,3)
  if(missing(f.data)) return('Missing model output out.list')
  f.data$sig <- ifelse(f.data$p < p.val,'dif','no_dif')
  if(!missing(parms)){
    if(nrow(f.data) == nrow(parms)){
      f.data <- cbind(f.data,parms)
    }else{
      r.warnings = "There number of item parameters doesn't match the number of items "
      r.warnings = paste(r.warnings,"given to get.dif.items. Item parameters omitted.")
    }
  }
  dif.items <- subset(f.data, sig == 'dif')
  no.dif.items <- subset(f.data, sig == 'no_dif')
  if(!missing(parms) && nrow(f.data) == nrow(parms)){
    if(nrow(no.dif.items)>1){
      no.dif.items <- no.dif.items[order(-no.dif.items$a1),]
    }
  }
  r.list <- list(dif_items = dif.items, no_dif = no.dif.items, warnings = r.warnings)
  return(r.list)
}

## Make data function
make.data<-function(N){ 
  set.seed(12345) 
  a <- matrix(abs(rnorm(15,1,.3)), ncol=1) 
  d <- matrix(rnorm(15,0,.7), ncol=1) 
  d1 <- d2 <- cbind(d) # b parameters for both groups 
  d2[c(2,5,7,9,15),] <- d1 [c(2,5,7,9,15),] + .52 # here is the DIF 
  itemtype <- rep('2PL',nrow(a)) 
  dataset1 <- simdata(a,d1,N,itemtype) 
  dataset2 <- simdata(a,d2,N,itemtype) 
  dat<-rbind(dataset1, dataset2) 
  return(dat) 
} 

## Data details
N<-1000 
dat <- make.data(N) 
group <- c(rep('Ref',N), rep('Foc',N)) 
foc.data <- dat[1:1000,] 
ref.data <- dat[1001:2000,] 
rm(N)

### Step 1: Get Item Parameters - Freely estimated model
model.free <- multipleGroup(dat, 1, group, verbose = FALSE)
coef(model.free, simplify = TRUE)          

### Step 2: Baseline Model 
model.constrained <- multipleGroup(dat, model = 1, group = group, verbose=FALSE,
                                   invariance = c(colnames(dat), 'free_means', 'free_var'))
(constrained.parameters <- coef(model.constrained, simplify = TRUE)[[1]][[1]])  

### Step 3: First round of DIF analyses (LRTs) - All Others As Anchors 
(dif.drop <- DIF(model.constrained, c('a1','d'), scheme = 'drop', seq_stat = .05))

## use the optional function to table the output
get.dif.items(f.data = dif.drop, p.val=.05, parms=constrained.parameters)

### Step 4: Specify a New Baseline Model Using Anchor Items (MaxA5 Approach) 
itemnames <- colnames(dat)
anc.items.names <- itemnames[c(12,13,1,11,4)] 
test.items <- c(2,3,5,7,8,9,10,14,15) 

model_anchor <- multipleGroup(dat, model = 1, group = group, verbode=FALSE, SE = TRUE,
                              invariance = c(anc.items.names, 'free_means', 'free_var'))
(anchor.parms <-coef(model_anchor,simplify = TRUE)[[1]][[1]])

### Step 6: Run Final Invariance Tests  
(dif.anchor <- DIF(model_anchor, c('a1','d'), items2test = test.items, 
                   p.adjust="BH"))

## use the optional function to table the output
get.dif.items(f.data=dif.anchor,p.val=.05)

### Step 5: Compute Item and Test Effect Sizes 
empirical_ES(model_anchor, plot=FALSE, as.table = TRUE) # item level effect sizes
empirical_ES(model_anchor, DIF=FALSE, plot=FALSE) # test level effect sizes

## Item and Test Plots
empirical_ES(model_anchor, plot=TRUE) # expected item score plots
empirical_ES(model_anchor, DIF=FALSE, plot=TRUE) # expected test score plot



