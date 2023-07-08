library(tidyverse)
library(vegan)
library(hillR)
library(directlabels)
library(scales)
library(grid)
library(gridExtra)
library(metR)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(gridtext)
library(here)
library(rmarkdown)
library(knitr)
library(patchwork)

#assumptions:
#sex ratio at birth=1:1
#groupsize is constant, if new arrives takes previous
#no breeding outside the unit
reproduce_fun<-function(dataorig,datares){
  #some default 0 
  fbin<-fbinall<-0
  main<-mbin<-mbinall<-mabinall<-0

  #if there are more than 1 female
  if(dataorig$fgs>1){
    #all female beta production
    #what f alpha did not get from each of them
    fbinall<-(1-dataorig$farepmono)*(dataorig$fgs-1)*dataorig$foff    
    #fbetas should get their production tamed down by farepmono, l.36
    fbin<-fbinall/(dataorig$fgs-1)
  }
  
  #what is female alpha taking from other females share?
  fatake<-(dataorig$fgs-1)*dataorig$foff*dataorig$farepmono
  #average f alpha produced offspring 
  fain<-ifelse(dataorig$coopbr==0,#cooperative breeding?
              dataorig$foff,#if no fa takes baseline f reprod
              #if there is cooperative breeding, add benefits from suppressing each additional female in the group
              dataorig$foff+(fatake*dataorig$fagainrepsup))
  
  #If there are no betas, all female offspring sired by alpha male no matter marepmono
  if(dataorig$mgs==1){
    main<-(fain+fbinall)
  }else{
    share<-(1/dataorig$mgs)*(fain+fbinall) #equal theoretical share among males from female reproduction
    #if there is no monopolization, 0 means everything is shared equally
    if(dataorig$marepmono==0){
      main<-share
      mbin<-share
      mbinall<-share*(dataorig$mgs-1)
    }else{#if monopolization
      #male alpha gets marepmono from each share of other males
      main<-share+(share*(dataorig$mgs-1))*dataorig$marepmono
      #other males what is left
      mbinall<-(fain+fbinall)-main
      mbin<-mbinall/(dataorig$mgs-1)
    }
    
  }
  
  #take the relevant portion of the reproduction table and add results
  datares<-datares[which(datares$fbd==dataorig$fbd & datares$mbd==dataorig$mbd),]
  datares$noffs[which(datares$sexpa=="f"&datares$rankpa=="a")]<-fain
  datares$noffs[which(datares$sexpa=="f"&datares$rankpa=="b")]<-fbin
  datares$noffs[which(datares$sexpa=="m"&datares$rankpa=="a")]<-main
  datares$noffs[which(datares$sexpa=="m"&datares$rankpa=="b")]<-mbin
  return(datares)
}

#makes a function to apply statistical tests per group and variable using rstatix package
make_test <- function (data,variable,grouping_variable,type){
  if(type=="dunn"){
    tst<-data %>% 
      rstatix::dunn_test(reformulate(grouping_variable, variable)) %>%
      rstatix::adjust_pvalue("holm",p.col="p")%>%
      rstatix::add_xy_position(x = grouping_variable) %>% 
      rstatix::add_significance(p.col="p")
  }
  if(type=="wilcox"){
    tst<-data %>% 
      rstatix::wilcox_test(reformulate(grouping_variable, variable)) %>%
      rstatix::add_xy_position(x = grouping_variable) %>% add_significance
  }
  return(tst)

}

#Answer from ytsaig copied from
#https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))


#END----------------