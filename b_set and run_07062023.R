source("scripts/a_fun_07062023.R")
#1.Generate data----------
#1.1 social variation----
#Group size----
varsoc<-c(1:25)


#Calculate population size-----
# allowing a number of units needed for all offspring to become alpha
#in the most produtive system with highest group size

#number of offspring produced by females, key variable
foffs<-10

nbaoa<-max(varsoc)*foffs/2
mgsp<-max(varsoc)+1#minimum group size that can afford most productivity all females + 1 male
maxgs<-max(varsoc)*2
totalpop<-maxgs*nbaoa#how much larger than max(varsoc) should the population be?


#Reproductive monopolization-----
varrepmono.fct<-c("none","medium","total")
params<-list(fgs=varsoc,# female group size
             mgs=varsoc,# male group size
             farepmono.fct=varrepmono.fct, #female alpha monopolization
             marepmono.fct=varrepmono.fct) #male alpha monopolization 
params<-expand.grid(params)#create all parameter combinations


#repmono values indicate
#the part of the equal share from each individual 
#that the alpha "steals" from the rest
params<-params %>% mutate(farepmono=case_when(farepmono.fct=="none"~0,
                                              farepmono.fct=="medium"~0.5,
                                              farepmono.fct=="total"~1),
                          marepmono=case_when(marepmono.fct=="none"~0,
                                              marepmono.fct=="medium"~0.5,
                                              marepmono.fct=="total"~1))
params

#Different taxa commonly studied
#mbd=male dispersal, fbd=female dispersal, fngenes=female nb genes, mngenes=idem male,
#ffpg=proportion of female genes to female offspring,fmpg=proportion of female genes to female offspring,
#mfpg,mmpg=idem from male parent,
#coopbr= falpha gain reprod from fbeta suppression (cooperative breeding)
mamm<-expand.grid(taxa="mammal",mbd=1,fbd=0,fngenes=c(10),mngenes=c(10),ffpg=c(0.5),fmpg=c(0.5),mfpg=c(0.5),mmpg=c(0.5),coopbr=c(0,1))
bird<-expand.grid(taxa="bird",  mbd=0,fbd=1,fngenes=c(10),mngenes=c(10),ffpg=c(0.5),fmpg=c(0.5),mfpg=c(0.5),mmpg=c(0.5),coopbr=c(0,1))
wasp<-expand.grid(taxa="wasp",  mbd=1,fbd=0,fngenes=c(10),mngenes=c(5),ffpg=c(0.5),fmpg=c(0.5),mfpg=c(1),mmpg=c(0),coopbr=c(0,1))
taxa<-bind_rows(list(mamm,bird,wasp)) %>% mutate(taxa=factor(taxa,levels=c("wasp","mammal","bird"),labels=c("wasp","mammal","bird")))

#get both taxa variation and social variation in single df
taxa.params<-expand.grid.df(taxa,params)
str(taxa.params)

#Further reproductive parameters
#fcare,mcare : who cares for offspring
 #not needed current investigation

#foffs: the baseline nb of offspring females can produce
#fgainrepsup: the proportion of the beta female reproduction that is added to the alpha female when eusocial
fagainrepsupr<-1
taxa.params<-taxa.params %>% mutate(
  fcare=1,
  mcare=0,
  foff=foffs,
  fagainrepsup=fagainrepsupr,#case_when(taxa=="mammal.eu"~0.2,taxa=="bird.eu"~0.2,taxa=="wasp.eu"~0.9,TRUE~0)
  tnbunits=totalpop/(fgs+mgs),#calculate number of units in the population
  gs=(fgs+mgs)
  )



#2. Anagenesis-----
#2.1 offspring production per unit-----
#2.1.1 Previous generation-------
#establish number of genes inherited from outside group or from inside to produce current generation
pa.rep.prev<-expand.grid(list(sex=c("f","m"),
                              rank=c("a","b"),
                              orig=c("in","out")))
pa.rep.prev<-pa.rep.prev %>% mutate(when="prev",
                                    nids=1,
                                    ngenes=0)

taxa.pa.rep.prev<-map(1:nrow(taxa.params),~{
  paramorig<-taxa.params[.x,]
  if(paramorig$mbd==1){
    #if males disperse
    #all genes in current adult males came from outside
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="m"& pa.rep.prev$orig=="in")]<-0
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="m" & pa.rep.prev$orig=="out")]<-paramorig$mngenes
  }else{
    #if males stay
    #origin "in" in males are genes coming from the father transmitted to sons
    #origin "out" in males are genes coming from the mother transmitted to sons
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="m" & pa.rep.prev$orig=="in")]<-paramorig$mngenes*paramorig$mmpg#transmitted from father
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="m" & pa.rep.prev$orig=="out")]<-paramorig$fngenes*paramorig$fmpg#transmitted from mother
  }

  if(paramorig$fbd==1){
    #if females disperse
    #all genes in current adult females came from outside
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="f" & pa.rep.prev$orig=="in")]<-0
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="f" & pa.rep.prev$orig=="out")]<-paramorig$fngenes
  }else{
    #if females stay
    #origin "in" in females are genes coming from the mother
    #origin "out" in females are genes coming from the father
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="f" & pa.rep.prev$orig=="in")]<-paramorig$fngenes*paramorig$ffpg
    pa.rep.prev$ngenes[which(pa.rep.prev$sex=="f" & pa.rep.prev$orig=="out")]<-paramorig$mngenes*paramorig$mfpg
  }
  #add number of individuals per rank
  pa.rep.prev$nids[which(pa.rep.prev$rank=="a")]<-1
  pa.rep.prev$nids[which(pa.rep.prev$sex=="f" & pa.rep.prev$rank=="b")]<-paramorig$fgs-1
  pa.rep.prev$nids[which(pa.rep.prev$sex=="m" & pa.rep.prev$rank=="b")]<-paramorig$mgs-1

  #actualize nb of genes according to the number of individuals
  pa.rep.prev$ngenes<-pa.rep.prev$ngenes*pa.rep.prev$nids
  return(pa.rep.prev)
})


#2.1.2 Current generation---------
#calculate the nb of offspring produced for each sex and rank class
pa.rep<-expand.grid(list(sexpa=c("f","m"),rankpa=c("a","b"),fbd=c(0,1),mbd=c(0,1)))
pa.rep<-pa.rep %>% mutate(origpa=case_when(sexpa=="f" & fbd==1~"out",
                                           sexpa=="f" & fbd==0~"in",
                                           sexpa=="m" & mbd==1~"out",
                                           sexpa=="m" & mbd==0~"in"),
                          when="curr",
                          noffs=0,
                          ngenes=0)
pa.rep

#same probabilities apply to all because all groups are assumed same
taxa.pa.rep<-map(1:nrow(taxa.params),~{reproduce_fun(taxa.params[.x,],pa.rep)})

#now divide offspring into respective sexes, assuming 1:1 sexratio at birth
taxa.pa.rep<-map(1:nrow(taxa.params),~{
  paramorig<-taxa.params[.x,]
  reprod<-expand.grid.df(taxa.pa.rep[[.x]],data.frame(sexoffs=c("f","m"),stringsAsFactors = F)) %>%
    transform(noffs=noffs/2)
  #count number of genes transmitted from father and mother, important for haplodiploidy diplodiplody difference
  reprod<-reprod %>% mutate(ngenesoffs=case_when(sexpa=="m"&sexoffs=="m"~noffs*paramorig$mngenes*paramorig$mmpg,
                                                 sexpa=="m"&sexoffs=="f"~noffs*paramorig$mngenes*paramorig$mfpg,
                                                 sexpa=="f"&sexoffs=="f"~noffs*paramorig$fngenes*paramorig$ffpg,
                                                 sexpa=="f"&sexoffs=="m"~noffs*paramorig$fngenes*paramorig$fmpg))
  reprod
})


#2.1.3 Next generation ----------
#carry forward the reproduction process,
#If you are the most successful individual in your group because you have an advantageous mutation
# and your offspring inherited it, they attain all breeding positions possible,
# how many grand offspring genes will you have produced?

taxa.pa.disp<-map(1:nrow(taxa.params),~{
  #this does not work for both sexes dispersing yet
  paramorig<-taxa.params[.x,]
  reprod<-taxa.pa.rep[[.x]]
  #how many alpha and beta positions are there in the population?
  #total number of social units for a fixed pop size
  tnbunits<-paramorig$tnbunits

  #number of male alpha, female alpha and betas in the population
  nma<-tnbunits*paramorig$mbd*1
  nfa<-tnbunits*paramorig$fbd*1
  nmb<-tnbunits*paramorig$mbd*(paramorig$mgs-1)
  nfb<-tnbunits*paramorig$fbd*(paramorig$fgs-1)

  #take the alpha individual from the dispersing sex (it has a mutation that is advantageous for becoming breeder)
  mam<-with(reprod,which(sexpa=="m" & rankpa=="a" & sexoffs=="m"))
  faf<-with(reprod,which(sexpa=="f" & rankpa=="a" & sexoffs=="f"))

  #number of male and female dispersing
  nmdisp<-reprod[mam,"noffs"]*paramorig$mbd
  nfdisp<-reprod[faf,"noffs"]*paramorig$fbd

  #for wasps, need to start with female line because no transmission through the male to male line
  if(paramorig$mbd==1 & paramorig$mmpg==0){
    fam<-with(reprod,which(sexpa=="f" & rankpa=="a" & sexoffs=="m"))
    nmdisp<-reprod[fam,"noffs"]*paramorig$mbd
  }

  #calculate the number of breeding positions obtained by offspring depending on availability in pop
  nmdispa<-ifelse(nma>=nmdisp,nmdisp,nma)#
  leftm<-nmdisp-nmdispa
  nmdispb<-ifelse(nmb>=leftm,leftm,nmb)

  nfdispa<-ifelse(nfa>=nfdisp,nfdisp,nfa)#
  leftf<-nfdisp-nfdispa
  nfdispb<-ifelse(nfb>=leftf,leftf,nfb)


  #add results to table
  #take the alpha individual from the dispersing sex (it has a mutation that is advantageous for becoming breeder)
  sextofind<-case_when(paramorig$mbd==1 & paramorig$fbd==1~c("m","f"),
                       paramorig$mbd==1 & paramorig$fbd==0~c("m",NA_character_),
                       paramorig$mbd==0 & paramorig$fbd==1~c("f",NA_character_),
                       TRUE~NA_character_)
  tf<-with(reprod,which(sexpa%in%sextofind & rankpa=="a" & sexoffs%in%sextofind))
  
  #for wasps, need to start with female line because no transmission through the male to male line
  if(paramorig$mbd==1 & paramorig$mmpg==0){
    sextofind1<-"f"
    sextofind2<-"m"
    tf<-with(reprod,which(sexpa%in%sextofind1 & rankpa=="a" & sexoffs%in%sextofind2))
  }

  #with number of dispersed offspring per sex and rank
  #calculate the number of genes from the alpha parent that were carried by them
  dispersed<-reprod[tf,]
  dispersed<-expand.grid.df(dispersed,data.frame(rankoffs=c("a","b"),stringsAsFactors = F))
  dispersed$noffs<-dispersed$ngenesoffs<-0
  dispersed<-dispersed %>% mutate(noffs=case_when(sexoffs=="m"&rankoffs=="a"~nmdispa,
                                                  sexoffs=="m"&rankoffs=="b"~nmdispb,
                                                  sexoffs=="f"&rankoffs=="a"~nfdispa,
                                                  sexoffs=="f"&rankoffs=="b"~nfdispb),
                                  ngenesoffs=case_when(sexpa=="m"&sexoffs=="m"~noffs*paramorig$mngenes*paramorig$mmpg,
                                                       sexpa=="m"&sexoffs=="f"~noffs*paramorig$mngenes*paramorig$mfpg,
                                                       sexpa=="f"&sexoffs=="f"~noffs*paramorig$fngenes*paramorig$ffpg,
                                                       sexpa=="f"&sexoffs=="m"~noffs*paramorig$fngenes*paramorig$fmpg))



  dispersed

  return(dispersed)
  })


#calculate the reproduction by dispersed offspring
taxa.pa.rep.disp<-map(1:nrow(taxa.params),~{
  #this does not work for both sexes dispersing yet
  paramorig<-taxa.params[.x,]
  reprod<-taxa.pa.rep[[.x]]
  disp<-taxa.pa.disp[[.x]]

  #first get to know reference values for reproduction depending on sex and rank
  reprod.ref<-reprod %>% group_by(sexpa,rankpa,sexoffs) %>%
    summarise(noffs=sum(noffs),.groups="keep") %>% ungroup

  #disp contains the offspring that dispersed successfully
  #now, how many grand offspring of each sex will they produce?
  disp2<-expand.grid.df(disp,data.frame(sexgrandoffs=c("m","f"),stringsAsFactors = F))
  disp2
  #calculate number of male and female grand offspring produced by emigrant offspring
  #and equivalence in number of genes
  disp.grandoffs<-disp2 %>% mutate(ngrandoffs=noffs*reprod.ref$noffs[match(paste0(sexoffs,rankoffs,sexgrandoffs),
                                                                           paste0(reprod.ref$sexpa,reprod.ref$rankpa,reprod.ref$sexoffs))],
                                   ngenesgrandoffs=case_when(sexoffs=="m"&sexgrandoffs=="m"~(ngenesoffs/noffs)*paramorig$mmpg*ngrandoffs,
                                                             sexoffs=="m"&sexgrandoffs=="f"~(ngenesoffs/noffs)*paramorig$mfpg*ngrandoffs,
                                                             sexoffs=="f"&sexgrandoffs=="f"~(ngenesoffs/noffs)*paramorig$ffpg*ngrandoffs,
                                                             sexoffs=="f"&sexgrandoffs=="m"~(ngenesoffs/noffs)*paramorig$fmpg*ngrandoffs))
  disp.grandoffs$ngenesgrandoffs[which(is.na(disp.grandoffs$ngenesgrandoffs))]<-0

  disp.grandoffs
  })


#2.2 Potential for gene spread-----------

#calculate total population gene production per reproductive bout
reprod.pop<-map_dbl(1:nrow(taxa.params),~{
  paramorig<-taxa.params[.x,]
  reprod<-taxa.pa.rep[[.x]]
  reprodunit<-reprod %>% filter(sexpa=="f")
  reprodunit$noffs[which(reprodunit$rankpa=="b")]<-reprodunit$noffs[which(reprodunit$rankpa=="b")]*(paramorig$fgs-1)
  #now translate the number of genes not to those transmitted by the mother but the total produced by the unit
  reprodunit<-reprodunit %>% mutate(ngenesoffs=case_when(sexoffs=="f"~paramorig$fngenes*noffs,
                                                         sexoffs=="m"~paramorig$mngenes*noffs))
  reprodunit<-sum(reprodunit$ngenesoffs)
  #calculate total population production
  reprodpop<-paramorig$tnbunits*reprodunit
  reprodpop
  })

#calculate
#ratio between genes in grand offspring descending from alpha grandparent/total number of genes produced in grand offspring
pspread<-map_dbl(1:nrow(taxa.params),~{
  reprod.disp<-taxa.pa.rep.disp[[.x]]
  reprodpop<-reprod.pop[.x]
  #total production by descendants of grandparent of dispersing sex
  reprod.disp<-sum(reprod.disp$ngenesgrandoffs)
  rprodrec<-reprod.disp/reprodpop
  rprodrec
})

#add this variable to the df
taxa.params$pspread<-pspread
#reorder
taxa.params<-taxa.params %>% select(taxa,fgs,mgs,farepmono,marepmono,coopbr,pspread,everything())


#3. Cladogenesis--------
#calculate the diversity of genetic lines in offspring of population

#the genetic lines will depend on which sex disperses
#all genes produced by the phylopatric sex as one line
#then what each immigrant contributed to offspring

mix<-map(1:nrow(taxa.params),~{
  paramorig<-taxa.params[.x,]
  reprod<-taxa.pa.rep[[.x]]

  #now cdivide genes produce into the different genetic lines
  #case only males disperse
  if(paramorig$mbd==1&paramorig$fbd==0){
    genes.diversity<-
      with(reprod,
           c(inma=ngenesoffs[which(sexpa=="m"&rankpa=="a")] %>% sum,#male alpha all
             inmb=rep(ngenesoffs[which(sexpa=="m"&rankpa=="b")]%>% sum,(paramorig$mgs-1)),#males beta
             infab=sum(ngenesoffs[which(sexpa=="f"&rankpa=="a")],#all females in and out as one genetic line
                     ngenesoffs[which(sexpa=="f"&rankpa=="b")]*(paramorig$fgs-1)))
           )
  }
  #only females
  if(paramorig$mbd==0&paramorig$fbd==1){
    genes.diversity<-
      with(reprod,
           c(infa=ngenesoffs[which(sexpa=="f"&rankpa=="a")] %>% sum,#female alpha all
             infb=rep(ngenesoffs[which(sexpa=="f"&rankpa=="b")]%>% sum,(paramorig$fgs-1)),#females beta
             inmab=sum(ngenesoffs[which(sexpa=="m"&rankpa=="a")],#all male production as one genetic line
                       ngenesoffs[which(sexpa=="m"&rankpa=="b")]*(paramorig$mgs-1))
             ))
  }
  #both sexes
  if(paramorig$mbd==1&paramorig$fbd==1){
    genes.diversity<-
      with(reprod,
           c(inma=ngenesoffs[which(sexpa=="m"&rankpa=="a")] %>% sum,#male alpha all
             inmb=rep(ngenesoffs[which(sexpa=="m"&rankpa=="b")]%>% sum,(paramorig$mgs-1)),#males beta
             infa=ngenesoffs[which(sexpa=="f"&rankpa=="a")] %>% sum,#male alpha all
             infb=rep(ngenesoffs[which(sexpa=="f"&rankpa=="b")]%>% sum,(paramorig$fgs-1)),#males beta
             ))

  }
  return(genes.diversity)
})

#rescale to a comparable population size
#how many genetic lines in pop
lines<-map(1:nrow(taxa.params),~{
  paramorig<-taxa.params[.x,]
  mixed<-mix[[.x]]
  lines<-rep(mixed,floor(paramorig$tnbunits))
  return(lines)
})
#how many reproducing in pop
repr.lines<-map(lines,~{
  r.lines<-.x[.x!=0&is.na(.x)==F]
  return(r.lines)
})
#and apply Shannon's diversity index
mix.hill<-map_dbl(repr.lines,~{
  #Diversity index function----
  #x is a vector with dispersing sex ids in a group and their respective abundance of genes produced
  mixed.hill.shannon<-hill_taxa(.x, q = 1, MARGIN = 1, base = exp(1))
  return(mixed.hill.shannon)
})

taxa.params$hill.shannon<-mix.hill
#standardize by population size which differ slightly because of floor(paramorig$tnbunits)
taxa.params$stnd.hill.shannon<-(mix.hill/(taxa.params$gs*floor(taxa.params$tnbunits)))*100

#Save data---------------
saveRDS(taxa.params,file=paste0(getwd(),"/objects/taxa_params",".RDS"))
save(taxa.params,mix.hill,repr.lines,lines,mix,pspread,reprod.pop,
     taxa.pa.rep.disp,taxa.pa.disp,taxa.pa.rep,pa.rep,taxa.pa.rep.prev,
     pa.rep.prev,fagainrepsupr,totalpop,maxgs,mgsp,nbaoa,foffs,varsoc,
     file=paste0(getwd(),"/objects/all_objects",".Rdata"))

#Turning the variation into factors for plotting results
taxa.params.pl<-taxa.params %>% 
  mutate(
    pspread=pspread*100,
    so=factor(case_when(mgs==0 & fgs ==1 ~"f",
                        mgs==1 & fgs ==0 ~"m",
                        mgs==1 & fgs ==1 ~"fm",
                        mgs==1 & fgs > 1 ~"mff",
                        mgs > 1 & fgs == 1 ~"fmm",
                        mgs==0 & fgs>1 ~"ff",
                        mgs>1 & fgs == 0 ~"mm",
                        mgs>1 & fgs>1 ~"ffmm",
                        TRUE~"unknown"),
              levels=c("f","ff","m","mm","fm","mff","fmm","ffmm","unknown")),
    ms1=factor(case_when(
      marepmono.fct=="total" & farepmono.fct=="total"~"monog",
      marepmono.fct!="total" & farepmono.fct!="total"~"promisc",
      marepmono.fct=="total" & farepmono.fct!="total"~"pgyn",
      marepmono.fct!="total" & farepmono.fct=="total"~"pandr",
      TRUE~"unknown"),levels=c("monog","pgyn","pandr","promisc")),
    ms2=factor(case_when(
      marepmono.fct=="total" & farepmono.fct=="total"~"monog",
      marepmono.fct!="total" & farepmono.fct!="total" & (fgs+mgs)>2~"promisc",
      marepmono.fct=="total" & farepmono.fct!="total" & fgs>1~"pgyn",
      marepmono.fct!="total" & farepmono.fct=="total" & mgs>1~"pandr",
      TRUE~"unknown"),levels=c("monog","pgyn","pandr","promisc")),
    coopbr=factor(coopbr,levels=c(0,1))
  )

saveRDS(taxa.params.pl,file=paste0(getwd(),"/objects/taxa_params_pl",".RDS"))

#END-------------
