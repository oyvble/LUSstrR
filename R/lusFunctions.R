
#General script to extract LUS_plus format
getLUSplus_general = function(motifReps,rules, lusSep="_") {
  lusTypes = c("LUS","Sec","Tert") #find longest motif rep for each type
  motifTypes = unlist(rules[lusTypes]) #obtain motif to extract LUS with
  nrepeats = NULL #matrix(NA, setNames(rep(NA,3),lusTypes)
  for(t in seq_along(motifTypes)) {
    motif = motifTypes[t] #get motif to look for
    if(motif=="") next 
    rank = sum(motifTypes[seq_len(t)]%in%motif) #obtain which ranked order to look for
    
    getLUS = function(vec) {
      checkInd = names(vec)%in%motif 
      if(!any(checkInd)) return(0) #return 0 if not fund
      sorted = sort(vec[checkInd],decreasing=TRUE)
      lus = sorted[rank] #extract correct rank
      ifelse(is.na(lus), 0, lus)
    } 
    lus = sapply(motifReps, getLUS)
    nrepeats = cbind(nrepeats,lus)
  }  
  colnames(nrepeats) = lusTypes[seq_len(ncol(nrepeats))]
  lusPlus = apply(nrepeats,1,function(x) paste0(x,collapse=lusSep)) #obtain LUS+ format
  return(lusPlus)
}

#From https://github.com/bioforensics/lusSTR/blob/master/lusSTR/marker.py:
#Special handling is required because the LUS repeat motif is the last 'TCTA' repeat set 
#and the secondary repeat motif is the first set of 'TCTA' repeats in the sequence.
getLUSplus_D21S11 = function(motifReps,rules, lusSep="_") {
  lusTypes = c("LUS","Sec","Tert") #find longest motif rep for each type
  motifTypes = unlist(rules[lusTypes]) #obtain motif to extract LUS with
  
  lusPlus = rep(NA,length(motifReps))
  for(s in seq_along(motifReps)) {
    motifs = motifReps[[s]] #obtain motifs
    TCATs = names(motifs)%in%motifTypes[1]
    midSearch = ceiling(sum(TCATs)/2) #split index of searching TCAT motif
    rangeFirst = head(which(TCATs),midSearch) #search range of SEC (first 2 blocks)
    rangeLast = setdiff(which(TCATs),rangeFirst)  #search range of LUS (last 2 blocks)
    rangeLast = tail(rangeLast,midSearch)
    #SEC = max(motifs[which(TCATs)[range]]) #First part is SEC
    SECind = rangeFirst[ motifs[rangeFirst]>=2 ] #require at least 2 repeats
    SEC = motifs[SECind[1]] #First repeated part is SEC
    LUS = max(motifs[rangeLast]) #Last part is LUS
    TER =  max(motifs[ names(motifs)%in%motifTypes[3]]) #Last type is other motif
    lusPlus[s] = paste0(LUS,lusSep,SEC,lusSep,TER) 
  }
  return(lusPlus)
}