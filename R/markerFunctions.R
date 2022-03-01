
#General function to convert sequences
conv_general = function(sequences, repeat_size, repeats, bylength) {
  if( length(bylength)!=length(sequences) ) stop("bylength argument must have same length as sequences!")
  #if(any(!bylength)) stop("asdasd")
  collapseseq = sequences #copy vector
  if(any(bylength)) collapseseq[bylength] = collapse_repeats_by_length(sequences[bylength], repeat_size)
  if(any(!bylength)) collapseseq[!bylength] = sequence_to_bracketed_form(sequences[!bylength], repeat_size, repeats)
  return(collapseseq)
}


conv_D1S1656 = function(sequences, repeats, repeat_size=4) {
  finals =  stringi::stri_sub(sequences,to = 2) #obtain two first letters
  sequence_filt = stringi::stri_sub(sequences,from=3) #split the remaining
  
  collapseseq = sequences
  #Traverse each sequence
  for(s in seq_along(sequences)) {
    final = finals[s] 
    
    splitted = split_sequence_into_two_strings(sequence=sequence_filt[s], repeat_for_split='CACA')
    
    first_string = splitted[1]
    second_string  = splitted[2]
    
    if(first_string == "") {
      final = append(final,'CACA')        
    } else {
      final = append(final, collapse_repeats_by_length(first_string, repeat_size))
    }
    
    if (nchar(second_string) %% repeat_size != 0) {
  #    tmp = sequence_to_bracketed_form(second_string, 4, repeats)
      #convertBracket2seq(tmp)==second_string
      final = append(final, sequence_to_bracketed_form(second_string, repeat_size, repeats))
    } else {
      final = append(final, collapse_repeats_by_length(second_string, repeat_size))
    }
    collapseseq[s] = paste0(final,collapse=" ") #insert final
  }
  return(collapseseq)  
}


conv_FGA = function(sequences, repeats, repeat_size=4) {
  isOK = nchar(sequences)%%repeat_size==0 | !grepl('GGAA',sequences)

  collapseseq = sequences
  if(any(isOK))  collapseseq[isOK] = collapse_repeats_by_length(sequences[isOK], repeat_size)
  if(all(isOK)) return(collapseseq) #return if done with all sequences

  for(s in which(!isOK)) {
# seq=sequences[!isOK][1]
    seq = sequences[s]
    
    #Note Find where GGAA repeat's ar broke 
    locates = stringi::stri_locate_all(seq,fixed='GGAA')[[1]] #locate positions of motif (possibly repeated)

    #obtain last index of repeated motif (prev)
    indlocate = 1
    if(nrow(locates)>1) {
      indlocate = which(diff(locates[,2])!=repeat_size) 
      if(length(indlocate)==0) indlocate = nrow(locates)
    } 
    prev = locates[indlocate[1],2]
    
    first_string =  stringi::stri_sub(seq,to = prev)
    second_string = stringi::stri_sub(seq,from = prev+1) #notice the added index
    
    locate2 = stringi::stri_locate_first(second_string,fixed='AAAA') #locate positions of motif
    checkString = stringi::stri_sub(second_string,from=locate2[1], to = locate2[1]+5) #obtain string to check

    #Performing special checks:
    if(checkString=='AAAAAA') {
      third_string = stringi::stri_sub(second_string, to=locate2[1]+1) #[:prev+2]
      fourth_string = stringi::stri_sub(second_string, from=locate2[1]+2) #[prev+2:]
    } else if(locate2[1]==1) { #if  prev == 0 (first position)
      third_string = stringi::stri_sub(second_string,to=nchar(second_string)-6) #remove last 6 letters
      fourth_string = stringi::stri_sub(second_string,from=nchar(second_string)-5) #keep last 6 letters      
    } else {
      third_string = stringi::stri_sub(second_string, to=locate2[1]-1) #[:prev]
      fourth_string = stringi::stri_sub(second_string, from=locate2[1]) #[prev:]
    }
    
    final = collapse_repeats_by_length(first_string, repeat_size)
    final = append(final, sequence_to_bracketed_form(third_string, repeat_size, repeats))
    
    parts = strsplit( fourth_string, 'GAAA')[[1]]
    tmp = NULL
    count = 0
    for(i in parts) {
      if(i == '') {
        count = count + 1
      }  else {
        if(count==1) {
          tmp = append(tmp,'GAAA')
        } else if( count>=2 ) {
          tmp = append(tmp, paste0('[GAAA]',count)) 
        }
        count = 1
    
        if(i == 'AAAAAA') {
          tmp = append(tmp, 'AA AAAA')
        } else if( nchar(i)>repeat_size ) {
          splitted = split_by_n(i, repeat_size )[[1]]
          tmp = append(tmp,splitted)
        } else  {
          tmp = append(tmp,i)
        }
      }
    }
    if(parts[length(parts)] == '') { #if last part is empty
     if(count>1) { #note slight modification here (-1)
       tmp = append(tmp, paste0('[GAAA]',count)) #note slight modification here (-1)
     } else {
       tmp = append(tmp, 'GAAA')  
     }
    }
    #convertBracket2seq( paste0(tmp,collapse=" "))==fourth_string
    collapseseq[s] = paste0(c(final,tmp),collapse=" ")
  } #end for each sequence
  return(collapseseq)
}

#conv_D7S820("AAACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCT",c("TATC","TGTC"),9.1)=="AAAC TATC AATC TGTC [TATC]9 T"
conv_D7S820 = function(sequences, repeats, canonicals, repeat_size=4) {
  isOK = round(canonicals)==canonicals
  
  collapseseq = sequences
  if(any(isOK))  collapseseq[isOK] = sequence_to_bracketed_form(sequences[isOK], repeat_size, repeats)
  if(all(isOK)) return(collapseseq) #return if done with all sequences
  
  for(s in which(!isOK)) {
    sequence = sequences[s]
    dec = round( (canonicals[s]%%1)*10 ) #obtain decimal number (whole integer)
    
    seqLen = nchar(sequence) #obtain sequence length
    if( dec==1 ) { 
      lastletter = stringi::stri_sub(sequence,from=seqLen) #obtaining last letter
      if(lastletter == "T") {
        forward_strand_brack_form = sequence_to_bracketed_form(sequence, repeat_size, repeats) 
      } else {
        firstletter = stringi::stri_sub(sequence,to=1) #obtain first letter
        sequence2 = stringi::stri_sub(sequence,from=2) #skip first letter
        bf = sequence_to_bracketed_form(sequence2, repeat_size, repeats)
        forward_strand_brack_form = paste0(firstletter,"", bf) #White space here??
      }
    } else if( dec==2 ) {
      new_repeat_list = c('TATC', 'TGTC', 'AATC')
      forward_strand_brack_form = sequence_to_bracketed_form( sequence, repeat_size, new_repeat_list )
    }  else {
      firstletters = stringi::stri_sub(sequence,to=3) #obtain first letters (3)
      sequence2 = stringi::stri_sub(sequence,from=4) #skip first letters (3)
      bf = sequence_to_bracketed_form( sequence2 , repeat_size, repeats)
      forward_strand_brack_form = paste0(firstletters," ", bf)
    }
    collapseseq[s] = forward_strand_brack_form #insert final
  }
  return(collapseseq)
}

conv_D18S51 = function(sequences, repeats, canonicals, repeat_size=4) {
  isOK = round(canonicals)==canonicals #check which are integers 
  collapseseq = sequences
  if(any(isOK))  collapseseq[isOK] = collapse_repeats_by_length(sequences[isOK], repeat_size) #In case Of Int
  if(all(isOK)) return(collapseseq) #return if done with all sequences
  collapseseq[!isOK] =  sequence_to_bracketed_form(sequences[!isOK], repeat_size, repeats)  #In case Of Int.Dec
  return(collapseseq)
}


#conv_TH01("AATGAATGAATGAATGAATGATGTTAATGAATGAATG","AATG")=="[AATG]5 ATG TT [AATG]3"
conv_TH01 = function(sequences, repeats, repeat_size=4) {
  
  #Obtain number of motifs per sequence:
  blockCounts = getBlockCounts( collapse_all_repeats(sequences, repeat_size, repeats))
  
  #check if motif length exceeds 4:
  isOK = sapply(blockCounts, function(x) all(nchar(names(x))<=repeat_size ))
  
  #Update block counts for those not proper
  for(s in which(!isOK)) {
    blockCount = blockCounts[[s]] #obtain number of blocks of different types
    blockCount_new = NULL
    for(i in seq_along(blockCount)) {
      block = blockCount[i]
      unit = names(block) #obtain unit/motif

      if(nchar(unit)>repeat_size) { #if motif length was very large
        group1  =  stringi::stri_sub(unit,to=3) #obtain first 3 letters
        group2  =  stringi::stri_sub(unit,from=4) #obtain first 3 letters
        if(group1 == 'ATG') {
          blockCount_new = append(blockCount_new, setNames(1,group1)) #add ATG
          
          counts = getBlockCounts(split_by_n(group2, n=repeat_size))[[1]] #obtain block counts of remainder motif seq
          blockCount_new = append(blockCount_new,counts) #add ATG
        } else {
          blockCount_new = append(blockCount_new,block)
        }
      } else {
        blockCount_new = append(blockCount_new,block)
      }
    }
    blockCounts[[s]] = blockCount_new
  }
  return( getBracketForm( blockCounts) )
}

#conv_D13S317("TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTC")=="[TATC]10 [AATC]2 [ATCT]3 TTCT [GTCT]2 TTTT GGGC TGCC TA [TATC]7 [AATC]2 [ATCT]3 TTCT GTCT GTC"
conv_D13S317 = function(sequences, repeat_size=4) {
  isOK = nchar(sequences) < 110
  collapseseq = sequences
  if(any(isOK))  collapseseq[isOK] = collapse_repeats_by_length(sequences[isOK], repeat_size)
  if(all(isOK)) return(collapseseq) #return if done with all sequences
  
  for(s in which(!isOK)) {
    seq = sequences[s]
    locate1 = stringi::stri_locate_last(seq,fixed='GGGCTGCCTA') #locate positions of motif

    first_string =  stringi::stri_sub(seq,to = locate1[2])
    second_string = stringi::stri_sub(seq,from = locate1[2]+1) #notice the added index
    
    collapseseq[s] = paste(collapse_repeats_by_length(first_string, repeat_size), 
      collapse_repeats_by_length(second_string , repeat_size)) #add together
  }
  return(collapseseq)
}

#sequences="CTCTCTTTCTTCCTCTCTCCTTCCTTCCTTCCTTCCTTCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTACCTTCTTTCCTT"
#repeats=c("CCTT","CCTA","TTCT")
#sequences="CTCTCTTTCTTCCTCTCTCCTTCCTTCCTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTACCTTCTTTCCTT"
conv_D19S433 = function(sequences, repeats, repeat_size=4) {
  nSeq = length(sequences)
  finals = rep(NA,nSeq)
  for(s in seq_len(nSeq)) {
    sequence = sequences[s]
    locates = stringi::stri_locate_all(sequence,fixed='CCTT')[[1]]#locate positions of motif

    #obtain last index of LAST repeated motif (locate1[2])
    if( is.na(locates[1,1] )) {
      prev=0
    } else {
      indlocate = 1
      if(nrow(locates)>1) {
        indlocate = which(diff(locates[,2])==repeat_size)  
        if(length(indlocate)==0) indlocate = nrow(locates)
      } 
      lastIndex = indlocate[length(indlocate)]
      if(nrow(locates)>lastIndex) lastIndex = lastIndex + 1 #correct to get last index
      prev = locates[lastIndex,2] #use last 
    }
    
    final = stringi::stri_sub(sequence,to=2) #obtain first 2 letters
    if(prev != 0) {
      first_string =  stringi::stri_sub(sequence,from=3, to = prev)
      second_string = stringi::stri_sub(sequence,from = prev+1) #notice the added index
      #paste0(final,first_string,second_string)==sequence
      
      if( nchar(first_string)%%4 != 0) {
        final = append(final, sequence_to_bracketed_form(first_string, repeat_size, repeats))
      } else {
        final = append(final, collapse_repeats_by_length(first_string, repeat_size))
      }
    } else {
      second_string = stringi::stri_sub(sequence,from=3)
    }
    if (second_string != "") {
      len = nchar(second_string) #obtain string length
      if( len %% 4 != 0 ) {
        if( len > 6) {
          third_string =  stringi::stri_sub(second_string, to = -7) #notice index change (to), string[:-6]
          final = append(final, collapse_repeats_by_length(third_string, 4))
        }
        fourth_string = stringi::stri_sub(second_string, from=-6, to = -5) #notice index change (to),  string[-6:-4]
        fifth_string =  stringi::stri_sub(second_string, from=-4) # string[-4:]
        final = append(final, fourth_string)
        final = append(final, fifth_string)
      } else {
        final = append(final, collapse_repeats_by_length(second_string, 4))
      }
    }
    finals[s] = paste0(final,collapse=" ")
  }   
  return(finals)
}

conv_PentaD = function(sequences, repeats, repeat_size=5) {
  isOK = nchar(sequences) < 18
  collapseseq = sequences
  if(any(isOK))  collapseseq[isOK] = sequence_to_bracketed_form(sequences[isOK], repeat_size, repeats)
  if(all(isOK)) return(collapseseq) #return if done with all sequences
  
  for(s in which(!isOK)) {
    sequence = sequences[s]
    prefix =  stringi::stri_sub(sequence, to=5) #[:5]
    suffix = stringi::stri_sub(sequence, from=6) # [5:]
    brack_form = sequence_to_bracketed_form(suffix, repeat_size, repeats)
    collapseseq[s] = paste(prefix, brack_form) #add together
  }
  return(collapseseq)
}

#conv_D21S11("CTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATATCTA",
# 4,c("TCTA","TCTG"))==
conv_D21S11 = function(sequences, repeats, repeat_size=4) {
  brack_forms = sequence_to_bracketed_form( sequences, repeat_size, repeats)
  finals = brack_forms #default is same bracket format
  for(s in seq_along(brack_forms)) {
    brack_form = brack_forms[s]
    len = nchar(brack_form)
    prev = stringi::stri_locate_last(brack_form,fixed=']')[2] #locate positions of motif
    if( prev %in%(len-1:5) ) next #skip loop

    first_string = stringi::stri_sub(brack_form, to=prev+2) #[:prev+2]
    second_string = stringi::stri_sub(brack_form, from=prev+3) #[prev+2:]
    
    #Special handling here: If last letter is ws (remove it)
    if( stringi::stri_sub(first_string,from=nchar(first_string))==" " ) {
      first_string = stringi::stri_sub(first_string,to=-2) #remove last string
    }
    second_string_final = gsub(' ', '',second_string)
    len = nchar(second_string_final) #obtain length of string
    
    final = first_string #append to this (different situations)
    if( len %% 4 == 0) {
      second_string = collapse_repeats_by_length(second_string_final, 4)
      final = c(final, second_string)
    } else if( len == 6) {
      third_string = stringi::stri_sub(second_string_final, from=-6, to = -5) #[-6:-4]
      fourth_string = stringi::stri_sub(second_string_final, from=-4) #[-4:]
      final = c(final, third_string, fourth_string)
    } else if( len %% 4 == 2) {
      third_string = stringi::stri_sub(second_string_final,to=-7) #[:-6]
      fourth_string = stringi::stri_sub(second_string_final, from=-6, to = -5)  #[-6:-4]
      last_string = stringi::stri_sub(second_string_final, from=-4) #[-4:]
      third_string_final = collapse_repeats_by_length(third_string, 4)
      final_string =  c(final, third_string_final, fourth_string, last_string)
    } else {
      third_string = collapse_repeats_by_length(second_string_final, 4)
      final = c(final, third_string)
    }  
    finals[s] = paste0(final,collapse=" ")
  }
  
  return(finals)
}

