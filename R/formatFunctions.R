#Function to split a sequence into two separate strings at a specified repeat unit.
split_sequence_into_two_strings = function(sequence,repeat_for_split) {
  #sequences = "ACACACACACCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA"
  #repeat_for_split='CACA'
  
  #Note Find where GGAA repeat's ar broke 
  locates = stringi::stri_locate_all(sequence,fixed=repeat_for_split)[[1]] #locate positions of motif (possibly repeated)
  
  #obtain last index of repeated motif (prev)
  prev <- 0
  if(!is.na(locates[1,1])) {
    indlocate <- 1
    if(nrow(locates)>1) {
      indlocate = which(diff(locates[,2])!=nchar(repeat_for_split)) 
      if(length(indlocate)==0) indlocate = nrow(locates)
    } 
    prev = locates[indlocate[1],2]
  }
  if(prev != 0) { #if motif was found
    first_string =  stringi::stri_sub(sequence,to = prev)
    second_string = stringi::stri_sub(sequence,from = prev+1) #notice the added index
  } else {
    first_string = ""
    second_string = sequence
  }
  return(c(first_string,second_string))
}


#Function converting block counts to bracket format
getBracketForm = function(block_counts, brksign=c("[","]")) {
  brack_form = rep(NA,length(block_counts)) #sequence string in bracket format
  for(s in seq_along(block_counts)) {
    counts = block_counts[[s]]
    bracks = names(counts) #obtain motifs
    indIns = counts>1
    bracks[indIns] = paste0(brksign[1],bracks[indIns],brksign[2])
    bracks[indIns] = paste0(bracks[indIns],counts[indIns]) #insert count numbers
    brack_form[s] = paste0(bracks, collapse=" ") #collapse to vector
  }  
  return(brack_form)
}


getBlockCounts = function(splittedList) {
  block_count = list()
  for(s in seq_along(splittedList)) {
    block_count[[s]] = as.integer()
    nBlocks = length(splittedList[[s]]) #number of blocks
    
    prevBlock = splittedList[[s]][1] 
    counter = 1
    if(nBlocks>1) {
      for(m in 2:nBlocks) { #for each motif
        motif = splittedList[[s]][m]
        isequal = motif==prevBlock
        if(isequal) {
          counter = counter + 1
        } else {
          block_count[[s]] = append(block_count[[s]], setNames(counter,prevBlock))
          prevBlock = motif #update what was previous block
          counter = 1 #reset counter (counted ones)
        }
      }
      if(!isequal) counter =  1 #reset counter
    } #LAST: ADD LAST BLOCK
    block_count[[s]] = append(block_count[[s]], setNames(counter,prevBlock))
  }
  return(block_count)
}

#Split a sequence into chunks of length n, and count adjacent repeated chunks.
get_blocks = function(sequences,n, rev=FALSE) {
  splittedList = split_by_n(sequences,n,rev)
  block_count = getBlockCounts(splittedList) #get counted blocks for each type
  return(block_count)
}



#Function to split string based on repeat
#getSplittedList(sequences="AATCTGTC",repeat_for_split="TGTC")

getSplittedList = function(sequences, repeat_for_split) {
  splitSym = " "
  #splittedList = strsplit(sequences,repeat_for_split,fixed=TRUE) #has problem with example above
  splittedList = list()
  for(s in seq_along(sequences)) {
    #splitted = splittedList[[s]]
    splitted = stringi::stri_split_fixed(sequences[s],repeat_for_split,simplify = TRUE)[1,] #one splitted string at the time
    splitted = paste0(splitted,collapse=paste0(splitSym,repeat_for_split,splitSym))  #collapse wrt ws
    splitted = strsplit(splitted," +")[[1]] #gsub("  "," ",split) #use only one whitespace
    splitted = splitted[splitted!=""]
    #if(splitted[length(splitted)]==repeat_for_split) splitted = append(splitted, repeat_for_split) #append last motif if using strsplit()
    splittedList[[s]] = splitted
  }
  return( splittedList )
}

#Convert to bracketed annotation form by splitting the sequence into blocks of size n.
collapse_repeats_by_length = function(sequences, n) {
  # sequence="ACACA";n=4
  return( getBracketForm( get_blocks(sequences, n, FALSE) ) )
 
}


# Collapse tandem stretches of the specified repeat sequence in a larger sequence.
# collapse_tandem_repeat('TAGATTATTATTTAGTAGATTTAGTAG', 'ATT') =='TAG [ATT]3 TAGTAG ATT TAGTAG'
# collapse_tandem_repeat('TAGATTATTATTTAGTAGATTTAGTAG', 'TAG') == 'TAG ATTATTATT [TAG]2 ATT [TAG]2'
collapse_tandem_repeat = function(sequence, repeat_for_split) {
#sequences='TAGATTATTATTTAGTAGATTTAGTAG'
#repeat_for_split = 'ATT'
  splittedList = getSplittedList(sequence,repeat_for_split)  #returing only 1 sequence
  block_counts = getBlockCounts(splittedList)
  return( getBracketForm(block_counts) ) 
}


#Convert a sequence to bracketed form by collapsing stretches of tandem repeats.
# getBracketForm(getBlockCounts(collapse_all_repeats('TAGATTATTATTTAGTAGATTTAGTAG', 3, c('ATT', 'TAG')))) =='TAG [ATT]3 [TAG]2 ATT [TAG]2'
# getBracketForm(getBlockCounts(collapse_all_repeats('AACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATC', 4, c("TATC","TGTC"))))

#Modified function to deal with combination of repeat and length
collapse_all_repeats = function(sequences, repeat_size, repeats) {
  collapsedList = list() #init with full seq
  for(s in seq_along(sequences)) { #for each sequence
    sequenceNew = sequences[s]
    for(r in seq_along(repeats)) { #for each repeats
      rep = repeats[r]
      nMotifSize = nchar(sequenceNew)
      traverseMore = nMotifSize > repeat_size #split more
      sequence2 = NULL
      for(t in seq_along(traverseMore)) {
        seq = sequenceNew[t] #sequence to add
        if(traverseMore[t]) { #need to split sequence further
          seq = getSplittedList(seq,rep)[[1]] #use this function directly
        }
        sequence2 = append(sequence2, seq)
      }
      sequenceNew = sequence2 #override with existing
    } #end for each repeats
    
    #Check again if long blocks must be splitted (split by repeat_size)
    nMotifSize = nchar(sequenceNew)
    splitMore = nMotifSize > repeat_size #split more
    if(any(splitMore)) {
      sequence2 = NULL
      for(r in seq_along(sequenceNew)) {
        seq = sequenceNew[r] #sequence to add
        if(splitMore[r]) { #must split sequence more
          seq = split_by_n(seq, repeat_size)[[1]] #use this function directly
        }
        sequence2 = append(sequence2, seq)
      }
      sequenceNew = sequence2 #override with existing
    }
    collapsedList[[s]] = sequenceNew #updated
  }
  return( collapsedList )
} 
    

#  Uses a combination of repeat-based and length-based methods to convert a sequence containing tandem repeats into a concise bracketed representation.
sequence_to_bracketed_form = function(sequences, repeat_size, repeats) {
    collapsedList = collapse_all_repeats(sequences, repeat_size, repeats)
    
    #paste0(collapsedList[[1]],collapse="")==sequences[1]
    return( getBracketForm( getBlockCounts(collapsedList) ) )
}




#Helpfunction to obtain number of motif repeats for each type in sequence
getMotifReps = function(bf, brksign=c("\\[","\\]")) {
  nreps = rep( 0,length(bf) )
  splittedList = strsplit(bf," ")
  nMotifReps = list()
  for(s in seq_along(splittedList)) {
    seq = splittedList[[s]] #look on specific bracket format sequence
    hasReps = grepl(brksign[1],seq ) #obtain which motifs which has repeats (must obtain motif and rep-number)
    seq = gsub(brksign[1],"",seq) #remove left bracket (Ex: "AGAT]3" "ATAT"   "AGAT]7")

    #Store in new format
    nmotifReps = rep(1,length(seq))
    names(nmotifReps)[!hasReps] = seq[!hasReps]

    if(any(hasReps)) {
      #obtain number of repeats per motif (stored in matrix)
      nmotifMat = matrix(unlist(strsplit(seq[hasReps],brksign[2])),nrow=2)
    
      nmotifReps[hasReps] = as.integer( nmotifMat[2,])
      names(nmotifReps)[hasReps] = nmotifMat[1,]
    }    
    nMotifReps[[s]] = nmotifReps
  }
  return(nMotifReps)
}
  


#Split a sequence into non-overlapping chunks of length n.
#split_by_n("ATGATG",4)[[1]]
split_by_n = function(sequences, n, rev=FALSE) {
  strsplit(sequences, paste0("(?<=.{",n,"})"), perl = TRUE)
}
convertBracket2seq = function(seqs, brksign=c("\\[","\\]")) {
  # All sequences are split with pattern " ". Output is a list of lists.Every sequence
  # is a list of tetranucleotides!!!
  blockSplitted = strsplit(seqs," ")
  fullseq = rep("",length(blockSplitted)) #full sequence to return
  for(s in seq_along(blockSplitted))  {
    #s=1
    seqSplitted = blockSplitted[[s]]
    for(b in seq_along(seqSplitted)) {
      #b=1
      block = seqSplitted[b]
      repMotif = strsplit(block,brksign[2])[[1]] #obtain repeat motif
      repMotif[1] = gsub(brksign[1],"",repMotif[1]) #remove possible [
      if(length(repMotif)==2) repMotif = paste0(rep(repMotif[1],as.integer(repMotif[2])),collapse="")
      fullseq[s] = paste0(fullseq[s],repMotif[1]) #add sequence motif
    }
  }
  return(fullseq)
}

#Helpfunction to get reverse complement of sequence
getComplement = function(seqs) {
  if(!all(is.character(seqs))) stop("Please only input character vectors")
  seqs = toupper(seqs) #make sure that letters are upper case
  fromLetters = c("A","T","C","G")
  toLetters = c("t","a","g","c")
  for(i in seq_along(fromLetters))  seqs = gsub(fromLetters[i],toLetters[i],seqs)
  return( stringi::stri_reverse(toupper(seqs)) )
}

#reverse_complement_bracketed("[ATCT]7 ACCT [ATCT]3")=="[AGAT]3 AGGT [AGAT]7"
reverse_complement_bracketed = function(brackets) {
  fromLetters = c("A","T","C","G")
  toLetters = c("t","a","g","c")
  for(i in seq_along(fromLetters))  brackets = gsub(fromLetters[i],toLetters[i],brackets)
  
  splittedList = strsplit(brackets, " ")
  for(s in seq_along(splittedList)) {
    tmp = strsplit(splittedList[[s]],"")
    for(i in seq_along(tmp)) {
      isLetter = tmp[[i]]%in%toLetters
      tmp[[i]][isLetter] = toupper(rev(tmp[[i]][isLetter])) #reverse these and make upperCase
      tmp[[i]] = paste0(tmp[[i]],collapse="")
    }
    brackets[s] = paste0(rev(unlist(tmp)),collapse=" ")
  }
  return(brackets)
}
