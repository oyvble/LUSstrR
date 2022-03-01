
#' @title convert
#' @author Oyvind Bleka 
#' @description Convert sequences to bracket format
#' @details The dataframe must contain the following columns: Locus, Sequence, others are kept
#' This is the main function for converting sequences into bracket format as done by lusSTR
#' @param df A dataframe which must contain the following columns: Locus, Sequence (upper case Strings), other columns are kept
#' @param hasFlanks Whether flanks are part of sequence (should be true if format=STRaitRazor)
#' @param addFlanks Whether flanks should be added to bracket format
#' @param lusSep Separator for LUS format
#' @param panel Assumed panel (ForenSeq). PowerPlex not yet implemented
#' @param format Format of data (UAS or STRaitRazor): Assumes that STRaitRazor (and others) are already in forward direction, while UAS must revComp for some markers
#' @return A dataframe with an additional column bracket format (converted sequences)
#' Traditional_STR_Allele: Traditional CE/RU format
#' Forward_Strand_Bracketed_form: Bracket format format in forward direction
#' UAS_Output_Bracketed_Form: Bracket format in forward/backward direction wrt UAS definitions
#' LUS_Plus: Provides LUS+ format of sequences
#' @export
#' @examples
#' \dontrun{ 
#' pgkPath <- path.package("LUSstrR") # Get package path.
#' df = readRDS(paste0(pgkPath,"/examples/exampleDataframe.RDS"))
#' df = convert(df)
#' }
convert = function(df=NULL, hasFlanks=FALSE, addFlanks=FALSE, panel="ForenSeq", format="UAS", lusSep="_") {
  if(is.null(df)) return(NULL)

  if(!panel%in%c("ForenSeq")) stop("Indicated panel is not supported!")

#  if(format=="STRaitRazor" && !hasFlanks)
    
  #Import marker rules
  pgkPath <- path.package("LUSstrR", quiet = TRUE) # Get package path.
  rules_str_source = paste0(pgkPath,"/ext/str_markers.json") #get data with rules (STR)
  rules_str = rjson::fromJSON(file=rules_str_source)#,simplify = TRUE)
  names(rules_str) = toupper(names(rules_str))
  
  if(! all(c("Locus","Sequence")%in%colnames(df)) ) stop("Missing column names in data")
  if(addFlanks && !hasFlanks) stop("addFlanks require flanks to be included.") 
  
  df$Locus = toupper(df$Locus) #force upper case
  locs = unique(df$Locus) #obtain columns to evaluate
  
  #INIT NEW VARIABLES:
  df$Traditional_STR_Allele <- df$Forward_Strand_Bracketed_form <- df$UAS_Output_Bracketed_Form <- df$LUS_Plus <- df$Sequence 
  if(hasFlanks && !addFlanks) { #ADD extra columns if flanks are not included as part of bracket form
    df[['5_Flank']] <- df[['3_Flank']] <- df$Sequence
  }  
  
  for(locus in locs) {
#   locus=locs[2]
    indUse = df$Locus==locus
  #View(df[indUse,])
    if(sum(indUse)==0) next #no observations
    sequences = df$Sequence[indUse] #obtain sequences to convert

    if(locus=="PENTAE") locus="PENTA E"
    if(locus=="PENTAD") locus="PENTA D"
    #Extract rules from JSON
    rules = rules_str[[locus]]
    
    if(is.null(rules)) {
      print(paste0("Marker ",locus," was not recognized. Skipping.."))
      next
    }
    
    #Should the string be reverse compleemented? If yes, do so.
    #ENSURES TO ALWAYS BE FORWARD COMPLEMENTED
    revComp = toupper(rules$ReverseCompNeeded)=="YES"
	  if( format!="UAS" ) revComp = FALSE #SEQUENCE IS ALREADY ASSUMED FORWARD IF NOT UAS
    if(revComp) sequences = getComplement(sequences)
        
    #If flanks sequence included
    if(hasFlanks) {
      seq_len = nchar(sequences) #obtain sequence length
      sequences_before <- sequences_after <- ""
      if(rules$Foren_5>0) sequences_before = substr(sequences,1,rules$Foren_5) #5' flanks
      if(rules$Foren_3>0) sequences_after = substr(sequences, seq_len-rules$Foren_3+1,seq_len) #3' flanks
      sequences = substr(sequences,rules$Foren_5+1, seq_len - rules$Foren_3)
    }
    
    repeats = rules$Repeats #obtain repeats from Rule-table
    repeat_size = nchar(repeats[1]) #rules$NumBasesToSeparate #obtain repeats from Rule-table
    
    #Obtain Canonical STR allele designation (CE/RU)  
    new_seqs <- sequences
    if(rules$BasesToSubtract > 0) {
      #stringi::stri_sub_all("abcde",to = -(3+1)) #[:-3]
      new_seqs = stringi::stri_sub_all(new_seqs,to = -(rules$BasesToSubtract+1)) #[:(-bases)]
    }
    new_seqs_len = nchar(new_seqs) #get length of sequences
    allele_int = as.integer(new_seqs_len/repeat_size)
    allele_dec = as.integer( new_seqs_len %% repeat_size)
    canon_alleles <- as.numeric( paste0(allele_int,".",allele_dec))
    
    #Recognize if sequence can be collapsed by length (or using priority on a particular repeat/motif)
    mustSplit = c('D13S317', 'D18S51', 'DYS643', 'DYS635', 'DYS635', 'DYS612', 'DYS576', 'DYS570','DYS549', 'DYS533', 'DYS505', 'DYS481', 'DYS460', 'DYS439', 'DYS438',
                  'DYS437', 'DYS392', 'DYS391', 'DYS390', 'DYS389II', 'DYS389I', 'DYS385A-B', 'DYS19','DYF387S1', 'DYS393', 'DYS456', 'HPRTB', 'DXS8378', 'DXS7423', 'DXS10103')
    issplit_compatible = locus%in%mustSplit | nchar(sequences)%%repeat_size==0 #check if compatible to split directly
    
    
    bylength  = issplit_compatible | locus%in%"D16S539"# c("D3S1358",)
    
    #Obtaining bracket format (coverted) for forward sequence
    sequences_bracket_forward = switch(locus,
       "D13S317"= conv_D13S317(sequences),
       "D1S1656"= conv_D1S1656(sequences, repeats),
       "FGA"= conv_FGA(sequences, repeats),
       "D7S820"= conv_D7S820(sequences, repeats, canon_alleles),
       "D18S51" = conv_D18S51(sequences, repeats, canon_alleles),
       "TH01"= conv_TH01(sequences, repeats),
       "D19S433"= conv_D19S433(sequences, repeats),
       "D21S11"= conv_D21S11(sequences, repeats),
       "PENTA D"= conv_PentaD(sequences, repeats),
       conv_general(sequences,  repeat_size, repeats, bylength) #default statement
    )
    
    #Test that sequence is correctly converted back  
    #if(checkString) {
    #  sequences2 = convertBracket2seq(sequences_bracket_forward) 
    #  if(!all(sequences==sequences2)) stop("WRONG SEQUENCE OBTAINED!")
    #} 
    
    #INSERT TO NEW VARIABLES:
    df$Traditional_STR_Allele[indUse] = canon_alleles
    df$Forward_Strand_Bracketed_form[indUse] <-  df$UAS_Output_Bracketed_Form[indUse] <-  sequences_bracket_forward #insert converted sequences (forward)
    
    #Need to reverse complement bracket format if REVERSE MARKER:
    if(revComp) df$UAS_Output_Bracketed_Form[indUse] <- reverse_complement_bracketed(brackets=sequences_bracket_forward)
    
    #LAST: OBTAINING LUS+ STRUCTURE:
    #obtain motif repeat structure (list) for each sample
    motifReps = getMotifReps(sequences_bracket_forward) 
    
    #Last oBtain LUS, LUS+
    lusPlus = switch(locus,
        "D21S11"= getLUSplus_D21S11(motifReps,rules,lusSep),
        getLUSplus_general(motifReps,rules,lusSep) #default statement
    )
    
    #Special handling of LUS (post-process:
    if(locus=="PENTA D") {
      indIns =  match(c("2.2","3.2"),canon_alleles)
      notNA = !is.na(indIns)
      lusPlus[indIns[notNA]] = c("5","6")[notNA]
    } 
    
    if(locus=="D7S820") {
      addedTert = rep(0,length(lusPlus))
      lastLetter = stringi::stri_sub(sequences_bracket_forward,from=-1) #obtain last letter
      indOne = lastLetter=="T" & canon_alleles!=round(canon_alleles) #add Ter=one for these
      addedTert[indOne] = 1 
      lusPlus = paste0(lusPlus,lusSep,addedTert)
    }
    
    #df$LUS[indUse] = paste0(canon_alleles,lusSep,nrepeats[,1]) #only LUS (always)
    df$LUS_Plus[indUse] = paste0(canon_alleles,lusSep,lusPlus) #reminind
    
    #LAST: ADD BRACKET flanks if wanted
    if(hasFlanks) {
      #Convert flank sequence to bracket (always if provided)
      #sequences_before
      #sequences_after
      if(addFlanks) {
        sequence_middle = df$UAS_Output_Bracketed_Form[indUse]
        df$UAS_Output_Bracketed_Form[indUse] = paste(sequences_before,sequence_middle,sequences_after)
      } else {
        df[['5_Flank']][indUse] <- sequences_before
        df[['3_Flank']][indUse] <- sequences_after
      }
    }
    
  } #end for each marker
  return(df)
}
