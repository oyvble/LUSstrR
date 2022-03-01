#' @title convertFromUASreport
#' @author Oyvind Bleka 
#' @description Function which converts sequences to bracket format for UAS report file input
#' @details This function is a wrapper of the convert function
#' @param file File name of report
#' @param AT A coverage threshold when extracting alleles from the report 
#' @param typedOnly Whether only typed alleles should be converted
#' @param hasFlanks Whether flanks are part of sequence
#' @return A dataframe with an additional column bracket format (converted sequences)
#' @export
#' @examples
#' \dontrun{ 
#' pgkPath <- path.package("LUSstrR", quiet = TRUE) # Get package path.
#' file = paste0(pgkPath,"/examples/UAS_Sample_Details_Report_test.xlsx")
#' df = convertFromUASreport(file)
#' }

convertFromUASreport = function(file, AT=11, typedOnly=FALSE, hasFlanks=FALSE) {
  sheet0="Autosomal STRs" #this is UAS report sheet to extract data from
  suppressMessages(
    dat <- readxl::read_excel(file, sheet = sheet0, col_names = FALSE, progress=FALSE) #no header
  )
    
  headrow = which(dat[,1]=="Coverage Information")+1 #recognize header
  header = dat[headrow,] #extract header
  dat = dat[-seq_len(headrow),,drop=FALSE]
  colnames(dat) = header
  
  #TRAVERSE EACH LOCUS AND EXTRACT RELEVANT SEQUENCES
  locs = unique(dat$Locus)
  addws = "PENTA" #adding whitespace after this string (PentaD/PentaE)
  df = NULL
  for(loc in locs) {
    # loc = locs[27]
    loc2 = toupper(loc)
    if( grepl("AM",loc2) ) next #skip AMEL
    if( grepl(addws,loc2) ) loc2 = gsub(addws,paste0(addws," "),loc2) 

    #Obtain AT for marker and Update table regarding AT
    AT0 = AT[loc2]
    if(is.na(AT0)) AT0 = AT
    sub = subset(dat,dat$Locus==loc & as.numeric(dat$Reads)>=AT0) #restrict the coverages here!
    
    if(typedOnly) {
      coluse = grep("typed",tolower(colnames(sub)))
      if(length(coluse)==0) {
        print("Typed column not found. Ignoring!")
      } else if(length(coluse)==1) {
        sub = subset(sub,tolower( sub[[coluse]] )=="yes")
      } else {
        print("Multiple typed column found. Ignoring!")
      }
    }
    
    df_new = data.frame(Locus=loc2,Reads=sub$Reads,Sequence = sub$`Repeat Sequence`)
    df = rbind(df, df_new)
  }
  df = convert(df,format="UAS",hasFlanks=hasFlanks)
  return(df) #returning dataframe
}

  
  
