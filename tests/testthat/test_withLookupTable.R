
#using the Lookup file to test the implementation 
#rm(list=ls())
#library(LUSstrR)

pgkPath <- path.package("LUSstrR", quiet = FALSE) # Get package path.
file = paste0(pgkPath,"/ext/LookupTable_081220.xlsx")

sheets = readxl::excel_sheets(file)

df = NULL
for(sheet in sheets) {
  #sheet = sheets[1]
  suppressMessages( dat <- readxl::read_excel(file, sheet = sheet, col_names = TRUE))
  
  df_new = data.frame(Locus=toupper(sheet), Sequence=dat$`UAS SEQUENCE`, ConvertedUAS=dat$`UAS REGION BRACKETED FORM`, lusUAS= dat$`LUS+ ALLELE`)
  df = rbind(df,df_new)
}

df = convert(df) #convert

locs = unique(df$Locus)
for(locus in locs) {
#locus=locs[22]
  sub = subset(df,df$Locus==locus)
#View(subset(sub,select=(-Sequence)))
  
  if(locus=="D21S11") sub$ConvertedUAS = gsub("TCCATA","TCCA TA", sub$ConvertedUAS)
  if(locus=="D7S820") sub$ConvertedUAS = gsub("GTTT T","GTTTT", sub$ConvertedUAS)
  
  test_that(paste0("Bracket format at ",locus), {
    expect_equal(sub$ConvertedUAS,sub$UAS_Output_Bracketed_Form)
  })
  test_that(paste0("LUS+ format at ",locus), {
    expect_equal(sub$LUS_Plus,sub$lusUAS)
  })
}

  #MANUAL CODE
if(0) {
  indDiff = sub$ConvertedUAS!=sub$UAS_Output_Bracketed_Form
  seqs = sub$Sequence[indDiff]
  View(sub[indDiff,])
  #indDiff = sub$lusUAS!=sub$LUS_Plus
  #if(sum(indDiff)>0) {
    #View(subset(sub,indDiff,select=(-Sequence))) 
    #stop(paste0("LusDiff=",which(indDiff),collapse="/" ))
  #  diffList = rbind(diffList, subset(sub,indDiff))
} 
  #if(length(diffList)==0) stop("No differences observed!")
#View(subset( diffList,select=(-Sequence)))