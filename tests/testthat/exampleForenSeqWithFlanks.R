#rm(list=ls())
#using the Lookup file to test the implementation 
library(LUSstrR)

pgkPath <- path.package("LUSstrR", quiet = FALSE) # Get package path.
file = paste0(pgkPath,"/examples/UASwithFlanksExample.csv")

#Obtain data with ForenSeq+flanks (UAS or STRait_Razor)
df = read.csv2(file) #read.table(file,header = T,sep="\t")
#View(df)

df1 = convert(df,hasFlanks = TRUE,addFlanks = FALSE)
df2 = convert(df,hasFlanks = TRUE,addFlanks = TRUE)

#flank3 = convertBracket2seq(df1$X3_Flank_Bracketed_Notation)
#flank5 = convertBracket2seq(df1$X5_Flank_Bracketed_Notation)

df1$UAS_Output_Bracketed_Form

df1$'3_Flank'
df1$'5_Flank'