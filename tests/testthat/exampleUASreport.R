
#using the Lookup file to test the implementation 
library(LUSstrR)

pgkPath <- path.package("LUSstrR", quiet = FALSE) # Get package path.
file = paste0(pgkPath,"/examples/UAS_Sample_Details_Report_test.xlsx")

#No flanks:
df = convertFromUASreport(file)
df$UAS_Output_Bracketed_Form

