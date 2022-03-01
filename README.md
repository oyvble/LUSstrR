## ABOUT
LUSstrR is a reimplementation of the lusSTR tool (Python):
- https://github.com/bioforensics/lusSTR

Current implementation is only supporting aSTR sequences typed with ForenSeq. \
Required sequence format can be either UAS or STRraitRazor.

Main function to use is convert(), which takes a dataframe which must include following column names:
- Sequence (a string with upper case letters)
- Locus (not case sensitive)
- Other columns are retained

Function convertFromUASreport can be used to import data from a UAS report file directly

## INSTALLATION
From GitHub (requires R-package devtools installed)
- devtools::install_github("oyvble/LUSstrR")

When released on CRAN:
- install.packages('LUSstrR')

## USAGE
- pgkPath = path.package("LUSstrR") # Get package path
- df = readRDS(paste0(pgkPath,"/examples/exampleDataframe.RDS")) #read a dataframe
- df = convert(df) #update with converted sequences

Output columns: \
- Traditional_STR_Allele: Traditional CE/RU format
- Forward_Strand_Bracketed_form: Bracket format format in forward direction
- UAS_Output_Bracketed_Form: Bracket format in forward/backward direction wrt UAS definitions
- LUS_Plus: The LUS+ format of sequences (can be used in EuroForMix)

## ACKNOWLEDGMENT
- Rebecca Just (nomenclature)
- Rebecca Mitchell (lusSTR developer)
- Daniel Standage (lusSTR developer)

