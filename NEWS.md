
=============================================
LUSstrR v0.2.1 (Release date: 2022-04-20)
=============================================
- Fixed bug in conv_FGA:L90 (markerFunctions.R):
 -- Use of strsplit instead of stri_sub caused missing '' in rare situations. 
 -- Also had to modify count variable to adjust for this change (code is now more similar to original).

=============================================
LUSstrR v0.2.0 (Release date: 2022-03-01)
=============================================
- Init release on GitHub

=============================================
LUSstrR v0.2 (Release date: 2022-01-15)
=============================================
- Sequences with flanks (from UAS) can now also be included. Be sure that convert function is called with argument "hasFlanks=TRUE". 
It's optional of whether to include the flank sequence as part of the allele or not (argument addFlanks).

=============================================
LUSstrR v0.1 (Release date: 2021-11-01)
=============================================
- Completed implementation of ForenSeq UAS format (non-flanks).