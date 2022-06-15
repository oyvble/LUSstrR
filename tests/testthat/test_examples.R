#rm(list=ls())
#using the Lookup file to test the implementation 
#library(LUSstrR);library(testthat);sessionInfo()

test_that("check that FGA convertion is correct", {
 df = data.frame(Locus="FGA",Sequence="CCAGCAAAAAAGAAAGGAAGAAAGGAAGGAAGGAGAAAGAAAGAAAGAAAGAAA")
 ret = convert(df, hasFlanks=TRUE, format="STRaitRazor")  
 expect_equal(ret$Forward_Strand_Bracketed_form,"[GGAA]2 GGAG [AAAG]3 A AA GAAA")
 
 df = data.frame(Locus="FGA",Sequence="TTTCTTTCTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTCCTTCCTTCC")
 ret = convert(df, hasFlanks=FALSE, format="UAS")  
 expect_equal(ret$UAS_Output_Bracketed_Form,"[TTTC]3 TTTT TT [CTTT]9 CTCC [TTCC]2")
})


test_that("check D21 bracket format conversion", {
  seq = "TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGATATCTA"
  df = data.frame(Locus="D21S11",Sequence=seq)
  ret = convert(df, hasFlanks=FALSE, format="UAS")
  seq2 = getComplement(convertBracket2seq(ret$Forward_Strand_Bracketed_form))
  expect_equal(seq,seq2)
  expect_equal(ret$Forward_Strand_Bracketed_form,"[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]11 TCGA TA TCTA")
})


test_that("check D7S820 bracket format conversion", {
  #seq = getComplement(convertBracket2seq("AAAAC TATC AATC TGTC [TATC]8"))
  seq = "GATAGATAGATAGATAGATAGATAGATAGATAGACAGATTGATAGTTTT"
  df = data.frame(Locus="D7S820",Sequence=seq)
  ret = convert(df, hasFlanks=FALSE, format="UAS")
  expect_equal(ret$Forward_Strand_Bracketed_form,"A AAAC TATC AATC TGTC [TATC]8")
})

test_that("check TH01 bracket format conversion (note slight difference to lusSTR)", {
  #split_by_n("ATGATG", 4)[[1]]
  #seq = convertBracket2seq("[AATG]5 ATG ATG [AATG]3")
  seq = "AATGAATGAATGAATGAATGATGATGAATGAATGAATG"
  df = data.frame(Locus="TH01",Sequence=seq)
  ret = convert(df, hasFlanks=FALSE, format="UAS")
  expect_equal(ret$Forward_Strand_Bracketed_form,"[AATG]5 ATGA TG [AATG]3")
})

test_that("check D21 LUS+ conversion", {
  brack = "[TCTA]2 ACTA TCTA [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]11"
  seq = convertBracket2seq(brack)
  df = data.frame(Locus="D21S11",Sequence=seq)
  ret = convert(df, hasFlanks=FALSE, format="UAS")
  expect_equal(ret$Forward_Strand_Bracketed_form,brack)
  expect_equal(ret$LUS_Plus,"29_11_2_6") 

  brack = "TCTA CCTA [TCTA]2 [TCTG]7 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]10"
  seq = convertBracket2seq(brack)
  df = data.frame(Locus="D21S11",Sequence=seq)
  ret = convert(df, hasFlanks=FALSE, format="UAS")
  expect_equal(ret$Forward_Strand_Bracketed_form,brack)
  expect_equal(ret$LUS_Plus,"29_10_2_7") 
})


test_that("check PentaD LUS+ conversion", {
  
  bracks = c( "GA [AAAGA]2 AACGA","GA [AAAGA]3","GA AAAGA ACAGA AAAGA",
              "GAAAA TA [AAAGA]2","GAAAA GC [AAAGA]2")
  
  for(brack in bracks) {
#    brack=bracks[1]
    df = data.frame(Locus="Penta D",Sequence=convertBracket2seq(brack))
    ret = convert(df, hasFlanks=FALSE, format="UAS")
    expect_equal(ret$Forward_Strand_Bracketed_form,brack)
    expect_equal(ret$LUS_Plus,"2.2_5") #
  }

  seq = "GATCACTTGAGCCTGGAAGGTCGAAGCTGAAGTGAGCCATGATCACACCACTACACTCCAGCCTAGGTGACAGAGCAAGACACCATCTCAAGAAAGAAAAGAAAAGAAAAGAAAAGAAAAAACGAA"
  df = data.frame(Locus="Penta D",Sequence=seq)
  ret = convert(df, hasFlanks=TRUE, format="STRaitRazor")
  ret$LUS_Plus
  
  seq = "GATCACTTGAGCCTGGAAGGTCGAAGCTGAAGTGAGCCATGATCACACCACTACACTCCAGCCTAGGTGACAGAGCAAGACACCATCTCAAGAAAGAAAAGAAAAGAAAAGAAAAGAAAAAACGAA"
  df = data.frame(Locus="Penta D",Sequence=seq)
  ret = convert(df, hasFlanks=TRUE, format="UAS")
  ret$LUS_Plus
    
})



if(0) {
  seq="TTTCTTTCTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTCCTTCCTTCC"
  #seq = LUSstrR::getComplement(seq)
  df = data.frame(Locus="FGA",Sequence=seq)

}

