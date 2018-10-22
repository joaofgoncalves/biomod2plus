
# An damn ugly hack to avoid notes in R CMD CHECK 
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".","algoName","avgResp","cnames","evalScoreAVG","evalScoreSTE",
                                                        "modAlgo","stdErr","varImpAVG","varImpSTE","varNames","targetVar",
                                                        "modEvalScore"))

