                                        #' This function runs the simulation for the settings in the paper.
                                        #'
                                        #' @param case Case number, the possible choices are "I", "II","III", "IV".
                                        #' @param beta The value of beta, which is between 0 and 1. In the paper, beta is chosen as 0.3 and 0.4 respectively.
                                        #' @export

CLAT_examplePaper <- function( case, beta, numSim=1000 )
  {
    if(case =="I")
      caseI(beta, numSim)
    if(case=="II")
      caseII(beta, numSim)
    if(case=="III")
      caseIII(beta, numSim)
    if(case=="IV")
      caseIV(beta, numSim)
  }
