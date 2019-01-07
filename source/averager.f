C--------------------------------------------------------------
C The program performs the linear averaging of the experimental
C data assuming the uncorrelated measurements form different
C experiments.
C--------------------------------------------------------------
      Program averager
      implicit none
      include 'common.inc'
      logical doCovconversion


C
C Reading data and preparing for averaging:
C
      Call hf_ErrLog(19031301,'I:Start HAverager')      
      
C
C Perhaps we want to convert covariance matrix -> nuisance representation ?
C
      call covartonui( doCovconversion) 
      if ( .not. doCovconversion) then

         Call InitAverager

C
C Perform the linear averaging:
C
         Call Averaging

C
C Write averaged data to disk:
C

         Call Output
 
      endif     
C
C Error logging:
C

      if (IDEBUG.gt.-1) then
          call HF_errsum(6)
      endif
C--------------------------------------------------------------
      end
