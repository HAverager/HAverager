C--------------------------------------------------------------
C The program performs the linear averaging of the experimental
C data assuming the uncorrelated measurements form different
C experiments.
C--------------------------------------------------------------
      Program averager
      implicit none
      include 'common.inc'
      logical doCovconversion
      
C--------------------------------------------------------------
      print *,'--------------------------------'
     $     //'---------------------------------------------------'
      print *,'    Initiating Iterative Linear Averager With'//
     $     ' Systematic Uncertainties'
      print *,' '
      print *,'    HAverager '
      print *,' '
      print *,'---------------------------------------------------'
     $     //'--------------------------------'
C--------------------------------------------------------------

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

 
      call HF_errsum(6)
C--------------------------------------------------------------
      end
