C     Averaging of data from different experiments
      Subroutine Averaging

      implicit none
      include 'common.inc'
C     Intermediate dimensions:
      real*8
     $     DIAG(NF2MAX),LAST(NF2MAX+NSYSTMAX),CORR(NSYSTMAX,NF2MAX)
      real*8  Box(NSYSTMAX,NSYSTMAX)

C     A copy:
      real*8 boxs(NSYSTMAX,NSYSTMAX)
      real*8 Diags(NF2MAX),Lasts(NF2MAX)
     $     ,CORRs(NSYSTMAX,NF2MAX)

      integer i,iP,idata,iSyst
C--------------------------------------------------------------
      print *,"Number of systematics",NSYSTOT,NSYSOTOT

C Allocate arrays:
      

C     Loop over offset systematics
      do i=1,(2*NSYSOTOT+1)

          print *,"Current systematic: ",i,"/",(2*NSYSOTOT+1)


C     Recalculate central values to estimate given offset systematics
C     Do nothing is case of nominal
          call ActivateSystematic(i)

C     Loop over iterations
          do iItr=0,NIteration

              print *,"Iteration",iItr,"/",NIteration

C     If there is a multiplicative treatment, recalculate stat errors and repeat the average:
              if (iItr.ne.0) then
                  Call StatRecalc
              endif

              Call FillArrays(diag,last,corr,box)      ! Prepare the system of equations

C     Find averaged value and systematics:
              Call ToBlockDiag(diag,last,corr,box)     ! Perform fast inversion
          enddo

          Call AnalyseShifts() ! Analyse systematic shifts over iterations
          
          call LastIteration(diag,last,corr,box)       ! Prepare output, rotate syst.
          call SaveAverageValue(i)
      enddo

C     Calculate offset systematics
      call CalcOffsetSyst()


      end
