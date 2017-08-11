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
      print *,"Run offset systematics"
      do i=1,(2*NSYSOTOT)

          print *,"Current offset systematic: ",i,"/",(2*NSYSOTOT)


C     Recalculate central values to estimate given offset systematics
C     Do nothing is case of nominal
C     Run nominal at the end
          call ActivateOffsetSystematic(i)

          call RunIterativeAveraging(diag,last,corr,box)

          call SaveAverageOValue(i)
      enddo
C     Calculate offset systematics
      call CalcOffsetSyst()

C     Check impact of systametics 
C     Loop over all non-offset systematics     
      print *,"Run syst. shifts" 
      do i=1,(2*NSYSTOT)

          print *,"Current systematic: ",i,"/",(2*NSYSTOT)

C     Recalculate central values to estimate given systematics
          call ActivateSystematic(i)

          call RunIterativeAveraging(diag,last,corr,box)

          call SaveAverageValue(i)
      enddo


C     Run ToyMC for statistical uncertainties
C     Loop over ToyMC shifts     
      print *,"Run stat ToyMC" 
      do i=1,(nToyMC)

          print *,"Current ToyMC: ",i,"/",(nToyMC)

C     Recalculate central values to estimate given systematics
          call ActivateStatToyMC(i)

          call RunIterativeAveraging(diag,last,corr,box)

          call SaveAverageToyMCValue(i)
      enddo
C     Calculate offset systematics
      call CalcToyMCStat()



C     Run nominal averaging
      print *,"Run nominal averaging"

C     Loop over all point and measurements
      do iP=1,NMeas
        do idata=1,NMeasF2(iP)
          F2TAB(ip,idata) = F2TABOrig(ip,idata)
        enddo
      enddo



C     Calculate offset systematics
      call CalcSystImpact()

      call RunIterativeAveraging(diag,last,corr,box)

      Call AnalyseShifts() ! Analyse systematic shifts over iterations

C     Prepare output, rotate syst.
C     Do not run for offset systematics
      call LastIteration(diag,last,corr,box)



      end



C     Averaging of data from different experiments
      Subroutine RunIterativeAveraging(diag,last,corr,box)

      implicit none
      include 'common.inc'

      real*8
     $     DIAG(NF2MAX),LAST(NF2MAX+NSYSTMAX),CORR(NSYSTMAX,NF2MAX)
      real*8  Box(NSYSTMAX,NSYSTMAX)



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

      end

