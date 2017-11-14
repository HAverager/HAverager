C     Averaging of data from different experiments
      Subroutine Averaging

      implicit none
      include 'common.inc'
C     Intermediate dimensions:
      real*8
     $     DIAG(NF2MAX),LAST(NF2MAX+NSYSTMAX),CORR(NSYSTMAX,NF2MAX)
      real*8  Box(NSYSTMAX,NSYSTMAX)

      integer i,iP,idata,iSyst
C--------------------------------------------------------------
      print *,"Number of systematics",NSYSTOT,NSYSOTOT

C Allocate arrays:

c      if (.false.) then
      

C     Loop over offset systematics
      print *,"Run offset systematics"
      do i=1,0! (2*NSYSOTOT)

          print *,"Current offset systematic: ",i,"/",(2*NSYSOTOT)


C     Recalculate central values to estimate given offset systematics
C     Do nothing is case of nominal
C     Run nominal at the end
          call ActivateOffsetSystematic(i)

          call RunIterativeAveraging(diag,last,corr,box,.true.,.false.)

          call SaveAverageOValue(i)
      enddo
C     Calculate offset systematics
      call CalcOffsetSyst()


C     Check impact of systametics 
C     Loop over all non-offset systematics     
      print *,"Run syst. shifts" 
      do i=1,0 ! (2*NSYSTOT)

          print *,"Current systematic: ",i,"/",(2*NSYSTOT)

C     Recalculate central values to estimate given systematics
          call ActivateSystematic(i)

          call RunIterativeAveraging(diag,last,corr,box,.true.,.false.)

          call SaveAverageValue(i)
      enddo

C     Run ToyMC for statistical uncertainties
C     Loop over ToyMC shifts     
      print *,"Run stat ToyMC" 
      do i=1,0 ! (nToyMC)

          print *,"Current ToyMC: ",i,"/",(nToyMC)

C     Recalculate central values to estimate given systematics
          call ActivateStatToyMC(i)

          call RunIterativeAveraging(diag,last,corr,box,.true.,.false.)

          call SaveAverageToyMCValue(i)
      enddo

C     Calculate ToyMC systematics
cc      call CalcToyMCStat()
c      endif


C     Run nominal averaging
      print *,"Run nominal averaging"

C     Loop over all point and measurements
      do iP=1,NMeas
        do idata=1,NMeasF2(iP)
          F2TAB(ip,idata) = F2TABOrig(ip,idata)
        enddo
      enddo

      call RunIterativeAveraging(diag,last,corr,box,.true.,.true.)


C     Analyse systematic shifts over iterations
      Call AnalyseShifts()

C     Calculate offset systematics
      call CalcSystImpact()

      end



C     Averaging of data from different experiments
      Subroutine RunIterativeAveraging(diag,last,corr,box,
     & onlyLast,lastItr)

      implicit none
      include 'common.inc'
      logical onlyLast
      logical lastItr
      integer i,j,k

      real*8
     $     DIAG(NF2MAX),LAST(NF2MAX+NSYSTMAX),CORR(NSYSTMAX,NF2MAX)
      real*8  Box(NSYSTMAX,NSYSTMAX)

C     A copy:
      real*8 boxs(NSYSTMAX,NSYSTMAX)
      real*8 Diags(NF2MAX),Lasts(NF2MAX+NSYSTMAX)
     $     ,CORRs(NSYSTMAX,NF2MAX)


C     Loop over iterations
      do iItr=0,NIteration

          print *,"Iteration",iItr,"/",NIteration,onlyLast

C     If there is a multiplicative treatment, recalculate stat errors and repeat the average:
          if (iItr.ne.0) then
c              Call StatRecalc
          endif

C         Prepare the system of equations
          Call FillArrays(diag,last,corr,box,onlyLast)

C         Copy all arrays
          diags(1:NMeas) = diag(1:NMeas)
          corrs(1:NSysTot,1:NMeas) = corr(1:NSysTot,1:NMeas) 
          boxs(1:NsysTot,1:NsysTot) = box(1:NsysTot,1:NsysTot)
          lasts(1:NMeas+NSysTot) = last(1:NMeas+NSysTot)

          print *,'BOX',box(1,1)
          
C         Find averaged value and systematics:
          Call ToBlockDiag(diags,lasts,corrs,boxs) 
      enddo
      stop
      if(lastItr)then
C       Prepare output, rotate syst.
C       Do not run for offset systematics
        call LastIteration(diags,lasts,corrs,boxs)
      endif
      end

