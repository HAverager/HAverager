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

C     Loop over offset systematics
      print *,"Run offset systematics"
      do i=1,(2*NSYSOTOT)

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
      if(doSystImpact)then
C     Loop over all non-offset systematics     
      print *,"Run syst. shifts" 
      do i=1,(2*NSYSTOT)

          print *,"Current systematic: ",i,"/",(2*NSYSTOT)

C     Recalculate central values to estimate given systematics
          call ActivateSystematic(i)

          call RunIterativeAveraging(diag,last,corr,box,.true.,.false.)

          call SaveAverageValue(i)
      enddo
      endif

C     Run ToyMC for statistical uncertainties
      if(nToyMC .gt. 0)then
C     Loop over ToyMC shifts     
      print *,"Run stat ToyMC" 
      do i=1,(nToyMC)

          print *,"Current ToyMC: ",i,"/",(nToyMC)

C     Recalculate central values to estimate given systematics
          call ActivateStatToyMC(i)

          call RunIterativeAveraging(diag,last,corr,box,.true.,.false.)

          call SaveAverageToyMCValue(i)
      enddo

C     Calculate ToyMC systematics
      call CalcToyMCStat()
      endif

C     Run nominal averaging
      print *,"Run nominal averaging"

C     Loop over all point and measurements
      do iP=1,NMeas
        do idata=1,NMeasF2(iP)
          F2TAB(ip,idata) = F2TABOrig(ip,idata)

          F2ETAB(ip,idata) = F2ETABOrig(ip,idata)
          do iSyst=1,NSYSOTOT

        if(SysForm(iSyst).eq.12 .or.
     &      SysForm(iSyst).eq.22)then
          SysTab(iSyst,ip,idata) =
     &     (SYSTABOrig(iSyst,ip,idata,1) -
     &     SYSTABOrig(iSyst,ip,idata,2))/2
        else if(SysForm(iSyst).eq. 11 .or. 
     &    SysForm(iSyst).eq. 21 ) then
          SysTab(iSyst,ip,idata) =
     $     SYSTABOrig(iSyst,ip,idata,1)
        endif

          enddo
        enddo
      enddo

      call RunIterativeAveraging(diag,last,corr,box,.true.,.true.)


C     Analyse systematic shifts over iterations
      Call AnalyseShifts()

C     Calculate offset systematics
      if(doSystImpact)then
          call CalcSystImpact()
      endif

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

      double precision time1,time2

C     Loop over iterations
      do iItr=0,NIteration

          print *,"Iteration",iItr,"/",NIteration,onlyLast

C     If there is a multiplicative treatment, recalculate stat errors and repeat the average:
          if (iItr.ne.0) then
              Call StatRecalc
          endif

C         Prepare the system of equations
          Call FillArrays(diag,last,corr,box,onlyLast)

C         Copy all arrays
          diags(1:NMeas) = diag(1:NMeas)
          corrs(1:NSysTot,1:NMeas) = corr(1:NSysTot,1:NMeas) 
          boxs(1:NsysTot,1:NsysTot) = box(1:NsysTot,1:NsysTot)
          lasts(1:NMeas+NSysTot) = last(1:NMeas+NSysTot)
          
C         Find averaged value and systematics:
          Call ToBlockDiag(diags,lasts,corrs,boxs) 
      enddo

      call cpu_time(time1)
      print *,'here3',time2,time1,time1-time2
      if(lastItr)then
C       Prepare output, rotate syst.
C       Do not run for offset systematics
        call LastIteration(diags,lasts,corrs,boxs,box)
      endif
      call cpu_time(time2)
      print *,'LastItr: ',time2,time1,time2-time1
      end

