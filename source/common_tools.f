C
C Helper functions:
C
         !> Get bin index from bin name for a bin-grid
      integer function GetBinIndexGrid(Name,IGrid)
      implicit none
      include 'common.inc'
      character *(*) name
      integer IGrid
      integer i
C---------------------------------------------------
      GetBinIndexGrid = 0
      
      if (IGrid.le.0 .or. IGrid .gt. NGrid ) then
         call hf_errlog(2041201,
     $        'S:GetBindIndexGrid Grid index outside range')
         Return
      endif

      do i=1,NDimensionGrid(IGrid)
C-----------------------
         if (trim(Name) .eq. trim(GridBinNames(i,IGrid))) then
            GetBinIndexGrid = i
            Return
         endif
      enddo
C--------------------------------------------------
      end


         !> Get bin index from bin name
      integer function GetBinIndex(Name,IDataSet)
      implicit none 
      character *(*) name
      integer IDataSet
      include 'common.inc'
      integer i
C---------------------------------------------
      GetBinIndex  = 0
      do i=1,DATASETBinningDimension(IDataSet)
         if (Name .eq. DATASETBinNames(i,IDataSet) ) then
            GetBinIndex = i
            Return
         endif
      enddo
C---------------------------------------------
      end

      !> Function to get data point index, given reaction, grid index, CME and info fields.
      integer Function GetDataIndex(idxGrid, idxReaction)
C---------------------------------------------------------------
      implicit none
      integer idxGrid, idxReaction
      include 'common.inc'
      integer iData
      integer iProcClass

      integer GetProcessClass
C----------------------------------------------------------------
      
      GetDataIndex = NMeas+1  !> By default, add extra point
      
      do iData = 1,NMeas
        iProcClass = idxProcessClass(iData)  !> Reference to process class store

        if ( idxGrid .eq. idxGridMeas(iData)
     $        .and. idxReaction .eq. idxReactionMeas(iProcClass)) then
          GetDataIndex = iData  ! candidate point
          return;
        endif
      enddo

      !> Check if we want to add a new entry:
      if ( GetDataIndex .eq. NMeas + 1) then
        !> Check for overflow:
        NMeas = NMeas + 1
        if ( NMeas .gt. NF2MAX) then
           Call hf_ErrLog(20031301,
     $     'F:Too many data points.'
     $     //' Increase NF2MAX in settings.inc and rebuild')
        endif

        !> Store:
        idxGridMeas(NMeas)     = idxGrid
        idxProcessClass(NMeas) = GetProcessClass(idxReaction)
      endif
   
      end

C============================================
      !> Function to determine "Process Class". Add new entry if not present.
      integer function GetProcessClass(idxReaction)

      implicit none
      include 'common.inc'
      integer idxReaction
      integer iClass
C-----------------------------------------------------------------------------

      GetProcessClass = NProcClass + 1
C     First check if already present:
      do iClass=1,NProcClass
        if ( idxReaction .eq. idxReactionMeas(iClass) ) then
          !> Candidate
          GetProcessClass = iClass
          return

        endif
      enddo

      if (GetProcessClass .eq. NProcClass+1) then
        ! Add new entry:
        NProcClass = NProcClass + 1
        if (NProcClass .gt. NProcClassMax) then
          call hf_errlog(1,
     $           'F:Exceeding number of process classes.'// 
     $           'Increase value of NProcClassMax')
        endif
        idxReactionMeas(NProcClass) = idxReaction
      endif

      end

      
      !> Function to store data point. Returns data entry
      integer Function StoreData(idxData, Value, StatUncor, Stat, 
     $     Uncor, Total, NSys, Syst, iSys, systype, iDataFile )

      implicit none
      include 'common.inc'
      integer idxData, NSys, iDataFile
      double precision Value, StatUncor, Stat, Uncor, Total, Syst(*)
      integer SysType(*), iSys(*)
      integer idxP, i
C-----------------------------

      NMEASF2(idxData) = NMEASF2(idxData) + 1
       if (NMEASF2(idxData).gt.NMEASMAX) then
         Call HF_ErrLog(20031303,
     $ 'F:Number of measurements exceed limit per grid point.'//
     $ ' Increase NMEASMAX in settings.inc')
      endif
      idxP =  NMEASF2(idxData)
      StoreData = idxP   

      F2TAB(idxData,idxP)          = Value
      F2TABOrig(idxData,idxP)      = Value

      F2ETAB(idxData,idxP)         = StatUncor
      F2ETABOrig(idxData,idxP)     = StatUncor

      F2ETAB_STA(idxData,idxP)     = Stat
      F2ETAB_STAORIG(idxData,idxP) = Stat
      F2ETAB_UNC(idxData,idxP)     = Uncor
      F2ETAB_TOT(idxData,idxP)     = Total
      F2DataFile(idxData,idxP)     = iDataFile

C     Define global variable for type of shift >2 - Up/Down shift
      do i = 1, NSys
        if(iSys(i).ne.0)then
          SYSTABOrig(iSys(i),idxData,idxP,SysType(i)) = Syst(i)
        endif
      enddo

C     Fill SysTab array using SYSTABOrig
      do i = 1, NSys
        if(iSys(i).eq.0)then
          continue
        endif
        if(iSys(i).le. NSYSTMAX)then
        if(SysForm(iSys(i)).eq.12 .or.
     &     SysForm(iSys(i)).eq.22)then
          SysTab(iSys(i),idxData,idxP) =
     &     (SYSTABOrig(iSys(i),idxData,idxP,1) -
     &     SYSTABOrig(iSys(i),idxData,idxP,2))/2
        else if(SysForm(iSys(i)).eq. 11 .or. 
     &          SysForm(iSys(i)).eq. 21 ) then
          SysTab(iSys(i),idxData,idxP) =
     $     SYSTABOrig(iSys(i),idxData,idxP,1)
        endif
        endif
      enddo

      end


      !> Update data point using stat. uncertainty, pre-average
      Subroutine PreAverageData(idxData, idxP, Value, StatUncor, Stat, 
     $     Uncor, Total, NSys, Syst, SystematicType )

      implicit none
      include 'common.inc'
      integer idxData, NSys
      double precision Value, StatUncor, Stat, Uncor, Total, Syst(NSys)
      character*(*) SystematicType(NSys)
      integer idxP
      integer jj, k
      double precision W1, W2, Wtemp
C-----------------------------
      W1 = 1.0D0 / F2ETAB(idxData,idxP)**2
      W2 = 1.0D0 / Stat**2
      Wtemp = 1.0D0 / ( W1 + W2 )
      W1 = W1*Wtemp
      W2 = W2*Wtemp
      
      F2TAB(idxData,idxP)          = W1*F2TAB(idxData,idxP) + W2*Value
      F2ETAB_STA(idxData,idxP)     = sqrt(WTemp)
      F2ETAB_STAORIG(idxData,idxP) = F2ETAB_STA(idxData,idxP)
      F2ETAB_UNC(idxData,idxP)     = sqrt(
     $     W1**2*F2ETAB_UNC(idxData,idxP)**2 +
     $     W2**2*Uncor**2)

      F2ETAB(idxData,idxP)         = sqrt(
     $    F2ETAB_STA(idxData,idxP)**2+ F2ETAB_UNC(idxData,idxP)**2)
      F2ETABOrig(idxData,idxP)     = F2ETAB(idxData,idxP)

      F2ETAB_TOT(idxData,idxP)     = sqrt(
     $     W1**2*F2ETAB_TOT(idxData,idxP)**2 +
     $     W2**2*Total**2)

C     FixME
      do jj = 1, NSys
         do k  = 1, NSysTot
            if ( SystematicName(k).eq.SystematicType(jj) ) then
               SysTab(k,idxData,idxP) =
     $              W1*SysTab(k,idxData,idxP) +
     $              W2*Syst(jj) 
               SYSTABOrig(k,idxData,idxP,1) =
     &          SysTab(k,idxData,idxP)
            endif
         enddo
      enddo

      end



      integer function getreactionidx(reaction)
      implicit none      
      character*(*) reaction
      integer i
      include 'common.inc'
C-----------------------------------------------------
      do i=1,NGrid
         write(*,*) 'gridreaction(i) = ', gridreaction(i)
         if (reaction.eq.gridreaction(i)) then
            getreactionidx = i
            return
         endif
      enddo
      print *,'reaction ', reaction
      print *,'not found. Stop'
      call hf_stop
C-----------------------------------------------------
      end



      integer function getgrididx(bins,idxReaction,iNCM)
      implicit none
      include 'common.inc'
      integer iNCM
      real bins(iNCM), temp
      integer idxReaction
      integer i,j
      logical match
      match_parameter = 1.0D-5
C------------------------------------------------------
      do i=1,NPointGrid(idxReaction)
         match = .true.
         do j=1,NDimensionGrid(idxReaction)
            match = match .and.
     $           abs((bins(j)-GridPoints(i,j,idxReaction))
C protection against exact zero:
     $           /(bins(j)+match_parameter))
     $           .lt. match_parameter
         enddo
         if (match) then
            getgrididx = i
            Return
         endif
      enddo
C------------------------------------------------------
      getgrididx = -1
      end


C-----------------------------------------------------------------------------
!> Add systematic source (26.10.2015)
C
!>  Detect "+" and "-" signs, at the end of source name, for asymmetric errors (x2)
!>  Symmectic systematic uncertainties marked as (x1)
!>
!>  Detect ":" modifiers
!>
!>   :M  - "multiplicative" (1x)
!>
!>   :A  - "additive" (2x)
!>
C-----------------------------------------------------------------------------
      subroutine AddSystematics(SName, idx, iSys, sysType)

      implicit none
      include 'common.inc'

      character*(*) SName
      integer sysType(*) ! type of systematic idx (output)
      integer iSys(*) ! index of systematic idx in globar array (output)
      integer idx     ! index of systematic in local array

      character *32 CurrentSysName
      character *1 ctmp
      integer  CurrentSysForm

      integer j

C-----------------------------------------

C---- Cut +/- in the sys name
      if(Index(SName, "+").gt.0) then
           CurrentSysName =
     &    SName(:Index(SName, "+")-1) //
     &    SName(Index(SName, "+")+1:)
          sysType(idx) = 1
          CurrentSysForm = 12
      else if(Index(SName, "-").gt.0) then
           CurrentSysName =
     &    SName(:Index(SName, "-")-1) //
     &    SName(Index(SName, "-")+1:)
          sysType(idx) = 2
          CurrentSysForm = 12
      else
          CurrentSysName = SName
          sysType(idx) = 1
          CurrentSysForm = 11
      endif

C---- Cut systematics form
      if(Index(CurrentSysName, ":").gt.0) then
          ctmp =
     &     CurrentSysName(Index(CurrentSysName, ":")+1:)
          if(ctmp.eq.'A')then
            CurrentSysForm = CurrentSysForm + 10
          endif
          CurrentSysName =
     &     CurrentSysName(:Index(CurrentSysName, ":")-1)
      endif

C---  Case of offset systematic
      if(ctmp=="O")then
C---  Check if the source already exists
          do j=1,NSYSOTOT
               if (SystematicName(NSYSTMAX+j).eq.CurrentSysName) then
                  iSys(idx) = NSYSTMAX+j
                  return
               endif
          enddo

C---  Add new source
          NSYSOTOT = NSYSOTOT + 1
          if (NSYSOTOT.gt.NSystOMax) then
               call hf_errlog(10,
     $        'F:AddSystematics Error: exceeding NSystMax'
     $         //'Increase velue of NSystMax in settings.inc')
          endif

          print *,"Add Offset Sys Source ",NSYSOTOT,CurrentSysName
          iSys(idx) = NSYSTMAX+NSYSOTOT
          SystematicName(NSYSTMAX+NSYSOTOT) = CurrentSysName
          SysForm(NSYSTMAX+NSYSTOT) = 10+sysType(idx)
          return;
      else
C---  Check if the source already exists
          do j=1,NSYSTOT
               if ( SystematicName(j).eq.CurrentSysName) then
                  iSys(idx) = j
                  return
               endif
          enddo

C---  Add new source
          NSYSTOT = NSYSTOT + 1
          if (NSYSTOT.gt.NSystMax) then
               call hf_errlog(10,
     $        'F:AddSystematics Error: exceeding NSystMax'
     $         //'Increase velue of NSystMax in settings.inc')
          endif

          print *,"Add Sys Source ",NSYSTOT,trim(CurrentSysName),
     &    ", form: ",CurrentSysForm
          iSys(idx) = NSYSTOT
          SystematicName(NSYSTOT) = CurrentSysName
          SysForm(NSYSTOT) = CurrentSysForm
          return;
      endif

      end

C-----------------------------------------------------------------------------
!>    Recalculate central values to estimate given offset systematics
!>    Do nothing is case of nominal
C-----------------------------------------------------------------------------
      Subroutine ActivateOffsetSystematic(i)
        implicit none
        include 'common.inc'

        integer i
        integer iSyst, iP, idata

        if(mod(i,2).eq.0)then        
            iSyst = i/2
            print *,"Do variation Down for ",
     & SystematicName(NSYSTMAX+iSyst)
C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C     Asymmetric down case
                if(SysForm(NSYSTMAX+iSyst).eq. 12) then
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &             SYSTABOrig(NSYSTMAX+iSyst,ip,idata,2)
C     Symmetric case
                else
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) -
     &               SYSTABOrig(NSYSTMAX+iSyst,ip,idata,1)
                endif

              enddo
            enddo
        else
            iSyst = (i+1)/2
            print *,"Do variation Up for ",
     & SystematicName(NSYSTMAX+iSyst)
C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C     Asymmetric up/symmetric case
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &               SYSTABOrig(NSYSTMAX+iSyst,ip,idata,1)
              enddo
            enddo
            return
       endif

      end


C-----------------------------------------------------------------------------
!>    Recalculate central values to estimate impact of given systematics
C-----------------------------------------------------------------------------
      Subroutine ActivateSystematic(i)
        implicit none
        include 'common.inc'

        integer i
        integer iSyst, iP, idata

        if(mod(i,2).eq.0)then        
            iSyst = i/2
            print *,"Do variation Down for ",SystematicName(iSyst)
C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C     Asymmetric down case
                if(SysForm(iSyst).eq. 12 .or.
     &             SysForm(iSyst).eq. 22 ) then
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &             SYSTABOrig(iSyst,ip,idata,2)
C     Symmetric case
                else
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) -
     &               SYSTABOrig(iSyst,ip,idata,1)
                endif

              enddo
            enddo
        else
            iSyst = (i+1)/2
            print *,"Do variation Up for ",SystematicName(iSyst)
C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C     Asymmetric up/symmetric case
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &               SYSTABOrig(iSyst,ip,idata,1)
              enddo
            enddo
            return
       endif

      end


C-----------------------------------------------------------------------------
!>    Recalculate central values to run stat. ToyMC
C-----------------------------------------------------------------------------
      Subroutine ActivateStatToyMC(i)
        implicit none
        include 'common.inc'
        EXTERNAL RANLUX


        integer i
        integer N, iP, idata, ierr
        real amu, E
        real gaus

C     Do Shift
        print *,"Do ToyMC shift ",i
C     Loop over all point and measurements
        do iP=1,NMeas
          do idata=1,NMeasF2(iP)
            E = F2TABOrig(ip,idata) /
     &            F2ETAB_STAORIG(ip,idata)**2
            amu = F2TABOrig(ip,idata)*E
            CALL RNPSSN(amu,N,ierr)
            CALL RNORMX(gaus,1,RANLUX)

C           Use puissonin errors
C           F2TAB(ip,idata) = N / E

C           Use gaussian errors
            F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &       F2ETAB_STAORIG(ip,idata)*gaus
          enddo
        enddo

      end

C-----------------------------------------------------------------------------
!>    Save average value after a certain offset systematic
C-----------------------------------------------------------------------------
      Subroutine SaveAverageOValue(i)
        implicit none
        include 'common.inc'
        integer i, iP

C     Loop over all point
        do iP=1,NMeas
            F2VaveOSyst(iP,i) = f2vave(iP)
        enddo

      end

C-----------------------------------------------------------------------------
!>    Save average value after a certain systematic shift of measured point
C-----------------------------------------------------------------------------
      Subroutine SaveAverageValue(i)
        implicit none
        include 'common.inc'
        integer i, iP

C     Loop over all point
        do iP=1,NMeas
            F2VaveSyst(iP,i) = f2vave(iP)
        enddo

      end

C-----------------------------------------------------------------------------
!>    Save average value after a certain systematic shift of measured point
C-----------------------------------------------------------------------------
      Subroutine SaveAverageToyMCValue(i)
        implicit none
        include 'common.inc'
        integer i, iP

C     Loop over all point
        do iP=1,NMeas
            F2VaveToyMCStat(iP,i) = f2vave(iP)
        enddo

      end


C-----------------------------------------------------------------------------
!>    Save chi2 after a certain iterations
C-----------------------------------------------------------------------------
      Subroutine SaveChi2()
        implicit none
        include 'common.inc'
        integer ndf
        real chi2

        call CalcChi2(chi2, ndf)
        Chi2Itr(iItr) = chi2
      end


C-----------------------------------------------------------------------------
!>    Calculate offset uncertainties
C-----------------------------------------------------------------------------
      Subroutine CalcOffsetSyst()
        implicit none
        include 'common.inc'
        integer if2, isys

C     Loop over bins
        do if2=1,NMeas
C     Loop over systematics
          do isys=1,NSYSOTOT
              OSyst(if2,isys)=
     $         ((F2VaveOSyst(if2,2*isys-1)-f2vave(if2)) -
     $        (F2VaveOSyst(if2,2*isys)-f2vave(if2)))/2
              OSystPercent(if2,isys)=OSyst(if2,isys)/f2vave(if2)*100
          enddo
        enddo
      end

C-----------------------------------------------------------------------------
!>    Calculate ToyMC statistical uncertainties
C-----------------------------------------------------------------------------
      Subroutine CalcToyMCStat()
        implicit none
        include 'common.inc'
        integer if2, isys
        real*8 b,b2

C     Loop over bins
        do if2=1,NMeas
C     Loop over systematics
          b = 0
          b2 = 0
          do isys=1,nToyMC
            b = b + F2VaveToyMCStat(if2,isys)
            b2 = b2 + F2VaveToyMCStat(if2,isys)**2
          enddo
          b = b/nToyMC
          b2 = b2/nToyMC 
          StatToyMC(if2)= SQRT(b2 - b**2)
        enddo

      end
C-----------------------------------------------------------------------------
!>    Calculate impact of systematics
C-----------------------------------------------------------------------------
      Subroutine CalcSystImpact()
        implicit none
        include 'common.inc'
        integer if2, isys

C     Loop over bins
        do if2=1,NMeas
C     Loop over systematics
          do isys=1,NSYSTOT
              SystImpact(if2,isys)=
     $         ((F2VaveSyst(if2,2*isys-1)-f2vave(if2)) -
     $        (F2VaveSyst(if2,2*isys)-f2vave(if2)))/2
          enddo
        enddo
      end

C-----------------------------------------------------------------------------
!>    Analysis of shifts of asymmetric systematics in order to understand if
!>    iterative procedure will converge.
C-----------------------------------------------------------------------------
      Subroutine AnalyseShifts()
        implicit none
        include 'common.inc'
        real SysShDiff(NIteration)
        integer isys,itr

C     Loop over systematics
        do isys=1,nsystot

            do itr=1,NIteration
                SysShDiff(itr)=
     &      (SYSSHItr(iSys,itr)-SYSSHItr(iSys,itr-1))/SYSSHItr(iSys,itr)
            enddo

            if(abs(SysShDiff(NIteration))
     &           .gt.SysShTolerance)then
               print *,"Warning: Systematic ",trim(SystematicName(isys))
               print *," does not converge over iterations"
               print *,"Tolerance parameter: SysShTolerance = ",
     &               SysShTolerance

            endif
        enddo

      end


