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
     $     Uncor, Total, NSys, Syst, SystematicType, iDataFile )

      implicit none
      include 'common.inc'
      integer idxData, NSys, iDataFile
      double precision Value, StatUncor, Stat, Uncor, Total, Syst(*)
      character*(*) SystematicType(*)
      integer idxP
      integer jj, k

      Integer nPlus, nMinus
      integer systype
      character *32 CurrentSysName
      character *1 CurrentSysForm
C-----------------------------

      !> Add an entry:
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
C      print *, Nsys, NSysTot

      do jj = 1, NSys
C     Check up (2) and down (3) systematics
         if(Index(SystematicType(jj), "+").gt.0) then
           CurrentSysName =
     &    SystematicType(jj)(:Index(SystematicType(jj), "+")-1) //
     &    SystematicType(jj)(Index(SystematicType(jj), "+")+1:)
          systype = 2
         else if(Index(SystematicType(jj), "-").gt.0) then
           CurrentSysName =
     &    SystematicType(jj)(:Index(SystematicType(jj), "-")-1) //
     &    SystematicType(jj)(Index(SystematicType(jj), "-")+1:)
          systype = 3
         else
            CurrentSysName = SystematicType(jj)
            systype = 1
         endif

C---- Cut systematics form
         if(Index(CurrentSysName, ":").gt.0) then
           CurrentSysForm =
     &     CurrentSysName(Index(CurrentSysName, ":")+1:)
           CurrentSysName =
     &     CurrentSysName(:Index(CurrentSysName, ":")-1)
         else
            CurrentSysForm = "M"
         endif

C---  If it is offset systematic
         if(CurrentSysForm .eq. "O")then
             do k  = 1, NSysOTot
                if ( SystematicOName(k).eq.CurrentSysName ) then
C     Define global variable for type of shift >2 - Up/Down shift
C            print *,systype
                  if(systype.gt.1)then
                    ShiftOType(k,idxData,idxP) =
     &              ShiftOType(k,idxData,idxP) + 1
C                   print *,"DoIncreaseShiftType",ShiftType(k)
                  endif
                  SYSTABOOrig(k,idxData,idxP,systype) = Syst(jj)
                  if (IDebug.gt.2) then
                     print *, "Offset",k, idxData, idxP, systype,
     &               Syst(jj),ShiftOType(k,idxData,idxP)
                  endif

                  if(ShiftOType(k,idxData,idxP).eq.1) then
                     continue
                  endif
C                 print *,"OType",ShiftOType(k),k,SystematicOName(k)
                endif
             enddo
         else

C---  non-offset systematic
         do k  = 1, NSysTot
            if ( SystematicName(k).eq.CurrentSysName ) then
C     Define global variable for type of shift >2 - Up/Down shift
C            print *,systype
              if(systype.gt.1)then
                ShiftType(k,idxData,idxP) =
     &           ShiftType(k,idxData,idxP) + 1
C                print *,"DoIncreaseShiftType",ShiftType(k)
              endif
               SYSTABOrig(k,idxData,idxP,systype) = Syst(jj)
               if (IDebug.gt.2) then
                 print *, k, idxData, idxP, systype, Syst(jj),
     &           ShiftType(k,idxData,idxP)
               endif

C     Fill SysTab array using SYSTABOrig
               if(ShiftType(k,idxData,idxP).eq.2) then
                SysTab(k,idxData,idxP) =
     &         (SYSTABOrig(k,idxData,idxP,2) -
     &          SYSTABOrig(k,idxData,idxP,3))/2
               else if(ShiftType(k,idxData,idxP).eq.0) then
                 SysTab(k,idxData,idxP) = SYSTABOrig(k,idxData,idxP,1)
               endif

               if(ShiftType(k,idxData,idxP).eq.1) then
                 continue
               endif
C               print *,"ShiftType",ShiftType(k),k,SystematicName(k)
            endif
         enddo

         endif
      enddo


C     Print error message if Up/Down part is missing
      do k  = 1, NSysTot
        if(ShiftType(k,idxData,idxP).eq.1) then
          call hf_errlog(11,
     $        'F:Missing Up or Down systematics for: '
     $         //SystematicName(k))
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
               SYSTABOrig(k,idxData,idxP,ShiftType(k,idxData,idxP)) =
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
!>  Detect "+" and "-" signs, at the end of source name, for asymmetric errors
!>
!>  Detect ":" modifiers
!>
!>   :M  - "multiplicative"
!>
!>   :A  - "additive"
!>
C-----------------------------------------------------------------------------
      subroutine AddSystematics(SName)

      implicit none
      include 'common.inc'

      character*(*) SName
      character *32 CurrentSysName
      character *1  CurrentSysForm

      integer j

C-----------------------------------------

      print *,"New syst",SName

C---- Cut +/- in the sys name
      if(Index(SName, "+").gt.0) then
           CurrentSysName =
     &    SName(:Index(SName, "+")-1) //
     &    SName(Index(SName, "+")+1:)
      else if(Index(SName, "-").gt.0) then
           CurrentSysName =
     &    SName(:Index(SName, "-")-1) //
     &    SName(Index(SName, "-")+1:)
      else
            CurrentSysName = SName
      endif

C---- Cut systematics form
      if(Index(CurrentSysName, ":").gt.0) then
          CurrentSysForm =
     &     CurrentSysName(Index(CurrentSysName, ":")+1:)
          CurrentSysName =
     &     CurrentSysName(:Index(CurrentSysName, ":")-1)
      else
          CurrentSysForm = "M"
      endif

      print *,"Form: ",CurrentSysForm

C---  Case of offset systematic
      if(CurrentSysForm=="O")then
C---  Check if the source already exists
          do j=1,NSYSOTOT
               if ( SystematicOName(j).eq.CurrentSysName) then
                  return
               endif
          enddo

C---  Add new source
          NSYSOTOT = NSYSOTOT + 1
          if (NSYSOTOT.gt.NSystMax) then
               call hf_errlog(10,
     $        'F:AddSystematics Error: exceeding NSystMax'
     $         //'Increase velue of NSystMax in settings.inc')
          endif

          print *,"Add Offset Sys Source ",NSYSOTOT,CurrentSysName

          SystematicOName(NSYSOTOT) = CurrentSysName
          return;
      else
C---  Check if the source already exists
          do j=1,NSYSTOT
               if ( SystematicName(j).eq.CurrentSysName) then
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

          print *,"Add Sys Source ",NSYSTOT,CurrentSysName,"form: ",
     &    CurrentSysForm

          SystematicName(NSYSTOT) = CurrentSysName
          SystematicForm(NSYSTOT) = CurrentSysForm
          return;
      endif

      end

C-----------------------------------------------------------------------------
!>    Recalculate central values to estimate given offset systematics
!>    Do nothing is case of nominal
C-----------------------------------------------------------------------------
      Subroutine ActivateSystematic(i)
        implicit none
        include 'common.inc'

        integer i
        integer iSyst, iP, idata
        print *,"Activate",NSYSOTOT
C     Loop over offset systeatics
        do iSyst=1,NSYSOTOT

C     Do Shift Up
          if(((2*iSyst)-1) .eq. i)then
            print *,"Do variation Up for ",SystematicOName(iSyst)
C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C     Asymmetric case
                if(ShiftOType(iSyst,ip,idata) .eq. 2) then
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &             SYSTABOOrig(iSyst,ip,idata,2)
C     Symmetric case
                else
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &               SYSTABOOrig(iSyst,ip,idata,1)
                endif
C                  print *,F2TAB(ip,idata),F2TABOrig(ip,idata)
C                  print *,SYSTABOOrig(iSyst,ip,idata,1)
              enddo
            enddo
            return
          endif

C     Do Shift Down
          if((2*iSyst) .eq. i)then
            print *,"Do variation Down for ",SystematicOName(iSyst)
C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C     Asymmetric case
                if(ShiftOType(iSyst,ip,idata) .eq. 2) then
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) +
     &             SYSTABOOrig(iSyst,ip,idata,3)
C     Symmetric case
                else
                  F2TAB(ip,idata) = F2TABOrig(ip,idata) -
     &               SYSTABOOrig(iSyst,ip,idata,1)
                endif
C                  print *,F2TAB(ip,idata),F2TABOrig(ip,idata)
C                  print *,SYSTABOOrig(iSyst,ip,idata,1)
              enddo
            enddo
            return
          endif

        enddo

C     Loop over all point and measurements
            do iP=1,NMeas
              do idata=1,NMeasF2(iP)
C                  print *,ip,idata,F2TAB(ip,idata),F2TABOrig(ip,idata)
                  F2TAB(ip,idata) = F2TABOrig(ip,idata)
              enddo
            enddo

        print *,"Do nominal average "

      end

C-----------------------------------------------------------------------------
!>    Save average wavue after a certain offset systematic
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
!>    Save chi2 after a certain iterations
C-----------------------------------------------------------------------------
      Subroutine SaveChi2()
        implicit none
        include 'common.inc'
        integer i, j, isys
        real sum, chi2loc, chi2

C       Fill Chi^2 and NDoF
        chi2 = 0.0
        do i = 1, NMeas
         do j = 1, NMeasF2(i)
            sum = F2TAB(i,j)
            do isys = 1, NSYSTOT
               sum = sum + SYSTAB(isys,i,j)*SYSSH(isys)
            enddo
            chi2loc = (F2VAVE(i)-sum)/ (F2ETAB(i,j))
            chi2 = chi2 + chi2loc**2
         enddo
        enddo
        do isys=1,NSYSTOT
          chi2 = chi2 + SYSSH(isys)**2
        enddo
 
        Chi2Itr(iItr) = chi2 + SYSSH(isys)**2
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
     $         ((F2VaveSyst(if2,2*isys-1)-f2vave(if2)) -
     $        (F2VaveSyst(if2,2*isys)-f2vave(if2)))/2
              OSyst(if2,isys)=OSyst(if2,isys)/f2vave(if2)*100
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


