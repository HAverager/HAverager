      subroutine ReadXSTable(cfile)
C----------------------------------------------------------------
C
C Read cross-section table, fill info
C
C----------------------------------------------------------------
      implicit none
      integer ifile
C-----------------------------------------------------------------

      include 'common.inc'

      character *(*) CFile

      character *80 Name
      integer  NData
      integer  NUncert
      integer  NBinDimension
      

      character *80 BinName(NDimensionMax)
      character *80 Reaction

      double precision buffer(ncolumnMax)

C
C Name and type of columns:
C      
      integer   NColumn 
      character *32 ColumnName(ncolumnMax)
      character *32 ColumnType(ncolumnMax)

C Systematics:
      character *32 SystematicType(nsystMax)
      logical Percent(1:nsystMax)

C Namelist definition:
      namelist/Data/Name,NData
     $     ,Reaction,Percent
     $     ,ColumnName, ColumnType, NColumn

      double precision XSections(ndataMax)
      double precision AllBins(10,ndataMax)
      double precision Syst(nsystmax)

      double precision StatError   ! stat
      double precision StatErrorConst ! stat. error to be treated as constant
      double precision UncorError  ! uncorrelated systematics
      double precision TotalError  ! total uncertainty
      double precision StatUncorError   ! the sqrt(stat^2+uncor^2)

      double precision TotalErrorRead ! total error, provided by the data file

      integer idxSigma

      integer i,j,jj,ii,ij,k,iBin,iError,itemp

C Temporary buffer to read the data (allows for comments starting with *)
      character *40960 CTmp

      integer idxreaction,idxGrid
C Functions:
      integer getreactionidx, getgrididx
      
      real ClosestBin(ncolumnMax)  ! The correspondent bin in the bingrid
      
      double precision Wtemp, W1, W2, ty

C Index for data point:

      integer idxData

      integer iDataLast(NF2MAX) 

      integer GetBinIndex  ! Function to get index of bin given by name
      integer GetDataIndex ! Function to get index of data point 
      integer StoreData    ! Function to store data
C-------------------------------------------------------

C Reset to default:
      NUncert = 0
      NData = 0
      NBinDimension = 0
      Reaction = ' '
      Name     = ' '
      idxSigma = 0

      do i=1,NF2MAX
         iDataLast(i) = 0
      enddo

C Reset scales to 1.0
      do i=1,nsystmax
         ColumnType(i) = ' '
         ColumnName(i) = ' '
      enddo

      open(51,file=CFile,status='old',err=99)

      print *,'Reading data file ...'
      print *,CFile
      read(51,NML=Data,err=98)

C
C Save the name of dataset:
C
      NName = NName + 1
      write(*,*) 'NName= ', NName

C 
C Reaction index:
C 
      idxReaction = GetReactionIdx(reaction)
      write(*,*) ' idxReaction =', idxReaction

C
C Check number of data per file:
C
      if (NData.gt.ndataMax) then
          call hf_errlog(10,
     $        'F:Number of data points per file exceed maximum'
     $         //'Increase ndataMax in settings.inc')
      endif

C
C Check dimensions:
C
      if (NColumn.gt. Ncolumnmax) then
          call hf_errlog(10,
     $        'F:Error in ReadDataFile for File'//cfile
     $         //'number of columns exceed maximum. '
     $         //'Increase velue of NColumnMax in settings.inc')
      endif
C
C Store
C
      NDATASETS = NDATASETS + 1

      DATASETLABEL(NDATASETS)    = Name

C Reaction info:
      DATASETREACTION(NDATASETS) = Reaction
C
C Parse ColumnType, count systematics:
C
      do i=1,NColumn
         if (ColumnType(i).eq.'Bin') then
            NBinDimension = NBinDimension + 1
            if(NBinDimension.gt.NDimensionMax)then
                call hf_errlog(10,
     $        'F:Number of dimensions exceed maximum'
     $         //'Increase NDimensionMax in settings.inc')
            endif
            BinName(NBinDimension) = ColumnName(i)
         elseif (ColumnType(i).eq.'Sigma') then
            idxSigma = i
         elseif (ColumnType(i).eq.'Error') then
            NUncert = NUncert + 1
            ! Special case: ignore 
            if (index(ColumnName(i),'ignore').gt.0) then
               SystematicType(NUncert) = 'ignore'
            else
               SystematicType(NUncert) = ColumnName(i)
            endif
         elseif (ColumnType(i).eq.'Dummy') then
C Ignore dummy column
         else
             call hf_errlog(10,
     $        'F:Unknown Column type for dataset '//CFile
     $         //',type '//ColumnType(i)
     $         //'. STOP in ReadDataFile')
            call HF_stop
         endif
      enddo

      write(*,*) '*** NUNCERT = ', NUncert

C Binning info:
      DATASETBinningDimension(NDATASETS) = NBinDimension
C Filling with 'dummy' first three names for proper formation of fittedresults.txt
      do i=1,3
         DATASETBinNames(i,NDATASETS) = 'dummy'
      enddo
      do i=1,NBinDimension
         DATASETBinNames(i,NDATASETS) = BinName(i)
      enddo

C
C Prepare systematics:
C
      do i=1,NUncert

C--- Statistical: special case
         if (SystematicType(i).eq.'stat' 
     $        .or. SystematicType(i).eq.'Stat') then
           continue
         else if (SystematicType(i).eq.'stat const') then
            Call HF_ERRLOG(16020001,'I: Stat Const Error type used')
C--- Uncorrelated: special case
         else if (SystematicType(i).eq.'uncor') then
           continue
C--- Total error: special case
         else if (SystematicType(i).eq.'total') then
C--- Ignore: special case - will not be counted
         else if (SystematicType(i).eq.'ignore') then
           continue
         else
           call AddSystematics(SystematicType(i))
           continue
         endif
      enddo

C Output for debug:
      write(*,*) ' Number of stored syst. sources: ', NSYSTOT

      ii = 0

C
C Read data info:
C

      write(*,*) '*** NData = ', NData
      do j=1,NData
C Allow for comments:
 89      read (51,'(A)',err=1017,end=1018) ctmp
         if (ctmp(1:1).eq.'*') then
C Comment line, read another one
            goto 89
         endif

C Check coherence of the table info
         if (idxSigma.eq.0) then
             call hf_errlog(10,
     $        'F:No column contains Sigma keyword for'
     $         //' the x-section info!!!')
         endif
         do i=1,NColumn
            if (ColumnName(i) .eq. ' ') then
               print *,'Undefined ColumnName !!!'
               print *,'Check name for column number = ',i
               call HF_stop
            endif
         enddo

C Read the columns
         read (ctmp,*,err=1019)(buffer(i),i=1,NColumn)

C Decode the columns
         iBin   = 0
         iError = 0
         do i=1,NColumn
            if ( ColumnType(i).eq.'Bin' ) then
               iBin = iBin + 1
               allbins(iBin,j) = buffer(i)
            elseif ( ColumnType(i).eq.'Sigma' ) then
               XSections(j) = buffer(i)
C               write(*,*) j,' CS = ', XSections(j)
            elseif ( ColumnType(i).eq.'Error' ) then
               iError = iError + 1
C               write(*,*) 'iError = ',iError, buffer(i)
               syst(iError) = buffer(i)
            endif
         enddo
C--------------------------------------

C
C Find the nearest bin on the grid and make shift into the found bin:
C

C
C Linear search of the nearest bin:
C
      idxGrid = 0

      call    GetNearestSwimmingBinLinear ! linear search
     $        (idxReaction
     $        ,NBinDimension
     $        ,allbins(:,j)
     $        ,BinName(:)
     $        ,ClosestBin(:)
     $        ,idxGrid)



C Translate errors from % to normal:
      TotalError = 0.0d0
      UncorError = 0.0d0
      StatError = 0.0d0
      StatErrorConst = 0.0d0
      TotalErrorRead = 0.0d0

C
C CALCULATE UNCERTAINTIES:
C

      do i=1,NUncert
C     Make relative uncertainties:
         if (.not.Percent(i)) then
            Syst(i) = Syst(i)
         else
            Syst(i) = Syst(i)*XSections(j)/100.0d0
         endif

         if (SystematicType(i).eq.'total') then
            TotalErrorRead = Syst(i)
         elseif (SystematicType(i).eq.'ignore') then
C     Ignore error source called 'ignore'
            continue
         else
            TotalError = TotalError + Syst(i)**2
         endif
         if (SystematicType(i).eq.'uncor') then
C Uncorrelated error:
            UncorError = UncorError +  Syst(i)**2
         endif
         if (SystematicType(i).eq.'stat' 
     $        .or. SystematicType(i).eq.'Stat') then
C Stat error:
            StatError = StatError +  Syst(i)**2
         endif

         if (SystematicType(i).eq.'stat const') then
C Stat error:
            StatErrorConst = StatErrorConst +  Syst(i)**2
         endif
      enddo


      StatUncorError = sqrt(StatError+UncorError) ! construct stat+uncorr uncertainty
      StatError = sqrt(StatError)
      StatErrorConst = sqrt(StatErrorConst)
      UncorError = sqrt(UncorError)
      TotalError = sqrt(TotalError)


C Check total error
      if (TotalErrorRead.ne.0) then
         if ( abs(TotalError -TotalErrorRead)/TotalErrorRead
     $        .gt.0.01) then
            print 
     $'(''WARRNING IN READDATA, LARGE DEVIATION FOR TOTAL ERROR'')'
            print '(''Total calculated='',G10.4,'' READ='',G10.4)',
     $           TotalError,TotalErrorRead
         endif
      endif


C
C     Locate on the grid
C     Find the grid point for the correspondent bin:
C


      if (idxGrid.eq.0) then
         idxGrid = GetGridIdx(ClosestBin, idxReaction, ncolumnMax)
      endif

      if (IDebug.gt.0) then
         write(*,*) ' allbins = ',allbins(:,j)
         write(*,*) ' Closest = ',ClosestBin(1),ClosestBin(2)
         write(*,*) ' ======= Store data ======= ',
     &        idxGrid, idxReaction
      endif


C
C Now we want to see if we already have a point stored for the same reaction, grid index
C If yes, add one more point to the list for averaging. If not, create new point 
C
      idxData = GetDataIndex( 
     $     idxGrid,                 ! Closest bin-grid point 
     $     idxReaction)             ! reaction index


C
C STORE DATA:
C
      if (  iDataLast(idxData).eq.0 ) then

C Store extra point:
         iDataLast(idxData) = StoreData(idxData,  
     $        XSections(j), 
     $        StatUncorError,
     $        StatError,
     $        UncorError,
     $        TotalError,
     $        NUncert,
     $        Syst,
     $        SystematicType,
     $        NName)

      else 
         Call HF_ErrLog(20031304,
     $ 'W:Two points from the same experiment in the same grid point')
         if ( AveSameExp ) then
C Store extra point:
            iDataLast(idxData) = StoreData(idxData,  
     $           XSections(j), 
     $           StatUncorError,
     $           StatError,
     $           UncorError,
     $           TotalError,
     $           NUncert,
     $           Syst,
     $           SystematicType,
     $           NName)

         else
C Pre-average data point:
            Call PreAverageData(idxData, iDataLast(idxData),
     $           XSections(j), 
     $           StatUncorError,
     $           StatError,
     $           UncorError,
     $           TotalError,
     $           NUncert,
     $           Syst,
     $           SystematicType)            
         endif
      endif



 101  continue

C      print *,idxData,StatError

      enddo ! NData


      close (51)

      return

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 98   continue
      print '(''Error reading namelist Data'')'
      print *,CFile
      call HF_stop

 99   continue
      print '(''Can not open file '')'
      print *,CFile
      call HF_stop

 1017 continue
      print '(''Error reading file'')'
      call HF_stop
 1018 continue
      print '(''End of file while reading file'')'
      call HF_stop
 1019 continue
      print '(''Problem interpreting data line='',i6)',j
      call HF_stop

      end

C======================================================


C======================================================
      subroutine GetNearestSwimmingBinLinear
     &(idxReaction,NBD,bins,BN,res,idxGrid)
C
C Search for the nearest grid-bin step by step
C according to the columns in the dataset.

      implicit none
      include 'common.inc'
      integer iNCM
      integer idxReaction,idxGrid
      integer NBD
      double precision bins(NDimensionMax)
      character *(*) BN(NDimensionMax)           ! Bin names
      real res(NDimensionMax)
      
      double precision distMin,dist  !> smallest distance

      integer i,j, idxMin
      integer idxBin(NDimensionMax)

 !> function:
      integer GetBinIndexGrid
C------------------------------------------------------

C Check if number of bins in the file and in the grid are the same:
      if (NBD .ne. NDimensionGrid(idxReaction)) then
         Call hf_errlog(1,
     $'F:Number of dimensions in the grid and in the file do not match')
      endif

C Match bin-names with bin-indicies
      do i=1,NBD
         idxBin(i) = GetBinIndexGrid(bn(i),idxReaction)
         if (idxBin(i).eq.0) then
            Call hf_errlog(2,
     $       'F:GetNearestSwimmingBinLinear: Can not match bin-name')
         endif
      enddo

      idxMin = 0
      distMin = 1.0E20
      do i=1,NPointGrid(idxReaction)
         dist = 0.
         do j=1,NBD
            dist = max ( dist, abs(bins(j)
     $           -GridPoints(i,idxBin(j),idxReaction)))
         enddo
         if (dist .lt. distMin) then
            distMin = dist
            idxMin = i
         endif
      enddo

      if (idxMin.eq.0) then
         Call HF_ErrLog(3,
     $ 'F:GetNearestSwimmingBinLinear: Failed to find nearby bin')
      endif
      do i=1,NBD
         res(i) = GridPoints(idxMin,idxBin(i),idxReaction)
      enddo
      idxGrid = idxMin !> return exact location too.

      end subroutine

