      Subroutine InitAverager
C------------------------------------------------------
C Reading the data and preparing for the averaging
C------------------------------------------------------

      implicit none
      include 'common.inc'
      integer sysstatus
C Main steering file
      character*256 SteeringFileName   !> Main steering file name
      
C Basic namelists:
      character*32 OutputMode          !> Output type: original syst or ortogonal
      character*32 AverageType         !> Control bias corrections
      integer Iteration

C The temporary file name:
      character*256 cfile

      namelist/averager/OutputMode,
     $     OutputPrefix, OutputFolder, IDebug,
     $     WriteSysTexTable, PostRotateSyst, WriteOriginal

C     $     match_parameter,

      namelist/InFiles/NInputFiles,InputFileNames

      namelist/CommonGrid/GridType,GridFiles,AveSameExp
      
      namelist/BiasCorrection/AverageType,
     $ RescaleStatSep, CorrectStatBias,
     $ Iteration, FixStat

      namelist/Systematics/NSysTot

C Temporary indices:
      integer i,j,k, idxreaction

C Functions:
      integer iargc
      integer getreactionidx
      integer itemp

C CME:
      double precision sqrtcme
      character*5 CInfoMes

C------------------------------------------------------
      if (iargc().ne.1) then
         print *,'Program usage: averager steering_file_name'
         stop
      endif
C
C OPEN STEERING FILE:
C
      Call GetArg(1,SteeringFileName)
      open (51,file=SteeringFileName,status='old',err=51)
C
C DEFINE DEFAULT VALUES OF THE STEERING PARAMETERS:
C
      GridType = 'External'
      match_parameter = 1.D-5
      do i=1,NREACTMAX
         GridFiles(i) = ''
      enddo
      AverageType = 'MIXED'
      RescaleStatSep = .false.
      CorrectStatBias = .false.
      FixStat = .false.
      OutputMode  = 'ORTH'
      OutputPrefix = 'Ave'
      OutputFolder = '../output/'
C
C DEFINE DEFAULT VALUES OF THE COMMON PARAMETERS:
C
      NName = 0
      NMeas = 0
      

      IDebug = 0
      WriteOriginal = .false.
      WriteSysTexTable = .false.

      PostRotateSyst = .true. ! by default keep syst. align to the original sources as much as possible.

C
C READ THE NUMBER AND NAMES OF DATAFILES FROM THE STEERING FILE:
C
      rewind (51)
      read (51,nml=InFiles,end=54,err=56)
      
      write(*,*) 'Number of input files: ', NInputFiles
      if(NInputFiles.gt.NMEASMAX)then
          call hf_errlog(10,
     $        'F:number of input files exceed maximum'
     $         //'Increase NMEASMAX in settings.inc and rebuild')
      endif
      
C
C Read average parameters:
C
      read (51,nml=CommonGrid,end=57,err=58)
      write(*,*) 'GridType = ',GridType

      read (51,nml=averager,end=52,err=53)

      read (51,nml=BiasCorrection,end=59,err=60)


C Method for stat. uncertainty correction:
      if (RescaleStatSep.and.FixStat) then
          call hf_errlog(10,
     $        'F:inconsistent steering: '
     $         //'Both RescaleStatSep and FixStat are TRUE')
      endif


C
C Define the type of systematic output
C
      Call DecodeOutputMode(OutputMode)

C
C READ GRIDS:
C
      NGrid = 0
      if (GridType.eq.'External') then
         Call ReadGrids
      elseif (GridType.eq.'Auto') then
         do i = 1, NInputFiles
           cfile = InputFileNames(i)
           Call ReadGridFromDataFile(cfile)
         enddo
      endif

C READ THE PART WITH SYSTEMATICES:
      rewind (51)
      read (51,nml=Systematics,end=74,err=55)
 74   continue     !> It is Ok if the namelist is not present


C
C read datafiles:
C
      do i = 1, NInputFiles
        cfile = InputFileNames(i)
        Call ReadXSTable(cfile)
      enddo

C
C Info. message:
C
      write (CInfoMes,'(I5)') NMeas
      Call hf_errlog(20031302,'I:Read '//CInfoMes//' unique points')



C
C Define the additive, multiplicative or other
C Call function when all systematics are already stored
C
      Call DecodeAverageType(AverageType)

C Store number of iterations
      NIteration = Iteration
      write(*,*) 'NIteration = ',NIteration
      if(NIteration.gt.NIterMax)then
         call hf_errlog(10081601,
     $      'F:Number of iterations exceed maximum.'//
     $      ' Increase NIterMax in settings.inc')
      endif

C     Create output directory
      CALL system("mkdir -p "//trim(OutputFolder))

C------------------------------------------------------
      Return
 51   continue     

      call hf_errlog(1,'F:Can not open steering file '
     $     //trim(SteeringFileName))

 52   continue
      call hf_errlog(2,'F:Can not find averager namelist')

 53   continue
      call hf_errlog(3,'F:Error reading averager namelist')

 55   continue
      call hf_errlog(4,'F:Error reading Systematics namelist')

 54   continue
      call hf_errlog(5,'F:Can not find InFiles namelist')

 56   continue
      call hf_errlog(6,'F:Error reading InFiles namelist')

 57   continue
      call hf_errlog(5,'F:Can not find CommonGrid namelist')

 58   continue
      call hf_errlog(6,'F:Error reading CommonGrid namelist')


 59   continue
      call hf_errlog(5,'F:Can not find BiasCorrection namelist')

 60   continue
      call hf_errlog(6,'F:Error reading BiasCorrection namelist')
      end

C-------------------------------------------------------------
C Decode AverageType string
C Rewrite systematics type with respect of steering file
C Add - add additive, MULT - all multiplicative
C offset systematics are not affected
C------------------------------------------------------------

      subroutine DecodeAverageType(AverageType)

      implicit none
      include 'common.inc'
      character*(*) AverageType
      integer k
C------------------------------------------------------------
      if (AverageType .eq. 'ADD') then
        write(*,*) 'Uncertainties treated as Additive'
        do k  = 1, NSysTot
          SystematicForm(k) = 'ADD'
        enddo
      elseif (AverageType .eq. 'MULT') then
        write(*,*) 'Uncertainties treated as Multiplicative'
        do k  = 1, NSysTot
          SystematicForm(k) = 'MULT'
        enddo
      elseif (AverageType .eq. 'MIXED') then
        write(*,*) 'Uncertainties treatment is read from the datasets'
      else
         call hf_errlog(1,
     $        'F:Unknown averaging type read from averager namelist '
     $        //AverageType)
      endif
C------------------------------------------------------------
      end


      subroutine DecodeOutputMode(OutputMode)
C-------------------------------------------------------------
C Determines the type of output
C-------------------------------------------------------------
      implicit none
      include 'common.inc'
      character*(*) OutputMode
C------------------------------------------------------------
      if (OutputMode .eq. 'ORIG') then
        iOutput = 0
      elseif (OutputMode .eq. 'ORTH') then
        iOutput = 1
      else
         call hf_errlog(1,'F:Unknown OutputMode read from namelist '
     $        //OutputMode)
      endif
C------------------------------------------------------------
      end


      subroutine ReadGridFromDataFile(CFile)
C------------------------------------------------------------
C Read grid-bins from datafiles
C------------------------------------------------------------
      implicit none
      integer ifile
C------------------------------------------------------------

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
      logical Percent(1:nsystMax)

      integer IndexDataset
      double precision SystScales(nsystMax)


C Namelist definition:
      namelist/Data/Name,NData
     $     ,Reaction,Percent
     $     ,SystScales, IndexDataset
     $     ,ColumnName, ColumnType, NColumn

      real AllBins(NDimensionMax,ndataMax)

      integer i,j,iBin,iError, iGrid, itemp, k

C Temporary buffer to read the data (allows for comments starting with *)
      character *4096 CTmp
      
      real ClosestBin(ncolumnMax)  ! The correspondent bin in the bingrid

      integer idxreaction,idxGrid
C Functions:
      integer getreactionidx, getgrididx
      
C------------------------------------------------------------

C Reset to default:
      NUncert = 0
      NData = 0
      NBinDimension = 0
      Reaction = ' '
      Name     = ' '
      IndexDataSet = 0


C Reset scales to 1.0
      do i=1,nsystmax
         SystScales(i) = 1.0
         ColumnType(i) = ' '
         ColumnName(i) = ' '
      enddo

      open(51,file=CFile,status='old',err=99)

      print *,'Reading data file for getting grid points...'
      print *,CFile
      read(51,NML=Data,err=98)
      
C     
C Check if the grid is given for existing reaction, if yes: update
C
      do j=1,NGrid
        if (Reaction.eq.GridReaction(j)) then
          iGrid = j
          goto 666
        endif
      enddo
C Add new grid
      NGrid = NGrid + 1
C Add a title of a new readed grid
      GridReaction(NGrid) = Reaction
      write(*,*) ' GRID subroutine NGRID: ', NGrid
C
      iGrid = NGrid
 666  continue

C 
C Reaction index:
C 
      idxReaction = GetReactionIdx(reaction)


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

C Parse ColumnTypes
      do i=1,NColumn
         if (ColumnType(i).eq.'Bin') then
            NBinDimension = NBinDimension + 1
            if(NBinDimension.gt.NDimensionMax)then
                call hf_errlog(10,
     $        'F:Number of dimensions exceed maximum'
     $         //'Increase NDimensionMax in settings.inc')
            endif
            BinName(NBinDimension) = ColumnName(i)
         endif
      enddo
      NDimensionGrid(iGrid) = NBinDimension

C Read data info:

      do j=1,NData
C Allow for comments:
 89      read (51,'(A)',err=1017,end=1018) ctmp
         if (ctmp(1:1).eq.'*') then
C Comment line, read another one
            goto 89
         endif

C Read the columns
         read (ctmp,*,err=1019)(buffer(i),i=1,NColumn)

C Decode the columns
         iBin   = 0
         iError = 0
         do i=1,NColumn
            if (ColumnType(i).eq.'Bin') then
               iBin = iBin + 1
C               write(*,*) 'iBin= ', iBin, ' buffer = ', buffer(i)
               allbins(iBin,j) = buffer(i)
            endif
         enddo

C
C STORE THE GRIDBINS
C
         idxGrid = 0
         idxGrid = GetGridIdx( allbins(:,j), idxReaction, ndimensionmax)

      if (idxGrid.gt.0) then ! check that this bin coincides with some previous one
c	  write(*,*) idxGrid
c	  write(*,*) 'Such bin already exist'
	  goto 555
      else ! if this bin is a new bin
	  NPointGrid(idxReaction) = NPointGrid(idxReaction) + 1
          if ( NPointGrid(idxReaction).gt.NGridPointMax) then
             call hf_errlog(2016053101,
     $            'F:Too many grid points,' 
     $            // ' increase NGridPointMax in steering.inc')
          endif
	  itemp = NPointGrid(idxReaction)
	  do k = 1, iBin
	    GridPoints(itemp,k,iGrid) = allbins(k,j) ! Store bins for given line of the file
C	    WRITE(*,*) ' GP: ',GridPoints(itemp,k,iGrid)
	    GridBinNames(k,iGrid) = BinName(k) ! Store bin names
	  enddo
      endif
 555    continue
      enddo

      close (51)

      
      return

 98   continue

 99   continue
      call hf_errlog(1,'F:Error reading namelist Data from file '
     $     //Trim(CFile))

 1017 continue
      call hf_errlog(2,'F:Error reading file '
     $     //Trim(CFile))

 1018 continue
      call hf_errlog(3,'F:End of file while reading file')

 1019 continue
      print '(''Problem interpreting data line='',i6)',j
      call HF_stop

C------------------------------------------------------------
      end subroutine


      subroutine ReadGrids
C------------------------------------------------------------
C Read grid files
C------------------------------------------------------------
      implicit none
      include 'common.inc'
      integer i,j,k,iGrid

      character*32 Reaction
      integer  NDimension
     $     ,   NPoints
      character*32 BinNames(NDimensionMax)
      namelist/Grid/Reaction,NDimension,NPoints, BinNames

C------------------------------------------------------------

C Initialize:

      NGrid = 0
      do i=1,NREACTMAX
         NPointGrid(i) = 0
      enddo

      do i=1,NREACTMAX
         if (GridFiles(i).ne.'') then
            print *,'Read grid file'
            print *,GridFiles(i)
            open (52,file=GridFiles(i),err=91)
            read (52,nml=Grid,err=92,end=93)
            
            write(*,*) ' GRID FUNCTION: ', Reaction
C Check if the grid is given for existing reaction, if yes: update
            do j=1,NGrid
               if (Reaction.eq.GridReaction(j)) then
                  iGrid = j
                  goto 61
               endif
            enddo
C Add new grid
            NGrid = NGrid + 1

            if (NGrid.gt.NGridMax) then
               print *,'Too many grids=',NGrid
               call hf_errlog(1,
     $          'F:Too many grid files. Increas NGridMax')
            endif

C Add a title of a new readed grid
            GridReaction(NGrid) = Reaction
            write(*,*) ' NGRID: ', NGrid
C
            iGrid = NGrid
 61         continue
C Store
            if (NPoints .gt. NGridPointMax) then
               print *,'Npoints = ',NPoints
               call hf_errlog(2,
     $   'F:Too many grid points. Increase NGridPointMax')
            endif

            NPointGrid(iGrid) = NPointGrid(iGrid) + NPoints
            write(*,*) 'NPointGrid = ', NPointGrid(iGrid)
            NDimensionGrid(iGrid) = NDimension
            do k=1,NDimension
               GridBinNames(k,iGrid) = BinNames(k)
               write(*,*) 'GridBinNames: ',k,iGrid,BinNames(k)
            enddo

C Read the file
            do j=1,NPoints
               read (52,*,err=94,end=95)
     $              (GridPoints(j,k,iGrid),k=1,NDimension)
            enddo
            close (52)
         else
            goto 51
         endif
      enddo
 51   continue
      return
 91   continue
      call hf_errlog(3,'F:Error opening grid file '//
     $     trim(GridFiles(i)))

 92   continue
      call hf_errlog(4,'F:Error reading Grid namelist in file '//
     $     trim(GridFiles(i)))

 93   continue
      call hf_errlog(5,'F:Did not find Grid namelist in file '//
     $     trim(GridFiles(i)))

 94   continue
      print *,'Error reading Grid file'
      print *,GridFiles(i)
      print *,'Line=',j
      print *,'Stop'
      stop
 95   continue
      print *,'Line=',j,' Expected length=',NPoints
      print *,'Stop'
      call hf_errlog(6,'F:Unexpected (early) end of grid file '//
     $     trim(GridFiles(i)))
C------------------------------------------------------------
      end


