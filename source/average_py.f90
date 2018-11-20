C     Module with input paratemetrs
      module AvIn
      integer init
      integer InDEBUG
      logical InWriteOriginal
      logical InWriteSysTexTable
      logical InPostRotateSyst

      integer InIteration
      integer InNToyMC
      logical InRescaleStatSep
      logical InCorrectStatBias
      logical InFixStat
      logical InUseBlas
      logical IndoSystImpact
      character(len=128), allocatable :: Insname(:)
      character(len=128), allocatable :: Infname(:)
      integer, allocatable ::  Inidxsys(:)
      integer, allocatable ::  Insystype(:)

      contains

C     Deallocate all vars
      subroutine cleanInVars
      if ( allocated( Insname )) Deallocate ( Insname )
      if ( allocated( Infname )) Deallocate ( Infname )
      if ( allocated( Inidxsys )) Deallocate ( Inidxsys )
      if ( allocated( Insystype )) Deallocate ( Insystype )
      end subroutine cleanInVars

C     Initialization of steerable parameters with default values
      subroutine initVariables
          implicit none
          include 'common.inc'
          InDebug = 0
          InWriteOriginal = .false.
          InWriteSysTexTable = .false.
          InPostRotateSyst = .false.

          iOutput = 1
          OutputPrefix = 'Ave'
          OutputFolder = 'output'

          InIteration = 10
          InNToyMC = 0
          InRescaleStatSep = .false.
          InCorrectStatBias = .false.
          InFixStat = .false.
          InUseBlas = .false.
          IndoSystImpact = .false.
          init = 777
      end subroutine initVariables

C     Set a name of the output folder
      subroutine SetOutputFolder(datain)
        character (*) datain
Cf2py intent(in) dataIn
        include 'common.inc'
        OutputFolder = datain
      end subroutine SetOutputFolder

C     Set a name of the output prefix
      subroutine SetOutputPrefix(datain)
        character (*) datain
Cf2py intent(in) dataIn
        include 'common.inc'
        OutputPrefix = datain
      end subroutine SetOutputPrefix

C     Set names of the systematic uncertainties
      subroutine SetSNames(snameIn,nSnameIn)
       character *(*) snameIn(nSnameIn)
Cf2py intent(in) snameIn
cf2py intent(in) nSnameIn
        include 'common.inc'
        allocate(Insname(nSnameIn))
        allocate(Inidxsys(nSnameIn))
        allocate(Insystype(nSnameIn))
        do k=1,nSnameIn
          Insname(k)=snameIn(k)
          call AddSystematics(snameIn(k), k, Inidxsys, Insystype)
        enddo
      end subroutine SetSNames

C     Set names of input files
      subroutine SetFNames(fnameIn,nFnameIn)
       character *(*) fnameIn(nFnameIn)
Cf2py intent(in) fnameIn
cf2py intent(in) nFnameIn
        include 'common.inc'
        allocate(Infname(nFnameIn))
        do k=1,nFnameIn
          Infname(k)=fnameIn(k)
        enddo
      end subroutine SetFNames

C     Set binning
      subroutine SetBinning(binIn,binNameIn,nBinIn,nDimIn)

          implicit none
          character *(*) binNameIn(nDimIn)
          real*8 binIn(nBinIn,nDimIn)
          integer nBinIn, nDimIn, j

C python helper:

Cf2py intent(in) binIn
Cf2py intent(in) binNameIn

Cf2py intent(in) nBinIn
cf2py intent(in) nDimIn

        include 'common.inc'

C     Fill dummy grid information (always 1D grid)
      NProcClass    = 1
      idxReactionMeas(1) = 1
      gridreaction(1) = 'Bla'
      NDimensionGrid(1) = nDimIn
      GridBinNames(:,1) = binNameIn(:)


C     Loop over data points and fill grid
      do j=1,nBinIn
         idxGridMeas(j) = j
         idxProcessClass(j) = 1
         GridPoints(j,:,1) = binIn(j,:)
      enddo
      end subroutine SetBinning

      end module AvIn

C     Mudule with output variables
      module AvOut
      real*8,allocatable :: pullsyst(:)
      real*8,allocatable :: shiftsyst(:)
      real*8,allocatable :: squeezesyst(:)
      real*8,allocatable :: pulldata(:,:)
      real*8,allocatable :: sysimpact(:,:)
      real*8,allocatable :: toystat(:)
      real*8 :: chi2
      integer :: ndof
      end module AvOut

      subroutine cleanOutVars
      use AvOut
      if (allocated(pullsyst)) deallocate(pullsyst)
      if (allocated(shiftsyst)) deallocate(shiftsyst)
      if (allocated(squeezesyst)) deallocate(squeezesyst)
      if (allocated(pulldata)) deallocate(pulldata)
      if (allocated(sysimpact)) deallocate(sysimpact)
      if (allocated(toystat)) deallocate(toystat)
      end subroutine cleanOutVars

C     Averaging
      subroutine average( dataIn, statIn, systIn,
     $     dataOut, statOut, systOut,
     $     nmeasIn, nDataIn,nSystIn)

      use AvIn
      use AvOut
      implicit none

      integer nmeasIn, nDataIn, nSystIn
      real*8 dataIn(ndataIn,nmeasIn), statin(ndataIn,nmeasIn)
      real*8 systIn(nSystIn,ndataIn,nmeasIn)
      real*8 dataOut(ndataIn), statOut(ndataIn)
     $     , systOut(nsystIn,ndataIn)
      character*8 ctmp
C python helper:

Cf2py intent(in) dataIn 
Cf2py intent(in) statIn 
Cf2py intent(in) systIn 

Cf2py intent(in) nmeasIn
cf2py intent(in) nDataIn 
cf2py intent(in) nSystIn

cf2py intent(out) dataOut
cf2py intent(out) statOut
cf2py intent(out) systOut


      include 'common.inc'
      integer i,j,k

      integer iFile,iF2,iexp, isys
      double precision sum

C     Print size of the input information
      print *,'Measured points          ',nDataIn
      print *,'Data samples             ',nmeasIn
      print *,'sources of uncertainties ',nSystIn
C-------------------------------------------------------------------

C     Initialize default values
      if(init.ne.777)then
         Call initVariables
      endif

C     This ensures that the code can be called twice
      call cleanOutVars

C     Fill input parameters
      IDebug = InDebug
      WriteOriginal = InWriteOriginal
      WriteSysTexTable = InWriteSysTexTable
      PostRotateSyst = InPostRotateSyst
      NIteration = InIteration

      RescaleStatSep = InRescaleStatSep
      CorrectStatBias = InCorrectStatBias
      FixStat = InFixStat
      UseBlas = InUseBlas
      doSystImpact = IndoSystImpact
      NToyMC = InNToyMC

C     Create output directory
      CALL system("mkdir -p "//trim(OutputFolder))

C     Print initial variables
      print *,'Debug:                   ',IDebug
      print *,'Number of iterations:    ',NIteration
      print *,'Output mode:             ',iOutput,' Orthogonal'
      print *,'WriteOriginal:           ',WriteOriginal
      print *,'WriteSysTexTable:        ',WriteSysTexTable
      print *,'PostRotateSyst:          ',PostRotateSyst
      print *,'CorrectStatBias:         ',CorrectStatBias
      print *,'RescaleStatSep:          ',RescaleStatSep
      print *,'FixStat:                 ',FixStat
      print *,'UseBlas:                 ',UseBlas
      print *,'Check syst impact:       ',doSystImpact
      print *,'Number of Toy MC:        ',NToyMC
      print *,'Output folder:           ',OutputFolder



C     Fill internal arrays and variables
      NInputFiles   = nmeasIn
      NMeas         = nDataIn

C     Check if the systematic names were already given
C     If not, give default names

      if(NSysTot.eq.0  .or. (.not. allocated(InsName)  ) ) then
         print *,'Set dummy names for systematics:'
         NSysTot = nSystIn
         do k=1,nsystot
            if(len(trim(SystematicName(k))).eq.132)then
               write (ctmp,'(''syst'',i0)') k
               SystematicName(k) = ctmp
               SysForm(k) = 11
            endif
         enddo
         call SetSNames(SystematicName,nsystot)
      endif


      if(NProcClass .ne. 1)then
C     Fill dummy grid information (always 1D grid)
      NProcClass    = 1
      idxReactionMeas(1) = 1
      gridreaction(1) = 'Bla'
      NDimensionGrid(1) = 1
      GridBinNames(1,1) = 'Y'

C     Loop over data points and fill grid
      do j=1,NMeas
         idxGridMeas(j) = j
         idxProcessClass(j) = 1
         GridPoints(j,1,1) = (j-1)
      enddo
      endif

C     Fill input filenames
      do j=1,nmeasIn
         if(allocated(InfName))then
             InputFileNames(j)= Infname(j)
         else
             write (ctmp,'(''File'',i0)') j
             InputFileNames(j)=ctmp
         endif
      enddo



      
C     Fill data and uncertainties
C     Loop over measurements
      do i=1,NMeas
         NMeasF2(i) = 0

C        Loop over data points
         do j=1,nmeasIn
            if(datain(i,j).ne.0) then
                call StoreData(i,
     $           datain(i,j),
     $           statin(i,j),
     $           statin(i,j),
     $           0.,
     $           0.,
     $           nSystIn,
     $           systin(:,i,j),
     $           Inidxsys, Insystype,
     $           j)

            endif
         enddo
      enddo
      
C     Perform the averaging:
      call Averaging

C     Print output
      Call Output

C     Fill output information
      allocate(pullsyst(nsysTot))
      allocate(shiftsyst(nsysTot))
      allocate(squeezesyst(nsysTot))
      allocate(pulldata(NInputFiles,NMeas))
      allocate(sysimpact(nsysTot,NMeas))
      allocate(toystat(NMeas))

      do i=1,NMeas
         dataOut(i) = f2vave(i)
         statOut(i) = F2EstAve(i)
         if(nToyMC .gt. 0)then
            toystat(i) = StatToyMC(i)
         else
            toystat(i) = 0
         endif
         do k=1,nsysTot
            systOut(k,i) =  SystDiag(k,i)
            pullsyst(k) = SYSSH(k) /
     $            sqrt(1 - errsyst(k)*errsyst(k))
            shiftsyst(k) = SYSSH(k)
            squeezesyst(k) = errsyst(k)
            if(doSystImpact)then
               sysimpact(k,i) = 
     $    (F2VaveSyst(i,2*k)-F2VaveSyst(i,2*k-1))*0.5
            else
               sysimpact(k,i) = 0
            endif
         enddo
      enddo

C     Fill data pulls
      do iFile=1,NInputFiles
         do if2=1,NMeas
            do iexp=1,NMeasF2(if2) !NMeas
               if (F2DataFile(if2,iexp).eq.IFile) then
                  sum = F2TAB(if2,iexp)
                  do isys=1,NSYSTOT
                     sum = sum + SYSTAB(isys,if2,iexp)*SYSSH(isys)
                  enddo

                  if (NMeasF2(if2).gt.1) then
                     pulldata(iFile,if2) = (F2VAVE(if2)-sum)/
     $                    sqrt(abs(F2ETAB(if2,iexp)**2
     $                    -F2EstAve(if2)**2))
                  else
                     pulldata(iFile,if2) = 0.
                  endif

               endif
            enddo
         enddo
      enddo

C     Fill Chi^2 and NDoF
      call CalcChi2(chi2, ndof)

C--------------------------------------------------------------------
      end subroutine average

