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

      contains

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


C     Set binning
      subroutine SetBinning(binIn,binNameIn,nBinIn,nDimIn)

          implicit none
          character *(*) binNameIn(nDimIn)
          real*4 binIn(nBinIn,nDimIn)
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

      if (nBinIn .gt. NGridPointMax) then
          print *,'nBinIn = ',nBinIn
          call hf_errlog(2,
     $   'F:Too many grid points. Increase NGridPointMax')
      endif

C     Loop over data points and fill grid
      do j=1,nBinIn
         idxGridMeas(j) = j
         idxProcessClass(j) = 1
         GridPoints(j,:,1) = binIn(j,:)
      enddo
      end subroutine SetBinning

      end module AvIn

      subroutine getpulls(pulls, nPulls)
      implicit none

      integer nPulls, k
      real*8 pulls(nPulls)

cf2py intent(in) nPulls
cf2py intent(out) pulls

      include 'common.inc'
      do k=1,nsysTot
         pulls(k) = SYSSH(k) /
     $         sqrt(1 - errsyst(k)*errsyst(k))
      enddo
      end subroutine getpulls

      subroutine getdatapulls(pulls, nPulls, nFiles)
      implicit none

      integer nPulls, nFiles
      integer iFile, if2, iexp, isys
      real*8 sum
      real*8 pulls(nFiles, nPulls)

cf2py intent(in) nPulls
cf2py intent(in) nFiles
cf2py intent(out) pulls

      include 'common.inc'
      do iFile=1,NInputFiles
         do if2=1,NMeas
            do iexp=1,NMeasF2(if2) !NMeas
               if (F2DataFile(if2,iexp).eq.IFile) then
                  sum = F2TAB(if2,iexp)
                  do isys=1,NSYSTOT
                     sum = sum + SYSTAB(isys,if2,iexp)*SYSSH(isys)
                  enddo

                  if (NMeasF2(if2).gt.1) then
                     pulls(iFile,if2) = (F2VAVE(if2)-sum)/
     $                    sqrt(abs(F2ETAB(if2,iexp)**2
     $                    -F2EstAve(if2)**2))
                  else
                     pulls(iFile,if2) = 0.
                  endif

               endif
            enddo
         enddo
      enddo

      end subroutine getdatapulls

      subroutine getchi2(chi2, ndof)

      real*8 chi2, ndof
cf2py intent(out) chi2
cf2py intent(out) ndof
      include 'common.inc'
      call CalcChi2(chi2, ndof)

      end subroutine getchi2

      subroutine getshiftsyst(shiftsyst, nshiftsyst)
      implicit none

      integer nshiftsyst, k
      real*8 shiftsyst(nshiftsyst)

cf2py intent(in) nshiftsyst
cf2py intent(out) shiftsyst

      include 'common.inc'
      do k=1,nsysTot
            shiftsyst(k) = SYSSH(k)
      enddo
      end subroutine getshiftsyst



      subroutine gettoystat(toystat, nData)
      implicit none

      integer nData
      integer i
      real*8 toystat(nData)

cf2py intent(in) nData
cf2py intent(out) toystat

      include 'common.inc'
      do i=1,NMeas
         if(nToyMC .gt. 0)then
            toystat(i) = StatToyMC(i)
         else
            toystat(i) = 0
         endif
      enddo

      end subroutine gettoystat


      subroutine getsysimpact(sysimpact, nData, nSyst)
      implicit none

      integer nData, nSyst
      integer i, k
      real*8 sysimpact(nSyst, nData)

cf2py intent(in) nData
cf2py intent(in) nSyst
cf2py intent(out) sysimpact

      include 'common.inc'
      do i=1,NMeas
         do k=1,nsysTot
            if(doSystImpact)then
               sysimpact(k,i) =
     $    (F2VaveSyst(i,2*k)-F2VaveSyst(i,2*k-1))*0.5
            else
               sysimpact(k,i) = 0
            endif
         enddo
      enddo

      end subroutine getsysimpact



C     Averaging
      subroutine average( dataIn, statIn, systIn, snameIn, fnameIn,
     $     dataOut, statOut, systOut,
     $     nmeasIn, nDataIn,nSystIn, nSnameIn, nFnameIn)

      use AvIn
      implicit none

      integer nmeasIn, nDataIn, nSystIn, nSnameIn, nFnameIn
      character *(*) snameIn(nSnameIn)
      character *(*) fnameIn(nFnameIn)
      real*8 dataIn(ndataIn,nmeasIn), statin(ndataIn,nmeasIn)
      real*8 systIn(nSystIn,ndataIn,nmeasIn)
      real*8 dataOut(ndataIn), statOut(ndataIn)
     $     , systOut(nsystIn,ndataIn)
C      character*8 ctmp
C python helper:

Cf2py intent(in) dataIn 
Cf2py intent(in) statIn 
Cf2py intent(in) systIn 
Cf2py intent(in) snameIn
Cf2py intent(in) fnameIn

Cf2py intent(in) nmeasIn
cf2py intent(in) nDataIn 
cf2py intent(in) nSystIn
cf2py intent(in) nSnameIn
cf2py intent(in) nFnameIn

cf2py intent(out) dataOut
cf2py intent(out) statOut
cf2py intent(out) systOut


      include 'common.inc'
      integer i,j,k

      integer ::  Inidxsys(nSystIn)
      integer ::  Insystype(nSystIn)
C     integer iFile,iF2,iexp, isys

C     Print size of the input information
      print *,'Measured points          ',nDataIn
      print *,'Data samples             ',nmeasIn
      print *,'sources of uncertainties ',nSystIn
C-------------------------------------------------------------------

C     Initialize default values
      if(init.ne.777)then
         Call initVariables
      endif

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
      NSysTot       = 0

C     Check if the systematic names were already given
C     If not, give default names


      if(NProcClass .ne. 1)then
C     Fill dummy grid information (always 1D grid)
      NProcClass    = 1
      idxReactionMeas(1) = 1
      gridreaction(1) = 'Bla'
      NDimensionGrid(1) = 1
      GridBinNames(1,1) = 'Y'

      if (NMeas .gt. NGridPointMax) then
          print *,'NMeas = ',NMeas
          call hf_errlog(2,
     $   'F:Too many grid points. Increase NGridPointMax')
      endif
C     Loop over data points and fill grid
      do j=1,NMeas
         idxGridMeas(j) = j
         idxProcessClass(j) = 1
         GridPoints(j,1,1) = (j-1)
      enddo
      endif

C     Add systematics

      do k=1,nSnameIn
          call AddSystematics(snameIn(k), k, Inidxsys, Insystype)
      enddo

C     Fill input filenames
      do j=1,nmeasIn
         InputFileNames(j)= fnameIn(j)
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
      do i=1,NMeas
         dataOut(i) = f2vave(i)
         statOut(i) = F2EstAve(i)
         do k=1,nsysTot
            systOut(k,i) =  SystDiag(k,i)
         enddo
      enddo

C--------------------------------------------------------------------
      end subroutine average

