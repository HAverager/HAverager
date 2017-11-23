C
C Include file for averager. Main parameters and scalar variables.
C
      double precision match_parameter ! grid zapping relative uncertainty
      integer NDimensionMax ! Max. number of bins
      integer nmeasmax  ! Max. number of measurements (H1 CC, NC etc)
      integer nf2max    ! max. number of F2 points
      integer nsystmax  ! max. number of systematic uncertainties
      integer nsystomax  ! max. number of systematic uncertainties
      integer NREACTMAX ! max. number of types of measurements
      integer NMaxShiftType ! Max. number of type of shifts: Symmetric,Up,Down

C Sizes:
      parameter (NMEASMAX=20,NF2MAX=1000,
     $     NSYSTMAX=1000,NSYSTOMAX=10,NREACTMAX=10)
      parameter (NDimensionMax=10)
      parameter (NMaxShiftType=2)   !> Max. number of type of shifts: Symmetric,Up,Down

!> Grid information
      integer NGridMax     !> Total number of grids
      parameter (NGridMax = NREACTMAX)

      integer NGridPointMax
      parameter (NGridPointMax=800) !> Max number of points per grid file
      integer NGrid        !> Number of grids
      common/CGridPars/NGrid

C     Max number of iterations
      integer NIterMax
      parameter (NIterMax=100)     

C     Number of datasets:      
      integer ndatasets
      real SysShTolerance
      parameter (SysShTolerance=0.05) 
      common/Cdatasetpar/NDatasets

      integer NInputFiles  ! number of input files
      common/CNInputFiles/NInputFiles

C     Number of ToyMC tests for stat. error
      integer nToyMC
      parameter (nToyMC=1000) 

C Namelist  variables:
      integer ndataMax,ncolumnMax
      parameter (ndataMax=1000)  !> Max. number of entries per data file
      parameter (ncolumnMax = nsystMax+NDimensionMax+1) !> Max. number of columns



C !> Process classes: 
      integer NProcClassMax
      parameter (NProcClassMax = 50)
      integer NProcClass

C     Variables for offset systematics
      integer NSYSOTOT ! Number of offset systematics
      Common/CHFParameters/NSYSOTOT



C Number of common Bins (=NDiag):
      integer NumCommonBins

      integer NName ! number of files
      integer NIteration, iItr ! number of iteration in case of multiplicative error treatment
      integer iOutput  ! OutputType


C-
      logical PostRotateSyst  ! rotate syst. after combination along the original sources
      Common/CPostRot/ PostRotateSyst


C Control:
      integer IDEBUG
      common/cdebug/idebug

      logical WriteOriginal
      logical WriteSysTexTable
      common/CWrite1/WriteOriginal
      common/CWrite2/WriteSysTexTable

C-----------------------------------------------

      integer NMEAS,NF2TOT,NSYSTOT, NSFSEC

C
C Sqrt scaling for stat. errors:
C       
      logical RescaleStatSep,CorrectStatBias,FixStat

C Average control:
      logical AveSameExp
      logical UseBlas

      Common/CF2Pars2/NMEAS,NF2TOT,NSYSTOT
     $     ,NSFSEC 
     $     ,RescaleStatSep
     $     ,CorrectStatBias
     $     ,FixStat
     $     ,AveSameExp
     $     ,NName
     $     ,NIteration
     $     , iItr
     $     ,NumCommonBins
     $     ,NMEASF2
     $     ,iOutput
     $     ,UseBlas

C Parameters:
      Common/CPARS/Match_parameter

