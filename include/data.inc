
C Input tables:

      character*132 SystematicName(NSYSTMAX+NSYSTOMAX)   ! name of systematic errors

C     Way of precessing systematic errors: 
C     1x - Multiplicative, 2x - additive
C     x1 - symmetric, x2 asymmetric    
      integer SysForm(NSYSTMAX+NSYSTOMAX)   ! 
      

      real*8 SYSTAB(NSYSTMAX,NF2MAX,NMEASMAX)  ! Table of all syst. dependences
      real*8 SYSTABOrig(NSYSTMAX+NSYSTOMAX,NF2MAX,
     &                 NMEASMAX,NMaxShiftType) ! Table of all syst. dependences
      real*8 F2TAB(NF2MAX,NMEASMAX)            ! Table of all meas. values 
      real*8 F2TABOrig(NF2MAX,NMEASMAX)        ! Table of all original meas. values      
      real*8 F2ETAB(NF2MAX,NMEASMAX)           ! Table of all errors.
      real*8 F2ETABOrig(NF2MAX,NMEASMAX)       ! Table of all original errors.


C !> Indicies to grid, reaction:
      integer idxGridMeas(NF2MAX)
      integer idxProcessClass(NF2MAX)  !> type of "process class, which is a reaction, CME and other properties"


      integer idxReactionMeas(NProcClassMax)

      Common/CIdxMeas/idxGridMeas,idxReactionMeas
     $     ,idxProcessClass, NProcClass


C Save unorrelated and stat. errors separately:
      real F2ETAB_STA(NF2MAX,NMEASMAX) 
      real F2ETAB_UNC(NF2MAX,NMEASMAX) 
      real F2ETAB_STAORIG(NF2MAX,NMEASMAX)
      real F2ETAB_TOT(NF2MAX,NMEASMAX)

C Reference back to the data file:
      integer F2DataFile(NF2Max,NMeasMax)
      integer NMEASF2(NF2MAX)    ! number of measurements per point
      

      common/CHFSyst/SystematicName,
     $     SysForm

      Common/CF2/F2TAB,F2ETAB,F2ETABOrig,SYSTAB,
     $     SYSTABOrig,F2TABOrig,
     $     F2ETAB_STA,F2ETAB_UNC,F2ETAB_STAORIG,
     $     F2DataFile,
     $     F2ETAB_TOT


      real*8 SystDiag(NSystMax,NF2MAX)  ! Diagonalised syst. erorrs
      real*8 SystOrig(NSystMax,NF2MAX)  ! Original syst. erorrs
      real*8 OSyst(NF2MAX,NSYSTMAX) ! offset syst. erorrs
      real*8 SystImpact(NF2MAX,NSYSTMAX) ! impact of syst. erorrs
      real*8 StatToyMC(NF2MAX) ! ToyMC for stat. erorrs
      real*8 SystDiagPercent(NSystMax,NF2MAX)  ! Diagonalised syst. erorrs in %
      real*8 SystOrigPercent(NSystMax,NF2MAX)  ! Original syst. erorrs in %
      real*8 OSystPercent(NF2MAX,NSYSTMAX) ! offset syst. erorrs in %
      Common/CDiagSyst/SystDiag,SystOrig,OSyst,SystImpact,
     $ SystDiagPercent,SystOrigPercent,OSystPercent,StatToyMC


C Output:
      real*8 F2Vave(NF2MAX)
     $     ,F2Eave(NF2Max)
     $     ,F2EstAve(NF2Max)
     $     ,F2EstAveTrue(NF2Max)
     $     ,F2Esyave(NF2Max)
      real*8 SYSSH(NSYSTMAX)       ! Shifts of the systematics
      real*8 SYSSHItr(NSYSTMAX,0:NIterMax)    ! Shifts of the systematics for all iterations
      real*8 F2VaveOSyst(NF2MAX,NSYSTMAX) ! Averaged values with all offset systematics
      real*8 F2VaveSyst(NF2MAX,NSYSTMAX) ! Averaged values with all syst. shifts
      real*8 F2VaveToyMCStat(NF2MAX,nToyMC) ! Averaged values with all ToyMC systematics
      real Chi2Itr(0:NIterMax) ! Chi2 value for all iterations

      Common/CF2out/F2Vave,F2Eave,F2Estave,F2EstAveTrue,F2Esyave,
     $ syssh,SYSSHItr,F2VaveOSyst,F2VaveSyst,F2VaveToyMCStat,Chi2Itr

      real*8 ERRSYST(NSYSTMAX)   ! Reduction of the systematics
      Common/CF2outd/ERRSYST


      integer NPointGrid(NGridMax)   
      integer NDimensionGrid(NGridMax)
      character*32 GridReaction(NGridMax)
      real GridPoints(NGridPointMax,NDimensionMax,NGridMax)   ! values of grid points
      character*32 GridBinNames(NDimensionMax,NGridMax)

      Common/CGrids/       !> Common block to store grids
     $     NPointGrid,NDimensionGrid,GridPoints,
     $     GridReaction,GridBinNames


      Character*132 datasetlabel(NmeasMax)
      Character*132 datasetreaction(NmeasMax)
      character*32 datasetbinnames(NDimensionMax,NmeasMax)
      integer DATASETBinningDimension(NmeasMax)

      Common/datasetinfo/datasetlabel
     $     ,Datasetreaction
     $     ,Datasetbinningdimension
     $     ,datasetbinnames

      character*256 InputFileNames(NMEASMAX)
      
      character*128 
     $     GridType,                   !> Grid type: External or Auto
     $     GridFiles(NREACTMAX)        !> Grid file names
      
      Common/CFilesGrids/InputFileNames,GridFiles

      character*16 OutputPrefix
      character*16 OutputFolder      
      Common/CHeraAverager/OutputPrefix,OutputFolder
      
