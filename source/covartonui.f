      Subroutine covartonui( doCovconversion )
C-----------------------------------------------
C
C Convert covariance matrix to nuisance param. representation
C
C-----------------------------------------------
      implicit none
C------------------------------------------------------------------
      logical doCovconversion
      integer narg
      integer command_argument_count
      character*80 command

C------------------------------------------------------------------      
      narg = command_argument_count()
      if (narg.gt.0) then
         call get_command_argument(1,command)
         if (index(command,'--convert-covar').ne.0) then
            doCovconversion = .true.
            if (index(command,'=').ne.0) then
               command = command(index(command,'=')+1:len(command))
            else
               command = 'covar.in'
            endif
            call CovMatrixConverter(command)
         else
            doCovconversion = .false.
         endif
      endif
      end

      subroutine CovMatrixConverter(fileName)
      implicit none
      character*(*) fileName
      integer ncovar
      double precision tolerance
      logical SeparateDiag
      logical Verbose
      namelist/Covar/NCovar,Tolerance, SeparateDiag, Verbose
      integer i,j,k, ncorr
      double precision test,devmax,dev
      double precision, allocatable :: cov(:,:)
      double precision, allocatable :: corr_syst(:,:)
      double precision, allocatable :: alpha(:)
      double precision, allocatable :: beta(:,:)

C----------------------------------------

C      call hf_errlog(2,'I:Read covariance matrix from file')

C init
      NCovar = 0
      Tolerance = 0.
      SeparateDiag = .true.
      Verbose = .false.

      print *,'Read namelist ...'

      open (51,file=trim(FileName),status='old',err=91)
      read (51,nml=Covar,err=92,end=93)

C --- allocate memory
      Allocate(cov(ncovar,ncovar))
      Allocate(corr_syst(ncovar,ncovar))
      Allocate(beta(ncovar,ncovar))
      Allocate(alpha(ncovar))

      print *,'Read covariance matrix ...'
      do i=1,NCovar
         read (51,*,err=99) ( Cov(i,j),j=1,NCovar)
      enddo
      do i=1,NCovar
         do j=1,NCovar
            corr_syst(i,j) = cov(i,j)
         enddo
      enddo

      close(51)
      print *,'Computing nuisance representation ...'

      Call GetNuisanceFromCovar(NCovar,NCovar,NCovar,Cov,Beta,
     $     Tolerance, NCorr, alpha, SeparateDiag)

      print *,' '
      print *,' RESULTS RESUTLS RESUTS'
      print *,' '
      print '('' N Nuisance parameters = '',i3)',NCorr
      write (*,'(A5,A12)',advance='no') ' Bin ', ' Diag  '
      print '(300(A8,i2,A2))',('param ',i,' ',i=1,NCorr)
      do i=1,NCovar
         print '(I4,300E12.4)',i,alpha(i),(beta(j,i),j=1,NCorr)
      enddo
      print *,' '
      print *,' END RESULTS '
      print *,' '

      
      dev = 0.
      if (Verbose) then
         print *,'  '
         print *,'Test Covariance matrix:'
         print '(2A4,3A12)','i','j',' nui ',' orig ','corr.diff%'
      endif
      do i=1,NCovar
         do j=1,NCovar
            test = 0.
            do k=1,NCorr
               test = test + beta(k,i)*beta(k,j)
            enddo
            if (i.eq.j) then
               test = test + alpha(i)*alpha(i)
            endif
            dev = (test-corr_syst(i,j))/sqrt(corr_syst(i,i)
     $           *corr_syst(j,j))*100.
            if (verbose) then
               print '(2i4,3F12.2)',i,j,test,corr_syst(i,j)
     $              ,dev
            endif
            devmax = max(abs(dev),devmax)
         enddo
      enddo
      print '(''Maximum deviation from the original correlation = ''
     $     ,F10.2,''%'')',devmax
      
      return
 91   call hf_errlog(1,'F:Can not open covar.in file')      
 92   call hf_errlog(2,'F:Can not read Covar namelist')
 93   call hf_errlog(3,'F:Can not find Covar namelist')
 99   call hf_errlog(3,'F:Error reading cov. matrix')
      end
      

C------------------------------------------------------------------------------
C
C> @brief Conversion from covariance matrix to nuisance parameter representation.
C
C> @param NDimCovar   -- Dimension of the covariance matrix
C> @param NDimSyst    -- Dimension of systematics matrix
C> @param NCovar -- Actual number of elements in the covariance matrix
C> @param Covar  -- Input covariance matrix. Output: nuisance parameters.
C> @param ANuisance -- Output nuisance parameter representation
C> @param Tolerance -- fractional sum of eigenvalues for the sourced treated as uncorrelated uncertainty. 0: NCorrelated = NCovar, 1: NCorrelated = 0.
C> @param Ncorrelated -- Output number of correlated nuisance parameters
C> @param Uncor      -- Output uncorrelated uncertainty 
C
C--------------------------------------------------------------------------------
      subroutine GetNuisanceFromCovar( NDimCovar, NDimSyst, NCovar,
     $     Covar, ANuisance, Tolerance, 
     $     Ncorrelated, Uncor, LSepDiag)
      implicit none
C--------------------------------------------------------------------------------
      integer NDimCovar, NDimSyst, NCovar
      double precision Covar(NDimCovar, NDimCovar)
      double precision ANuisance(NDimSyst, NDimCovar)
      double precision Tolerance
      integer Ncorrelated
      double precision Uncor(NDimCovar)
      logical LSepDiag
      
      double precision Eigenvalues(NDimCovar)
      integer ifail

      double precision factor, facMax, facMin

      double precision, allocatable :: testm(:,:),diag(:)
      double precision Sum,Run
      integer i,j,k

C--------------------------------------------------------------------------------
      
C Try to remove diagonal term first:

      if ( LSepDiag ) then

         allocate(testm(NCovar,NCovar))
         allocate(diag(NCovar))

         facMax = 1.0D0
         facMin = 0.0D0

         do while ((facMax-facMin.gt.0.01).or.(Eigenvalues(1).lt.0))
            factor = 0.5*(facMax + facMin)
            do i=1,NCovar
               do j=1,NCovar
                  testm(i,j) = Covar(i,j)
               enddo
            enddo
            do j=1,NCovar
               diag(j) = sqrt(factor*covar(j,j))
               testm(j,j) = Covar(j,j) - diag(j)*diag(j)
            enddo
            Call MyDSYEVD(NCovar,testm,NCovar, EigenValues,IFail)

c            print *,EigenValues(1)
            if (EigenValues(1).lt.0) then
               facMax = factor
            else
               facMin = factor
            endif
c            print *,'ha',factor,facMax,facMin
         enddo
         ! Ok, subtract diagonal:
         do j=1,NCovar
            Covar(j,j) = Covar(j,j) - diag(j)*diag(j)
         enddo
         DeAllocate(testm)

      endif

      
      Call MyDSYEVD(NCovar,Covar,NDimCovar, EigenValues,IFail)
      
      Sum = 0
      do i=1,NCovar
c         print *,i,EigenValues(i)
         Sum = Sum + EigenValues(i)
      enddo

      Ncorrelated = NCovar
      Run = 0.0D0
      do i=NCovar,1,-1
         Run = Run + EigenValues(i)/Sum
C         print *,i,EigenValues(i)/Sum,Run, 1.D0-Tolerance
         if (Run  .gt. 1.D0 - Tolerance) then
            Ncorrelated = NCovar - i + 1
            goto 17
         endif
      enddo
 17   continue
c      print *,NCorrelated

      do i=1,NCovar
         do j=1,NCovar
            Covar(j,i) = Covar(j,i)*sqrt(max(0.0D0,EigenValues(i)))
            ANuisance(i,j) = Covar(j,i)
C            print *,j,i, ANuisance(i,j),Covar(i,j)
         enddo
      enddo

      do j=1, NCovar
         do i=1,NCorrelated
            ANuisance(i,j) = Covar(j,NCovar-i+1)
         enddo
         do i=NCorrelated+1,NCovar
            Uncor(j) = Uncor(j) + Covar(j,NCovar-i+1)**2
         enddo
         Uncor(j) = sqrt(Uncor(j))
      enddo

      if (LSepDiag) then
       ! Add diag back to uncor:
         do j=1, NCovar
            Uncor(j) = sqrt( Uncor(j)**2 + diag(j)**2 )
         enddo
         DeAllocate(diag)
      endif

      end

C--------------------------------------------------------  	 
C> @Brief Interface to lapack, to dynamically allocate work arrays 
      subroutine MyDSYEVD(NCovar,Covar,NDimCovar, EigenValues,ifail)
      implicit none
      integer NCovar, NDimCovar
      double precision Covar(NDimCovar,NDimCovar), EigenValues(NCovar)
      integer IFail
      double precision Work
      integer IWork
      Character*80 msg
C---------------------------------------------------------------
C Determine optimal size of the work array:                                                                                                             
      Call DSYEVD('V','U',NCovar,Covar,NDimCovar, EigenValues, Work,
     $     -1, IWork, -1, IFail)


      write (msg,
     $ '(''I:MyDSYEVD: optimal dimensions for work arrays:'',2i6)')
     $     int(work)+1, iwork
      call HF_ERRLOG(14121701,msg)
      call MyDSYEVD2(NCovar,Covar,NDimCovar, EigenValues,
     $     int(work)+1,iwork,ifail)

      end

      subroutine MyDSYEVD2(NCovar,Covar,NDimCovar, EigenValues, nwork,
     $     nlwork,ifail)
      implicit none
      integer NCovar, NDimCovar
      double precision Covar(NDimCovar,NDimCovar), EigenValues(NCovar)
      integer nwork, nlwork
      double precision Work(nwork)  ! Dynamic array                                                                                                     
      integer IWork(nlwork)         ! Dynamic array                                                                                                     
      integer IFail
C---------------------------------------------------------------------                                                                                  
      Call DSYEVD('V','U',NCovar,Covar,NDimCovar, EigenValues, Work,
     $     nwork, IWork, nlwork, IFail)


      end
