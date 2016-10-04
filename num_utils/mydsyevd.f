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
