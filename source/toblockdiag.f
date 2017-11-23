      Subroutine ToBlockDiag(diag,last,corr,box)
C---------------------------------------------------------------
C Transform problem to block-diagonal matrix. Matrices "box" and "last"
C are modified
C---------------------------------------------------------------
      implicit none
      include 'common.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 Box(NSYSTMAX,NSYSTMAX) ! Original syst --> covariance matrix
C Local:
      integer if2,isys,j,iexp
      real*8  SUM

      real*8 work(NSYSTMAX)
      integer ifail

      integer NDiag

      double precision time1,time2
C---------------------------------------------------------------------


C Prepare vaiables:
      NDiag = NMeas 

      call cpu_time(time1)

C Invert box matrix corresponding in order to get systematic uncertainties:
      if (NSYSTOT.ne.0) then
        if (iItr.ne.NIteration) then
C Solve without inversion (faster):
          Call DEQN(nsystot,box,nsystmax,work,ifail,1,last(NDiag+1))
        else
C Calc As'^-1 (box) and (last) Bave = As'^-1 (Cs − Asm^T Am^-1 Cm)
          Call DEQINV(nsystot,box,nsystmax,work,ifail,1,last(NDiag+1))
        endif
      endif

      call cpu_time(time2)
      print *,'Inversion = ',time1,time2,time2-time1,NIteration

      if (IFail.ne.0) then
         call hf_errlog(1,'F:Failed to invert syst. matrix !!!') 
      endif

      
      do j=1,NSYSTOT
         SYSSHItr(j,iItr) = last(NDiag+j)
         SYSSH(j) = last(NDiag+j)
      enddo


C Calculate averaged value
      do if2=1,NDiag
         Sum = Last(if2)
         do isys=1,NSYSTOT
            Sum = Sum - SYSSH(isys)*Corr(isys,if2)
         enddo
         F2VAVE(if2) = Sum/diag(if2)
      enddo

      call SaveChi2()
      end 



C---------------------------------------------------

      subroutine LastIteration(diag,last,corr,boxinv,box)
      implicit none
      include 'common.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 Box(NSYSTMAX,NSYSTMAX) ! Covariance matrix As'
      real*8 Boxinv(NSYSTMAX,NSYSTMAX) ! Inverted Covariance matrix As'^-1

      real*8 Corr2(NSYSTMAX,NF2MAX)  !  Rotated systematics:
      real*8 chi2
      integer i,j,if2, ifail, imeas, ipc, isys1, isys2, isys
      real*8 work(NSYSTMAX)
      real*8 WWW(NSYSTMAX)  ! Eigenvalues
      real*8 sum
      integer ndf
      real sterr,uncerr,toterr

      double precision time1,time2

      call cpu_time(time1)

C Calc error reduction
      do j=1,nsystot
         errsyst(j) = sqrt(boxinv(j,j))
      enddo

C Print the central values of the sys uncert and the errors of the systematics
      call PrintSystSummary()

      call cpu_time(time2)
      print *,'LI1 = ',time1,time2,time2-time1
C Get NEW SYST Errors
      do if2=1,NMeas
         SUM = 0
         do isys1=1,NSYSTOT
            do isys2=1,NSYSTOT
               SUM = SUM + Corr(isys1,if2)*
     $              Corr(isys2,if2)
     $              * boxinv(isys1,isys2)
            enddo
         enddo
         F2EsyAve(if2) = sqrt(Sum)/diag(if2)
         if (IDEBUG.gt.3) then
            print *,if2,sqrt(Sum)/diag(if2)
         endif
      enddo

      call cpu_time(time1)
      print *,'LI2 = ',time1,time2,time1-time2

C Get stat and total uncertainty
      do if2=1,NMeas
         STERR = 0.0d0
         UNCERR = 0.0d0
C Stat+Uncorr uncertainties:
C Delta_ave (F2EstAve) = sqrt(Am^-1)
         F2EstAve(if2) = sqrt(1.0d0/DIAG(if2))
         
C True stat. errors:
C (F2EstAveTrue) = Delta_ave^4 * Delta_stat^2/Delta^4
         do imeas = 1, NMeasF2(if2)

            STERR = STERR  +
     $           F2ETAB_STA(if2,imeas)**2/F2ETAB(if2,imeas)**4
            UNCERR = UNCERR +
     $           (F2ETAB(if2,imeas)**2-F2ETAB_STA(if2,imeas)**2)
     $           /F2ETAB(if2,imeas)**4
            TOTERR = TOTERR +
     $           (F2ETAB_TOT(if2,imeas)**2)
     $           /F2ETAB(if2,imeas)**4
         enddo
         STERR = sqrt(STERR/(DIAG(if2)**2))
         UNCERR = sqrt(UNCERR/(DIAG(if2)**2))
         TOTERR = SQRT(TOTERR/(DIAG(if2)**2))

         F2EstAveTrue(if2) = STERR 
         F2EAVE(if2)      = sqrt(F2EsyAve(if2)**2
     $     +F2EstAve(if2)**2)

      enddo

C Calculate CHI2 according to the formula for chi2:
      call CalcChi2(chi2, ndf)

C Write Averaged values and chi2
      call PrintAveSummary(chi2, ndf)


C Zero close to zero correlations to avoid strange results
C for a single dataset
      do isys1=1,nsystot
         do isys2=1,nsystot
            if (abs(box(isys1,isys2)).lt.1.D-10) then
               box(isys1,isys2) = 0
            endif
         enddo
         if (abs(box(isys1,isys1)-1.0).lt.1.D-10) then
            box(isys1,isys1) = 1.0D0
         endif
      enddo

      call cpu_time(time2)
      print *,'LI5 = ',time1,time2,time2-time1

C     Print variance matrix
C      Call PrintMatrix(box)

C Get eigenvectors and eigenvalues:
C WWW - eigenvalues D^2 in ascending order of As'
C (box) = eigenvectors U of the matrix As' (box)
      Call MyDSYEVD(NSysTot, Box, NSystMax,WWW,ifail)

      call cpu_time(time1)
      print *,'Decomposition = ',time1,time2,time1-time2

C Scale Box to take into account error reduction
C (box) = UD
      do i=1,NSysTot
         do j=1,NSysTot
            Box(j,i) = Box(j,i)/dsqrt(WWW(i))
         enddo
      enddo

      call cpu_time(time2)
      print *,'LI7 = ',time1,time2,time2-time1

C Post-process Box3, make it triangular and positive along the diagonal
      if (PostRotateSyst) then
         call HF_ERRLOG(15030701,
     $ 'I:Rotate systematic error sources along original directions')
         Call PostRotate(nsystot,nsystmax,box)
      endif

C     Print the eigenvalues and eigenvectors
      call PrintEigInfo(WWW,box)

      call cpu_time(time1)
      print *,'Rotation = ',time1,time2,time1-time2

C Get rotated systematic matrix:
      if(UseBlas)then
        call dgemm('N','N',nsystot,NMeas, nsystot, 1.D0,
     $           Box(1:nsystot,1:nsystot), nsystot, 
     $   Corr(1:nsystot,1:NMeas), nsystot, 0.D0, 
     $   Corr2(1:nsystot,1:NMeas), nsystot)
      endif

      call cpu_time(time2)
      print *,'LI8 = ',time1,time2,time2-time1

      do if2=1,NMeas
         do isys1=1,NSYSTOT
             if(.not. UseBlas) then
              do isys=1,NSYSTOT
                Corr2(isys1,if2) = Corr2(isys1,if2) +
     $              Corr(isys,if2)*Box(isys,isys1)
              enddo
             endif

C (SystDiag) = Asm Am^-1 D^{-1} U^T
            SystDiag(isys1,if2) = 
     $           Corr2(isys1,if2)/diag(if2)

            SystDiagPercent(isys1,if2) = 
     $           -100.D0*SystDiag(isys1,if2)
     $           /f2vave(if2)

C (SystOrig) = Asm Am^-1 sqrt(Diag(As'^-1))
            SystOrig(isys1,if2) = CORR(isys1,if2)/DIAG(if2)
     $           *ErrSyst(isys1)

            SystOrigPercent(isys1,if2) = -100.*SystOrig(isys1,if2)
     $           /f2vave(if2)  
     

         enddo
      enddo

      call cpu_time(time1)
      print *,'LI9 = ',time1,time2,time1-time2


C---------------------------------------------------------------------
      end



      Subroutine PostRotate( Nsystot, nsystmax,box3)
      implicit none
      integer NSysTot, NSystMax
      double precision Box3(NSystMax,NSystMax)
C Dynamic:
      double precision RR(NSysTot,Nsystot)
      double precision AA2(NSysTot,Nsystot)
      double precision C(NSysTot)

      integer k,i,j,l
      integer ifail
C---------------------------------------------------------------

C Starting with last source, up to the first:
      do k=Nsystot,1,-1

         do i=1,Nsystot
            do j=1,Nsystot
               RR(i,j) = 0.
               if (i.eq.j) then
                  RR(i,j) = 1.
               endif

               RR(i,j) = RR(i,j) + box3(k,i)*box3(k,j)

            enddo
         enddo

        if (NSYSTOT.ne.0) then
          Call MyDSYEVD(k,RR,NSysTot,C,ifail)
        endif
C rotate rotation matrix:

         do i=1,k
            do j=1,k
               AA2(i,j) = 0.
               do l=1,k
                  AA2(i,j) = AA2(i,j) + box3(i,l)*RR(l,j)
               enddo
            enddo
         enddo


         do i=1,k
            do j=1,k
               box3(i,j) = AA2(i,j)
            enddo
         enddo
      enddo                     ! loop over k.  

C   Last loop to keep the direction of the original vectors                                                                                                                                                         
      do i=1,NSysTot
         if (box3(i,i).lt.0) then
            do j=1,NSysTot
               box3(j,i) = -box3(j,i)
            enddo
         endif
      enddo


C---------------------------------------------------------------
      end

C Print symmary for systematic uncertainty
      Subroutine PrintSystSummary()
      implicit none
      include 'common.inc'
      integer j

      print *,' ' 
      print *,' ' 
      print *,' ===== RESULTS RESULTS RESULTS RESULTS RESULTS ======'
      print *,' ' 
      print *,' '
      print *,' ============== Fitted systematics =================='
      print *,' '
      print '((A3),(A32),''  '',(A8),(A8))',
     $     'Ind',' Name','Shift','Error'
      print *,' ----------------------------------------------------'
      do j=1,nsystot
         print '(i3,(A32),''  '',1F8.2,2F8.2)',
     $        j,TRIM(ADJUSTL(SystematicName(j))),
     $        SYSSH(j),errsyst(j)
      enddo
      print *,' ----------------------------------------------------'

      end

C Print symmary for averaged values 
      Subroutine PrintAveSummary(chi2, ndf)
      implicit none
      include 'common.inc'
      integer i,if2, ipc, isys
      real*8 chi2
      integer ndf

      print *,' ' 
      print *,' ' 
      print *,' ' 
      write(*,*) ' =============== Averaged Data '//
     $     ' ====================== '
      print *,' ' 
      print *,' ' 
      print *,' ' 
         
       !> Loop over process classes, to make neat printout
      do iPC= 1,NProcClass
!> Do little header:
         print *,'======================================='
         print 
     $        '(''Reaction IDX='',i3,'' Process name='',A)',
     $        idxReactionMeas(iPC),gridreaction(idxReactionMeas(iPC))
         
         print *,'======================================='
         write (*,'(10A12)',advance='no')
     $        (GridBinNames(i,idxReactionMeas(iPC))
     $        ,i=1,NDimensionGrid(idxReactionMeas(iPC)))
         write (*,'(4A12)') 'Average ', 'StatErr ', 
     $        'SystErr', 'TotalErr' 
         
         do if2=1,NMeas
            if (idxProcessClass(if2).eq.iPC) then
               
               write (*,'(10ES12.4)',advance='no') 
     $         ( GridPoints(idxGridMeas(if2),i,idxReactionMeas(iPC))
     $              ,i=1,NDimensionGrid(idxReactionMeas(iPC)))
               
               write(*,'(4ES12.4)') 
     $              F2VAVE(if2),
     $              F2EstAve(if2), F2EsyAve(if2), F2EAVE(if2)
            endif
         enddo
      enddo

      print *,' '
      print '(''TOTAL Chi2/ndf='',F10.4,''/'',i4)',chi2,ndf
      print *,' '

      end

C Calculate CHI2 according to the formula for chi2:
      Subroutine CalcChi2(chi2, ndf)
      implicit none
      include 'common.inc'
      real*8 chi2,chi2loc, sum
      integer ndf, i, j, isys

      chi2 = 0.0
      ndf  =  0
      do i = 1, NMeas
         do j = 1, NMeasF2(i)
            ndf = ndf + 1
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
      ndf = ndf-NMeas
      end
