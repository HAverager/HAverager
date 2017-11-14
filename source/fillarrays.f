      Subroutine FillArrays(diag,last,corr,box,onlyLast)
C---------------------------------------------------------------
C Fill axillary arrays for the matrix inversion
C---------------------------------------------------------------
      implicit none
      include 'common.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 Box(NSystmax,NSystmax)
      logical onlyLast

      real*8 oneOverDiag(NF2Max)
      
C Local:
      integer isys,if2,iexp,i,j,k,isys1,isys2
      integer noff
      integer NDiag, NMatr
      real*8 coef, erro

      real*8 work(nf2max+nsystmax)
      integer ifail
      integer ip(nf2max+nsystmax)

      real*8, allocatable :: err2tab(:,:)

      double precision time1,time2
      logical useBlas

      real*8, allocatable :: sys1(:,:),sys2(:,:)
      real*8, allocatable :: sys1b(:,:),sys2b(:,:)
      logical lfirst
      data lfirst/.true./
      integer NdataT,idataT

      double precision box2(nsystmax,nsystmax)
      
C---------------------------------------------------------------------



      if (lfirst)  then
         NdataT = 0
         do i=1,NMeas
            do j=1,NMeasF2(i)
               NdataT = NdataT + 1
            enddo
         enddo
         lfirst = .false.
      endif
      
C Count number of SF / X-section points to fit 

      NDiag = NMeas !min(NMeas,NPointGrid(1)) !idxReaction <= IMPROVE      !NF2TOT
      NMatr = NDiag+NSysTot

      NF2TOT = NDiag
      NSFSEC = NF2TOT

      !NSysTot = NSysTot - 1  ! because the stat uncertainty is not included


      Allocate (err2tab(NDiag,NMEASMAX)) ! !

      do i=1,NDiag
         do j=1, NMeasF2(i)
            err2tab(i,j) = 1./F2ETAB(i,j)**2
         enddo
      enddo


C
C Construct the matrix:
C

      if(onlyLast)then

C Set values to zero:
         do i = 1,NDiag
            diag(i)   = 0.0         
            do j = 1,NSysTot
               corr(j,i) = 0.0
            enddo
         enddo


         do i=1,NSysTot
            do j=1,NSysTot
               box(i,j) = 0.0
               box2(i,j) = 0.0
            enddo
         enddo

C Fill in the diagonal elements (matrix Am)
         do i = 1,NDiag
            do j = 1,NMeasF2(i)
               diag(i)   = diag(i) + err2tab(i,j)
            enddo

            oneOverDiag(i) = 1/diag(i)
         enddo



C Fill in the rectangle correlation matrix (Asm):
         do k = 1,NSysTot
            do i = 1,NDiag
               do j = 1,NMEASF2(i)
                  CORR(k,i) = CORR(k,i) 
     $                 -SYSTAB(k,i,j)
     $                 *err2tab(i,j)
               enddo
            enddo
         enddo


C Fill in the box matrix (As):

         call cpu_time(time1)
         print *,'here',time1

         useBlas = .true.
c       useBlas = .false.

         if (useBlas) then
         
! perpare arrays:
            allocate(sys1(nsystot,NDataT))
            allocate(sys2(NDataT,nsystot))
            do isys=1,NSysTot
               iDataT = 0
               do i=1,NDiag
                  do j=1,NMeasF2(i)
                     iDataT = iDataT + 1
                     sys1(isys,iDataT) = SYSTAB(isys,i,j)
                     sys2(iDataT,isys) = SYSTAB(isys,i,j)*err2tab(i,j)
                  enddo
               enddo
            enddo
         
         ! BLAS:
            call dgemm('N','N',nsystot,nsystot, nDataT, 1.0D0,
c         call cublas_dgemm('N','N',nsystot,nsystot, nDataT, 1.0D0,
     $           sys1, nsystot, sys2, NDataT, 0.D0, box, nsystmax)

            deallocate(sys1)
            deallocate(sys2)        
         else
            do i     = 1,NDiag
               do j     = 1, NMeasF2(i)
                  do isys1 = 1,NSysTot
                     do isys2 = isys1,NSysTot
                        box(isys1,isys2) = box(isys1,isys2)
     $                       + SYSTAB(isys1,i,j)
     $                       * SYSTAB(isys2,i,j) 
     $                       * err2tab(i,j)
                     enddo
                  enddo
               enddo
            enddo
            do isys1 = 1,NSysTot
               do isys2 = isys1+1,NSysTot
                  box(isys2,isys1) = box(isys1,isys2)
               enddo
            enddo
         endif
         call cpu_time(time2)
         print *,'here2',time1,time2,time2-time1

      

C     Add 1 to box diagonal:
         do j=1,NSYSTOT
            BOX(j,j) = BOX(j,j) + 1.0
         enddo
      else
         print *,'NOT ONLY LAST'
      endif  ! box part only

C Fill in the right side vector:
      do i = 1,NMatr
         last(i) = 0.0
      enddo

       do i = 1,NDiag
         do j = 1,NMeasF2(i)
            last(i) = last(i) + 
     $           F2TAB(i,j)*err2tab(i,j)
         enddo
       enddo

C Last column correlation:
       do k = 1,NSysTot
         do i = 1,NDiag
            do j = 1,NMeasF2(i)
               last(NDiag+k) = last(NDiag+k) 
     $              - SYSTAB(k,i,j)
     $              * F2TAB(i,j)*err2tab(i,j)
            enddo
         enddo
       enddo


C    Fill matrix A' = As - Asm^T Am^-1 Asm
      call cpu_time(time1)
      print *,'here3',time2,time1,time1-time2

      do if2=1,NDiag
         do isys=1,NSYSTOT
            Coef = - Corr(isys,if2) * oneOverDiag(if2)
            last(NDiag+isys) = last(NDiag+isys) + last(if2)*Coef
         enddo
      enddo

      if (onlylast) then


c         if (.false.) then      ! needs debuging.
         if (useBlas) then
            allocate(sys1b(nsystot,NDiag))
            allocate(sys2b(NDiag,nsystot))
            do isys=1,nsystot
               do if2=1,ndiag
                  sys1b(if2,isys) =  -corr(isys,if2)
                  sys2b(isys,if2) =  corr(isys,if2)*oneOverDiag(if2)
               enddo
            enddo

         ! BLAS:
            call dgemm('N','N',nsystot,nsystot, ndiag, 1.D0,
c            call cublas_dgemm('N','N',nsystot,nsystot, nDataT, 1.0D0,
     $           sys1b, nsystot, sys2b, Ndiag, 1.D0, box, nsystmax)

c            do if2=1,NDiag
c               do isys=1,NSYSTOT
c                  do j=1,NSYSTOT
c                     box2(isys,j) = box2(isys,j)
c     $                    + sys1b(isys,if2)*sys2b(if2,j)
c                  enddo
c               enddo
c            enddo
                        
            deallocate(sys1b)
            deallocate(sys2b)

         else
            print *,box(1,2),box(2,1)
            do if2=1,NDiag
               do isys=1,NSYSTOT
                  do j=isys,NSYSTOT
                     box(isys,j) = box(isys,j)
     $                    - corr(j,if2)  
     $                    * corr(isys,if2) * oneOverDiag(if2)
                  enddo
               enddo
            enddo
            
            do isys1 = 1,NSysTot
               do isys2 = isys1+1,NSysTot
                  box(isys2,isys1) = box(isys1,isys2)
               enddo
            enddo         
         endif
      endif

      print *,box(1,2),box(2,1)
      
      call cpu_time(time2)
      print *,'here4',time1,time2,time2-time1
      
      if (IDebug. ge. 4) then
         write(*,*) 'debug:'
         print *,'cor'
         do if2=1,NDiag
            print '(20F10.4)',(corr(j,if2),j=1,NSYSTOT)
         enddo
         print *,'box'

         do j=1,NSYSTOT
            print '(20F10.4)',(box(i,j),i=1,NSYSTOT)
         enddo
         print *,'last'
         print '(6F12.4)',(last(j),j=1,NDiag+NSYSTOT)
      endif

      call cpu_time(time1)
      print *,'here5',time2,time1,time1-time2

      DeAllocate( err2tab)
C-------------------------------------
      end

