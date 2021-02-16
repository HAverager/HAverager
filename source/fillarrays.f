      Subroutine FillArrays(diag,last,corr,box,fillSyst)
C---------------------------------------------------------------
C Fill axillary arrays for the matrix inversion
C---------------------------------------------------------------
      implicit none
      include 'common.inc'
      real*8 diag(NMeas),Last(NMeas+NSystot),Corr(NSystot,NMeas)
      real*8 Box(NSystot,NSystot)
      logical fillSyst

C     Local:
      integer isys,if2,iexp,i,j,k,isys1,isys2
      integer noff
      integer NDiag, NMatr
      real*8 coef, erro

      real*8, allocatable :: err2tab(:,:)
      real*8, allocatable :: oneOverDiag(:)

      double precision time1,time2

      real*8, allocatable :: sys1(:,:),sys2(:,:)
      real*8, allocatable :: sys1b(:,:),sys2b(:,:)
      integer NdataT,idataT
   
C---------------------------------------------------------------------
      
C     Count number of SF / X-section points to fit
      NDiag = NMeas !min(NMeas,NPointGrid(1)) !idxReaction <= IMPROVE
      NMatr = NDiag+NSysTot

      Allocate (err2tab(NMeas,NMEASMAX))
      Allocate (oneOverDiag(NMeas))

      do i=1,NDiag
         do j=1, NMeasF2(i)
            err2tab(i,j) = 1./F2ETAB(i,j)**2
         enddo
      enddo


C Construct the matrix with uncertainties:
      if(fillSyst)then

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
         if (useBlas) then
           NdataT = 0
           do i=1,NMeas
             NdataT = NdataT + NMeasF2(i)
           enddo
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
     $           sys1, nsystot, sys2, NDataT, 0.D0, box, nsystot)
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
         if (IDEBUG.gt.0) then
           print '(" Time Fill As" 3(e9.2))',time1,time2,time2-time1
         endif


C     Add 1 to box diagonal:
         do j=1,NSYSTOT
            BOX(j,j) = BOX(j,j) + 1.0
         enddo
      else
         print *,'NOT ONLY LAST'
      endif

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
      do if2=1,NDiag
         do isys=1,NSYSTOT
            Coef = - Corr(isys,if2) * oneOverDiag(if2)
            last(NDiag+isys) = last(NDiag+isys) + last(if2)*Coef
         enddo
      enddo

      if (fillSyst) then

         if (useBlas) then
            allocate(sys1b(NDiag,nsystot))
            allocate(sys2b(nsystot,NDiag))
            do isys=1,nsystot
               do if2=1,ndiag
                  sys1b(if2,isys) =  -corr(isys,if2)
                  sys2b(isys,if2) =  corr(isys,if2)*oneOverDiag(if2)
               enddo
            enddo

         ! BLAS:
            call dgemm('N','N',nsystot,nsystot, ndiag, 1.D0,
     $           sys2b, nsystot, sys1b, Ndiag, 1.D0, box, nsystot)

            deallocate(sys1b)
            deallocate(sys2b)

         else
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
      
      call cpu_time(time2)
      if (IDEBUG.gt.0) then
          print '(" Time Fill A^" 3(e9.2))',time1,time2,time2-time1
      endif

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

      DeAllocate(err2tab)
      deAllocate(oneOverDiag)
C-------------------------------------
      end

