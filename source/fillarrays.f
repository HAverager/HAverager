      Subroutine FillArrays(diag,last,corr,box,onlyLast)
C---------------------------------------------------------------
C Fill axillary arrays for the matrix inversion
C---------------------------------------------------------------
      implicit none
      include 'common.inc'
      real*8 diag(NF2MAX),Last(NF2MAX+NSYSTMAX),Corr(NSYSTMAX,NF2MAX)
      real*8 Box(NSystmax,NSystmax)
      logical onlyLast

C Local:
      integer isys,if2,iexp,i,j,k,isys1,isys2
      integer noff
      integer NDiag, NMatr
      real*8 coef, erro

      real*8 work(nf2max+nsystmax)
      integer ifail
      integer ip(nf2max+nsystmax)

      real, allocatable :: err2tab(:,:)

C---------------------------------------------------------------------




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
         enddo
      enddo

C Fill in the diagonal elements (matrix Am)
      do i = 1,NDiag
         do j = 1,NMeasF2(i)
            diag(i)   = diag(i) + err2tab(i,j)
         enddo
      enddo



C Fill in the rectangle correlation matrix (Asm):
      do k = 1,NSysTot
         do i = 1,NDiag
            do j = 1,NMEASF2(i)
               CORR(k,i) = CORR(k,i) 
     $              -SYSTAB(k,i,j)
     $              *err2tab(i,j)
            enddo
         enddo
      enddo


C Fill in the box matrix (As):
      do isys1 = 1,NSysTot
         do isys2 = isys1,NSysTot
            do i     = 1,NDiag
               do j     = 1, NMeasF2(i)
                  box(isys1,isys2) = box(isys1,isys2)
     $                 + SYSTAB(isys1,i,j)
     $                 * SYSTAB(isys2,i,j) 
     $                 * err2tab(i,j)
               enddo
            enddo
         enddo
      enddo

      do isys1 = 1,NSysTot
         do isys2 = isys1+1,NSysTot
            box(isys2,isys1) = box(isys1,isys2)
         enddo
      enddo

C Add 1 to box diagonal:
      do j=1,NSYSTOT
         BOX(j,j) = BOX(j,j) + 1.0
      enddo

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

C Fill matrix A' = As - Asm^T Am^-1 Asm
      do if2=1,NDiag
         do isys=1,NSYSTOT
            Coef = - Corr(isys,if2)/diag(if2)
            if (Coef.ne.0) then
               if (onlyLast) then
               do j=1,NSYSTOT
                  box(j,isys) = box(j,isys) + corr(j,if2)*Coef
               enddo
               endif
               last(NDiag+isys) = last(NDiag+isys) + last(if2)*Coef
            endif
         enddo
      enddo

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


      DeAllocate( err2tab)
C-------------------------------------
      end

