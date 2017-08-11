      subroutine StatRecalc
C--------------------------------------------
C Rescale stat. error of the measurement assuming that relative 
C errors do not change.
C--------------------------------------------
      implicit none
      include 'common.inc'
      integer if2,iexp,isys, NDiag
      real f2orig, f2new, erorig, scal
      real uncor, stat
C 9 jan 09: add also expected value - systematic shift:
      real*8 f2newSys
C--------------------------------------------
      integer NMatr
C--------------------------------------------
C
C Prepare vaiables:
C
      NDiag = NMeas ! min(NMeas,NPointGrid(1))
      NMatr = NDiag + NSysTot  ! Range of the matrix

      do if2 = 1,NDiag
         do iexp = 1,NMeasF2(if2)

            f2orig = F2TAB(if2,iexp)
            erorig = F2ETABOrig(if2,iexp)

            if (erorig.eq.0) goto 1717
               f2new = F2VAVE(if2)
               
               f2newSys = f2new

               do isys=1,NSYSTOT
                  f2newSys = f2newSys -
     $                 SYSTAB(isys,if2,iexp)*SYSSH(isys)
               enddo
            if (f2orig.eq.0) then
               scal = 0.
            else
               if (CorrectStatBias) then
                  scal = f2newSys/f2orig
               else
                  scal = f2new/f2orig
               endif
            endif

            if (scal.lt.0) then
               scal = 1.0
            endif

            print *,f2new,f2orig,f2newSys,scal

C New error:
            if (RescaleStatSep.or.FixStat) then

C Sqrt scaling for stat:
               if (RescaleStatSep) then !LSqrtStat
                  stat  = F2ETAB_STAORIG(if2,iexp) * sqrt(scal)
                  if (IDebug.gt.2) then
                     write(*,*) 'scal = ', if2, f2orig, sqrt(scal)
                  endif
               else
                  stat =  F2ETAB_STAORIG(if2,iexp)
               endif
               F2ETAB_STA(if2,iexp) = stat

C Rescale uncor syst. part:
               uncor = F2ETAB_UNC(if2,iexp)
               uncor = uncor * scal

C Total uncor error:
               F2ETAB(if2,iexp) = sqrt(stat**2+uncor**2)
            else
               F2ETAB(if2,iexp) = erorig * scal
            endif

C New systematic sensitivity:
C Next iterator for asymmetric systematic
            do isys=1,NSYSTOT
C Multiplicative uncertainty
                if (SystematicForm(isys).eq.'M') then
C Asymmetric
                  if(ShiftType(isys,if2,iexp).eq.2) then
                   SYSTAB(isys,if2,iexp) =
     &              ((SYSTABOrig(isys,if2,iexp,2) -
     &              SYSTABOrig(isys,if2,iexp,3))/2) +
     &              ((SYSTABOrig(isys,if2,iexp,2) +
     &              SYSTABOrig(isys,if2,iexp,3))/2) * SYSSH(isys)
                   SYSTAB(isys,if2,iexp) =
     &              SYSTAB(isys,if2,iexp) * scal
                  endif
C Symmetric
                  if(ShiftType(isys,if2,iexp).eq.0) then
                    SYSTAB(isys,if2,iexp) =
     $              SYSTABOrig(isys,if2,iexp,1)
     $               * scal
                  endif
C Additive uncertainty
                elseif (SystematicForm(isys).eq.'A') then
C Asymmetric
                  if(ShiftType(isys,if2,iexp).eq.2) then
                   SYSTAB(isys,if2,iexp) =
     &              ((SYSTABOrig(isys,if2,iexp,2) -
     &              SYSTABOrig(isys,if2,iexp,3))/2) +
     &              ((SYSTABOrig(isys,if2,iexp,2) +
     &              SYSTABOrig(isys,if2,iexp,3))/2) * SYSSH(isys)
                  endif
                else
C Other type of asymmetric   uncertainty
                  if(ShiftType(isys,if2,iexp).eq.2) then
                   SYSTAB(isys,if2,iexp) =
     &              ((SYSTABOrig(isys,if2,iexp,2) -
     &              SYSTABOrig(isys,if2,iexp,3))/2) +
     &              ((SYSTABOrig(isys,if2,iexp,2) +
     &              SYSTABOrig(isys,if2,iexp,3))/2) * SYSSH(isys)

                  endif
                  write(*,*) 'StatRecalc: WARNING: Treated as additive'
                endif
            enddo
 1717    enddo
      enddo

      end
