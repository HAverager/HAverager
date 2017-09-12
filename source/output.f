      Subroutine Output
      implicit none
      include 'common.inc'



C     Write default output in user-frendly format
      call writeSimpleOut()

C     Write default output in simple format (easy to parce)
      call outtable()

C     Write output in X-Fitter format
      call writeXFitterOut()

C     Write summary for off-set systematics
      call outOffsetSyst()

C     Write summary for impact of systematics
      call outSystImpact()

C     Write summary for stat. ToyMC
      call outToyMCStat

C     Write systematics shifts and pulls
      call WriteSystShifts()

C     Write pulls of combined values
      call WritePulls()

C     Write information about iterations
      call WriteItrInfo()

      end


      subroutine WriteItrInfo()
        implicit none
        include 'common.inc'
        integer isys,itr

        print *,'Print: ',98,Chi2Itr(98)

        open (55,file=trim(OutputFolder)//'/ItrInfo.dat',
     $     status='unknown')

        write (55, '(200F10.2)'),
     $       (Chi2Itr(itr),itr=1,NIteration)

C       Loop over systematics
        do isys=1,nsystot
          write (55, '(200F10.4)'),
     $       (SYSSHItr(isys,itr),itr=1,NIteration)
        enddo

      end

      subroutine writeXFitterOut()
      implicit none
      include 'common.inc'

      character*80 ctmp
      integer iPC, iF2,i,i2, NPoint, NCol, NBin, isys, isysOff
      double precision uncor_syst
      character*8 date
C---------------------------------------------------------
      Call Date_and_time(DATE=date)

            do iPC=1, NProcClass
         ! Write out the header:

         NBin = NDimensionGrid(idxReactionMeas(iPC))

        ! grid  ! sigma, stat, uncor     All syst.
         NCol = NBin         + 3        + NSysTot + NSysOTot

         write (ctmp,'(A,''/'',A,i3,''.dat'')')
     $        trim(OutputFolder),trim(OutputPrefix), 100+iPC
         open (53,file=ctmp,status='unknown')

C Prepare header:
         write (53,*)
         write (53,'('' !*  File produced by averager'')')
         write (53,'('' !*  Created on '',A)') date
         write (53,*)
         write (53,'(''&Data '')')
         write (53,'(''  Name = "Average '',a,'' '',i3,''" '')')
     $        Trim(OutputPrefix), 100+iPC
         write (53,'(''  Reaction = "'',A,''"'')')
     $        Trim(gridreaction(idxReactionMeas(iPC)))
         write (53,'(''  NData = '',i5)') NPoint
         write (53,'(''  NColumn = '',i5)') NCol
         write (53,
     $'(''  ColumnType = '',i1,''*"Bin","Sigma", '',i3,''*"Error"'')')
     $        NBin,2+NSysTot+NSysOTot
         write (53,'(''  ColumnName = '',10(''"'',A,''",''))'
     $        , advance='no')
     $        (Trim(GridBinNames(i,idxReactionMeas(iPC)))
     $        ,i=1,NDimensionGrid(idxReactionMeas(iPC)))


         write (53,'(''Sigma", "stat", "uncor"'')',advance='no')
         write (53,'('','',500(''"'',A,i4''",''))', advance='no')
     $  ( 'sys'//Trim(OutputPrefix),i,i=1000+1,1000+NSysTot+NSysOTot-1)
         write (53,'(A,i4,A,''"'')') 'sys'//Trim(OutputPrefix)
     $        ,1000+NSysTot+NSysOTot,':O'

         write (53,'(''  IndexDataset = '', i5)') 100+iPC

         write (53,'(''  Percent = '',500(A,'',''))')
     $        ( 'false',i=1,2+NSysTot)

         write (53,'(''&End '')')

C        Header is ready, time to write the cross sections:

C ... but a comment line first ...
         write (53,'(''*'' A10,9A12)',advance='no')
     $        (trim(GridBinNames(i,idxReactionMeas(iPC)))
     $        ,i=1,NDimensionGrid(idxReactionMeas(iPC)))

         write (53,'(5A12,100A32)') 'Average', 'Stat','Uncor  ',
     $        (SystematicName(isys),isys=1,NSYSTOT),
     $        (SystematicOName(isysOff),isysOff=1,NSYSOTOT)


         do if2=1,NMeas
            if (idxProcessClass(if2).eq.iPC) then

               write (53,'(10ES12.4)',advance='no')
     $            ( GridPoints(idxGridMeas(if2),i,idxReactionMeas(iPC))
     $              ,i=1,NDimensionGrid(idxReactionMeas(iPC)))


               uncor_syst = sqrt(abs(f2estave(if2)**2
     $              -f2estaveTrue(if2)**2))

               if (iOutput.eq.1) then
                  write(53,'(3ES12.4,700F8.3)')
     $                 F2VAVE(if2),
     $                 F2EstAveTrue(if2)
     $                 ,uncor_syst
     $                 ,(SystDiag(isys,if2),isys=1,NSYSTOT)
     $                 ,(OSyst(if2,isysOff),isysOff=1,NSYSOTOT)
               else
C Write original sources, neglecting correlations (!)
                  write(53,'(3ES12.4,700F8.3)')
     $                 F2VAVE(if2),
     $                 F2EstAveTrue(if2)
     $                 ,uncor_syst
     $                 ,(SystOrig(isys,if2),isys=1,NSYSTOT)
     $                 ,(OSyst(if2,isysOff),isysOff=1,NSYSOTOT)
               endif
            endif
         enddo
         close (53)
      enddo

      end

      subroutine writeSimpleOut()
      implicit none
      include 'common.inc'

      character*80 ctmp
      integer iPC, iF2,i,i2, NPoint, isys, isysOff
      double precision sum
      logical found


      ctmp = trim(OutputFolder)//'/'//trim(OutputPrefix)//'.dat'
      open (52, file=ctmp, status='unknown')

      do iPC=1, NProcClass
         ! Write out the header:

         write (52,
     $   '(''Reaction IDX='',i3,'' Process name='',A)')
     $    idxReactionMeas(iPC),gridreaction(idxReactionMeas(iPC))

         write (52,'(10A12)',advance='no')
     $        (GridBinNames(i,idxReactionMeas(iPC))
     $        ,i=1,NDimensionGrid(idxReactionMeas(iPC)))
         write (52,'(4A12)', advance='no') 'Average ', 'StatErr ',
     $        'SystErr', 'TotalErr'

         if ( WriteOriginal ) then
            do i2=1,NInputFiles
               write (52,'(A10,I2,3A12)',advance='no')
     $              'Data',i2,'TotalErr','DataShift'
     $           ,'UncorErr'
            enddo
         endif

         write (52,*)


         NPoint = 0

         do if2=1,NMeas
            if (idxProcessClass(if2).eq.iPC) then

               NPoint = NPoint + 1

               write (52,'(10ES12.4)',advance='no')
     $            ( GridPoints(idxGridMeas(if2),i,idxReactionMeas(iPC))
     $              ,i=1,NDimensionGrid(idxReactionMeas(iPC)))

               write(52,'(4ES12.4)', advance='no')
     $              F2VAVE(if2),
     $              F2EstAve(if2), F2EsyAve(if2), F2EAVE(if2)

               if (WriteOriginal) then
                  do i=1,NInputFiles
                     found = .false.
                     do i2=1,NMeasF2(if2)
                        if (F2DataFile(if2,i2).eq.i) then
                           found = .true.
C Get the shifted data:
                           sum = F2TAB(if2,i2)
                           do isys = 1, NSYSTOT
                              sum = sum + SYSTAB(isys,if2,i2)
     $                             *SYSSH(isys)
                           enddo


                           write (52,'(4ES12.4)',advance='no')
     $                          F2Tab(if2,i2),
     $                          F2ETAB_TOT(if2,i2),
     $                          sum,
     $                          F2ETAB(if2,i2)
                        endif
                     enddo
                     if ( .not. found ) then
                        write (52,'(4ES12.4)',advance='no')
     $                       0., 0., 0., 0.
                     endif
                  enddo
               endif
               write (52,*)
            endif
         enddo
      enddo

      close(52)

      end subroutine writeSimpleOut


      subroutine outTable()
C------------------------------------------------------------
      implicit none
      include 'common.inc'
          logical flag
      real*8 uncor_syst
      integer isys, if2, isysOff

C
C Write in output using original syst. sources:
C
      open (51,file=trim(OutputFolder)//'/tab.dat',status='unknown')
      if (iOutput.eq.0) then
         do if2=1,NMeas
 1718       format (2F14.7,200F9.4)

C For now Y uses CM for the first experiment:
            uncor_syst = sqrt(f2estave(if2)**2
     $           -f2estaveTrue(if2)**2)

            write (51,1718) 
     $           ,f2vave(if2)
     $           ,f2estaveTrue(if2)
     $           ,uncor_syst
     $           ,(SystOrig(isys,if2),isys=1,NSYSTOT)
     $           ,(OSyst(if2,isysOff),isysOff=1,NSYSOTOT)
         enddo
C
C Write in output using using orthogonal syst. sources:
C

      elseif (iOutput.eq.1) then
         do if2=1,NMeas
                        
C For now Y uses CM for the first experiment:
            uncor_syst = sqrt(abs(f2estave(if2)**2
     $           -f2estaveTrue(if2)**2))

            write (51,1718) 
     $           ,f2vave(if2)
     $           ,f2estaveTrue(if2)
     $           ,uncor_syst
     $           ,(SystDiag(isys,if2),isys=1,NSYSTOT)
     $           ,(OSyst(if2,isysOff),isysOff=1,NSYSOTOT)
         enddo
      endif

      close(51)

C------------------------------------------------------------
      end subroutine


      subroutine  outOffsetSyst()
        implicit none
        include 'common.inc'
        integer isys, if2

        open (55,file=trim(OutputFolder)//'/offsettab.dat',
     $     status='unknown')

         do if2=1,NMeas
 1818       format (7F14.7,200F8.4)

            write (55,1818)
     $           ,f2vave(if2)
     $           ,(F2VaveOSyst(if2,isys),isys=1,2*NSYSOTOT)
         enddo

      end subroutine


      subroutine  outSystImpact()
        implicit none
        include 'common.inc'
        integer isys, if2

        open (55,file=trim(OutputFolder)//'/systimpact.dat',
     $     status='unknown')

         do if2=1,NMeas
 1818      format (200ES12.4)

           write (55,1818)
     $     ,f2vave(if2)
     $     ,((F2VaveSyst(if2,2*isys)-F2VaveSyst(if2,2*isys-1))*0.5,
     $     isys=1,NSYSTOT)
         enddo

      end subroutine

      subroutine  outToyMCStat()
        implicit none
        include 'common.inc'
        integer isys, if2

        open (55,file=trim(OutputFolder)//'/toymcstat.dat',
     $     status='unknown')

         do if2=1,NMeas
 1818       format (200ES12.4)

            write (55,1818)
     $           ,f2vave(if2)
     $           ,StatToyMC(if2)
         enddo

      end subroutine




C     Calculate and write pulls for all input data vs average
      Subroutine WritePulls()
      implicit none
      include 'common.inc'
      integer iFile,iF2,iexp, isys, ndf, i, iPC
      double precision pull, sum
C---------------------------------------------------------

      open(58,file=trim(OutputFolder)//'/chi2map.dat')
C Order using input data files
      do iFile=1,NInputFiles
         ndf  =  0   !> Count points in each file
         do if2=1,NMeas
            do iexp=1,NMeasF2(if2) !NMeas
               if (F2DataFile(if2,iexp).eq.IFile) then
                  ndf = ndf + 1
                  sum = F2TAB(if2,iexp)
                  do isys=1,NSYSTOT
                     sum = sum + SYSTAB(isys,if2,iexp)*SYSSH(isys)
                  enddo

C For pulls, divide by err^2 - err_ave^2 instead of err^2:

                  if (NMeasF2(if2).gt.1) then
                     pull = (F2VAVE(if2)-sum)/
     $                    sqrt(abs(F2ETAB(if2,iexp)**2
     $                    -F2EstAve(if2)**2))
                  else
                     pull = 0.
                  endif

                  iPC = idxProcessClass(if2)

                  write (58,'(I5,I5,10ES12.4)',advance='no')
     $                 ,ndf,NMeasF2(if2),
     $          ( GridPoints(idxGridMeas(if2),i,idxReactionMeas(iPC))
     $                 ,i=1,NDimensionGrid(idxReactionMeas(iPC)))
                  write (58,'(ES12.4,2X,(A))')
     $                 pull, Trim(InputFileNames(iFile))
               endif
            enddo
         enddo
      enddo

      close (58)
      end


C     Write Fitted Systematics, their shifts and pulls to separate file as tex or txt file
      Subroutine WriteSystShifts()
      implicit none
      include 'common.inc'
      integer j,renstat
      double precision pull
C     ---------------------------------------------------------

      if ( WriteSysTexTable ) then
         open (55,file=trim(OutputFolder)//'/sys.tex',status='unknown')

         if (modulo(nsystot,2).gt.0) then
            renstat = (nsystot+1)/2
            do j=1,renstat-1
               write (55,
     $'(i3'' &'',(A32)'' &'',F8.2'' &'',
     $F8.2'' &'',i3'' &'',(A32)'' &'',F8.2'' &'',F8.2'' \\'')')
     $              j,TRIM(ADJUSTL(SystematicName(j))),
     $              SYSSH(j),ErrSyst(j),
     $              j+renstat,TRIM(ADJUSTL(SystematicName(j+renstat))),
     $              SYSSH(j+renstat),ErrSyst(j+renstat)
            enddo
            do j=renstat,renstat
               write (55,
     $'(i3'' &'',(A32)'' &'',F8.2'' &'',
     $F8.2'' \\'',i3'' &'',(A32)'' &'',F8.2'' &'',F8.2'' \\'')')
     $              j,TRIM(ADJUSTL(SystematicName(j))),
     $              SYSSH(j),ErrSyst(j)
            enddo
         elseif (modulo(nsystot,2).eq.0) then
            renstat = nsystot/2
            do j=1,renstat
               write (55,
     $'(i3'' &'',(A32)'' &'',F8.2'' &'',
     $F8.2'' &'',i3'' &'',(A32)'' &'',F8.2'' &'',F8.2'' \\'')')
     $           j,TRIM(ADJUSTL(SystematicName(j))),
     $           SYSSH(j),ErrSyst(j),
     $           j+renstat,TRIM(ADJUSTL(SystematicName(j+renstat))),
     $           SYSSH(j+renstat),ErrSyst(j+renstat)
            enddo
         endif
         close(55)
      else
          open (55,file=trim(OutputFolder)//'/sys.txt',status='unknown')
          do j=1,nsystot
              pull = 0
              if(errsyst(j) .lt. 0.9999999 )then
                  pull = SYSSH(j) /
     $            sqrt(1 - errsyst(j)*errsyst(j))
              endif
              write (55,'(i3,(A32),F8.2,F8.2,F8.2,F8.2)')
     $        j,TRIM(ADJUSTL(SystematicName(j))),
     $        SYSSH(j),ErrSyst(j),pull
          enddo
          close(55)
      endif

      end

C     Print the eigenvalues and eigenvectors
      Subroutine PrintEigInfo(WWW,box3)
        implicit none
        include 'common.inc'
        integer isys1,isys2
        real*8 WWW(NSYSTMAX)
        real*8 Box3(NSYSTMAX,NSYSTMAX)

         open (51,file=trim(OutputFolder)//'/eigvalues.dat',
     $       status='unknown')
         do isys1=1,nsystot
            write (51,'(E16.8)') WWW(isys1)
         enddo
         close (51)
         open (51,file=trim(OutputFolder)//'/eigvectors.dat',
     $        status='unknown')
         do isys1=1,nsystot
            write (51,'(100E16.8)') (box3(isys2,isys1),isys2=1,nsystot)
         enddo
         close (51)
      end

C     Print veriance matrix
      Subroutine PrintMatrix(box3)
        implicit none
        include 'common.inc'
        integer isys1,isys2
        real*8 Box3(NSYSTMAX,NSYSTMAX)

         open (51,file=trim(OutputFolder)//'/matrix.dat',
     $       status='unknown')
         do isys1=1,nsystot
            write (51,'(100E16.8)') (box3(isys2,isys1),isys2=1,nsystot)
         enddo
         close (51)
      end
