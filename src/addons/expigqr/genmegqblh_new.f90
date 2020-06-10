subroutine genmegqblh_new(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
integer, intent(in) :: iq
integer, intent(in) :: ikloc
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvmt1(lmmaxapw*nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxapw,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)

integer wfsize
integer ivg1(3)
integer i,j,ik,jk,igkq,n1,ispn1,ispn2,ist1,ist2,ic
integer ig,ig1,ig2,ias,ifg,ir
logical l1
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: wfir1(:)

!--begin Convert to true ZGEMM
  INTEGER :: nmt                    ! Number of muffin-tin elements
  INTEGER, PARAMETER :: nb = 64     ! Block size for ZGEMM batching
  INTEGER :: k1, k2, ki, nsize      ! Dummy variables for batching
  COMPLEX(KIND((0.D0,1.D0))), &
    DIMENSION( lmmaxapw*nufrmax, nb, natmtot, ngq(iq) ) :: b1, b2
!--end Convert to true ZGEMM

#ifdef _DEBUG_bmegqblh_
  INTEGER :: dbgcnt, dbgunit
#endif /* _DEBUG_bmegqblh_ */

INTEGER :: idxhiband, iband, ntran, idxtran
EXTERNAL :: zcopy

wfsize=lmmaxapw*nufrmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngq(iq))) ! TODO: Change dimensions appropriately
allocate(wftmp2(wfsize,nstsv))   ! TODO: Check contiguity of ZCOPY transfers
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)

#ifdef _DEBUG_bmegqblh_
  dbgunit = 1000+iproc ! Make sure this matches the definition in mod_expigqr
  dbgcnt=1
  WRITE( dbgunit, '(A,I3,A,I5)' ) 'nmegqblh(ikloc=', ikloc, ') = ', nmegqblh(ikloc)
#endif /* _DEBUG_bmegqblh_ */

do ispn1=1,nspinor
  if (expigqr22.eq.1) ispn2=ispn1

  ! Load number of matching |ist2=n'> ket states for each <ist1=n| bra
  IF( ltranconst ) ntran = ntranblhloc(ikloc)

  ! Convert to the corresponding ist1 loop in getmeidx() line 56-149
  ! as stored into idxhibandblh(ikloc=1:nkptnr) at getmeidx() line 155,
  ! skipping as necessary (with a warning message... it should NOT happen!)
  idxhiband = idxhibandblhloc(ikloc)
  IF( idxhiband == 0 ) THEN
     ! Unit 151 is either 'CRPA.OUT' or 'RESPONSE.OUT'i
     WRITE(151, '( "Warning[genmegqblh]: highest band is zero for iq=", &
                 &I6, " ikloc=", I6, " ispn1=", I1 )' ) iq, ikloc, ispn1
     CYCLE
  END IF

!--begin Convert to true ZGEMM

  call timer_start(3)
  call papi_timer_start(pt_megqblh_mt)

  ! Note that the loop order has been switched
  ! such that iband loop is now the innermost loop

  ! Number of muffin-tin elements
  nmt = lmmaxapw*nufrmax

  ! Batching by block size nb
  DO k1 = 1, idxhiband, nb
     k2 = MIN( iband, k1+nb-1 )
     nsize = k2 - k1 + 1

     !$OMP PARALLEL PRIVATE(b1,b2)
     b1(:,:,:,:) = zzero
     b2(:,:,:,:) = zzero

     !$OMP DO COLLAPSE(3) PRIVATE(iband,i,ist1,l1)
     do ig=1,ngq(iq)
        do ias=1,natmtot

           ! Loop for a single batch
           DO ki = 1, nsize

              iband = k1 + ki - 1
              i = idxtranblhloc( iband, ikloc )
              ist1 = bmegqblh(1,i,ikloc)

              ! TODO: Dump bmegqblh and inspect, then get rid of the l1 check
              l1=.true.
              if (spinpol) then
                 if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
              endif
              if (l1) then
                 ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
                 b1(:,ki,ias,ig) = DCONJG( wfsvmt1(:,ias,ispn1,ist1) * &
                                           sfacgq(ig,ias) )
              END IF ! l1

           END DO ! ki

        END DO ! ias
     END DO ! ig
     !$OMP END DO

     ! Original code for historical purpose
     !do j=1,ngntuju(ic,ig)
     !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
     !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
     !enddo

     !$OMP DO COLLAPSE(2) PRIVATE(ic)
     do ig=1,ngq(iq)
        do ias=1,natmtot
           ic = ias2ic(ias)

           ! Perform ZGEMM by batch
           CALL zgemm( 'N', 'N', nmt, nsize, nmt, &
                       zone,  gntuju(1:nmt,1:nmt,ic,ig), nmt, &
                              b1(1:nmt,1:nsize,ias,ig),  nmt, &
                       zzero, b2(1:nmt,1:nsize,ias,ig),  nmt )

        enddo !ias
     enddo !ig
     !$OMP END DO
     !$OMP END PARALLEL

  END DO ! k1

  call timer_stop(3)
  call papi_timer_stop(pt_megqblh_mt)

  ! Start the bounded do loop
  DO iband = 1, idxhiband

! left <bra| state 
     wftmp1=zzero

     ! The starting point of the index "i" for accessing bmegqblh(:,i,:)
     ! for each iband and ikloc was stored as idxtranblhloc
     i = idxtranblhloc( iband, ikloc )
     ist1 = bmegqblh(1,i,ikloc)

     ! Fetch muffin-tin part
     ! TODO: Split wftmp1 into muffin-tin and interstitial parts
     wftmp1( (ias-1)*nmt+1:ias*nmt, ig ) = b2( :, iband, ias, ig )

     ! Same as above
     l1=.true.
     if (spinpol) then
        if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
     endif

     if (l1) then

!--end Convert to true ZGEMM

! interstitial part
      call papi_timer_start(pt_megqblh_it)
      call timer_start(4)
      wfir1=zzero
      do ig1=1,ngknr1
        ifg=igfft(igkignr1(ig1))
        wfir1(ifg)=wfsvit1(ig1,ispn1,ist1)
      enddo
      call zfftifc(3,ngrid,1,wfir1)
      do ir=1,ngrtot
        wfir1(ir)=wfir1(ir)*cfunir(ir)
      enddo
      call zfftifc(3,ngrid,-1,wfir1)
      do ig=1,ngq(iq)
        do ig2=1,ngknr2
! G1=G2-G-Gkq
          ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
          ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
          wftmp1(lmmaxapw*nufrmax*natmtot+ig2,ig)=dconjg(wfir1(ifg))
        enddo
      enddo
      call timer_stop(4)      
      call papi_timer_stop(pt_megqblh_it)
    endif !l1
    call timer_start(5)

#ifdef _DEBUG_bmegqblh_
IF( ntran > 0 ) THEN
  WRITE( dbgunit, '(7(1X,I5))' ) dbgcnt, ikloc, iq, iband, i, ntran, i+ntran-1
  dbgcnt = dbgcnt + 1
END IF
#endif // _DEBUG_bmegqblh_

    ! Note: seeing the pattern, this shouldn't happen, but just in case
    IF( .NOT. ltranconst ) THEN
       IF( iband == idxhiband ) THEN
          ntran = nmegqblh(ikloc) - idxtranblhloc(idxhiband,ikloc) + 1
       ELSE
          ntran = idxtranblhloc(iband+1,ikloc) - idxtranblhloc(iband,ikloc)
       END IF
    END IF ! ltranconst

! collect right |ket> states into matrix wftmp2
    ! Note: ntran can be zero (zero-trip loop)
    DO n1 = 1, ntran

       ist2 = bmegqblh(2,i+n1-1,ikloc) ! Now n1 starts from 1 instead of 0

       ! Following Ed's advice, use ZCOPY() from BLAS instead of memcopy
       ! TODO: check whether it's better to transfer all in one go
       !       or overlap computation & data movement

       ! Muffin tin
       CALL zcopy( lmmaxapw*nufrmax*natmtot, &
                   wfsvmt2(1,1,1,ispn2,ist2), 1, &
                   wftmp2(1,n1), 1 )
       ! Interstitial
       CALL zcopy( ngknr2, &
                   wfsvit2(1,ispn2,ist2), 1, &
                   wftmp2(lmmaxapw*nufrmax*natmtot+1,n1), 1 )

    END DO ! n1; replaced do while loop (i+n1) <= nmegqblh(ikloc)

! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)

    ! This particular ZGEMM() call corresponds with line 9 of Algorithm 2
    ! in the Gordon Bell paper
    CALL zgemm( 'T', 'N', ntran, ngq(iq), wfsize, zone, &
                wftmp2, wfsize, wftmp1, wfsize, zone, &
                megqblh_new(i,1,ikloc), nstsv*nstsv )

    ! No need to add n1 to i anymore to move on to the next <nk| bra
    ! since it is already stored as ntranblhloc

    CALL timer_stop(5) ! Same as before

 END DO ! iband; replaces do while loop i <= nmegqblh(ikloc)

#ifdef _DEBUG_bmegqblh_
     WRITE( dbgunit, '(A,I3)') 'highest band = ', idxhiband
#endif // _DEBUG_bmegqblh_

!--end Convert do while into bounded do loop

enddo !ispn
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)

call papi_timer_stop(pt_megqblh)

return
end subroutine genmegqblh_new
