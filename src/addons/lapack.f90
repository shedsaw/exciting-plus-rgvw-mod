subroutine invzge(mtrx,ndim)
implicit none
! passed var
integer, intent(in) :: ndim
complex(8), intent(out) :: mtrx(ndim,ndim)
! local var
integer lwork,nb,info, dummy
real*8 ,allocatable :: work(:)
integer ,allocatable :: ipiv(:)

!integer, external :: ilaenv
!nb=ilaenv(1,'zgetri','U',ndim,-1,-1,-1)
!lwork=ndim*nb

!-- IBM ESSL fix
CALL zgetri( ndim, mtrx, ndim, dummy, lwork, -1, info)
!--

allocate(work(2*lwork),ipiv(ndim))
call zgetrf(ndim,ndim,mtrx,ndim,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invzge) zgetrf returned ",I4)')info
  call pstop
endif
call zgetri(ndim,mtrx,ndim,ipiv,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invzge): zgetri returned ",I4)')info
  call pstop
endif
deallocate(work,ipiv)
end subroutine invzge

subroutine invdsy(n,mtrx)
implicit none
integer, intent(in) :: n
real(8), intent(inout) :: mtrx(n,n)

real(8) t1
integer lwork,info
integer, allocatable :: ipiv(:)
real(8), allocatable :: work(:)

!-- IBM ESSL fix
INTEGER :: dummy
REAL(KIND=KIND(1.D0)) :: query
!--

!allocate(ipiv(n))
!lwork=-1
!call dsytrf('U',n,mtrx,n,ipiv,t1,lwork,info)
!lwork=int(t1)+1
!allocate(work(lwork))
!call dsytrf('U',n,mtrx,n,ipiv,work,lwork,info)
!if (info.ne.0) then
!  write(*,*)
!  write(*,'("Error(invdsy) dsytrf returned ",I4)')info
!  call pstop
!endif
!call dsytri('U',n,mtrx,n,ipiv,work,info)
!if (info.ne.0) then
!  write(*,*)
!  write(*,'("Error(invdsy) dsytri returned ",I4)')info
!  call pstop
!endif
!deallocate(ipiv,work)

!-- IBM ESSL fix
! For some reason, ESSL doesn't provide *sytri
! Use *getrf and *getri instead

! Query workspace
CALL DGETRI( n, mtrx, n, dummy, query, -1, info )
lwork = CEILING(query)

! LU factorization using DGETRF
ALLOCATE( ipiv(n) )
CALL DGETRF( n, n, mtrx, n, ipiv, info )
IF ( info /= 0 ) THEN
  WRITE(*,*)
  WRITE(*,'("Error(invdsy) dgetrf returned ",I4)') info
  DEALLOCATE( ipiv )
  CALL pstop
END IF

! Allocate workspace
ALLOCATE( work(lwork) )

! Invert matrix using DGETRI
CALL DGETRI( n, mtrx, n, ipiv, work, lwork, info )
! Regardless of results, ipiv and work are no longer needed at this point
DEALLOCATE( ipiv )
DEALLOCATE( work )
IF ( info /= 0 ) THEN
  WRITE(*,*)
  WRITE(*,'("Error(invdsy) dgetri returned ",I4)') info
  CALL pstop
END IF

!--
return
end subroutine invdsy

subroutine diagzhe(ndim,mtrx,evalue)
implicit none
integer ,intent(in) :: ndim
complex(8), intent(out) :: mtrx(ndim,ndim)
real(8), intent(out) :: evalue(ndim)
integer nb,lwork,inf, dummy
real*8, allocatable :: work(:),rwork(:)

!integer, external :: ilaenv
!nb=ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
!lwork=(nb+1)*ndim
!-- IBM ESSL fix
CALL zheev( 'V', 'U', ndim, mtrx, evalue, lwork, -1, dummy, inf )
!--

allocate(work(lwork*2))
allocate(rwork(3*ndim+2))
call zheev('V','U',ndim,mtrx,ndim,evalue,work,lwork,rwork,inf)
if (inf.ne.0) then
  write(*,*)
  write(*,'("Error(diagzhe) zheev returned ",I4)')inf
  call pstop
endif
deallocate(work,rwork)
end

subroutine diagdsy(n,mtrx,eval)
implicit none
integer, intent(in) :: n
real(8), intent(inout) :: mtrx(n,n)
real(8), intent(out) :: eval(n)
integer lwork,info
real(8), allocatable :: work(:)
real(8) t1
lwork=-1
call dsyev('V','U',n,mtrx,n,eval,t1,lwork,info)
lwork=int(t1)+1
allocate(work(lwork))
call dsyev('V','U',n,mtrx,n,eval,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Warning(diagdsy) dsyev returned ",I4)')info
  write(*,*)
endif
deallocate(work)
return
end subroutine diagdsy

subroutine isqrtzhe(ndim,mtrx,ierr)
use mod_mpi_grid
implicit   none
! arguments
integer, intent(in) :: ndim
integer, intent(out) :: ierr
complex(8), intent(inout) :: mtrx(ndim,ndim)

integer nb,lwork,info,i,j,n, dummy
real(8), allocatable :: work(:),rwork(:),ev(:),ev1(:)
complex(8), allocatable :: z1(:,:)

!integer, external :: ilaenv
!nb=ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
!lwork=(nb+1)*ndim

!-- IBM ESSL fix
CALL zheev( 'V', 'U', ndim, mtrx, ev1, lwork, -1, dummy, info )
!--

allocate(z1(ndim,ndim))
allocate(work(lwork*2))
allocate(rwork(3*ndim+2))
allocate(ev(ndim))
allocate(ev1(ndim))

!if (mpi_grid_root()) write(*,*)mtrx

call zheev('V','U',ndim,mtrx,ndim,ev1,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(isqrtzhe): zheev returned ",I4)')info
  write(*,*)
  call pstop
endif

!if (mpi_grid_root()) write(*,*)ev1
!if (mpi_grid_root()) write(*,*)mtrx

ierr=0
do i=1,ndim
  if (ev1(i).lt.1.d-12) then
    ierr=i
    ev(i)=0.d0
  else
    ev(i)=1.d0/dsqrt(ev1(i))
  endif
enddo

z1(:,:)=mtrx(:,:)
mtrx=dcmplx(0.d0,0.d0)
do i=1,ndim
  do j=1,ndim
    do n=1,ndim
      mtrx(i,j)=mtrx(i,j)+z1(i,n)*dconjg(z1(j,n))*ev(n) 
    enddo
  enddo
enddo
deallocate(z1,work,rwork,ev,ev1)
return
end subroutine isqrtzhe

subroutine isqrtdsy(n,mtrx,ierr)
use mod_mpi_grid
implicit   none
! arguments
integer, intent(in) :: n
integer, intent(out) :: ierr
real(8), intent(inout) :: mtrx(n,n)
!
integer i,j,k
real(8), allocatable :: eval(:),eval_isq(:)
real(8), allocatable :: evec(:,:)

allocate(evec(n,n))
allocate(eval(n),eval_isq(n))
call diagdsy(n,mtrx,eval)

ierr=0
do i=1,n
  if (eval(i).lt.1.d-12) then
    ierr=i
    eval_isq(i)=0.d0
  else
    eval_isq(i)=1.d0/dsqrt(eval(i))
  endif
enddo
evec(:,:)=mtrx(:,:)
mtrx=0.d0
do i=1,n
  do j=1,n
    do k=1,n
      mtrx(i,j)=mtrx(i,j)+evec(i,k)*evec(j,k)*eval_isq(k) 
    enddo
  enddo
enddo
deallocate(evec,eval,eval_isq)
return
end subroutine isqrtdsy

subroutine diagzge(ndim,mtrx,evalue)
implicit   none
integer, intent(in) :: ndim
complex*16, intent(inout) :: mtrx(ndim,ndim)
complex*16, intent(out) :: evalue(ndim)
integer lwork,inf
complex(8) zt1
real*8, allocatable :: rwork(:)
complex(8), allocatable :: work(:)
complex(8), allocatable :: evec(:,:)

!-- IBM ESSL fix
!integer, external :: ilaenv
!--

lwork=-1
call zgeev('N','V',ndim,mtrx,ndim,evalue,evec,ndim,evec,ndim,zt1,lwork,rwork,inf)
lwork=dble(zt1)+1
allocate(work(lwork))
allocate(rwork(2*ndim))
allocate(evec(ndim,ndim))
call zgeev('N','V',ndim,mtrx,ndim,evalue,evec,ndim,evec,ndim,work,lwork,rwork,inf)
if (inf.ne.0) then
  write(*,'("Error(diagzge) zgeev returned ",I4)')inf
  call pstop
endif
mtrx=evec
deallocate(work,rwork,evec)
end subroutine diagzge

subroutine diagzheg(n,nv,ld,etol,a,b,eval,evec)
implicit none
integer, intent(in) :: n
integer, intent(in) :: nv
integer, intent(in) :: ld
real(8), intent(in) :: etol
complex(8), intent(inout) :: a(n,n)
complex(8), intent(inout) :: b(n,n)
real(8), intent(out) :: eval(nv)
complex(8), intent(out) :: evec(ld,nv)
!
integer m,info,i,nb,lwork, dummy1, dummy2
real(8) vl,vu
integer, allocatable :: iwork(:)
integer, allocatable :: ifail(:)
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)

!integer, external :: ilaenv 
!nb=ilaenv(1,'ZHETRD','U',n,-1,-1,-1)
!lwork=(nb+1)*n

!-- IBM ESSL fix
CALL zhegvx( 1, 'V', 'I', 'U', n, a, n, b, n, vl, vu, 1, nv, etol, m, w, evec,&
             ld, lwork, -1, dummy1, dummy2, ifail, info )
!--

allocate(iwork(5*n))
allocate(ifail(n))
allocate(w(n))
allocate(rwork(7*n))
allocate(work(lwork))
call zhegvx(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld,&
  work,lwork,rwork,iwork,ifail,info)
eval(1:nv)=w(1:nv)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(diagzheg): diagonalisation failed")')
  write(*,'("  zhegvx returned info = ",I8)') info
  if (info.gt.n) then
    i=info-n
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') n
    write(*,*)
  end if
  call pstop
end if
deallocate(iwork,ifail,w,rwork,work)
return
end subroutine diagzheg
