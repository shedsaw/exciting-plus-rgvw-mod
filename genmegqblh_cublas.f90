subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr
use ISO_C_BINDING
use cublas_f
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
complex(8) b1(lmmaxapw*nufrmax),b2(lmmaxapw*nufrmax)

complex(8),dimension(:,:,:),pointer :: b1Batch, b2Batch
complex(8),dimension(:,:,:),pointer :: gntujuBatch
type(C_PTR) :: d_b1, d_b2, d_gntuju
type(C_PTR), dimension(:), pointer :: h_d_b1, h_d_b2, h_d_gntuju
type(C_PTR) :: handle,stream
integer :: sizeof_complex,sizeof_ptr,idx,stat,batch_count
integer(8) :: b1Size,b2Size,gntujuSize,bytes
parameter (sizeof_complex=16,sizeof_ptr=8)
wfsize=lmmaxapw*nufrmax*natmtot+ngknr2

allocate(wftmp1(wfsize,ngq(iq)))
allocate(wftmp2(wfsize,nstsv))
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)

batch_count = ngq(iq)*natmtot

!allocate(gntujuBatch(lmmaxapw*nufrmax, lmmaxapw*nufrmax, batch_count))
allocate(b1Batch(lmmaxapw*nufrmax, 1, batch_count))
allocate(b2Batch(lmmaxapw*nufrmax, 1, batch_count))
 
allocate(h_d_b1(batch_count))
allocate(h_d_b2(batch_count))
allocate(h_d_gntuju(batch_count))
b1Size = lmmaxapw*nufrmax*sizeof_complex
b2Size = b1Size
gntujuSize = lmmaxapw*nufrmax*b1Size

do i=1,batch_count
  stat = cudaMalloc(h_d_b1(i),b1Size)
  stat = cudaMalloc(h_d_b2(i),b2Size)
  stat = cudaMalloc(h_d_gntuju(i),gntujuSize)
enddo

bytes = batch_count*sizeof_ptr

stat = cudaMalloc(d_b1, bytes)
stat = cudaMalloc(d_b2, bytes)
stat = cudaMalloc(d_gntuju, bytes)
stat = cudaMemcpy(d_b1, C_LOC(h_d_b1(1)),bytes,cudaMemcpyHostToDevice); 
stat = cudaMemcpy(d_b2, C_LOC(h_d_b2(1)),bytes,cudaMemcpyHostToDevice); 
stat = cudaMemcpy(d_gntuju, C_LOC(h_d_gntuju(1)),bytes,cudaMemcpyHostToDevice); 

stat = cublasCreate(handle)
stat = cudaStreamCreate(stream)

write(*,*) 'ngkmax=',ngkmax
write(*,*) 'nstsv=',nstsv
write(*,*) 'iq=',iq
write(*,*) 'ngq(iq)=',ngq(iq)
write(*,*) 'ngrtot=',ngrtot
write(*,*) 'lmmaxapw=',lmmaxapw
write(*,*) 'nufrmax=',nufrmax
write(*,*) 'natmtot=',natmtot
write(*,*) 'wfsize=',wfsize
write(*,*) 'batch_count=',batch_count
write(*,*) 'b1Size=',b1Size
write(*,*) 'gntujuSize=',gntujuSize
write(*,*) 'bytes=',bytes
write(*,*) 'nspinor=',nspinor
write(*,*) 'ig=1..',ngq(iq)
write(*,*) 'ias=1,..',natmtot

do ispn1=1,nspinor
  if (expigqr22.eq.1) ispn2=ispn1
! index of the interband transitions
  i=1
! go through the interband transitions    
  do while (i.le.nmegqblh(ikloc))
! left <bra| state 
    ist1=bmegqblh(1,i,ikloc)
    wftmp1=zzero
    l1=.true.
    if (spinpol) then
      if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
    endif
    if (l1) then
      call timer_start(3)
      call papi_timer_start(pt_megqblh_mt)
      do ig=1,ngq(iq)
! precompute muffint-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
        do ias=1,natmtot
          idx = (ig-1)*natmtot + ias 
          b1Batch(:,1,idx)=dconjg(wfsvmt1(:,ias,ispn1,ist1)*sfacgq(ig,ias))
          ic=ias2ic(ias)
          b2Batch(:,1,idx)=zzero
          !b2=zzero
          !do j=1,ngntuju(ic,ig)
          !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
          !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
          !enddo

          !call memcopy(b1(1), b1Batch(1,1,idx), b1Size) 
          !call memcopy(gntuju(1,1,ic,ig), gntujuBatch(1,1,idx),gntujuSize) 
          stat  = cublasSetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
           &C_LOC(b1Batch(1,1,idx)),lmmaxapw*nufrmax, h_d_b1(idx), &
           &lmmaxapw*nufrmax, stream)
          stat = cublasSetMatrixAsync(lmmaxapw*nufrmax,lmmaxapw*nufrmax,&
           &sizeof_complex,C_LOC(gntuju(1,1,ic,ig)),lmmaxapw*nufrmax,&
           &h_d_gntuju(idx),lmmaxapw*nufrmax, stream)

          stat  = cublasSetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
           &C_LOC(b2Batch(1,1,idx)),lmmaxapw*nufrmax, h_d_b2(idx), &
           &lmmaxapw*nufrmax, stream)
          !!call zgemm('N','N',lmmaxapw*nufrmax,1,lmmaxapw*nufrmax,&
          !!  &zone,gntuju(1,1,ic,ig),lmmaxapw*nufrmax,b1,lmmaxapw*nufrmax,&
          !!  &zzero,b2,lmmaxapw*nufrmax)
          !!wftmp1((ias-1)*lmmaxapw*nufrmax+1:ias*lmmaxapw*nufrmax,ig)=b2(:)
        enddo !ias
      enddo !ig  
      !do ig=1, batch_count
      !  cublasSetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
      !   &C_LOC(b1Batch(1,1,ig)),lmmaxapw*nufrmax, h_d_b1, &
      !   &lmmaxapw*nufrmax, stream)
      !  cublasSetMatrixAsync(lmmaxapw*nufrmax,lmmaxapw*nufrmax,sizeof_complex,&
      !   &C_LOC(gntujuBatch(1,1,ig)),lmmaxapw*nufrmax, h_d_gntuju, &
      !   &lmmaxapw*nufrmax, stream)
      !enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timer_stop(3)
      call papi_timer_stop(pt_megqblh_mt)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   call cudaStreamSynchronize(stream)

    stat = cublasZgemmBatched(handle,CUBLAS_OP_N,CUBLAS_OP_N,lmmaxapw*nufrmax,&
       &1,lmmaxapw*nufrmax,zone,d_gntuju,lmmaxapw*nufrmax,d_b1,&
       &lmmaxapw*nufrmax,zzero,d_b2,lmmaxapw*nufrmax,batch_count)
    

    do ig=1,batch_count
      stat = cublasGetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
        &C_LOC(b2Batch(1,1,ig)),lmmaxapw*nufrmax, h_d_b2(ig), &
        &lmmaxapw*nufrmax, stream)
    enddo
 
    call timer_start(5)
    n1=0
! collect right |ket> states into matrix wftmp2
    do while ((i+n1).le.nmegqblh(ikloc))
      if (bmegqblh(1,i+n1,ikloc).ne.bmegqblh(1,i,ikloc)) exit
      ist2=bmegqblh(2,i+n1,ikloc)
      n1=n1+1
      call memcopy(wfsvmt2(1,1,1,ispn2,ist2),wftmp2(1,n1),16*lmmaxapw*nufrmax*natmtot)
      call memcopy(wfsvit2(1,ispn2,ist2),wftmp2(lmmaxapw*nufrmax*natmtot+1,n1),16*ngknr2)
    enddo !while

    call cudaStreamSynchronize(stream)
      do ig=1,ngq(iq)
        do ias=1,natmtot
           idx = (ig-1)*natmtot + ias
           wftmp1((ias-1)*lmmaxapw*nufrmax+1:ias*lmmaxapw*nufrmax,ig)=b2Batch(:,1,idx)
        enddo
      enddo

! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)
    call zgemm('T','N',n1,ngq(iq),wfsize,zone,wftmp2,wfsize,wftmp1,wfsize,&
      &zone,megqblh(i,1,ikloc),nstsv*nstsv)
    i=i+n1
    call timer_stop(5)
  enddo !while
enddo !ispn
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)

do i=1,batch_count
  call cudaFree(h_d_b1(i))
  call cudaFree(h_d_b2(i))
  call cudaFree(h_d_gntuju(i))
enddo
deallocate(b1Batch)
deallocate(b2Batch)
!deallocate(gntujuBatch)
call cublasDestroy(handle)
stat =  cudaStreamDestroy(stream)
call papi_timer_stop(pt_megqblh)

return
end
