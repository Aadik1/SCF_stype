module GreensFunctions
  use DefineHamiltonian
  use OMP_LIB
  implicit none
  integer :: INFO

  real*8, parameter :: epsilon = 1e-9
  real*8, allocatable, dimension(:) :: Hub, omega, Ev
  real*8 :: pulay, dw, up
  
  complex*16, allocatable, dimension(:,:) :: GammaL, GammaR, Eigenvec, G_nil
!  complex*16, allocatable, dimension(:,:) ::  SigmaL, Sigma1, SigmaR
  complex*16, allocatable, dimension(:,:) :: work1, work2, work3, work4

  type :: GF
     complex*16, allocatable, dimension(:,:,:) :: r, a, L, G
  end type GF
  type(GF) :: GF0
  
  type :: GF_full
     complex*16, allocatable, dimension(:,:,:) :: R, A, L, G
  end type GF_full
  type(GF_full) :: GFf

contains 
  
   !.....fermi-dirac distribution
    real*8 function fermi_dist(w, V)
    implicit none
    real*8 :: arg, w, V

    arg = (w - mu + V/hbar)*beta
    fermi_dist =  1.d0/(exp(hbar*arg) +1.d0)
    
  end function fermi_dist

!====================================================
!========== Self consistency field calculations =====
!====================================================

!.............need G0f%L to calculate GL_of_0
!.............Calculates GFs at every omega and Voltage simultaneously 
subroutine SCF_GFs(Volt,first)
  implicit none
  integer :: iw, iteration,i, Vname, wheel, it
  real*8 :: Volt, err, diff
  character(len=30) :: fn, fn1
  logical :: first

   iteration = 0

  write(22,*) '........SCF Calculations at Voltage:', Volt, '..........'

  if (first) then 
     call G0_R_A()
     call G0_L_G(Volt)
  end if

  print *, '>>>>>>>>>>VOLTAGE:', Volt

  Vname = abs(Volt)
  write(fn1,'(i0)') Vname
  if (Volt .ge. 0) then 
     open(17,file='err_V_'//trim(fn1)//'.dat',status='unknown')
  else
     open(17,file='err_V_n'//trim(fn1)//'.dat',status='unknown')
  end if

  DO
     iteration = iteration + 1
     write(*,*) '.... ITERATION = ',iteration,' ....'

!.......real variable interactions turns off the Interaction component of the sigmas 
!.................full Gr and Ga, Eq. (5) and (6)

     call GL_of_0()

     !$OMP PARALLEL DO &
     !$OMP& PRIVATE(iw, INFO)
     
     do iw = 1, N_of_w
        call G_full(iw, Volt)
     end do
     !$OMP END PARALLEL DO
     
     err=0.0d0
     do iw = 1, N_of_w
        do i=1,Natoms
           diff=2.d0*hbar*(AIMAG(GFf%R(i,i,iw))-AIMAG(GF0%R(i,i,iw)))
           err=err +diff*diff
        end do
     end do
     write(*,*) 'err = ',sqrt(err)

     GF0%R = pulay*GFf%R + (1.0d0-pulay)*GF0%R
     GF0%L = pulay*GFf%L + (1.0d0-pulay)*GF0%L
     
     !...... do the advanced and greater components
     
     do iw=1,N_of_w
        work1=GF0%r(:,:,iw) 
        call Hermitian_Conjg(work1, Natoms, work2)
        GF0%a(:,:,iw)=work2
     end do
     GF0%G = GF0%L + GF0%R - GF0%A

  !____________ useful if one would like to check the convergence     
     !     write(13,*) iteration, current(Volt)

     call save_GFs()
 
     if (sqrt(err) .lt. epsilon .or. order .eq. 0) then
        write(*,*)'... REACHED REQUIRED ACCURACY ...'
        exit
     end if
    
     
  END DO
  close(17)
end subroutine SCF_GFs


!=====================================================
!================== Full GFs =========================
!===================================================== 
subroutine G_full(iw, Volt) !... Full Greens function, leaves Retarded and Advanced in the work arrays, application of Eq. (16) and (17), but with the full Sigmas, Eq. (3), (7) and (8) in CHE
  implicit none
  integer :: i, j, iw, sp, sp1, ii, jj, N
  real*8 :: Volt, w 
  complex*16 :: SigL, SigG, Omr
  complex*16, allocatable, dimension(:,:) ::  SigmaL, Sigma1, SigmaR, SigmaG,  w1, w2

  allocate(SigmaL(Natoms, Natoms), SigmaR(Natoms, Natoms), SigmaG(Natoms, Natoms), Sigma1(Natoms, Natoms))
  allocate(w1(Natoms, Natoms), w2(Natoms, NAtoms))
  !............full SigmaR due to interactions Eq. (7)
  SigmaR = (0.d0, 0.d0); SigmaL = (0.d0, 0.d0); SigmaG = (0.d0, 0.d0)
  
  if (order .eq. 1) then
     call first_order_sigma(Sigma1)
     SigmaR = Sigma1
  else if (order .eq. 2) then
     call first_order_sigma(Sigma1)

     do i = 1, Natoms, 2
        do sp = 0, 1
           ii = i +sp
           
           do j = 1, Natoms, 2
              do sp1 = 0, 1
                 jj= j + sp1
                 
                 Omr = Omega_r(i,j,sp,sp1,iw)
                 SigmaR(ii,jj) =  Sigma1(ii,jj)  + (Hub(j)*Hub(i)*OmR)*hbar**2
                 
              end do
           end do
           
        end do
     end do
 
     !..............full SigmaL, Eq. (3) and (4)     
     !.....Interaction contribution of both Sigmas     

     do i = 1, Natoms, 2
        do sp = 0, 1
           ii = i +sp
           
           do j = 1, Natoms, 2
              do sp1 = 0, 1
                 jj= j + sp1                 

                 call int_SigLnG(i,j, sp, sp1, iw, SigL, SigG)
                 SigmaL(ii,jj) = Hub(j)*Hub(i)*SigL*hbar**2
                ! SigmaG(ii,jj) = Hub(j)*Hub(i)*SigG*hbar**2
                 
              end do
           end do
           
        end do
     end do
  end if
  
  !............full Gr and Ga, Eq. (5) and (6)
  
  w = omega(iw)
  w1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) - SigmaR ! LK <== must be + for emb and minus for interaction sigma
  do i = 1 , Natoms
     w1(i,i) = w1(i,i) + hbar*w !(w+im*0.01)
  end do
  
  call Inverse_complex(Natoms, w1, info)
  call Hermitian_Conjg(w1, Natoms, w2)
  
  GFf%R(:,:,iw) = w1; GFf%A(:,:,iw) = w2

!.....Embedding contribution of both Sigmas

  SigmaL =  im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)/hbar + SigmaL 
  !SigmaG =  im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)/hbar + SigmaG

  !.............full GL and GG, Eq. (16) and (17)
  
  GFf%L(:,:,iw) = matmul(matmul(w1, SigmaL), w2) !.. GL = Gr * SigmaL * Ga
  !GFf%G(:,:,iw) = matmul(matmul(w1, SigmaG), w2) !.. GG = Gr * SigmaG * Ga
  GFf%G(:,:,iw) = GFf%L(:,:,iw) + GFf%R(:,:,iw) - GFf%A(:,:,iw)

  deallocate(SigmaL, Sigma1, SigmaR, SigmaG,  w1, w2)
end subroutine G_full

!=====================================================
!======== Non-interacting GFs ========================
!=====================================================  
  
!......Calculates lesser Greens function at time = 0
subroutine GL_of_0()
  !.....Lesser Greens function at time = 0
  !.....ei - eigenvalues of the Hamiltonian 
  implicit none
  integer :: i, j, k1, s1, s2
  complex*16 :: s
  real*8 :: pp
  
  pp=delta/(2.d0*pi)
  G_nil = (0.d0, 0.d0)
  do i = 1, Natoms,2
     do s1=0,1
        do s2=0,1
           
           s =(0.d0, 0.d0)
           do k1 = 1, N_of_w
              s = s + GF0%L(i+s1,i+s2,k1)
           end do
           G_nil(i+s1,i+s2)=s*pp
           
        end do
     end do
  end do

end subroutine GL_of_0

subroutine G0_R_A()
  !............non-interacting Greens functions: GR and Ga,  Eq. (5) and Eq. (6) in 'Current_Hubbard_Equations' document (CHE)
  implicit none
  integer :: j, i
  real*8 :: w
  
    do j = 1, N_of_w
       work1 = -H + 0.5d0*(im/hbar) * (GammaL + GammaR) !LK <========= must be +
       w = omega(j)
       do i = 1 , Natoms
          work1(i,i) = work1(i,i) + hbar*w !(w +im*0.01)
        end do
       
       call Inverse_complex(Natoms, work1, info)
       call Hermitian_Conjg(work1, Natoms, work2)
       
       GF0%r(:,:,j) = work1
       GF0%a(:,:,j) = work2
    end do

end subroutine G0_R_A

 subroutine G0_L_G(Volt)
!............non-interacting Greens functions: G> and G< for all omega on the grid, Eq. (16) and (17) in CHE
    implicit none
    real*8 :: Volt, w
    integer :: j, i
    
    work1 = (0.d0, 0.d0) ; work2 =(0.d0, 0.d0) ; work3 = (0.d0, 0.d0); work4 = (0.d0, 0.d0)
    do j = 1 , N_of_w
       w = omega(j)
       work1 = GF0%r(:,:,j) 
       work2 = GF0%a(:,:,j)

       work3 = matmul(matmul(work1, im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)), work2) 
      ! work4 = matmul(matmul(work1, im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)), work2)
       GF0%L(:,:,j) = work3
      ! GF0%G(:,:,j) = work4
       GF0%G(:,:,j) =  GF0%L(:,:,j) + GF0%R(:,:,j) - GF0%A(:,:,j)
    end do
    
  end subroutine G0_L_G 


  !=====================================================
!========Calcualtions needed for full GFs=============
!=====================================================

  subroutine first_order_sigma(Sigma1)
    implicit none
    integer :: i, s, s1, N 
    complex*16 :: Hartree
    complex*16 :: Sigma1(:,:)
    
    Sigma1 = (0.d0, 0.d0); Hartree = (0.d0, 0.d0)
    
    do i = 1, Natoms, 2 !..orbitals 
       
       Hartree = G_nil(i,i)+G_nil(i+1,i+1)
       !.. spins 
       do s = 0, 1
          do s1 = 0, 1
             if (s .eq. s1) then 
                Sigma1(i+s, i+s1) = im*hbar*Hub(i)*(G_nil(i+s, i+s1) - Hartree)
             else
                Sigma1(i+s, i+s1) = im*hbar*Hub(i)*G_nil(i+s, i+s1)
             end if
          end do
       end do
       
    end do
    
  end subroutine first_order_sigma
  
  !......................Calculation of Omega terms for the self-energies, Eq. (9) in CHE
complex*16 function Omega_r(i, j, sp, sp1, iw)
  implicit none
  integer :: i, j, iw, k_1, k_2, k_3, m, n, s, s1, sp, sp1, ii, jj
  complex*16 :: Omr
  real*8 :: pp
  
  ii = i +sp; jj= j + sp1
  
  Omr = (0.d0, 0.d0)
  do k_1 = 1, N_of_w
     do k_2 = 1, N_of_w
        k_3 = iw- k_1 +k_2

        if (k_3 .ge. 1 .and. k_3 .le. N_of_w) then       
           do s = 0, 1 !...Sum over orbitals
              do s1 = 0, 1
                 m = i+s1; n= j+s
                 !.....both second order diagram contributions
                 !_________ 3rd diagram
                 Omr = Omr + GF0%r(ii,jj,k_1)*GF0%L(n,m,k_2)*GF0%G(m,n,k_3) & 
                      + GF0%L(ii,jj,k_1)*GF0%L(n,m,k_2)*GF0%r(m,n,k_3) & 
                      + GF0%L(ii,jj,k_1)*GF0%a(n,m,k_2)*GF0%L(m,n,k_3) &
                      
                      - GF0%r(ii,n,k_1)*GF0%L(n,m,k_2)*GF0%G(m,jj,k_3) &
                      - GF0%L(ii,n,k_1)*GF0%L(n,m,k_2)*GF0%r(m,jj,k_3) &
                      - GF0%L(ii,n,k_1)*GF0%a(n,m,k_2)*GF0%L(m,jj,k_3) 
                      !_________ 4th diagram
                   
              end do
           end do
        end if
        
     end do
  end do
  
  pp = delta/(2.d0*pi)
  Omega_R = Omr*pp*pp
end function Omega_r

subroutine int_SigLnG(i,j,sp,sp1,iw,SigL,SigG) !... interaction contributions of Eq. (23) + Eq. (26) 
  implicit none
  integer :: i, j, k1, k2, k3, iw, m, n, sp, sp1, s, s1, ii, jj
  complex*16 ::  SigL, SigG
  real*8 :: pp

  ii = i +sp; jj= j + sp1

  SigL=(0.0d0, 0.0d0) ! SigG=(0.0d0, 0.0d0)
  
  do k1 = 1, N_of_w
     do k2 = 1, N_of_w
        k3 = iw- k1 +k2

        if (k3 .ge. 1 .and. k3 .le. N_of_w) then   
           do s = 0, 1 !...Sum over spin
              do s1 = 0, 1
                 m = j+s1; n= i+s
                 
                 SigL = SigL + GF0%L(ii,jj,k1)*GF0%G(n,m,k2)*GF0%L(m,n,k3) &
                      - GF0%L(ii,n,k1)*GF0%G(n,m,k2)*GF0%L(m,jj,k3) 
                      
                
               !  SigG = SigG + GF0%G(m,n,k1)*GF0%L(n,m,k2)*GF0%G(ii,jj,k3) &
                      !- GF0%G(ii,m,k1)*GF0%L(m,n,k2)*GF0%G(n,jj,k3) 
                 
              end do
           end do
        end if
        
     end do
  end do
  pp = delta/(2.d0*pi)
  SigL = SigL*pp*pp !; SigG = SigG*pp*pp
end subroutine int_SigLnG


!====================================================
!========== restart routines  =======================
!====================================================

subroutine save_GFs()
  implicit none
  open(1,file='Greens_functions.dat',form='unformatted',status='unknown')
  write(1) GF0%r,GF0%a,GF0%l,GF0%g
  close(1)
end subroutine save_GFs

subroutine read_saved_GFs()
  implicit none
  open(1,file='Greens_functions.dat',form='unformatted',status='old')
  read(1) GF0%r,GF0%a,GF0%l,GF0%g
  close(1)
end subroutine read_saved_GFs

end module GreensFunctions
 
