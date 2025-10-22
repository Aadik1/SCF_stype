module DefineHamiltonian
  implicit none
  integer :: Natoms, N_of_w, Volt_range, N_turns, N_ions, order

  real*8, dimension(:,:) :: C(3,3), pauli_z(2,2)
  real*8, dimension(3) :: Rij, w0
  real*8 :: T, mu, beta, V, delta, w_init, w_fin, delv, Vf
  real*8 :: E_CC, t_hop
  real*8 :: hel_radius, hel_length, lamb, hand, Hubbard
  real*8 :: Gamma, del_Gamma
  
  real*8, parameter :: hbar = 1.d0 !6.582119569e-16 !ev s
  real*8, parameter :: kb = 8.6173303e-5 !ev/K
  real*8, parameter :: pi = 4.d0*atan(1.d0)
  real*8, parameter :: el = 1.d0 !.1.602176634e-19 !C
  
  complex*16, parameter :: im = (0.d0, 1.d0)
  complex*16, allocatable, dimension(:,:) :: H

  logical :: restart
contains 
  
  subroutine helix_coordinate(r, i) !....hand chooses helicity- RH helix = 1, LH helix = -1
    implicit none
    integer :: i
    real*8, dimension(3) :: r
    r = 0.d0
    r(1) = hel_radius*cos(hand*(i-1)*2*pi/N_ions)
    r(2) = hel_radius*sin(hand*(i-1)*2*pi/N_ions)
    r(3) = hel_length*(i-1)/(N_turns*(N_ions - 1))
  end subroutine helix_coordinate
  
  subroutine SOC_Hamiltonian() !...Eq(1)
    implicit none    
    integer :: i,j,ii,jj,Nat
    real*8 :: r1(3),r2(3),r(3),d1(3),d2(3),dm1,dm2,v(3),mod_vec
    complex*16 :: im=(0.0d0,1.0d0)

    H=(0.0d0,0.0d0)

    Nat=Natoms/2 
    do i= 1, Nat
       ii=2*i
!..................diagonal        
       H(ii-1,ii-1)= E_CC
       H(ii,ii)    = E_CC
       
!.................. hopping       
       if(i.le.Nat-1) then
          j=i+1 ; jj=2*j
          H(ii-1,jj-1) = t_hop
          H(ii,jj)     = t_hop
       end if
       
!...................spin-orbit          
      if(i.le.Nat-2) then
          j=i+2 ; jj=2*j 
          call helix_coordinate(r,i)
          call helix_coordinate(r1,i+1)
          call helix_coordinate(r2,i+2)
          dm1=mod_vec(r-r1) ; d1=(r-r1)/dm1
          dm2=mod_vec(r-r2) ; d2=(r-r2)/dm2
          call vec_cross_product(d1,d2,v) 
          H(ii-1,jj-1) = im*lamb * v(3)
          H(ii-1,jj)   = im*lamb * (v(1)-im*v(2))
          H(ii,jj-1)   = im*lamb * (v(1)+im*v(2))
          H(ii,jj)     = im*lamb * (-v(3))
          
      end if
    end do
    
    !........ make Hermitian
    
    do i= 1,Natoms
       do j=i,Natoms
          H(j,i)=conjg(H(i,j))
       end do
    end do
    
  end subroutine SOC_Hamiltonian
 
end module DefineHamiltonian
