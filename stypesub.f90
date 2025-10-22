subroutine input_SOC()
  use GreensFunctions 
  implicit none
  character(len=256) :: line, key
  integer :: ios, eq_pos
  
  open(22, file='inputSOC.dat', status='old')

  do
    read(22,'(A)', iostat=ios) line
    if (ios .ne. 0) exit   ! end of file
    if (trim(line) .eq. '' .or. line(1:1) .eq. '#') cycle  ! skip blank or comment

    eq_pos = index(line, "=")
    if (eq_pos .eq. 0) cycle  ! skip malformed lines

    key = adjustl(trim(line(:eq_pos-1)))

    select case (trim(key))
    case("T");            read(line(eq_pos+1:),*) T
    case("V");            read(line(eq_pos+1:),*) V
    case("Vf");           read(line(eq_pos+1:),*) Vf
    case("delv");         read(line(eq_pos+1:),*) delv

    case("N_ions");       read(line(eq_pos+1:),*) N_ions
    case("N_turns");      read(line(eq_pos+1:),*) N_turns
    case("hel_length");   read(line(eq_pos+1:),*) hel_length
    case("hel_radius");   read(line(eq_pos+1:),*) hel_radius
    case("hand");         read(line(eq_pos+1:),*) hand
       
    case("mu");           read(line(eq_pos+1:),*) mu
    case("order");        read(line(eq_pos+1:),*) order
    case("dw");           read(line(eq_pos+1:),*) dw
    case("up");           read(line(eq_pos+1:),*) up
    case("delta");        read(line(eq_pos+1:),*) delta

    case("pulay");        read(line(eq_pos+1:),*) pulay
    case("restart");      read(line(eq_pos+1:),*) restart
       
    case("E_CC");      read(line(eq_pos+1:),*) E_CC
    case("t_hop");        read(line(eq_pos+1:),*) t_hop
    case("lamb");         read(line(eq_pos+1:),*) lamb
    case("Hubbard");      read(line(eq_pos+1:),*) Hubbard
       
    case("Gamma");        read(line(eq_pos+1:),*) Gamma
    case("del_Gamma");    read(line(eq_pos+1:),*) del_Gamma

    end select
 end do

  close(22) 
  
  beta = 1.d0/(kb*T)
  Natoms = 2*N_ions*N_turns !...# of ions per turn of the helix  and the # of turns multiplied by 2 for spins, then add the leads on either end

  Volt_range = (Vf - V)/delv

  pauli_z = 0.d0; pauli_z(1,1) = 1.d0; pauli_z(2,2) = -1.d0
  
  write(*,*) 'T:', T, 'V:', V, 'mu:', mu, 'Volt_range:', Volt_range, 'Hubbard:', Hubbard
  write(*,*) 'Order:', order, 'Natoms:', Natoms
  write(*,*) 'delta:', delta
  write(*,*) 'pulay:', pulay

end subroutine input_SOC

subroutine PrintFunctions()
  use GreensFunctions
  implicit none
  integer :: i, j
  complex*16 :: diff

  open(12, file='info_Hamiltonian.dat', status='replace')
  
  write(12,'(/a)') '... Hamiltonian:'
  do i = 1, Natoms
     write(12, '(10(a,2f10.5,a))') (' [',H(i,j),'] ', j=1,Natoms)
  end do
  
  !.... checking if Hermitian  
  
  write(12,'(/a)') '... Checking the Hamiltonian is Hermitian:'
  do i=1,Natoms
     do j=i,Natoms
        diff=H(i,j)-conjg(H(j,i))
        write(12,*) i,j,diff
     end do
  end do
  
  write(12,'(/a)') '... Eigenvalues:'
  do i = 1, Natoms
     write(12, *) i,'.', Ev(i)
  end do
  
  write(12,'(/a)') '... Eigenvectors:'
  do j = 1, Natoms
     write(12, '(i3,10(a,2f10.5,a))') j,(' [',Eigenvec(i,j),'] ', i = 1, Natoms)
  end do

  write(12,*) 'Natoms:', Natoms

  close(12)
end subroutine PrintFunctions


subroutine trans(iw, Volt, trans_up, trans_down) !....square bracket terms of Eq. (2) in CHE
  use GreensFunctions
  implicit none
  integer :: iw,i,ii
  complex*16 :: trace1, trace2
  real*8 :: Volt, w, trans_up, trans_down

  w = omega(iw)
  
  work1 = GF0%L(:,:,iw)
  work2 = GF0%G(:,:,iw)
  
  work3 = im*matmul(GammaL, (fermi_dist(w, Volt)-1.d0)*work1 - fermi_dist(w, Volt)*work2)

  trace1 = (0.d0,0.d0);   trace2 = (0.d0,0.d0)
  do i = 1, Natoms, 2
     ii = i + 1 

     trace1 = trace1 + work3(i,i)
     trace2 = trace2 + work3(ii,ii)
     
  end do

  trans_up = real(trace1)/(2.d0*pi)
  trans_down = real(trace2)/(2.d0*pi)
end subroutine trans

subroutine Current(Volt, J_up, J_down)
  use GreensFunctions
  implicit none
  real*8 :: Volt, J_up, J_down, trans_up, trans_down
  integer :: iw
  
  J_up = 0.d0; J_down = 0.d0
  do iw = 1, N_of_w
     call trans(iw, Volt, trans_up, trans_down)
     J_up = J_up + trans_up
     J_down = J_down + trans_down
  end do
  
  J_up = J_up*(delta/hbar)
  J_down = J_down*(delta/hbar)
  
end subroutine Current
