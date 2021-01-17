module wavefunction
    implicit none

    type :: psi
        double precision :: r(2),p(2)       ! mean coordinate and mean momentum for each Gaussian wavepacket
        complex(kind=8)  :: A(4)            ! coefficient in front of electronic wavefunction
    endtype
contains
    
end module wavefunction

!---------------------------------------------------------------------------------------------------------
! for a wavefunction with index j(configuration), s(electronic state)  M configuration, N electronic state
! the save sequence of A_j^(s) in an array A(M*N) is defined to be:
! 
! do s=1,N
!     do j=1,M
!         A((s-1)*N+j)=A_j^(s)
!     enddo
! enddo
!
!----------------------------------------------------------------------------------------------------------
