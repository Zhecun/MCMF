module wavefunction
    implicit none

    type :: psi
        double precision :: r(2),p(2)       ! mean coordinate and mean momentum for each Gaussian wavepacket
        complex(kind=8)  :: A(4)            ! coefficient in front of electronic wavefunction
    endtype
contains
    
    subroutine DensityMatrix(Phi,Rho)
        implicit none
        integer :: j,l,s,sprime
        type(psi) :: Phi 
        complex(kind=8) :: Rho(2,2,2,2)
        ! The first and second indexes stand for configuration, the third and fourth indexes stand for electronic states.

        do s=1,2
            do sprime=1,2
                do j=1,2
                    do l=1,2
                        Rho(j,l,s,sprime)=conjg(Phi%A((s-1)*2+j))*Phi%A((sprime-1)*2+l)
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine
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
