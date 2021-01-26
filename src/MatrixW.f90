!------------------------------------------------------------------------------------------------------------------
!
! This module contains programs for the calculation of overlap matrix W and its derivative (some Gaussian wavepacket 
! are replaced by their derivative of R or P ). Every W is a (M * N) * (M * N) matrix, where M is the number of 
! Gaussian configurations and N is the number of electronic states.
!
!--------------------------------------------------------------------------------------------------------------------
module Wmat
    implicit none
    
contains

    complex(kind=8) function OverlapNuclear(Rj,Pj,Rl,Pl)
    !-----------------------------------------------------
    ! This function calculate the Overlap integral between
    ! two Gaussian wavepackets, < \chi_j | \chi_l > ,
    ! j is bra, l is ket.
    !-----------------------------------------------------
        implicit none
        double precision :: Rj,Pj,Rl,Pl
        complex(kind=8) :: zj,zl

        zj=sqrt(gamma/2)*Rj+(0,1)*sqrt(1/(2*gamma))*Pj
        zl=sqrt(gamma/2)*Rl+(0,1)*sqrt(1/(2*gamma))*Pl

        OverlapNuclear=exp(conjg(zj)*zl-conjg(zj)*zj/2-conjg(zl)*zl/2)

        return
    endfunction

    double precision function OverlapElectronic(Rj,Rl,s,sprime)
    !---------------------------------------------------------------------
    ! This function calculate the Overlap integral between
    ! two adiabatic electronic wavefunction located at different
    ! centers Rj and Rl, < \phi_j_s | \phi_l_sprime >, j is bra, l is ket.
    !---------------------------------------------------------------------
        implicit none
        integer :: s,sprime
        double precision :: Rj,Rl
        double precision :: cj(2,2),cl(2,2)

        call Diagonal(cj,Rj)
        call Diagonal(cl,Rl)

        OverlapElectronic=conjg(cj(s,1))*cl(sprime,1)+conjg(cj(s,2))*cl(sprime,2)

        return 
    endfunction


    subroutine W_0000(Phi,W0000)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W0000(4,4)      ! the general dimension formula should be M*N, for a 2-configuration and 
                                            ! 2-electronic-state wavefunction, the number is 2*2=4
        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W0000(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNuclear(R(j),P(j),R(l),P(l))*OverlapElectronic(R(j),R(l),s,sprime)
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_R000(r,p,A)
        implicit none
    endsubroutine

    subroutine W_00R0(r,p,A)
        implicit none
    endsubroutine

    subroutine W_RR00(r,p,A)
        implicit none
    endsubroutine

    subroutine W_R00R(r,p,A)
        implicit none
    endsubroutine

    subroutine W_0RR0(r,p,A)
        implicit none
    endsubroutine

    subroutine W_0R00(r,p,A)
        implicit none
    endsubroutine

    subroutine W_000R(r,p,A)
        implicit none
    endsubroutine

    subroutine W_RP00(r,p,A)
        implicit none
    endsubroutine

    subroutine W_0PR0(r,p,A)
        implicit none
    endsubroutine

    subroutine W_0P00(r,p,A)
        implicit none
    endsubroutine

    subroutine W_P000(r,p,A)
        implicit none
    endsubroutine

    subroutine W_P00R(r,p,A)
        implicit none
    endsubroutine

    subroutine W_PR00(r,p,A)
        implicit none
    endsubroutine

    subroutine W_PP00(r,p,A)
        implicit none
    endsubroutine
    
end module Wmat