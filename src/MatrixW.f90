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

    double precision function OverlapNAC(Rj,Rl,s,sprime)
    !-----------------------------------------------------------------------
    ! This function calculate the "Overlap" integral between an adiabatic
    ! wavefunction and the derivative of another adiabatic wavefunction
    ! centers at Rj and Rl, < \phi_j_s | \pdv{\phi_l_sprime}{R_l} >, j is bra,
    ! l is ket
    !-----------------------------------------------------------------------
        implicit none
        integer :: s,sprime
        double precision :: Rj,Rl 
        integer :: spp

        OverlapNAC=0
        call NAC(d,Rl)

        do spp=1,2
            OverlapNAC=OverlapNAC+OverlapElectronic(Rj,Rl,s,spp)*d(spp,sprime)
        enddo

        return
    endfunction

    double precision function OverlapNAC_dagger(Rj,Rl,s,sprime)
    !-----------------------------------------------------------
    ! just like OverlapNAC above, this function calculate 
    ! < \pdv{\phi_j_s}{R_s} | \phi_l_sprime >
    !-----------------------------------------------------------
        implicit none
        integer :: s,sprime
        double precision :: Rj,Rl 
        integer :: spp

        OverlapNAC_dagger=0
        call NAC(d,Rj)

        do spp=1,2
            OverlapNAC_dagger=OverlapNAC_dagger-d(s,spp)*OverlapElectronic(Rj,Rl,spp,s)
        enddo

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
                        OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_R000(Phi,WR000)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WR000(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        WR000(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*sqrt(gamma/2)-(gamma/2)*Phi%R(j))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_00R0(Phi,W00R0)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W00R0(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W00R0(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNAC_dagger(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_RR00(Phi,WRR00)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WRR00(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        if (j/=l) then
                            WRR00(j+(s-1)*2,l+(sprime-1)*2)=&
                            OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                            (gamma/2&
                            +((sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*sqrt(gamma/2)-(gamma/2)*Phi%R(j))*&
                            ((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*sqrt(gamma/2)-(gamma/2)*Phi%R(l)))
                        elseif (j==l) then
                            if (s/=sprime) then
                                WRR00(j+(s-1)*2,l+(sprime-1)*2)=0
                            elseif (s==sprime) then
                                WRR00(j+(s-1)*2,l+(sprime-1)*2)=(Phi%P(j))**2/4+gamma/2
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_R00R(Phi,WR00R)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WR00R(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        WR00R(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNAC(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*sqrt(gamma/2)-(gamma/2)*Phi%R(j))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_0RR0(Phi,W0RR0)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W0RR0(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W0RR0(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNAC_dagger(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*sqrt(gamma/2)-(gamma/2)*Phi%R(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_0R00(Phi,W0R00)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W0R00(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W0R00(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*sqrt(gamma/2)-(gamma/2)*Phi%R(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_000R(Phi,W000R)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W000R(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W000R(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNAC(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_RP00(Phi,WRP00)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WRP00(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        if (j/=l) then
                            WRP00(j+(s-1)*2,l+(sprime-1)*2)=&
                            OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                            (((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*(0,1)*sqrt(1/(2*gamma))-&
                            (1/(2*gamma))*Phi%P(l))*&
                            ((sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*sqrt(gamma/2)-(gamma/2)*Phi%R(j))&
                            +(0,1)*0.5)
                        elseif (j==l) then
                            WRP00(j+(s-1)*2,l+(sprime-1)*2)=0
                        endif
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_0PR0(Phi,W0PR0)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W0PR0(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W0PR0(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNAC_dagger(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*(0,1)*sqrt(1/(2*gamma))-&
                        (1/(2*gamma))*Phi%P(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_0P00(Phi,W0P00)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: W0P00(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W0P00(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*(0,1)*sqrt(1/(2*gamma))-&
                        (1/(2*gamma))*Phi%P(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_P000(Phi,WP000)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WP000(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        WP000(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        (-(sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*(0,1)*sqrt(1/(2*gamma))-&
                        (1/(2*gamma))*Phi%P(j))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_P00R(Phi,WP00R)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WP00R(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        WP00R(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapNAC(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        (-(sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*(0,1)*sqrt(1/(2*gamma))-&
                        (1/(2*gamma))*Phi%P(j))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_PR00(Phi,WPR00)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WPR00(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        if (j/=l) then
                            WPR00(j+(s-1)*2,l+(sprime-1)*2)=&
                            OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                            ((-(sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*(0,1)*sqrt(1/(2*gamma))-&
                            (1/(2*gamma))*Phi%P(j))*&
                            ((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*sqrt(gamma/2)-(gamma/2)*Phi%R(l))&
                            -(0,1)*0.5)
                        elseif (j==l) then
                            WPR00(j+(s-1)*2,l+(sprime-1)*2)=0
                        endif
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_PP00(Phi,WPP00)
        implicit none
        integer :: j,l
        integer :: s,sprime
        type(psi) :: Phi
        double precision :: WPP00(4,4)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        if (j/=l) then
                            WPP00(j+(s-1)*2,l+(sprime-1)*2)=&
                            OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                            (1/(2*gamma)+
                            (-(sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*(0,1)sqrt(1/(2*gamma))&
                            -(1/(2*gamma))*Phi%P(j))&
                            *((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*(0,1)sqrt(1/(2*gamma))&
                            -(1/(2*gamma))*Phi%P(l)))
                        elseif (j==l) then
                            if (s/=sprime) then
                                WPP00(j+(s-1)*2,l+(sprime-1)*2)=0
                            elseif (s==sprime) then
                                WPP00(j+(s-1)*2,l+(sprime-1)*2)=1/(2*gamma)+(Phi%R(j))**2/4
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine
    
end module Wmat