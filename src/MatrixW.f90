!---------------------------------------------------------------------------------------------------------------------
!                                                                                                                    !
! This module contains programs for the calculation of overlap matrix W and its derivative (some Gaussian wavepacket !
! are replaced by their derivative of R or P ). Every W is a (M * N) * (M * N) matrix, where M is the number of      !
! Gaussian configurations and N is the number of electronic states.                                                  !
!                                                                                                                    !
!---------------------------------------------------------------------------------------------------------------------
module Wmat
    use wavefunction
    use Model
    implicit none
    
    complex(kind=8) :: W0000(4,4)
    complex(kind=8) :: WR000(4,4)
    complex(kind=8) :: W00R0(4,4)
    complex(kind=8) :: W00RR(4,4)
    complex(kind=8) :: WRR00(4,4)
    complex(kind=8) :: WR00R(4,4)
    complex(kind=8) :: W0RR0(4,4)
    complex(kind=8) :: W0R00(4,4)
    complex(kind=8) :: W000R(4,4)
    complex(kind=8) :: WRP00(4,4)
    complex(kind=8) :: W0PR0(4,4)
    complex(kind=8) :: W0P00(4,4)
    complex(kind=8) :: WP000(4,4)
    complex(kind=8) :: WP00R(4,4)
    complex(kind=8) :: WPR00(4,4)
    complex(kind=8) :: WPP00(4,4)
contains

    complex(kind=8) function OverlapNuclear(Rj,Pj,Rl,Pl)
    !-----------------------------------------------------!
    ! This function calculate the Overlap integral between!
    ! two Gaussian wavepackets, < \chi_j | \chi_l > ,     !
    ! j is bra, l is ket.                                 !
    !-----------------------------------------------------!
        implicit none

        double precision,intent(in) :: Rj,Pj,Rl,Pl
        
        complex(kind=8)             :: zj,zl

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
        
        integer,         intent(in) :: s,sprime
        double precision,intent(in) :: Rj,Rl
        
        double precision            :: cj(2,2),cl(2,2)

        call Diagonal(cj,Rj)
        call Diagonal(cl,Rl)

        OverlapElectronic=cj(s,1)*cl(sprime,1)+cj(s,2)*cl(sprime,2)    ! since in the model system the electronic wavefunction 
                                                                       ! is always real, so there is no need to use conjg() 

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
        
        integer,         intent(in) :: s,sprime
        double precision,intent(in) :: Rj,Rl 
        
        integer                     :: spp

        OverlapNAC=0
        call NAC(d,Rl)               ! maybe there will be some problem with d?

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
        
        integer,         intent(in) :: s,sprime
        double precision,intent(in) :: Rj,Rl 
        
        integer                     :: spp

        OverlapNAC_dagger=0
        call NAC(d,Rj)

        do spp=1,2
            OverlapNAC_dagger=OverlapNAC_dagger-d(s,spp)*OverlapElectronic(Rj,Rl,spp,sprime)
        enddo

        return
    endfunction

    double precision function OverlapddPhi(Rj,Rl,s,sprime)
    !-----------------------------------------------------------
    ! This function calculate 
    ! < \pdv{\phi_j_s}{R_j} | \pdv{\phi_l_sprime}{R_l} >
    !-----------------------------------------------------------
        implicit none
        
        integer,         intent(in) :: s,sprime
        double precision,intent(in) :: Rj,Rl 
        
        integer                     :: spp,sppp
        double precision            :: d_aux(2,2)

        call NAC(d,Rj)
        call NAC(d_aux,Rl)

        do spp=1,2
            do sppp=1,2
                OverlapddPhi=OverlapddPhi-d_aux(s,spp)*OverlapElectronic(Rj,Rl,spp,sppp)*d(sppp,sprime)
            enddo
        enddo

        return
    endfunction
    
    subroutine W_0000(Phi,W0000)
        implicit none
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W0000(4,4)    ! the general dimension formula should be M*N, for a 2-configuration and 
                                                        ! 2-electronic-state wavefunction, the number is 2*2=4

        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WR000(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W00R0(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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

    subroutine W_00RR(Phi,W00RR)
        implicit none
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W00RR(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        W00RR(j+(s-1)*2,l+(sprime-1)*2)=&
                        OverlapddPhi(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine W_RR00(Phi,WRR00)
        implicit none
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WRR00(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

        ! Be careful of the derivative when they are to the same index!
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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WR00R(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W0RR0(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W0R00(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W000R(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WRP00(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W0PR0(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: W0P00(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WP000(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi)                     :: Phi
        complex(kind=8),intent(inout) :: WP00R(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WPR00(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

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
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: WPP00(4,4)
        
        integer                       :: j,l
        integer                       :: s,sprime

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        if (j/=l) then
                            WPP00(j+(s-1)*2,l+(sprime-1)*2)=&
                            OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                            (1/(2*gamma)+&
                            (-(sqrt(gamma/2)*Phi%R(l)+(0,1)*sqrt(1/(2*gamma))*Phi%P(l))*(0,1)*sqrt(1/(2*gamma))&
                            -(1/(2*gamma))*Phi%P(j))&
                            *((sqrt(gamma/2)*Phi%R(j)-(0,1)*sqrt(1/(2*gamma))*Phi%P(j))*(0,1)*sqrt(1/(2*gamma))&
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