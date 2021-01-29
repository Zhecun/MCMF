module Hmat
    use wavefunction
    use Wmat
    use constant
    use Model
    use Wmat
    implicit none
contains
    subroutine H_0000(Phi,H0000)
        implicit none
        integer :: j,l,s,sprime
        type(psi) :: Phi
        double precision :: H0000(4,4)
        complex(kind=8) :: z(2)

        z(:)=sqrt(gamma/2)*Phi%R(:)+(0,1)*sqrt(1/(2*gamma))*Phi%P(:)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        H0000((s-1)*2+j,(sprime-1)*2+l)=(-1/(2*M))*&                                                              ! This part belongs to the kinetic operator
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&       ! kinetic
                        (-gamma-Phi%P(l)+&                                                                                        ! kinetic
                        gamma**2*(1/(2*gamma)-conjg(z(j))*(0,1)*sqrt(1/(2*gamma))+Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2))+&         ! kinetic
                        2*gamma*Phi%P(l)*(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-conjg(z(j))*(0,1)*sqrt(1/(2*gamma)))&                ! kinetic
                        +0.5*OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&  ! This part belongs to the electronic Hamiltonian
                        (E(s,Phi%R(j))+E(sprime,Phi%R(l))+&                                                                       ! electronic
                        D_E(sprime,Phi%R(l))*(0,1)*(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-conjg(z(j))*(0,1)*sqrt(1/(2*gamma)))+&    ! electronic
                        D_E(sprime,Phi%R(j))*&                                                                                    ! electronic
                        conjg((0,1)*(Phi%P(j)/(2*gamma)+(0,1)*Phi%R(j)/2-conjg(z(l))*(0,1)*sqrt(1/(2*gamma)))))                   ! electronic
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine H_R000(Phi,HR000)
        implicit none
        integer :: j,l,s,sprime
        type(psi) :: Phi
        double precision :: HR000(4,4)
        complex(kind=8) :: z(2)

        z(:)=sqrt(gamma/2)*Phi%R(:)+(0,1)*sqrt(1/(2*gamma))*Phi%P(:)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        HR000((s-1)*2+j,(sprime-1)*2+l)=(-1/(2*M))*&                                                              ! This part belongs to the kinetic operator
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&       ! kinetic
                        (-gamma**2*(0,1)/2-(0,1)*gamma*Phi%P(l)+(-gamma-(Phi%P(l))**2+gamma**2*(1/(2*gamma)-conjg(z(j))*(0,1)*&   ! kinetic
                        sqrt(1/(2*gamma))+Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2)+2*gamma*Phi%P(l)*(Phi%P(l)/(2*gamma)+&             ! kinetic
                        (0,1)*Phi%R(l)/2-conjg(z(j))*(0,1)*sqrt(1/(2*gamma))))*(z(l)*sqrt(gamma/2)-gamma/2*Phi%R(j)))+&            ! kinetic
                        0.5*OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&   ! This part belongs to the electronic Hamiltonian
                        ((z(l)*sqrt(gamma/2)-gamma/2*Phi%R(j))*(E(s,Phi%R(j))+E(sprime,Phi%R(l)))+D_E(sprime,Phi%R(l))*(0,1)*&    ! electronic
                        (-0.5*(0,1)+(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-conjg(z(j))*(0,1)*sqrt(1/(2*gamma)))*&                   ! electronic
                        (z(l)*sqrt(gamma/2)-gamma/2*Phi%R(j)))+D_E(s,Phi%R(j))*(0.5-conjg(z(j))*(0,1)*sqrt(gamma/2)+Phi%P(l)/2+&  ! electronic
                        (0,1)*gamma*Phi%R(l)/2)-&                                                                                 ! electronic
                        Phi%P(j)/2*(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-conjg(z(j))*(0,1)*sqrt(1/(2*gamma))))                     ! electronic
                    enddo
                enddo
            enddo
        enddo

    endsubroutine

    subroutine H_P000(Phi,HP000)
        implicit none
        integer :: j,l,s,sprime
        type(psi) :: Phi
        double precision :: HP000(4,4)
        complex(kind=8) :: z(2)

        z(:)=sqrt(gamma/2)*Phi%R(:)+(0,1)*sqrt(1/(2*gamma))*Phi%P(:)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        HP000((s-1)*2+j,(sprime-1)*2+l)=(-1/(2*M))*&                                                              ! This part belongs to the kinetic operator
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&       ! kinetic
                        (-gamma/2-Phi%P(l)+(-gamma-(Phi%P(l))**2+gamma**2*(1/(2*gamma)-conjg(z(j))*(0,1)*sqrt(1/(2*gamma))+&      ! kinetic
                        Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2)+2*gamma*Phi%P(l)*(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-conjg(z(j))*&  ! kinetic
                        (0,1)*sqrt(1/(2*gamma))))*(-z(l)*(0,1)*sqrt(1/(2*gamma))-Phi%P(j)/(2*gamma)))+&                           ! kinetic
                        0.5*OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&
                        ((E(s,Phi%R(j))+E(sprime,Phi%R(l)))*(-z(l)*(0,1)*sqrt(1/(2*gamma)))+&
                        D_E(s,Phi%R(j))*(-(0,1)/(2*gamma)+(0,1)*(Phi%P(j)/(2*gamma)-(0,1)*Phi%R(j)/2+z(l)*(0,1)*sqrt(1/(2*gamma)))&
                        *(z(l)*(0,1)*sqrt(1/(2*gamma))+Phi%P(j)/(2*gamma)))-&
                        D_E(sprime,Phi%R(l))*((0,1)/(2*gamma)+(0,1)*(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-&
                        conjg(z(j))*(0,1)*sqrt(1/(2*gamma)))*(z(l)*(0,1)*sqrt(1/(2*gamma))+Phi%P(j)/(2*gamma))))                    
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine H_00R0(Phi,H00R0)
        implicit none
        integer :: j,l,s,sprime,spp
        type(psi) :: Phi
        double precision :: H00R0(4,4),H0000(4,4)
        complex(kind=8) :: z(2)

        z(:)=sqrt(gamma/2)*Phi%R(:)+(0,1)*sqrt(1/(2*gamma))*Phi%P(:)
        call H_0000(Phi,H0000)

        do j=1,2
            call NAC(d,Phi%R(j))
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        H00R0((s-1)*2+j,(sprime-1)*2+l)=OverlapNAC_dagger(Phi%R(j),Phi%R(l),s,sprime)*&
                        (-1/(2*M))*&                                                                                              ! This part belongs to the kinetic operator
                        OverlapElectronic(Phi%R(j),Phi%R(l),s,sprime)*OverlapNuclear(Phi%R(j),Phi%P(j),Phi%R(l),Phi%P(l))*&       ! kinetic
                        (-gamma-Phi%P(l)+&                                                                                        ! kinetic
                        gamma**2*(1/(2*gamma)-conjg(z(j))*(0,1)*sqrt(1/(2*gamma))+Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2))+&         ! kinetic
                        2*gamma*Phi%P(l)*(Phi%P(l)/(2*gamma)+(0,1)*Phi%R(l)/2-conjg(z(j))*(0,1)*sqrt(1/(2*gamma)))                ! kinetic
                        do spp=1,2
                            H00R0=H00R0-d(s,spp)*H0000((spp-1)*2+j,(sprime-1)*2+l)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    
end module Hmat