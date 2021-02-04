!-----------------------------------------------------------------------------------------
!
! To propagate R,P and to calculte the matrix elements of involving matrices C, M, 
! and Y, it will be better to use a high-rank tensor style of W. This module provides 
! some definition and subroutines to transform W from matrix to tensor. Then the subroutines
! which generate C, M, Y are provided. 
!
!-------------------------------------------------------------------------------

module MatCMY
    use wavefunction
    use Wmat
    use Hmat
    use LAmethod
    implicit none

    complex(kind=8) :: YR(4)
    complex(kind=8) :: YP(4)
    complex(kind=8) :: CRR(4,4)
    complex(kind=8) :: CPR(4,4)
    complex(kind=8) :: MPR(4,4)
    complex(kind=8) :: MPP(4,4)
contains
    subroutine Transform(W,W_tensor)
        implicit none
        integer          :: j,l,s,sprime
        double precision :: W(4,4)
        double precision :: W_tensor(2,2,2,2)

        do s=1,2
            do sprime=1,2
                do j=1,2
                    do l=1,2
                        W_tensor(j,l,s,sprime)=W(j+(s-1)*2,l+(sprime-1)*2)
                    enddo
                enddo
            enddo
        enddo

        return 
    endsubroutine

    subroutine Y_R(Phi,YR)
        implicit none
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: YR(2)
        
        
        integer                       :: s,sprime,j,l
        complex(kind=8)               :: SS(4,4)
        complex(kind=8)               :: W0000inv(4,4)
        complex(kind=8)               :: Rho(4,4)

        YR(:)=0

        call DensityMatrix(Phi,Rho)
        call MKL_inv(4,W0000,W0000inv)
        SS=matmul(matmul(WR000+W00R0,W0000inv),H0000)    

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        YR(j)=YR(j)+Rho(2*(s-1)+j,2*(sprime-1)+l)*(HR000(2*(s-1)+j,2*(sprime-1)+l)+&
                        H00R0(2*(s-1)+j,2*(sprime-1)+l)-SS(2*(s-1)+j,2*(sprime-1)+l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine Y_P(Phi,YP)
        implicit none
        
        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: YP(2)
        
        integer                       :: s,sprime,j,l
        complex(kind=8)               :: SS(4,4)
        complex(kind=8)               :: W0000inv(4,4)
        complex(kind=8)               :: Rho(4,4)

        YP(:)=0

        call DensityMatrix(Phi,Rho)
        call MKL_inv(4,W0000,W0000inv)
        SS=matmul(matmul(WP000,W0000inv),H0000)    ! -- is the inverse of W0000 

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        YP(j)=YP(j)+Rho(2*(s-1)+j,2*(sprime-1)+l)*(HP000(2*(s-1)+j,2*(sprime-1)+l)-&
                        SS(2*(s-1)+j,2*(sprime-1)+l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine C_RR(Phi,CRR)
        implicit none

        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: CRR(4,4)

        integer                       :: s,sprime,j,l
        complex(kind=8)               :: SS0R00(4,4)
        complex(kind=8)               :: SS000R(4,4)
        complex(kind=8)               :: W0000inv(4,4)
        complex(kind=8)               :: Rho(4,4)

        CRR=0

        call DensityMatrix(Phi,Rho)
        call MKL_inv(4,W0000,W0000inv)
        SS0R00=matmul(matmul(WR000+W00R0,W0000inv),W0R00)
        SS000R=matmul(matmul(WR000+W00R0,W0000inv),W000R)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        CRR(j,l)=CRR(j,l)+Rho(2*(s-1)+j,2*(sprime-1)+l)*(WRR00(2*(s-1)+j,2*(sprime-1)+l)+&
                        WR00R(2*(s-1)+j,2*(sprime-1)+l)+W0RR0(2*(s-1)+j,2*(sprime-1)+l)+&
                        W00R0(2*(s-1)+j,2*(sprime-1)+l)-SS0R00(2*(s-1)+j,2*(sprime-1)+l)-SS000R(2*(s-1)+j,2*(sprime-1)+l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine C_PR(Phi,CPR)
        implicit none

        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: CPR(4,4)

        integer                       :: s,sprime,j,l
        complex(kind=8)               :: SS0R00(4,4)
        complex(kind=8)               :: SS000R(4,4)
        complex(kind=8)               :: W0000inv(4,4)
        complex(kind=8)               :: Rho(4,4)

        CPR=0

        call DensityMatrix(Phi,Rho)
        call MKL_inv(4,W0000,W0000inv)
        SS0R00=matmul(matmul(WP000,W0000inv),W0R00)
        SS000R=matmul(matmul(WP000,W0000inv),W000R)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        CPR(j,l)=CPR(j,l)+Rho(2*(s-1)+j,2*(sprime-1)+l)*(WP00R(2*(s-1)+j,2*(sprime-1)+l)+&
                        WPR00(2*(s-1)+j,2*(sprime-1)+l)&
                        -SS0R00(2*(s-1)+j,2*(sprime-1)+l)-SS000R(2*(s-1)+j,2*(sprime-1)+l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine M_RP(Phi,MRP)
        implicit none

        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: MRP(4,4)

        integer                       :: s,sprime,j,l
        complex(kind=8)               :: SS(4,4)
        complex(kind=8)               :: W0000inv(4,4)
        complex(kind=8)               :: Rho(4,4)

        MRP=0

        call DensityMatrix(Phi,Rho)
        call MKL_inv(4,W0000,W0000inv)
        SS=matmul(matmul(WR000+W00R0,W0000inv),W0P00)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        MRP(j,l)=MRP(j,l)+Rho(2*(s-1)+j,2*(sprime-1)+l)*(WRP00(2*(s-1)+j,2*(sprime-1)+l)+&
                        W0PR0(2*(s-1)+j,2*(sprime-1)+l)-SS(2*(s-1)+j,2*(sprime-1)+l))
                    enddo
                enddo
            enddo
        enddo

        return
    endsubroutine

    subroutine M_PP(Phi,MPP)
        implicit none

        type(psi),      intent(in)    :: Phi
        complex(kind=8),intent(inout) :: MPP(4,4)

        integer                       :: s,sprime,j,l
        complex(kind=8)               :: SS(4,4)
        complex(kind=8)               :: W0000inv(4,4)
        complex(kind=8)               :: Rho(4,4)

        MPP=0

        call DensityMatrix(Phi,Rho)
        call MKL_inv(4,W0000,W0000inv)
        SS=matmul(matmul(WP000,W0000inv),W0P00)

        do j=1,2
            do l=1,2
                do s=1,2
                    do sprime=1,2
                        MPP(j,l)=MPP(j,l)+Rho(2*(s-1)+j,2*(sprime-1)+l)*(WPP00(2*(s-1)+j,2*(sprime-1)+l)&
                        -SS(2*(s-1)+j,2*(sprime-1)+l))
                    enddo
                enddo
            enddo
        enddo

        return

    endsubroutine

end module MatCMY