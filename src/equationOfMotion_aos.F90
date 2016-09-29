!> This module calculates the yprime terms for the RK interpolation routine for the particle push
!> mini-app. Mostly copied and cleaned up from XGC1/kernels/pushe/vec/derivs_elec_vec.F90
!> @author Tuomas Koskela
!> @date 2 Sep 2016
!<
module eom

  use params
  
  implicit none

contains

  !> Evaluates the LHS of the equations of motion. See equation (1) in Chang et al, PoP, 2009
  !<
  function eom_eval(y,Bvec,jacB,Evec,dt,yprime,mu,charge,mass) result(err)

    double precision,    intent(in) :: dt  !> Time step length
    double precision,    intent(in),  dimension(veclen) :: charge  !> Charge of particle
    double precision,    intent(in),  dimension(veclen) :: mass    !> Mass of particle
    double precision,    intent(in),  dimension(veclen) :: mu !> Magnetic moment of particle
    double precision,    intent(in),  dimension(veclen,4) :: y !> Particle phase information
    double precision,    intent(in),  dimension(veclen,3) :: Evec !> Electric field vector (Er,Ephi,Ez)
    double precision,    intent(in),  dimension(veclen,3) :: Bvec !> Magnetic field vector (Br,Bphi,Bz)
    !> Magnetic field jacobian
    !> |  dBRdR,   dBRdphi,   dBRdz  |
    !> | dBphidR, dBphidphi, dBphidz |
    !> |  dBzdR,   dBzdphi,   dBzdz  |
    !<
    double precision,    intent(in),  dimension(veclen,3,3)  :: jacB

    !> LHS of 1st (1:3) and 2nd (4) equations of motion
    !<
    double precision,    intent(out), dimension(veclen,4)     :: yprime 

    integer :: err !> error flag
    integer :: iv  !> loop index
    double precision :: over_r  !> 1/R, stored to avoid divides
    double precision :: B       !> |B|
    double precision :: B2      !> |B|^2
    double precision :: over_B  !> |B|^{-1}
    double precision :: over_B2 !> |B|^{-2}
    double precision :: c_m     !> charge * mass
    double precision :: cmrho   !> charge * mass * rho_parallel
    double precision :: cmrho2  !> charge * mass * rho_parallel^2
    double precision :: D       !> LHS of 3rd eqn of motion
    double precision :: Fr      !> Force term, stored for brevity, R-component
    double precision :: Fp      !> Force term, stored for brevity, phi-component
    double precision :: Fz      !> Force term, stored for brevity, z-component

    err = 0

    do iv = 1,veclen

       ! Norms
       B2 = dot_product(bvec(iv,:),bvec(iv,:))
#ifdef EXPERIMENT
       B  = B2
       
       c_m = charge(iv) * mass(iv)
       cmrho = c_m * y(iv,4)
       cmrho2 = cmrho * y(iv,4)
       over_r = 1.0D0 * y(iv,1)
       over_B  = 1.0D0 * B
       over_B2 = 1.0D0 * B2       

       D = 1.0D0 * (1.0D0 + y(iv,4) * (&
            -1.D0 * over_B2 * ( jacb(iv,2,3) * bvec(iv,1) &
            + (jacb(iv,3,1) - jacb(iv,1,3)) * bvec(iv,2) &
            - (bvec(iv,2) * over_r + jacb(iv,2,1)) * bvec(iv,3) )))
#else
       B  = sqrt(B2)
       
       c_m = charge(iv) * mass(iv)
       cmrho = c_m * y(iv,4)
       cmrho2 = cmrho * y(iv,4)
       over_r = 1.0D0 / y(iv,1)
       over_B  = 1.0D0 / B
       over_B2 = 1.0D0 / B2       

       D = 1.0D0 / (1.0D0 + y(iv,4) * (&
            -1.D0 * over_B2 * ( jacb(iv,2,3) * bvec(iv,1) &
            + (jacb(iv,3,1) - jacb(iv,1,3)) * bvec(iv,2) &
            - (bvec(iv,2) * over_r + jacb(iv,2,1)) * bvec(iv,3) )))
#endif
       Fr = evec(iv,1) - &
            ( (mu(iv) + charge(iv) * cmrho2 * B) & ! murho2b
            * (bvec(iv,1) * jacb(iv,1,1) + bvec(iv,2) * jacb(iv,2,1) &
            + bvec(iv,3) * jacb(iv,3,1)) * over_B & ! dbdr
            ) 
       Fp = evec(iv,2) ! dbdphi is 0 by definition
       Fz = evec(iv,3) - &
            ( (mu(iv) + charge(iv) * cmrho2 * B) & ! murho2b
            * (bvec(iv,1) * jacb(iv,1,3) + bvec(iv,2) * jacb(iv,2,3) &
            + bvec(iv,3) * jacb(iv,3,3)) * over_B & ! dbdz
            )
       
       yprime(iv,1) = D*( (bvec(iv,3) * Fp - bvec(iv,2) * Fz) * over_b2         &
            + cmrho * bvec(iv,1) &
            + cmrho2 * (jacb(iv,3,2) * over_r - jacb(iv,2,3) ) )
       yprime(iv,2) = D*( (bvec(iv,2) * Fr - bvec(iv,1) * Fp ) * over_b2      &
            + cmrho * bvec(iv,3) &
            + cmrho2 * (bvec(iv,2)*over_r + jacb(iv,2,1)-jacb(iv,1,2)*over_r) )
       yprime(iv,3) = D*( (bvec(iv,1) * Fz - bvec(iv,3) * Fr) * over_b2         &
            + cmrho * bvec(iv,2)                     &
            + cmrho2 * ( jacb(iv,1,3) - jacb(iv,3,1)) ) * over_r
              
       yprime(iv,4) = D * over_b2 *( &
            bvec(iv,1) * Fr + bvec(iv,3) * Fz + bvec(iv,2) * Fp &
            + y(iv,4) * ( Fr * ( jacb(iv,3,2) * over_r - jacb(iv,2,3)) &
            + Fz * (bvec(iv,2) * over_r + jacb(iv,2,1) - jacb(iv,1,2) * over_r)  &
            + Fp * (jacb(iv,1,3) - jacb(iv,3,1))) )

       !yp_exb(iv,1)= D*(bvec(iv,3) * evec(iv,2) - Bvec(iv,2) * evec(iv,3)) * over_b2
       !yp_exb(iv,2)= D*(bvec(iv,2) * evec(iv,1) - bvec(iv,1) * evec(iv,2)) * over_b2
       !yp_exb(iv,3)= D*(bvec(iv,1) * evec(iv,3) - bvec(iv,3) * evec(iv,1)) * over_b2 * over_r
       
       
    end do
    
    
  end function eom_eval


end module eom
