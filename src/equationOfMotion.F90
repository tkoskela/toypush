!> This module calculates the yprime terms for the RK interpolation routine for the particle push
!> mini-app. Mostly copied and cleaned up from XGC1/kernels/pushe/vec/derivs_elec_vec.F90
!> @author Tuomas Koskela
!> @date 2 Sep 2016
!<
module eom

  use params
  
  implicit none

contains

#ifndef NODIVSQRT
  !> Evaluates the LHS of the equations of motion. See equation (1) in Chang et al, PoP, 2009
  !<
  function eom_eval(y,Bvec,jacB,Evec,dt,yprime,mu,charge,mass) result(err)

    double precision,    intent(in) :: dt  !> Time step length
    double precision,    intent(in),  dimension(veclen) :: charge  !> Charge of particle
    double precision,    intent(in),  dimension(veclen) :: mass    !> Mass of particle
    double precision,    intent(in),  dimension(veclen) :: mu !> Magnetic moment of particle
    double precision,    intent(in),  dimension(4,veclen) :: y !> Particle phase information
    double precision,    intent(in),  dimension(3,veclen) :: Evec !> Electric field vector (Er,Ephi,Ez)
    double precision,    intent(in),  dimension(3,veclen) :: Bvec !> Magnetic field vector (Br,Bphi,Bz)
    !> Magnetic field jacobian
    !> |  dBRdR,   dBRdphi,   dBRdz  |
    !> | dBphidR, dBphidphi, dBphidz |
    !> |  dBzdR,   dBzdphi,   dBzdz  |
    !<
    double precision,    intent(in),  dimension(3,3,veclen)  :: jacB

    !> LHS of 1st (1:3) and 2nd (4) equations of motion
    !<
    double precision,    intent(out), dimension(4,veclen)     :: yprime 

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
       B2 = dot_product(bvec(:,iv),bvec(:,iv))
       B  = sqrt(B2)
       
       ! Divisions
       c_m = charge(iv) / mass(iv)
       cmrho = c_m * y(4,iv)
       cmrho2 = cmrho * y(4,iv)
       over_r = 1.0D0 / y(1,iv)
       over_B  = 1.0D0 / B
       over_B2 = 1.0D0 / B2       

       D = 1.0D0 / (1.0D0 + y(4,iv) * (&
            -1.D0 * over_B2 * ( jacb(2,3,iv) * bvec(1,iv) &
            + (jacb(3,1,iv) - jacb(1,3,iv)) * bvec(2,iv) &
            - (bvec(2,iv) * over_r + jacb(2,1,iv)) * bvec(3,iv) )))

       Fr = evec(1,iv) - &
            ( (mu(iv) + charge(iv) * cmrho2 * B) & ! murho2b
            * (bvec(1,iv) * jacb(1,1,iv) + bvec(2,iv) * jacb(2,1,iv) &
            + bvec(3,iv) * jacb(3,1,iv)) * over_B & ! dbdr
            ) 
       Fp = evec(2,iv) ! dbdphi is 0 by definition
       Fz = evec(3,iv) - &
            ( (mu(iv) + charge(iv) * cmrho2 * B) & ! murho2b
            * (bvec(1,iv) * jacb(1,3,iv) + bvec(2,iv) * jacb(2,3,iv) &
            + bvec(3,iv) * jacb(3,3,iv)) * over_B & ! dbdz
            )
       
       yprime(1,iv) = D*( (bvec(3,iv) * Fp - bvec(2,iv) * Fz) * over_b2         &
            + cmrho * bvec(1,iv) &
            + cmrho2 * (jacb(3,2,iv) * over_r - jacb(2,3,iv) ) )
       yprime(2,iv) = D*( (bvec(2,iv) * Fr - bvec(1,iv) * Fp ) * over_b2      &
            + cmrho * bvec(3,iv) &
            + cmrho2 * (bvec(2,iv)*over_r + jacb(2,1,iv)-jacb(1,2,iv)*over_r) )
       yprime(3,iv) = D*( (bvec(1,iv) * Fz - bvec(3,iv) * Fr) * over_b2         &
            + cmrho * bvec(2,iv)                     &
            + cmrho2 * ( jacb(1,3,iv) - jacb(3,1,iv)) ) * over_r
              
       yprime(4,iv) = D * over_b2 *( &
            bvec(1,iv) * Fr + bvec(3,iv) * Fz + bvec(2,iv) * Fp &
            + y(4,iv) * ( Fr * ( jacb(3,2,iv) * over_r - jacb(2,3,iv)) &
            + Fz * (bvec(2,iv) * over_r + jacb(2,1,iv) - jacb(1,2,iv) * over_r)  &
            + Fp * (jacb(1,3,iv) - jacb(3,1,iv))) )

       !yp_exb(1,iv)= D*(bvec(3,iv) * evec(2,iv) - Bvec(2,iv) * evec(3,iv)) * over_b2
       !yp_exb(2,iv)= D*(bvec(2,iv) * evec(1,iv) - bvec(1,iv) * evec(2,iv)) * over_b2
       !yp_exb(3,iv)= D*(bvec(1,iv) * evec(3,iv) - bvec(3,iv) * evec(1,iv)) * over_b2 * over_r
       
       
    end do
    
    
  end function eom_eval

#else
  !> "Optimized" version of eom_eval to push for peak performance. NOTE: the calculation result is incorrect
  !<
  function eom_eval(y,Bvec,jacB,Evec,dt,yprime,mu,charge,mass) result(err)

    double precision,    intent(in) :: dt  !> Time step length
    double precision,    intent(in),  dimension(veclen) :: charge  !> Charge of particle
    double precision,    intent(in),  dimension(veclen) :: mass    !> Mass of particle
    double precision,    intent(in),  dimension(veclen) :: mu !> Magnetic moment of particle
    double precision,    intent(in),  dimension(4,veclen) :: y !> Particle phase information
    double precision,    intent(in),  dimension(3,veclen) :: Evec !> Electric field vector (Er,Ephi,Ez)
    double precision,    intent(in),  dimension(3,veclen) :: Bvec !> Magnetic field vector (Br,Bphi,Bz)
    !> Magnetic field jacobian
    !> |  dBRdR,   dBRdphi,   dBRdz  |
    !> | dBphidR, dBphidphi, dBphidz |
    !> |  dBzdR,   dBzdphi,   dBzdz  |
    !<
    double precision,    intent(in),  dimension(3,3,veclen)  :: jacB

    !> LHS of 1st (1:3) and 2nd (4) equations of motion
    !<
    double precision,    intent(out), dimension(4,veclen)     :: yprime 

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
       B2 = dot_product(bvec(:,iv),bvec(:,iv))
       B  = B2
       
       ! No Divisions
       c_m = charge(iv) * mass(iv)
       cmrho = c_m * y(4,iv)
       cmrho2 = cmrho * y(4,iv)
       over_r = 1.0D0 * y(1,iv)
       over_B  = 1.0D0 * B
       over_B2 = 1.0D0 * B2       

       D = 1.0D0 * (1.0D0 + y(4,iv) * (&
            -1.D0 * over_B2 * ( jacb(2,3,iv) * bvec(1,iv) &
            + (jacb(3,1,iv) - jacb(1,3,iv)) * bvec(2,iv) &
            - (bvec(2,iv) * over_r + jacb(2,1,iv)) * bvec(3,iv) )))

       Fr = evec(1,iv) - &
            ( (mu(iv) + charge(iv) * cmrho2 * B) & ! murho2b
            * (bvec(1,iv) * jacb(1,1,iv) + bvec(2,iv) * jacb(2,1,iv) &
            + bvec(3,iv) * jacb(3,1,iv)) * over_B & ! dbdr
            ) 
       Fp = evec(2,iv) ! dbdphi is 0 by definition
       Fz = evec(3,iv) - &
            ( (mu(iv) + charge(iv) * cmrho2 * B) & ! murho2b
            * (bvec(1,iv) * jacb(1,3,iv) + bvec(2,iv) * jacb(2,3,iv) &
            + bvec(3,iv) * jacb(3,3,iv)) * over_B & ! dbdz
            )
       
       yprime(1,iv) = D*( (bvec(3,iv) * Fp - bvec(2,iv) * Fz) * over_b2         &
            + cmrho * bvec(1,iv) &
            + cmrho2 * (jacb(3,2,iv) * over_r - jacb(2,3,iv) ) )
       yprime(2,iv) = D*( (bvec(2,iv) * Fr - bvec(1,iv) * Fp ) * over_b2      &
            + cmrho * bvec(3,iv) &
            + cmrho2 * (bvec(2,iv)*over_r + jacb(2,1,iv)-jacb(1,2,iv)*over_r) )
       yprime(3,iv) = D*( (bvec(1,iv) * Fz - bvec(3,iv) * Fr) * over_b2         &
            + cmrho * bvec(2,iv)                     &
            + cmrho2 * ( jacb(1,3,iv) - jacb(3,1,iv)) ) * over_r
              
       yprime(4,iv) = D * over_b2 *( &
            bvec(1,iv) * Fr + bvec(3,iv) * Fz + bvec(2,iv) * Fp &
            + y(4,iv) * ( Fr * ( jacb(3,2,iv) * over_r - jacb(2,3,iv)) &
            + Fz * (bvec(2,iv) * over_r + jacb(2,1,iv) - jacb(1,2,iv) * over_r)  &
            + Fp * (jacb(1,3,iv) - jacb(3,1,iv))) )

       !yp_exb(1,iv)= D*(bvec(3,iv) * evec(2,iv) - Bvec(2,iv) * evec(3,iv)) * over_b2
       !yp_exb(2,iv)= D*(bvec(2,iv) * evec(1,iv) - bvec(1,iv) * evec(2,iv)) * over_b2
       !yp_exb(3,iv)= D*(bvec(1,iv) * evec(3,iv) - bvec(3,iv) * evec(1,iv)) * over_b2 * over_r
       
       
    end do
    
    
  end function eom_eval
#endif  

end module eom
