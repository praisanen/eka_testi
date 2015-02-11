PROGRAM snow_ssp
!
! This Fortran program demonstrates the use of the snow single-scattering
! property parameterization presented in
!
! P. Räisänen, A. Kokhanovsky, G. Guyot, O. Jourdan and T. Nousiainen:
! Parameterization of single-scattering properties of snow.
! The Cryosphere Discuss., 9, 2015.
!
! If you have questions, comments, bugs to report etc., please
! send e-mail to P. Räisänen (petri.raisanen@fmi.fi).
!
!************************************************************************
!
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307) 

  INTEGER, PARAMETER :: nstr=32     ! Number os streams for DISORT
  INTEGER, PARAMETER :: ntheta=181  ! Number of scattering angles 
                                    ! for phase function output

! Input for PARAMETERIZATION
 
  REAL(dp) :: refre ! Real part of refractive index   
  REAL(dp) :: refim ! Imaginary part of refractive index
  REAL(dp) :: sizep ! size parameter defined with respect to 
                    ! the volume-to-projected area equivalent radius

  REAL(dp), DIMENSION(ntheta) :: theta ! Scattering angles for phase
                                       ! function output  

  LOGICAL :: lfull  ! Logical switch for using the full (.TRUE.) or
                    ! the simpler phase function parameterization (.FALSE.)

! Output from PARAMETERIZATION

  REAL(dp) :: qext  ! extinction efficiency
  REAL(dp) :: coalb ! single-scattering co-albedo (= 1-single-scattering albedo)
  REAL(dp) :: g     ! asymmetry parameter

  REAL(dp), DIMENSION(0:nstr) :: p ! phase function Legendre coefficients
  REAL(dp), DIMENSION(ntheta) :: p11 ! phase function computed using the 
                                     ! Legendre expansion for P_resid
  REAL(dp), DIMENSION(ntheta) :: p11_2 ! phase function computed using the 
                                       ! the ordinary polynomial expansion
                                       ! (just for a check, should equal p11!)

! Miscellaneous
  INTEGER:: itheta,i
 
! Scattering angles
  DO itheta=1,ntheta
    theta(itheta)=180._dp*REAL(itheta-1)/REAL(ntheta-1)
  ENDDO   

  refre = 1.3942999839782715     
  refim = 9.5650001807889851E-011
  sizep = 63.147590085381225  

  lfull=.TRUE.

  CALL parameterization(nstr,ntheta,refre,refim,sizep,lfull,theta, &
                        qext,coalb,g,p,p11,p11_2)

  write(*,*) qext
  write(*,*) coalb
  write(*,*) g

  WRITE(*,*)
  DO i=0,nstr
    write(*,'(I4,ES12.4)') i,p(i)
  ENDDO

  DO itheta=1,ntheta
    write(*,'(I4,F8.2,2ES12.4,ES12.4)') itheta,theta(itheta),&
           p11(itheta), p11_2(itheta), (p11_2(itheta)/p11(itheta)-1._dp)
!   write(*,*) p11_2(itheta)/p11(itheta)
  ENDDO

!*********
CONTAINS
!*********

  SUBROUTINE parameterization(nstr,ntheta,refre,refim,sizep,lfull,theta, &
                              qext,coalb,g,p,p11,p11_2)
!
! This subroutine contains the parameterization equations for
! single-scattering properties of snow
!
! INPUT:
!---------
!
! NSTR = length of phase function Legendre series for DISORT 
!        ("number of streams")
! NTHETA = number of scattering angles at which phase function values
!          are computed
! THETA(NTHETA) = scattering angles (in degrees) for the computation
!                 of phase function values
! LFULL = switch for using the "full" phase function patameterization
!         (LFULL=.TRUE.) or the simpler parameterization without the
!         P_resid part in phase function (LFULL=.FALSE.) (see Eqs. 14,30)
! SIZEP = size parameter defined with respect to the volume-to-projected
!         area equivalent radius
! REFRE = real part of refractive index
! REFIM = imaginary part of refractive index   
!
! OUTPUT:
!--------
! QEXT  = extinction efficiency
! COALB = single-scattering co-albedo (= 1-single-scattering albedo)
! G     = asymmetry parameter
! P(0:NSTR) = phase function Legendre coefficients
! P11(NTHETA) = values of phase function at the NTHETA scattering angles,
!                  computed using the Legendre expansion for P_resid
! P11_2(NTHETA) = values of phase function at the NTHETA scattering angles,
!                   computed using the ordinary polynomial expansion
!                   (just for checking that it works correctly!)
!**************************************************************************

    IMPLICIT NONE

    REAL, PARAMETER :: pi=3.1415926536_dp
          
! Input parameters
    INTEGER, INTENT(in) :: nstr, ntheta
    REAL(dp), INTENT(in) :: refre, refim, sizep
    REAL(dp), DIMENSION(ntheta), INTENT(in) :: theta
    LOGICAL, INTENT(in) :: lfull

! Output parameters
    REAL(dp), INTENT(out) :: qext, coalb, g
    REAL(dp), DIMENSION(0:nstr), INTENT(OUT) :: p
    REAL(dp), DIMENSION(ntheta), INTENT(OUT) :: p11, p11_2   

! Local variables
    INTEGER :: n
    REAL(dp), DIMENSION(0:nstr,ntheta) :: lpoly !Legendre polynomials
    REAL(dp), DIMENSION(0:6,4) :: c ! Coefficients needed in Eq. (27)
    REAL(dp), DIMENSION(0:5,4) :: d ! Coefficients needed in Eq. (29)

    REAL(dp), DIMENSION(ntheta) :: costh

! VAATISIVATKO NÄMÄ KOMMENTTEJA?
    REAL(dp), DIMENSION(0:nstr) :: a
    REAL(dp), DIMENSION(0:5) :: b

    REAL(dp), DIMENSION(ntheta) :: pdiff_plus_ray, presid, presid2

    REAL(dp) :: xabs, ssa, gdiff, gray, g1, wdiff, wray, w1

!**************************************************************************
! Define the coefficients needed in parameterizing phase function Legendre
! moments
    CALL define_coeffs(c,d)

! Compute Legdendre polynomials (0..NSTR) at the user-defined scattering 
! angles THETA
    costh(1:ntheta)=COS(PI/180._dp*theta(1:ntheta))
    CALL legendre_poly (nstr,ntheta,costh,lpoly)

!***************************************************
! And here starts the actual parameterization ...
!**************************************************
!
! Extinction efficiency (Section 6.1)
    qext = 2._dp

! Single-scattering co-albedo (Section 6.2)
! Size parameter for absorption (Eq. 12)
    xabs = sizep * refim * refre**2
! Eq. (11)
    coalb = 0.470_dp*(1._dp-exp(-2.69_dp*xabs* &
                     (1._dp-0.31_dp*MIN(xabs,2.)**0.67_dp)))

! Asymmetry parameter (Section 6.3, Eq. 13)
    g = 1._dp-1.146_dp*(refre-1._dp)**0.8_dp*(0.52_dp-coalb)**1.05_dp &
             *(1._dp+8._dp*sizep**(-1.5_dp))

!**************************************************
! Parameterization of phase function (Section 6.4)
!**************************************************
! Asymmetry parameter for diffraction (Eq. 19)
    gdiff=1._dp-0.60_dp/sizep

! Fractional weights for the diffraction and ray optics parts (Eqs. 15-16)
    ssa = 1._dp-coalb
    wdiff = 1._dp / (qext * ssa)
    wray = 1._dp-wdiff

! Asymmetry parameter for the ray optics part (Eq. 23)
    gray  = (g-wdiff*gdiff)/wray

! Weight factor for the Henyey-Greenstein part of the "ray optics"
! phase function (Eq. 22)
    w1 = 1._dp-1.53_dp*MAX(0.77_dp-gray,0._dp)**1.2_dp
! The corresponding asymmetry parameter, Eq. (24)
    g1= gray / w1

! Legendre expansion of the phase function residual 
! (used for the "full" parameterization only; LFULL=.TRUE.)
    IF (lfull) THEN
      DO n=0,6  ! Eq. (27)
        a(n) = c(n,1) + c(n,2)*coalb + c(n,3)*g + c(n,4)*coalb*g 
      ENDDO     
! For higher-order indices, use the coefficient A(6).
! Note that these become zero after delta-M-scaling in DISORT.
      a(7:nstr) = a(6) 
    ELSE
! For the simpler parametrization, set these coefficients to zero
      a(0:nstr) = 0._dp
    END IF

! Phase function Legendre coefficients for use in DISORT.
! This assumes delta-M-scaling is used in DISORT (i.e., DELTAM=.TRUE.)
! Eq. (31), taking into account that A(7:nstr) = A(6).

    p(0) = 1._dp
    DO n=1,nstr   
      p(n) = wdiff * gdiff**n + wray*w1*g1**n + a(n)
    ENDDO

!*********************************************************************
! Compute phase function values at the user-defined scattering angles
!********************************************************************* 
! "piff_plus_ray" covers the first three terms in Eq. (30)
!
    pdiff_plus_ray(:) = &
       wdiff*(1.-gdiff*gdiff)/(1+gdiff*gdiff-2*gdiff*costh(:))**1.5 &
     + wray*w1*(1.-g1*g1)/(1+g1*g1-2*g1*costh(:))**1.5 + wray*(1.-w1)  
    
! Computation of the residual part of the phase function 
! "presid" for the Legendre expansion
! "presid2" for the ordindary polynomial expansion (just a check)

    IF (lfull) THEN
      presid(:) = 0._dp
! Eq. (26) (without Dirac's delta function at theta=0):
      DO n=0,5
        presid(:)=presid(:)+(2*n+1)*(a(n)-a(6))*lpoly(n,:)
      ENDDO 
! Eqs. (28) and (29): (without Dirac's delta function at theta=0): 
      presid2(:) = 0._dp        
      DO n=0,5
        b(n) = d(n,1) + d(n,2)*coalb + d(n,3)*g + d(n,4)*coalb*g
        presid2(:) = presid2(:) + b(n)*costh(:)**n 
      ENDDO
    ELSE
! For the simpler parameterization, set the residuals to zero
      presid(:) = 0._dp 
      presid2(:) = 0._dp 
    END IF

! Eq. (30), using either the Legendre expansion or the ordinary
! poynomial expansion
    p11(:) = pdiff_plus_ray(:) + presid(:)    
    p11_2(:) = pdiff_plus_ray(:) + presid2(:)    

    RETURN
  END SUBROUTINE parameterization

!************************************************************************
  SUBROUTINE define_coeffs(c,d)
!*************************************************************************
! OUTPUT:
! - COEFFS_C = coefficients used to compute the Legendre coefficients
!              for use in the phase function parameterization (Eq. 27)
!          
! - COEFFS_D = corresponding coefficients for computing the coefficients
!              in an ordinary polynomial fit (Eq. 29)
!***********************************************************************
    IMPLICIT NONE
    REAL(dp), DIMENSION(0:6,4), INTENT(out) :: c
    REAL(dp), DIMENSION(0:5,4), INTENT(out) :: d

! Coefficients for computing Legendre coefficients in Eq. (27)

    c(0,1:4)=(/  0.00000,  0.00000,  0.00000,  0.00000/)
    c(1,1:4)=(/  0.00000,  0.00000,  0.00000,  0.00000/)
    c(2,1:4)=(/ -0.01400, -0.10367,  0.02144,  0.08903/)
    c(3,1:4)=(/ -0.13184, -0.01741,  0.16890, -0.06365/)
    c(4,1:4)=(/ -0.20878, -0.03438,  0.27353, -0.10418/)
    c(5,1:4)=(/ -0.29763, -0.06931,  0.38501, -0.11329/)
    c(6,1:4)=(/ -0.32153, -0.10691,  0.41282, -0.07934/)

! Coefficients used for computing ordinary polynomial coefficients
! in Eq. (29)

    d(0,1:4)=(/-0.06679,  0.34357,  0.09553, -0.42542/)
    d(1,1:4)=(/-0.53413,  0.15642,  0.74905, -0.62700/)
    d(2,1:4)=(/-1.49866, -2.42334,  1.76580,  2.10118/)
    d(3,1:4)=(/ 1.01884, -2.05239, -1.59160,  3.54237/)
    d(4,1:4)=(/ 4.43936,  2.85558, -5.48475, -0.97817/)
    d(5,1:4)=(/ 2.07065,  3.25673, -2.40933, -2.94094/)

  END SUBROUTINE define_coeffs
!*********************************************************************
  SUBROUTINE legendre_poly (npoly,nx,x,lpoly)
! Compute the values of Legendre polynomials 0...NPOLY
! for argument values X(1) ... X(NX).

    IMPLICIT NONE
    INTEGER, INTENT(in) :: npoly,nx
    REAL(dp), DIMENSION(nx), INTENT(in) :: x 
    REAL(dp), DIMENSION(0:npoly,nx), INTENT(out) :: lpoly 

    INTEGER :: i,l
    REAL(dp) :: sum

    DO i=1,nx
      lpoly(0,i) = 1.
      lpoly(1,i) = x(i)
      DO l=1,npoly-1
! Arfken (1970, p. 541) 
        lpoly(l+1,i) = 2*x(i)*lpoly(l,i)-lpoly(l-1,i) &
                      -(x(i)*lpoly(l,i)-lpoly(l-1,i))/(l+1._dp)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE legendre_poly
!********************************************************************  
END PROGRAM snow_ssp
