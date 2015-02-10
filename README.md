# eka_testi
First attempt at using GitHub
by P. Räisänen, FMI

PROGRAM koko_testi

  IMPLICIT NONE
  INTEGER, PARAMETER :: ndim=100

  REAL, DIMENSION(ndim,ndim,ndim) :: f3d

  INTEGER :: i,j,k

  DO i=1,ndim
    DO j=1,ndim
      DO k=1,ndim
       f3d(i,j,k) = 1.5*i+j+2*k  
       IF (k>1) f3d(i,j,k) = f3d(i,j,k)+0.001*f3d(i,j,k-1)       
      ENDDO
      write(*,*) i,j,f3d(i,j,ndim)
    ENDDO
  ENDDO

END PROGRAM koko_testi
