module my_modules
  implicit none
  contains
    FUNCTION transform(r,s,t,Bx,By,Bz,Px,Py,Pz) 
      real(8), DIMENSION(3) :: transform
      real(8) :: r,s,t,Bx,By,Bz,Px,Py,Pz
      real(8) :: alpha, beta, gamma
      !-------------------------------------------------
      ! Transformations:
      alpha = sqrt((Bz*Py-By*Pz)**2+(Bx*Pz-Bz*Px)**2+(By*Px-Bx*Py)**2)
      beta  = sqrt(Bx**2+By**2+Bz**2)
      gamma = sqrt((Bx**2+By**2+Bz**2)*((Bz*Py-By*Pz)**2+(Bx*Pz-Bz*Px)**2+(By*Px-Bx*Py)**2))
      transform(1) = (r*(-Bz*Py+By*Pz)+s*(-Bx*Pz+Bz*Px)+t*(-By*Px+Bx*Py)) / alpha
      transform(2) = (r*Bx+s*By+t*Bz) / beta
      transform(3) = -(r*(-By**2*Px+Bx*By*Py+Bx*Bz*Pz-Bz**2*Px)+ &
           s*(-Bz**2*Py+By*Bz*Pz+Bx*By*Px-Bx**2*Py)+t*(-Bx**2*Pz+Bx*Bz*Px+By*Bz*Py-By**2*Pz)) / gamma

    END FUNCTION transform

    FUNCTION inv_transform(r,s,t,Bx,By,Bz,Px,Py,Pz) 
      real(8), DIMENSION(3) :: inv_transform
      real(8) :: r,s,t,Bx,By,Bz,Px,Py,Pz
      real(8) :: a11,a12,a13,a21,a22,a23,a31,a32,a33
      real(8) :: inv11,inv12,inv13,inv21,inv22,inv23,inv31,inv32,inv33,det
      real(8) :: alpha, beta, gamma
      !-------------------------------------------------
      ! CHECKED !
      ! Transformations:
      alpha = sqrt((Bz*Py-By*Pz)**2+(Bx*Pz-Bz*Px)**2+(By*Px-Bx*Py)**2)
      beta  = sqrt(Bx**2+By**2+Bz**2)
      gamma = sqrt((Bx**2+By**2+Bz**2)*((Bz*Py-By*Pz)**2+(Bx*Pz-Bz*Px)**2+(By*Px-Bx*Py)**2))

      a11 = (-Bz*Py+By*Pz)/alpha
      a12 = (-Bx*Pz+Bz*Px)/alpha
      a13 = (-By*Px+Bx*Py)/alpha
      
      a21 = Bx/beta
      a22 = By/beta
      a23 = Bz/beta

      a31 = -(-By**2*Px+Bx*By*Py+Bx*Bz*Pz-Bz**2*Px)/gamma
      a32 = -(-Bz**2*Py+By*Bz*Pz+Bx*By*Px-Bx**2*Py)/gamma
      a33 = -(-Bx**2*Pz+Bx*Bz*Px+By*Bz*Py-By**2*Pz)/gamma

      inv11 =  (a22*a33-a23*a32)
      inv21 = -(a21*a33-a23*a31)
      inv31 =  (a21*a32-a31*a22)

      inv12 = -(a12*a33-a32*a13)
      inv22 =  (a11*a33-a31*a13)
      inv32 = -(a11*a32-a12*a31)
      
      inv13 =  (a12*a23-a13*a22)
      inv23 = -(a11*a23-a13*a21)
      inv33 =  (a11*a22-a12*a21)

      det = a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a31*a22)
      
      inv_transform(1) = (inv11*r+inv12*s+inv13*t)/det
      inv_transform(2) = (inv21*r+inv22*s+inv23*t)/det
      inv_transform(3) = (inv31*r+inv32*s+inv33*t)/det
    END FUNCTION inv_transform
  end module my_modules

program PTUMF_3D
  use my_modules
  implicit none
!Initialise:
  real(8) :: x,y,z,x0,y0,z0,step,u,v,w,theta,bin_dim
  real(8) :: Bx,By,Bz,B0x,B0y,B0z,Px,Py,Pz,P0x,P0y,P0z
  real(8) :: dx,dy,dz,du,dv,dw
  real(8) :: Pu,Pv,Pw,Bu,Bv,Bw
  real(8), DIMENSION(3) :: dummy
  real(8) :: size, rho, ppr
  real(8), parameter :: PI = atan2(0.0D0,-1.0D0)
  integer :: a,d,h,l,i,iter,j,ln
  character(2) :: file_no
  !character(11) :: shape
  
!--------------------------------------------------
!Parsing script:
  INTEGER :: nlines, io
  CHARACTER (len=8) :: file 
  TYPE :: field_array
     REAL :: X_COORD, Y_COORD, Z_COORD
     REAL(8) :: B_X, B_Y, B_Z
  END TYPE field_array
  TYPE(field_array), DIMENSION(127756) :: field_container
  TYPE(field_array) :: temp
  real(8) :: tmp_avg_Bx,tmp_avg_By,tmp_avg_Bz

  print *,'# Enter filename (hint: name it file.txt and place it in the same dir): '
  read *,file
  
  nlines = 0 
  OPEN (10, file = file)
  DO
     read(10,*,iostat=io)
     IF (io/=0) EXIT
     nlines = nlines + 1
  END DO
  CLOSE(10)
  OPEN (10, file = file)
  DO i=1,nlines
     READ(10,*) field_container(i)
     !print *,'#',i,' field_array vectors -->', field_container(i)
  END DO
  CLOSE(10)

!-------------------------------------------------- 
  !shape = 'square' 
!User input:
  !print *,'# Enter ',shape, ' 3D grid size, a [in units]='
  !read *,a
  print *,'# Enter initial position of particle, (x0,y0,z0) [in units]='
  read *,x0,y0,z0
  print *,'# Enter initial moment of particle, (P0x,P0y,P0z) [in units]='
  read *,P0x,P0y,P0z
  print *,'# Enter bin dimension (hint: consult the input .txt file for this value) = '
  read *,bin_dim
  !print *,'# Enter points per revolution = '
  !read *,ppr
  !print *,'# Enter (B0x,B0y,B0z) [in units], points per revolution = '
  !read *,B0x,B0y,B0z,ppr
  !print *,'# Enter number of iterations ='
  !read *,iter
  !print *,'# Dimensions of ',shape,': a [in units]= ',a
  !print *,'# points per revolution = ',ppr
  print *,'# bin dimension = ',bin_dim
  !print *,'# B [in units], points per revolution = ',B0x,B0y,B0z,', ',ppr
  print *,'# (x0,y0,z0) [in units]= ','(',x0,',',y0,',',z0,')'
  print *,'# (P0x,P0y,P0z) [in units]= ','(',P0x,',',P0y,',',P0z,')'

  !NAMING CONVENTION VALID FROM 10 to 99 GeV
  write(file_no,'(I2)') INT(sqrt(P0x**2+P0y**2+P0z**2))
  open(unit=11, file='PTUMF_3D_'// file_no // 'GeV.dat')
  !file='PTUMF_3D_'//CHAR(INT(sqrt(P0x**2+P0y**2+P0z**2)))//'GeV.dat')
  !-------------------------------------------------
  
  ! Initialise position:
  x = x0
  y = y0
  z = z0
  
  i = 0
  open(unit=12, file='GRID.dat')
  do while(.TRUE.)
     i=i+1
     temp = field_container(i)

     u = temp%X_COORD
     v = temp%Y_COORD
     w = temp%Z_COORD
     !print *,temp%X_COORD,temp%Y_COORD,temp%Z_COORD
     !print *,u,v,w,x,y,z
    
     write(12,*)temp%X_COORD,temp%Y_COORD,temp%Z_COORD

     Bu = temp%B_X
     Bv = temp%B_Y
     Bw = temp%B_Z

     if (x>=u .and. x<u+0.2D0 .and. &
         y>=v .and. y<v+0.2D0 .and. &
         z>=w .and. z<w+0.2D0) then

        tmp_avg_Bx = 0.0D0
        tmp_avg_By = 0.0D0
        tmp_avg_Bz = 0.0D0

        ! index location
        temp = field_container(i)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = Bu
        tmp_avg_By = Bv
        tmp_avg_Bz = Bw
        ! index location for point with increment in y direction
        temp = field_container(i + 1)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in z direction
        temp = field_container(i + 76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in x direction
        temp = field_container(i + 41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in y&z direction
        temp = field_container(i + 1+41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in y&x direction
        temp = field_container(i + 1 + 76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in z&x direction
        temp = field_container(i + 41+76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in y&z&x direction
        temp = field_container(i + 1+41+76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw

        tmp_avg_Bx = tmp_avg_Bx / 8.0D0
        tmp_avg_By = tmp_avg_By / 8.0D0
        tmp_avg_Bz = tmp_avg_Bz / 8.0D0

        Bx = tmp_avg_Bx
        By = tmp_avg_By
        Bz = tmp_avg_Bz
     else
        print *,'Iteration NO. ',i,'OUT OF BOUNDS OR BIN BOUNDARY ERROR!'
     end if
     if (i == 127756) then
        exit
     end if
  end do
  close(12)

  ! Initialise magnetic field:
  !Bx = 0.0D0
  !By = 1.0D0
  !Bz = 0.0D0
  
  ! Initialise magnetic field:
  !B0x = Bx
  !B0y = By
  !B0z = Bz

! Initialise moment & position:
  Px = P0x
  Py = P0y 
  Pz = P0z

!-------------------------------------------------
     
!-------------------------------------------------
  
  dummy = transform(Bx,By,Bz,Bx,By,Bz,Px,Py,Pz)
  Bu = dummy(1)
  Bv = dummy(2)
  Bw = dummy(3)

  ! Bending angle (theta) in radians:
  dummy = transform(Px,Py,Pz,Bx,By,Bz,Px,Py,Pz)
  Pu = dummy(1)
  Pv = dummy(2)
  Pw = dummy(3)

  dummy = transform(x,y,z,Bx,By,Bz,Px,Py,Pz)
  u = dummy(1)
  v = dummy(2)
  w = dummy(3)
  
  rho   = 3.33D0*sqrt(Pu**2+Pv**2+Pw**2)/sqrt(Bu**2+Bv**2+Bw**2)
  !theta = step / (rho * sqrt(1+(Pv/Pw)**2))

  !SELF-ADJUSTABLE PPR PARAMETER (1 tracepoint per bin objective)
  !perform stepwise scaling
  ppr =  max(abs(Pv/Pw*rho*2.0D0*PI),sqrt(2.0D0*rho**2+(Pv/Pw*rho*PI/2.0D0)**2), &
             sqrt(4.0D0*rho**2+(Pv/Pw*rho*PI)**2)) / bin_dim 
  theta =  PI/(ppr/2)

  du  =  rho*(cos(theta)-1.0D0)
  dw  =  rho*sin(theta)
  dv  =  Pv/Pw*rho*theta

  i = 0
  
!-------------------------------------------------
! Computations:
  do while (.TRUE.)
     !write(11,*)i,u,v,w,Pu,Pv,Pw,Bu,Bv,Bw,rho,theta
     ! Debug mode
     !write(11,*)i,u,v,w,x,y,z,Pu,Pv,Pw,Px,Py,Pz,Bu,Bv,Bw,Bx,By,Bz,du,dv,dw,rho,step
     write(11,*)x,y,z,Bx,By,Bz,Px,Py,Pz
     !write(11,*)x,y,z,Bx,By,Bz,step,rho
     ! Plot mode
     !write(11,*)i,u,v,w,x,y,z,rho
     i = i + 1


     !SELF-ADJUSTABLE PPR PARAMETER (1 tracepoint per bin objective)
     !perform stepwise scaling
     ppr =  max(abs(Pv/Pw*rho*2.0D0*PI),sqrt(2.0D0*rho**2+(Pv/Pw*rho*PI/2.0D0)**2), &
             sqrt(4.0D0*rho**2+(Pv/Pw*rho*PI)**2)) / bin_dim 
     theta =  PI/(ppr/2)
     
     !step  = sqrt((rho*theta)**2+dv**2)
     ! For constant B, equation of motion can be solved exactly
     du  =  rho*(cos(theta)-1.0D0)
     !dv doesn't describe motion as expected [ATTENTION HERE]!
     dv  =  Pv/Pw*rho*theta
     dw  =  rho*sin(theta)
     !dv = sqrt(step**2-du**2-dw**2)
     step = sqrt(du**2+dv**2+dw**2)
     
     ! Bending radius (rho):
     rho   = 3.33D0*sqrt(Pu**2+Pv**2+Pw**2)/sqrt(Bu**2+Bv**2+Bw**2)
     !step  = sqrt(du**2+dv**2+dw**2)
     !step = sqrt((rho*theta)**2+dv**2)
     !theta = sqrt(step**2 - dv**2) / rho

     !theta = step / (rho * sqrt(1+(Pv/Pw)**2))

     u = u+du
     v = v+dv
     w = w+dw

     !dummy = inv_transform(u,v,w,Bx,By,Bz,P0x,P0y,P0z)
     dummy = inv_transform(u,v,w,Bx,By,Bz,Px,Py,Pz)
     x = dummy(1)
     y = dummy(2)
     z = dummy(3)

!-------------------------------------------------
     
!-------------------------------------------------
     
     Pu = -Pw*sin(theta)
     Pv =  Pv
     Pw =  Pw*cos(theta)

     dummy = inv_transform(Pu,Pv,Pw,Bx,By,Bz,Px,Py,Pz)
     Px = dummy(1)
     Py = dummy(2)
     Pz = dummy(3)

     j = 0
     do while(.TRUE.)
     j=j+1
     temp = field_container(j)

     u = temp%X_COORD
     v = temp%Y_COORD
     w = temp%Z_COORD
     !print *,temp%X_COORD,temp%Y_COORD,temp%Z_COORD
     !print *,u,v,w
     
     Bu = temp%B_X
     Bv = temp%B_Y
     Bw = temp%B_Z

     if (x>=u .and. x<u+0.2D0 .and. &
         y>=v .and. y<v+0.2D0 .and. &
         z>=w .and. z<w+0.2D0) then

        tmp_avg_Bx = 0.0D0
        tmp_avg_By = 0.0D0
        tmp_avg_Bz = 0.0D0
        
        ! index location
        temp = field_container(j)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = Bu
        tmp_avg_By = Bv
        tmp_avg_Bz = Bw
        ! index location for point with increment in y direction
        temp = field_container(j + 1)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z
        
        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in z direction
        temp = field_container(j + 76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in x direction
        temp = field_container(j + 41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in y&z direction
        temp = field_container(j + 1+41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in y&x direction
        temp = field_container(j + 1 + 76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in z&x direction
        temp = field_container(j + 41+76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        ! index location for point with increment in y&z&x direction
        temp = field_container(j + 1+41+76*41)
        Bu = temp%B_X
        Bv = temp%B_Y
        Bw = temp%B_Z

        tmp_avg_Bx = tmp_avg_Bx + Bu
        tmp_avg_By = tmp_avg_By + Bv
        tmp_avg_Bz = tmp_avg_Bz + Bw
        
        tmp_avg_Bx = tmp_avg_Bx / 8.0D0
        tmp_avg_By = tmp_avg_By / 8.0D0
        tmp_avg_Bz = tmp_avg_Bz / 8.0D0

        Bx = tmp_avg_Bx
        By = tmp_avg_By
        Bz = tmp_avg_Bz

        print *,Bx,By,Bz
     else
        print *,'Iteration NO. ',i,j,' OUT OF BOUNDS OR BIN BOUNDARY ERROR!'
     end if
     if (j == 127756) then
        exit
     end if
  end do

!-------------------------------------------------
     
!-------------------------------------------------

     
     !dummy = inv_transform(Pu,Pv,Pw,Bx,By,Bz,P0x,P0y,P0z)
    !dummy = inv_transform(Pu,Pv,Pw,Bx,By,Bz,Px,Py,Pz)
    ! Px = dummy(1)
    ! Py = dummy(2)
    ! Pz = dummy(3)
     
     dummy = transform(x,y,z,Bx,By,Bz,Px,Py,Pz)
     u = dummy(1)
     v = dummy(2)
     w = dummy(3)

     dummy = transform(Px,Py,Pz,Bx,By,Bz,Px,Py,Pz)
     Pu = dummy(1)
     Pv = dummy(2)
     Pw = dummy(3)

     dummy = transform(Bx,By,Bz,Bx,By,Bz,Px,Py,Pz)
     Bu = dummy(1)
     Bv = dummy(2)
     Bw = dummy(3)

     write(11,*)x,y,z,Bx,By,Bz,Px,Py,Pz
     
     if (abs(x)>10.0D0 .or. abs(y)>4.0D0 .or. abs(z)>4.0D0) then
        exit ! exit the loop
     end if
  end do
  close(11)
  print *,'# *Recap*: ITERATIONS= ',i,' POINTS PER REVOLUTION= ',ppr 
  print *,'# The .dat file has been written!'
end program PTUMF_3D
