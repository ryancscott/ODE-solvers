program euler
implicit none
!====================================================================
!PROGRAM NAME: euler.f90
!PROGRAM AUTHOR: Ryan C. Scott
!--------------------------------------------------------------------
!This code computes a numerical approximation to the ODE:
!                      y' = f(t,y) = RHS
!using the first order (local truncation error) Euler method
!====================================================================
! VARIABLE/FUNCTION DECLARATIONS  
real :: a,b                 ! solution interval
real :: y_0                 ! initial value, step size       
real :: h                   ! time step            
integer :: N,i              ! number of time steps, loop index
     
!create vector of arbitrary size TBD 
real, dimension(:), allocatable :: y 
real, dimension(:), allocatable :: t

!prompt for...
!IC to being computation
write(*,*) "Enter the initial condition y_0:"
read(*,*) y_0
!step size
print* , "Enter the left domain point a:"
read(*,*) a
!decay rate
print* , "Enter the right domain point b:"
read(*,*) b
!number of steps
print* , "Enter the number of steps N:"
read(*,*) N

!preliminary input based calculations
h = (b-a)/N                          ! step size
write(*,*) "This corresponds to a step size dt:"
write(*,*) h

allocate(y(N+1))                     ! set solution array size
allocate(t(N+1))                     ! set time array size

!compute
do  i = 1,N,1
   y(1) = y_0                        ! initial value
   t(1) = a                          ! initial domain point
   t(i+1) = a+i*h                    ! partition of domain
   y(i+1)=y(i)+h*f(t(i),y(i))        ! successive approximations making use of function
end do

!print description of program output 
print* , "===================================================="
print* , "This program approximates the solution to the ODE..."
print* , "                 dy/dt = f(t,y)                     "
print* , "using the explicit Euler integration scheme. The    "
print* , "approximate solution value y_i at step t_i ...      " 
print* , "for i = 0,1,...,N is:                               " 
print* , "===================================================="
write(*,*) y

contains
!function definition
function f(t,y)
real :: f
real :: t,y
f = t*y+t**3
end function

end program euler


