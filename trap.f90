!============================================================================
!PROGRAM NAME: trap.f90
!PROGRAM AUTHOR: Ryan C. Scott
!----------------------------------------------------------------------------
!This code computes a numerical approximation to the ODE:
!                      y' = f(t,y) = RHS
!using the second order trapezoid method. The function on the
!righthand side, a user-defined FORTRAN function called f, is modified 
!before use. The user will be prompted for necessary input parameters.
!============================================================================
program trap
implicit none
!variable declarations 
real :: a,b                 ! solution interval
real :: y_0                 ! initial value, step size       
real :: h                   ! time step            
integer :: N,i              ! number of time steps, loop index     

!create vector of arbitrary size TBD 
real, dimension(:), allocatable :: y 
real, dimension(:), allocatable :: t

!Should add a prompt for the user to modify code prior to use

!prompt for...
!IC to being computation
write(*,*) "Enter the initial condition y_0:"
read(*,*) y_0
!leftmost solution domain pt
print* , "Enter the left domain point a:"
read(*,*) a
!rightmost solution domain pt
print* , "Enter the right domain point b:"
read(*,*) b
!number of steps
print* , "Enter the number of steps N:"
read(*,*) N

!preliminary input-based calculations
h = (b-a)/N                          ! step size
allocate(y(N+1))                     ! set solution array size
allocate(t(N+1))                     ! set time array size

write(*,*) "The corresponding step size is:"
write(*,*) h

!compute
!need to rewrite the form of the term multiplying h
do  i = 1,N,1
   y(1) = y_0                                        !initial value
   t(1) = a                                          !initial domain point
   t(i+1) = a+i*h                                    !partition of domain
   y(i+1)=y(i)+(h/2)*(f(t(i),y(i))+f(t(i+1),y(i)+h*f(t(i),y(i)))) !successive approximations
end do

!print description of program output 
print* , "===================================================="
print* , "This program approximates the solution to the ODE:  "
print* , "                 dy/dt = f(t,y)                     "
print* , "using the explicit trapezoid integration scheme.    "
print* , "This method is second order accurate (local         "
print* , "truncation error). The approximate solution y_i at  "
print* , "step t_i, for i = 0,1,...,N is:                     " 
print* , "===================================================="
write(*,*) y

!function corresponding to RHS of ODE to be solved
contains
function f(t,y)
  real :: f
  real :: t,y
 f = t+y          !EDIT HERE FOR EACH USE
end function

end program trap

 
