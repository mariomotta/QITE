!----------------------------------------------------------
!Program: Block error estimator
!VERSION: 2-Jan-2014
!NOTICE: Give n data, estimate the error bar from
!different block, the block number is scaned from the 
!all the factors of n.
!TYPE: Serial code
!COMMENT: 
!REFERENCE:subroutine Factors is from 
!http://rosettacode.org/wiki/Factors_of_an_integer#Fortran
!----------------------------------------------------------
program main
implicit none
integer::n
real(kind=8),allocatable::x(:),y(:)
integer::m,pn !pn is part n, m is number of factors
integer,allocatable::fact(:)
integer::i,j
character*90 :: filname
real(kind=8)::mean,err
!read the parameter
print*,"Input the total number of data:"
read(*,*) n
allocate(x(n),y(n))
print*,"Input data file name:"
read(*,*) filname
open(unit=111,file=trim(filname))
do i=1,n,1
   read(111,*) x(i)
end do

!Get the factor number
call factor_num(n,m)
allocate(fact(m))
call factors(n,m,fact)

!write(*,*) "The number of factors:",m
!do i=1,m,1
!   write(*,*) fact(i)
!end do

do i=1,m-1,1
   pn=n/fact(i)
   if(pn*fact(i).NE.n) then
     write(*,*) "Something is wrong with pn and fact(i):",pn,fact(i)
   end if
   do j=1,pn,1
      call avg(fact(i),x((j-1)*fact(i)+1),y(j)) 
   end do

   call err_anal(y(1),pn,mean,err)

   write(*,*) fact(i),mean,err
end do

deallocate(x,y,fact)
stop
end program main


!--------------------------------
!get the number of factors from n
!--------------------------------
subroutine Factor_num(n,m)
implicit none
integer,intent(IN)::n
integer,intent(OUT)::m
integer :: i

m=0
do i = 1, int(sqrt(real(n))) - 1
   !if (mod(n, i) == 0) write (*,*) i, n/i
   if (mod(n, i) == 0) m=m+2
end do

! Check to see if n is a square
i = int(sqrt(real(n)))
if (i*i == n) then
   !write (*,*) i
   m=m+1
else if (mod(n, i) == 0) then
   !write (*,*) i, n/i
   m=m+2
end if

end subroutine Factor_num


!---------------------------------------------
!get all the factors after we get the number m
!---------------------------------------------
subroutine Factors(n,m,fact)
implicit none
integer,intent(IN)::n,m
integer,intent(OUT)::fact(m)
integer ::i,j

j=0
do i = 1, int(sqrt(real(n))) - 1
   if (mod(n, i) == 0) then
      j=j+1
      fact(j)=i
      fact(m+1-j)=n/i
   end if
end do

! Check to see if n is a square
i = int(sqrt(real(n)))
if (i*i == n) then
   j=j+1
   fact(j)=i
else if (mod(n, i) == 0) then
   j=j+1
   fact(j)=i
   fact(m+1-j)=n/i
end if

end subroutine Factors




!---------------------------
!Get the error bar of dat(N)
!---------------------------
subroutine err_anal(dat,N,m,er)
implicit none
integer,intent(IN)::N
real(kind=8),intent(IN)::dat(N)
real(kind=8),intent(OUT)::m,er
integer::i,j,k

if(N.LE.0) then
  write(*,*) "N should not be smaller than or EQ 0", N
  !call mystop
  stop
else if(N.EQ.1) then
  !write(*,*) "N eq 1 warning", N
  m=dat(1)
  er=0.d0
  return
end if

!Get the mean
m=0.d0
do i=1,N,1
   m=m+dat(i)
end do
m=m/dble(N)

er=0.d0
do i=1,N,1
   er=er+dat(i)**2
end do
er=er/dble(N)

er=er-m**2

if(ABS(er).GT.1.d-10) then
  if(er.lt.0.d0) then
    write(*,*) "Something is wrong in err_anal",er
  end if
else
  er=0.d0
end if

er=sqrt(er/dble(N-1))
end subroutine err_anal


!-------------------------------
!get the average value of a(1:n)
!-------------------------------
subroutine avg(n,a,ea)
implicit none
integer,intent(IN)::n
real(kind=8),intent(IN)::a(n)
real(kind=8),intent(OUT)::ea
integer::i
ea=0.d0
do i=1,n,1
   ea=ea+a(i)
end do
ea=ea/dble(n)
end subroutine avg
