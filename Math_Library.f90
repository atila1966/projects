!
!  Math Library
!
      subroutine MATINV(a,ai,blow)
!
      implicit none
      real *8 a(3,3),ai(3,3),deta,blow
!
      deta = a(1,1)*( a(2,2)*a(3,3)-a(2,3)*a(3,2) ) &
           - a(1,2)*( a(2,1)*a(3,3)-a(2,3)*a(3,1) ) &
           + a(1,3)*( a(2,1)*a(3,2)-a(2,2)*a(3,1) )
!
!	write(40,*)"deta=",deta
	if(abs(deta).lt.1000) then
	 blow=1.0
	 goto 100
	endif
!
	deta=1.0/deta
	ai(1,1) = (a(2,2)*a(3,3) - a(2,3)*a(3,2))*deta
	ai(1,2) =-(a(2,1)*a(3,3) - a(3,1)*a(2,3))*deta
	ai(1,3) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))*deta
!	ai(2,1) =-(a(1,2)*a(3,3) - a(1,3)*a(3,2))*deta
	ai(2,2) = (a(1,1)*a(3,3) - a(1,3)*a(3,1))*deta
	ai(2,3) =-(a(1,1)*a(3,2) - a(1,2)*a(3,1))*deta
!	ai(3,1) = (a(1,2)*a(2,3) - a(2,2)*a(1,3))*deta
!	ai(3,2) =-(a(1,1)*a(2,3) - a(2,1)*a(1,3))*deta
	ai(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))*deta
	ai(2,1) = ai(1,2)
	ai(3,1) = ai(1,3)
	ai(3,2) = ai(2,3)
!
 100 return
!
	end subroutine MATINV
!
!
!######################################################################################################################################
!
!
      subroutine invers( n, a, y, np )
      implicit none
      integer n, np, indx(n), i, j
      real *8 a(np,np), y(np,np), d
      do i = 1 , n
         do j = 1 , n
            y(i,j) = 0.
         enddo
         y(i,i) = 1.
      enddo
      call ludcmp(a,n,np,indx,d)
      do j = 1 , n
         call lubksb(a,n,np,indx,y(1,j))
      enddo
      return
      end subroutine invers
!
!
!#################################################################################################################################
!
!
      subroutine solve( n, a, y, np )      
      implicit none
      integer n, np, indx(n), i, j
      real *8 a(np,np), y(np), d
      call ludcmp(a,n,np,indx,d)
      call lubksb(a,n,np,indx,y)
      return
      end subroutine solve
!
!
!#################################################################################################################################
!
!
      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit none
      INTEGER n,np,indx(n)
      REAL *8 d,a(np,np),TINY
      PARAMETER (TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL *8 aamax,dum,sum,vv(n)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if( aamax.eq.0.)then
           print *, ''
           print *, 'E R R O R: singular matrix in ludcmp'
           print *, ''
           stop
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END SUBROUTINE ludcmp
!
!
!#################################################################################################################################
!
!
      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit none
      INTEGER n,np,indx(n)
      REAL *8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL *8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END SUBROUTINE lubksb
!
!
!##########################################################################################################################
!
!
subroutine skew_symmetric(vec,mat)
implicit none
real *8  vec(3), mat(3,3)
!
mat(1,1) = 0.0d0
mat(1,2) = -vec(3)
mat(1,3) = vec(2)
!
mat(2,1) = vec(3)
mat(2,2) = 0.0d0
mat(2,3) = -vec(1)
!
mat(3,1) = -vec(2)
mat(3,2) = vec(1)
mat(3,3) = 0.0d0
!
return
!
end subroutine skew_symmetric
!
!
!##########################################################################################################################
!
!
function sign_check( x )
implicit none
real *8  x, sign_check
if( dabs(x) .lt. 1.0d-9 )then
    sign_check = 0.0d0
elseif( x .gt. 0.0d0 )then
    sign_check = 1.0d0
else
    sign_check = -1.0d0
endif
return
end function sign_check
!
!
!####################################################################################################################################
!
!
function heaviside( x )
implicit none
real *8  x, heaviside
if( x .lt. 0.0d0 )then
    heaviside = 0.0d0
else
    heaviside = 1.0d0
endif
return
end function heaviside
!
!
!####################################################################################################################################
!
!
subroutine isort(n,x,y)
implicit none

integer, intent(in) :: n
integer, intent(inout) :: x(n)
real *8, intent(inout) :: y(n)

integer i, j, itemp
integer jvec(1)
real *8 temp

do i=1,n-1
   jvec=minloc(x(i:n))
   j=jvec(1)+i-1
   if( j .ne. i )then
       itemp=x(i)
       x(i)=x(j)
       x(j)=itemp
       temp=y(i)
       y(i)=y(j)
       y(j)=temp       
   endif
enddo

return

end subroutine isort  
!
!
!#####################################################################################################################################
!
!
subroutine dsort(n,x,y)
implicit none

integer, intent(in) :: n
integer, intent(inout) :: x(n)
real *8, intent(inout) :: y(n)

integer i, j, itemp
integer jvec(1)
real *8 temp

do i=1,n-1
   jvec=minloc(y(i:n))
   j=jvec(1)+i-1
   if( j .ne. i )then
       itemp=x(i)
       x(i)=x(j)
       x(j)=itemp
       temp=y(i)
       y(i)=y(j)
       y(j)=temp       
   endif
enddo

return

end subroutine dsort  
!
!!#####################################################################################################################################