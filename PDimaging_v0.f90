!
!   Peridynamics
!
module dbase
implicit none
integer nmat, nlay, totnode, tot_avail_node, tot_missing_node, ntarget, numTimeSteps
integer output_frequency, restart_frequency, time_history_frequency     

Type Target_Point
     real *8 x, y, z
     integer ID
end type Target_Point
Type(Target_Point) targetPoints(100)    
!
real *8, allocatable :: coord(:,:)
integer, allocatable :: numfam(:)
!
    Type layer
       real *8 dens
       real *8 emod1
       real *8 emod2
       real *8 emod3
       real *8 pratio12
       real *8 pratio21
       real *8 kappa
       real *8 mu
       real *8 g12
       real *8 g13
       real *8 theta
       real *8 thick
       real *8 vol
       real *8 zbot
       real *8 length
       real *8 width
       real *8 hole_radius
       integer ndivx
       integer ndivy
       integer imat
    end type layer
    type(layer) ply(20)
!
    Type material_library
       real *8 dens
       real *8 emod1
       real *8 emod2
       real *8 emod3
       real *8 pratio12
       real *8 pratio21
       real *8 g12
       real *8 g13
    end type material_library
    Type(material_library) matlib(20)
!
    Type ImpactorSection
         real *8 x, y, z
         real *8 vol
         integer mask
    End Type ImpactorSection
!
    Type ImpactorSurface
         real *8 x, y, z, nx, ny, nz
    End Type ImpactorSurface
!
    Type DiffOperator_
         integer n1, n2, n3
         real *8 val, coef
    end type DiffOperator_
!
    Type BoundaryCondition
         integer num_diff_ops
         Type(DiffOperator_), allocatable :: diffs(:)
         real *8 val, xmin, xmax, ymin, ymax, zmin, zmax
    end type BoundaryCondition
!    
    Type PDoperator_
         Type(DiffOperator_) Diffs(100)
         Type(DiffOperator_) out(100)
         integer num_diff_ops, norder, n1order, n2order, n3order, nsize, morder, asymFlag
         integer num_bc, nteqs, ncons, atype, num_out, nwk
         Type(BoundaryCondition), allocatable :: bc(:)
         Type(DiffOperator_), allocatable :: order(:)
         character *80 fname
         integer nx, ny, n
    end type PDoperator_
    Type(PDoperator_) PDoperator
!
    Type nodefam_
         integer , allocatable :: node(:)
         integer i
         real *8 deltax, deltay, deltaz, delta
    end type nodefam_     
!
integer, ALLOCATABLE :: plynum(:,:)
!
real *8, ALLOCATABLE :: disp(:,:)

real *8, ALLOCATABLE :: vel(:,:) 
real *8, ALLOCATABLE :: massvec(:,:)
integer atype
!
real *8 length
real *8 width
real *8 hole_radius
!
Type(nodefam_), allocatable :: family(:)
real *8, allocatable :: DiffAmat3D(:,:), DiffAmatInv3D(:,:), DiffBvec3D(:), DiffAvec(:), plist(:), blist(:), weight(:)
real *8, allocatable :: DiffAmat2D(:,:), DiffAmatInv2D(:,:), DiffBvec2D(:)
real *8, allocatable :: DiffAmat1D(:,:), DiffAmatInv1D(:,:), DiffBvec1D(:)
real *8, allocatable :: rcvec(:), dvolume(:)
real *8, allocatable :: redvec(:), greenvec(:), bluevec(:), alfavec(:)
real *8, allocatable :: redvec_org(:), greenvec_org(:), bluevec_org(:), alfavec_org(:)
real *8, allocatable :: SpSysMat(:), SysMat(:,:), SysVec(:), Coefs(:), sysnorm(:), SpColMax(:)
integer, allocatable :: icol(:), irow(:), mask(:), miss_order(:), avail_order(:)
real *8, allocatable :: xsivec(:)
integer, allocatable :: xsi_order(:), xsi_order_org(:)

end module dbase
!
!
!
module GlobalVariables
!
integer, parameter :: node_above=1, node_below=2, ndof=3, nload=4, ndiv_impactor=30
!
real *8 dx, dy, dz, area
real *8, parameter :: pi=dacos(-1.0d0)
!
integer ilay, ijply
real *8 pratio12
!
real *8 emod1, emod2, emod3, g12, g13, ply_theta, ply_thick, dens
real *8 zloc
integer ndvxy
real *8 scale
character *6 strNum
character *80 restart_file_name
character *50 cskip
real *8 wall1, wall0, twall
end module GlobalVariables
!  
!
!  
program main
use dbase
use GlobalVariables
implicit none
integer i, j, k, ii, n_impactors, ij, iorder
real *8 TIME
real *8 OMP_GET_WTIME
integer START_CLOCK,END_CLOCK,ClOCK_RATE
real *8, allocatable :: matF(:,:,:)
real *8 delta_pass, delta_out
!
!
call system('rd results /s /q' )
call system('mkdir results' )
!
open(23,file='./results/summary.dat')
!
!call checkThreads()
!
call input()
!
    do i = 1 , totnode, 8
       mask(i) = 1
    enddo
!
    call SortAvailableAndMissing()
    !call GenerateNodeFamilies()
    call Derivative_Operator_2D_PASS(0,0)
    !call ResetNodeFamilies()
    !call GenerateNodeFamilies(delta_out,delta_out,delta_out)
    !call Derivative_Operator_2D()
!
close(23)
close(421)
!
end program main
!
!
!
      subroutine input()
      use dbase
      use GlobalVariables
      implicit none
      character line
      integer i,j, k, ii, imat, idiff, ibc
      real *8 radius, psi1, psi2, theta1, theta2
      integer mR, mPsi, mTheta, nPsi, nTheta
      real *8 Lx, Ly, Lz
      integer ndivx, ndivy, mdivx, mdivy
      integer ndivR, ndivT, ndivZ, mdivR, mdivT, mdivZ
      real *8 critical_stretch_fiber, critical_stretch_matrix, critical_stretch_fiber_comp, critical_stretch_matrix_comp
      real *8 critical_stretch_peel, critical_stretch_shear
      integer itarget, iout, ired, igreen, iblue, ialfa
      integer ij
      character *80 fname
!
      open(1,file='PDopr.inp')
!
! 
      read(1,'(A)') line
      read(1,*) PDoperator%fname
      open(2,file=trim(adjustl(PDoperator%fname)))
      read(1,'(A)') line
      read(1,*) PDoperator%n1order, PDoperator%n2order, PDoperator%n 
!
      PDoperator%atype = 0
      PDoperator%morder = 2
      PDoperator%n3order = 0
      PDoperator%asymFlag = 0
      
      read(2,*) PDoperator.nx , PDoperator.ny
      totnode = PDoperator.nx*PDoperator.ny
      allocate(redvec(totnode))
      allocate(greenvec(totnode))
      allocate(bluevec(totnode))
      allocate(alfavec(totnode))
      allocate(redvec_org(totnode))
      allocate(greenvec_org(totnode))
      allocate(bluevec_org(totnode))
      allocate(alfavec_org(totnode))
      allocate(mask(totnode))
      allocate(coord(totnode,3))
      allocate(dvolume(totnode))
      allocate(numfam(totnode))
      allocate(family(totnode))
      do i = 1 , totnode
         mask(i) = 0
      enddo
      ij = 0
      do i = 1 , PDoperator.nx
         do j = 1 , PDoperator.ny
            ij = ij + 1
            read(2,*) ired, igreen, iblue, ialfa
            redvec_org(ij) = ired
            bluevec_org(ij) = iblue
            greenvec_org(ij) = igreen
            alfavec_org(ij) = ialfa
            coord(ij,1) = i
            coord(ij,2) = j
            coord(ij,3) = 0
            dvolume(ij) = 1.0d0
         enddo
      enddo
      close(2)
      call getSize2D(PDoperator%n1order, PDoperator%n2order, PDoperator%nsize)
!
      close(1)
!
      return
      end
!
!
!
      subroutine SortAvailableAndMissing()
      use dbase
      use GlobalVariables
      implicit none
      integer, parameter :: maxFamNodes=10000
      real *8, parameter :: tol=1.0d-5
      integer i, j, k, ii, n, jj, jmem, imiss, iavail
      real *8 ijdist, xdist, ydist, zdist, delta_mag, xyz_mag, deltaxk,deltayk,deltazk
      real *8 deltax, deltay, deltaz, delta
!
      PRINT *,'GenerateNodeFamilies'
!
      n = PDoperator%n
      tot_missing_node = 0
      do k = 1,totnode
         if(mask(k) == 0) then
            tot_missing_node = tot_missing_node + 1
            redvec(k) =-1
            bluevec(k) =-1
            greenvec(k) = -1
            alfavec(k) = -1
         else
            redvec(k) = redvec_org(k)
            bluevec(k) = bluevec_org(k)
            greenvec(k) = greenvec_org(k)
            alfavec(k) = alfavec_org(k)         
         endif
      enddo
      tot_avail_node = totnode - tot_missing_node
      allocate(miss_order(tot_missing_node))
      allocate(avail_order(tot_avail_node))
      allocate(xsivec(tot_avail_node))
      allocate(xsi_order(tot_avail_node))
      allocate(xsi_order_org(tot_avail_node))
!      
      imiss = 0
      iavail = 0
      do k = 1 , totnode
         if(mask(k) == 0) then
            imiss = imiss + 1
            miss_order(imiss) = k
         else
            iavail = iavail + 1
            avail_order(iavail) = k
         endif
      enddo
!
	  return
      end
!
!
!
      subroutine GenerateNodeFamily(k)
      use dbase
      use GlobalVariables
      implicit none
      integer, parameter :: maxFamNodes=10000
      real *8, parameter :: tol=1.0d-5
      integer i, j, k, ii, n, jj, jmem, imiss, iavail
      real *8 ijdist, xdist, ydist, zdist, delta_mag, xyz_mag, deltaxk,deltayk,deltazk
      integer :: nodefam(maxFamNodes)
      real *8 deltax, deltay, deltaz, delta
!
      PRINT *,'GenerateNodeFamilies'
!
      n = PDoperator%n
!      
	     do jj = 1 , tot_avail_node
	        j = avail_order(jj)
	        xsi_order(jj) = j
	        xsivec(jj) = dsqrt( (coord(j,1)-coord(k,1))**2 + (coord(j,2) - coord(k,2))**2 )
	     enddo
	     call dsort(n,tot_avail_node,xsi_order,xsivec)
	     allocate(family(k)%node(n))
	     deltax = 0.0d0
	     deltay = 0.0d0
	     deltaz = 0.0d0
	     delta = 0.0d0
	     do jmem = 1 , n
	        j = xsi_order(jmem)
	        family(k)%node(jmem) = j
	        if( dabs(coord(j,1)-coord(k,1)) > deltax ) deltax = dabs(coord(j,1)-coord(k,1))
	        if( dabs(coord(j,2)-coord(k,2)) > deltay ) deltay = dabs(coord(j,1)-coord(k,1))
	     enddo
         family(k)%deltax = deltax
         family(k)%deltay = deltay
         family(k)%delta = xsivec(n)
         family(k)%i = i
         numfam(k) = n
  !
	  return
      end
!
!
!
      subroutine GenerateNodeFamilies()
      use dbase
      use GlobalVariables
      implicit none
      integer, parameter :: maxFamNodes=10000
      real *8, parameter :: tol=1.0d-5
      integer i, j, k, ii, n, jj, jmem, imiss, iavail
      real *8 ijdist, xdist, ydist, zdist, delta_mag, xyz_mag, deltaxk,deltayk,deltazk
      integer :: nodefam(maxFamNodes)
      real *8 deltax, deltay, deltaz, delta
!
      PRINT *,'GenerateNodeFamilies'
!
      n = PDoperator%n
!      
      do i = 1,tot_missing_node
         k = miss_order(i)
	     do jj = 1 , tot_avail_node
	        j = avail_order(jj)
	        xsi_order(jj) = j
	        xsivec(jj) = dsqrt( (coord(j,1)-coord(k,1))**2 + (coord(j,2) - coord(k,2))**2 )
	     enddo
	     call dsort(n,tot_avail_node,xsi_order,xsivec)
	     allocate(family(k)%node(n))
	     deltax = 0.0d0
	     deltay = 0.0d0
	     deltaz = 0.0d0
	     delta = 0.0d0
	     do jmem = 1 , n
	        j = xsi_order(jmem)
	        family(k)%node(jmem) = j
	        if( dabs(coord(j,1)-coord(k,1)) > deltax ) deltax = dabs(coord(j,1)-coord(k,1))
	        if( dabs(coord(j,2)-coord(k,2)) > deltay ) deltay = dabs(coord(j,1)-coord(k,1))
	     enddo
         family(k)%deltax = deltax
         family(k)%deltay = deltay
         family(k)%delta = xsivec(n)
         family(k)%i = i
         numfam(k) = n
	  enddo
  !
	  return
      end
!
!
!
      subroutine ResetNodeFamilies()
	  use dbase
	  implicit none
	  integer i
!	  
	  do i = 1 , totnode
	     deallocate(family(i)%node)
	  enddo
	  return
	  end
     
!
!////////////////////////////// 2-D ////////////////////////////////////////
!
!
!
!
    subroutine inverse2(nz, mat, invmat, i, rcond)
    implicit none
    real *8 mat(nz,nz), invmat(nz,nz)
    integer, allocatable :: ipiv(:), iwork(:)
    real *8, allocatable :: work(:)
    integer info, kz, k2, k1, k, i, nz
    real *8 anorm, a, rcond
!    
    allocate (ipiv(nz), work(4*nz), iwork(nz))
    call dgetrf( nz, nz, mat, nz, ipiv, info )   
    anorm=0.0d0
    do k2=1, nz
       a=0.0d0
       do k1=1, nz
          a=a+dabs(mat(k1,k2))
       enddo
       anorm=max(anorm,a)
    enddo
    call dgecon( '1', nz, mat, nz, anorm, rcond, work, iwork, info )
    rcond=1.0d0/rcond
    write(23,*) 'rcond(i)=', i, rcond
    if( rcond > 1.0d8 .or. rcond .ne. rcond )then
        print *, 'E R R O R: ill-conditioned for k=', i, 'rcond=', rcond
        !read(*,*)
        !stop
    endif
    call dgetri( nz, mat, nz, ipiv, work, 4*nz, info ) 
    do k=1, nz
       !invmat(k,:) = fact(k) * ds**sum(order(k,:)) * mat(k,:)
       invmat(k,:) = mat(k,:)
    enddo
    deallocate (ipiv)
    deallocate(work)
    deallocate(iwork)
    return
    end
!
!
!
    subroutine Derivative_Operator_2D_PASS(n1,n2)
    use dbase
    use GlobalVariables
    implicit none
    integer n1order, n2order, nsize, k, idiff, n1, n2, n3, m, i, j, jmem, l, morder
    integer ii, jj
    integer rval, gval, bval, aval
    real *8 fval_k, fval_j, xsi1, xsi2, xsi3, gfunval, delta_mag
    real *8  rcond, max_rcond
    real *8 redval_k, greenval_k, blueval_k, alfaval_k
    real *8 redval_j, greenval_j, blueval_j, alfaval_j
!    
    allocate(DiffAmat2D(PDoperator%nsize,PDoperator%nsize))
    allocate(DiffAmatInv2D(PDoperator%nsize,PDoperator%nsize))
    allocate(DiffBvec2D(PDoperator%nsize))
    allocate(DiffAvec(PDoperator%nsize))
    allocate(plist(PDoperator%nsize))
    allocate(blist(PDoperator%nsize))
    allocate(weight(PDoperator%nsize))
    allocate(rcvec(PDoperator%nsize))
!
    n1order = PDoperator%n1order
    n2order = PDoperator%n2order
    nsize = PDoperator%nsize
    morder = 2
    do i = 1 , tot_avail_node
       k = avail_order(i)
       redvec(k) = redvec_org(k)
       bluevec(k) = bluevec_org(k)
       greenvec(k) = greenvec_org(k)
       alfavec(k) = alfavec_org(k)
    enddo
    do i = 1 , tot_missing_node
       k = miss_order(i)
       call GenerateNodeFamily(k)
       call FormDiffAmat2D(n1order, n2order, nsize, k)
       !call recondition(DiffAmat2D,rcvec,nsize)
       !call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)
       call Inverse(nsize,DiffAmat2D,DiffAmatInv2D,nsize)
       !delta_mag = family(k)%delta
       delta_mag = 1.0d0
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do ii = 1 , nsize
             DiffAvec(ii) = 0.0d0
             do jj = 1 , nsize
                DiffAvec(ii) = DiffAvec(ii) + DiffAmatInv2D(ii,jj)*DiffBvec2D(jj)
             enddo
          enddo
          !call recover(DiffAvec,rcvec,nsize)
          redval_k = 0.0d0
          blueval_k = 0.0d0
          greenval_k = 0.0d0
          alfaval_k = 0.0d0
          do jmem = 1 , PDoperator%n
             j = family(k)%node(jmem)
             redval_j = redvec(j)
             blueval_j = bluevec(j)
             greenval_j = greenvec(j)
             alfaval_j = alfavec(j)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             redval_k = redval_k + redval_j*gfunval*dvolume(j)
             blueval_k = blueval_k + blueval_j*gfunval*dvolume(j)
             greenval_k = greenval_k + greenval_j*gfunval*dvolume(j)
             alfaval_k = alfaval_k + alfaval_j*gfunval*dvolume(j)
          enddo
          redval_k = (1.0d0/delta_mag**(n1+n2))*redval_k
          blueval_k = (1.0d0/delta_mag**(n1+n2))*blueval_k
          greenval_k = (1.0d0/delta_mag**(n1+n2))*greenval_k
          alfaval_k = (1.0d0/delta_mag**(n1+n2))*alfaval_k
          redvec(k) = redval_k
          bluevec(k) = blueval_k
          greenvec(k) = greenval_k
          alfavec(k) = alfaval_k
    enddo
!
    do i = 1 , totnode
       rval = idnint(redvec(i))
       gval = idnint(greenvec(i))
       bval = idnint(bluevec(i))
       aval = idnint(alfavec(i))
       write(12,"(5(i6,3x))") rval, gval, bval, aval, mask(i)    
    enddo
!    
    deallocate(DiffAmat2D)
    deallocate(DiffAmatInv2D)
    deallocate(DiffBvec2D)
    deallocate(DiffAvec)
    deallocate(plist)
    deallocate(blist)
    deallocate(weight)
    deallocate(rcvec)
!       
    return
    end
!
!
!
    subroutine Derivative_Operator_2D()
    use dbase
    use GlobalVariables
    implicit none
    integer n1order, n2order, nsize, k, idiff, n1, n2, n3, m, i, j, jmem, l, morder
    real *8 fval_k, fval_j, xsi1, xsi2, xsi3, gfunval, delta_mag
    real *8  rcond, max_rcond
    real *8 redval_k, greenval_k, blueval_k, alfaval_k
    real *8 redval_j, greenval_j, blueval_j, alfaval_j
    real *8 redval, greenval, blueval, alfaval
!    
    allocate(DiffAmat2D(PDoperator%nsize,PDoperator%nsize))
    allocate(DiffAmatInv2D(PDoperator%nsize,PDoperator%nsize))
    allocate(DiffBvec2D(PDoperator%nsize))
    allocate(DiffAvec(PDoperator%nsize))
    allocate(plist(PDoperator%nsize))
    allocate(blist(PDoperator%nsize))
    allocate(weight(PDoperator%nsize))
    allocate(rcvec(PDoperator%nsize))
!
    n1order = PDoperator%n1order
    n2order = PDoperator%n2order
    nsize = PDoperator%nsize
    morder = 2
    open(28,file='./results/derivatives.dat')
    do k = 1 , totnode
       redval_k = redvec(k)
       blueval_k = bluevec(k)
       greenval_k = greenvec(k)
       alfaval_k = alfavec(k)
       call FormDiffAmat2D(n1order, n2order, nsize, k)
       call recondition(DiffAmat2D,rcvec,nsize)
       !call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)
       call Inverse(nsize,DiffAmat2D,DiffAmatInv2D,nsize)
       !delta_mag = family(k).delta
       delta_mag = 1.0d0
       do idiff = 1 , PDoperator%num_out
          n1 = PDoperator%out(idiff)%n1
          n2 = PDoperator%out(idiff)%n2
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          redval = 0.0d0
          greenval = 0.0d0
          blueval = 0.0d0
          alfaval = 0.0d0
          do jmem = 1 , numfam(k)
             j = family(k)%node(jmem)
             redval_j = redvec(j)
             blueval_j = bluevec(j)
             greenval_j = greenvec(j)
             alfaval_j = alfavec(j)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             redval = redval + redval_j*gfunval*dvolume(j)
             greenval = greenval + greenval_j*gfunval*dvolume(j)
             blueval = blueval + blueval_j*gfunval*dvolume(j)
             alfaval = alfaval + redval_j*gfunval*dvolume(j)
          enddo
          PDoperator%out(idiff)%val = (1.0d0/delta_mag**(n1+n2))*PDoperator%out(idiff)%val
       enddo
       write(28,"(20(e16.9,2x))") (coord(k,j), j=1,3), fval_k, (PDoperator%out(idiff)%val,idiff=1,PDoperator%num_out)
    enddo
    close(28)
!    
    deallocate(DiffAmat2D)
    deallocate(DiffAmatInv2D)
    deallocate(DiffBvec2D)
    deallocate(DiffAvec)
    deallocate(plist)
    deallocate(blist)
    deallocate(weight)
    deallocate(rcvec)
!       
    return
    end
!
!
!
    subroutine FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
    use dbase
    implicit none
    integer morder, i, nsize, n1order, n2order, n1, n2, m
!    
    morder = 2
    call b_operator_2d( n1order, n2order, nsize, n1, n2, m )    
    do i = 1 , nsize
       DiffBvec2D(i) = blist(i)
    enddo
    return
    end
!
!
!    
    subroutine FormDiffAmat2D(n1order, n2order, nsize, k )
    use dbase
    use GlobalVariables
    implicit none
    integer morder, n1order, n2order, nsize, k, jmem, j, ii, jj
    real *8 xsi1, xsi2, delta_mag, tol
!    
    morder = 2
    !delta_mag = family(k)%delta
    delta_mag = 1.0d0
!
    do ii = 1 , nsize
       do jj = 1 , nsize
          DiffAmat2D(ii,jj) = 0.0d0
       enddo
    enddo
    
    do jmem = 1 , numfam(k)
       j = family(k)%node(jmem)
       xsi1 = coord(j,1) - coord(k,1)
       xsi2 = coord(j,2) - coord(k,2)
       call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
       call weights_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
       do ii = 1 , nsize 
          do jj = 1 , nsize
             DiffAmat2D(ii,jj) = DiffAmat2D(ii,jj) + weight(jj)*plist(ii)*plist(jj)*dvolume(j)
          enddo
       enddo
    enddo
    return
    end
!
!
!
    subroutine getSize2D(n1order, n2order, nsize)
    use dbase
    implicit none
    integer n1order, n2order, nsize, iterm
!        
        iterm = 0
        if( n1order >=0 .and. n2order >= 0 ) then
            iterm = iterm + 1
        endif
!        
        if( n1order >= 1 ) then
            iterm = iterm + 1
        endif
        if( n2order >= 1 ) then
            iterm = iterm + 1
        endif
!        
        if( n1order >= 2 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 2 .and. n2order >= 2 ) then
            iterm = iterm + 1
        endif
        if( n2order >= 2 ) then
            iterm = iterm + 1
        endif
!        
        if( n1order >= 3 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 3 .and. n2order >= 3) then
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 3 ) then
            iterm = iterm + 1
        endif
!
        if( n1order >= 4 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 4 .and. n2order >= 4 ) then
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 4 ) then
            iterm = iterm + 1
        endif
!
        if( n1order >= 5 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 5 .and. n1order >= 5) then
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 5 ) then
            iterm = iterm + 1
        endif
!
        if( n1order >= 6 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 6 .and. n2order >= 6 ) then
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 6 ) then
            iterm = iterm + 1
        endif
        nsize = iterm
    return
    end        
!
!
!
    subroutine SetDiffOperators2D(n1order, n2order, nsize)
    use dbase
    implicit none
    integer n1order, n2order, nsize, iterm
!        
    allocate(PDoperator.order(nsize))
        iterm = 0
        if( PDoperator%atype==0 .and. n1order >=0 .and. n2order >= 0 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
!        
        if( n1order >= 1 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 1 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
        endif
!        
        if( n1order >= 2 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 2 .and. n2order >= 2 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 2 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
        endif
!        
        if( n1order >= 3 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 3 .and. n2order >= 3) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 3 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=3
            PDoperator.order(iterm).n3=0
        endif
!
        if( n1order >= 4 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=4
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 4 .and. n2order >= 4 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 4 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=4
            PDoperator.order(iterm).n3=0
        endif
!
        if( n1order >= 5 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=5
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 5 .and. n1order >= 5) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=4
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=3
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=4
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 5 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=5
            PDoperator.order(iterm).n3=0
        endif
!
        if( n1order >= 6 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=6
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 6 .and. n2order >= 6 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=5
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=4
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=3
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=4
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=5
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 6 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=6
            PDoperator.order(iterm).n3=0
        endif
        nsize = iterm
    return
    end        
!
!
!
    subroutine p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
    use dbase
    implicit none
    integer morder, n1order, n2order, nsize, iterm
    real *8 xsi1, xsi2, xsi1p, xsi2p, delta_mag
!    
    xsi1p = xsi1/delta_mag
    xsi2p = xsi2/delta_mag
    iterm = 0
    if( n1order >= 0 .and. n2order >= 0 ) then
        iterm = iterm + 1
        plist(iterm) = 1
    endif
!    
    if( n1order >= 1 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p
    endif
    if( n2order >= 1 ) then
        iterm = iterm + 1
        plist(iterm) = xsi2p
    endif
!    
    if( n1order >= 2 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi1p
    endif
    if( n1order >= 2 .and. n2order >= 2 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi2p
    endif
    if( n2order >= 2 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi2p*xsi2p
    endif
!    
    if( n1order >= 3 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi1p*xsi1p
    endif
    if( n1order >= 3 .and. n2order >= 3 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi2p*xsi2p
    endif
    if( n2order >= 3 ) then
        iterm = iterm + 1
        plist(iterm) = xsi2p*xsi2p*xsi2p
    endif
!
    if( n1order >= 4 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p
    endif
    if( n1order >= 4 .and. n2order >= 4 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi2p*xsi2p*xsi2p
    endif
    if( n2order >= 4 ) then
        iterm = iterm + 1
        plist(iterm) = xsi2p*xsi2p*xsi2p*xsi2p
    endif
!
    if( n1order >= 5 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi1p
    endif
    if( n1order >= 5 .and. n2order >= 5 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi2p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
    endif
    if( n2order >= 5 ) then
        plist(iterm) = xsi2p*xsi2p*xsi2p*xsi2p*xsi2p
    endif
!    
    if( n1order >= 6 ) then
        iterm = iterm + 1
    endif
    if( n1order >= 6 .and. n2order >= 6 ) then
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi1p*xsi1p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi2p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi2p*xsi2p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
    endif
    if( n1order >= 6 ) then
        plist(iterm) = xsi2p*xsi2p*xsi2p*xsi2p*xsi2p*xsi2p
    endif
!    
    return
    end
!
!
!
    subroutine weights_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
    use dbase
    implicit none
    integer morder, n1order, n2order, nsize, iterm
    real *8 xsi1, xsi2, xsi_mag, wt, delta_mag
!    
    morder = 2
    xsi_mag = dsqrt(xsi1*xsi1+xsi2*xsi2)
!
    !wt = (delta_mag/xsi_mag)**2
    !wt = exp(-4*(xsi_mag/delta_mag)**2)
    wt = 1.0d0
    !print *, 'wt=',wt
    
    iterm = 0
    if( n1order >= 0 .and. n2order >= 0 ) then
        !wt = 1.0d0
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!
    if( n1order >= 1 ) then
        !wt = (delta_mag/xsi_mag)**2
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 1 ) then
        !wt = (delta_mag/xsi_mag)**2
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    if( n1order >= 2 ) then
        !wt = (delta_mag/xsi_mag)**3
        iterm = iterm + 1
        weight(iterm)  = wt
    endif
    if( n1order >= 2 .and. n2order >= 2 ) then
        !wt = (delta_mag/xsi_mag)**3
        iterm = iterm + 1
        weight(iterm)  = wt
    endif
    if( n2order >= 2 ) then
        !wt = (delta_mag/xsi_mag)**3
        iterm = iterm + 1
        weight(iterm)  = wt
    endif
!
    if( n1order >= 3 ) then
!        wt = (delta_mag/xsi_mag)**4
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 3 .and. n2order >= 3 ) then
!        wt = (delta_mag/xsi_mag)**4
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 3 ) then
!        wt = (delta_mag/xsi_mag)**4
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    if( n1order >= 4 ) then
!        wt = (delta_mag/xsi_mag)**5
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 4 .and. n2order >= 4) then
!        wt = (delta_mag/xsi_mag)**5
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 4 ) then
!        wt = (delta_mag/xsi_mag)**5
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!
    if( n1order >= 5 ) then
!        wt = (delta_mag/xsi_mag)**6
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 5 .and. n2order >= 5 ) then
!        wt = (delta_mag/xsi_mag)**6
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 5 ) then
!        wt = (delta_mag/xsi_mag)**6
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    if( n1order >= 6 ) then
!        wt = (delta_mag/xsi_mag)**7
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 6 .and. n2order >= 6 ) then
!        wt = (delta_mag/xsi_mag)**7
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 6 ) then
!        wt = (delta_mag/xsi_mag)**7
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    return
    end
!
!
!    
    subroutine b_operator_2d( n1order, n2order, nsize, n1, n2, m )
    use dbase
    implicit none
    integer n1, n2, m, morder, n1order, n2order, nsize, iterm, i
    real *8 fn1, fn2, coef
!    
    morder = 2
    if( n1 > n1order .or. n2 > n2order ) stop
    do i = 1 , nsize
       blist(i) = 0.0d0
    enddo
    fn1=1
    fn2=1
    do i = 1 , n1
       fn1 = fn1*i
    enddo
    do i = 1 , n2
       fn2 = fn2*i
    enddo
    coef = fn1*fn2
        iterm = 0
        if( n1order >=0 .and. n2order >= 0 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
!        
        if( n1order >= 1 ) then
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n2order >= 1 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 1 ) then 
                m = iterm
            endif
        endif
!        
        if( n1order >= 2 ) then
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 2 .and. n2order >= 2 ) then
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 1 ) then
                m = iterm
            endif
        endif
        if( n2order >= 2 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 2 ) then
                m = iterm
            endif
        endif
!        
        if( n1order >= 3 ) then
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 3 .and. n2order >= 3) then
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 2 ) then
                m = iterm
            endif
        endif
        if( n2order >= 3 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 3 ) then
                m = iterm
            endif
        endif
!
        if( n1order >= 4 ) then
            iterm = iterm + 1
            if( n1 == 4 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 4 .and. n2order >= 4 ) then
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 2 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 3 ) then
                m = iterm
            endif
        endif
        if( n2order >= 4 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 4 ) then
                m = 15
            endif
        endif
!
        if( n1order >= 5 ) then
            iterm = iterm + 1
            if( n1 == 5 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 5 .and. n1order >= 5) then
            iterm = iterm + 1
            if( n1 == 4 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 2 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 3 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 4 ) then
                m = iterm
            endif
        endif
        if( n2order >= 5 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 5 ) then
                m = iterm
            endif
        endif
!
        if( n1order >= 6 ) then
            iterm = iterm + 1
            if( n1 == 6 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 6 .and. n2order >= 6 ) then
            iterm = iterm + 1
            if( n1 == 5 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 4 .and. n2 == 2 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 3 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 4 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 5 ) then
                m = iterm
            endif
        endif
        if( n2order >= 6 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 6 ) then
               m = iterm
            endif
        endif
    blist(m) = coef                                       
    return                                      
    end
!
!
!////////////////////////////////UTILITY SUBROUTINES /////////////////////////
!
!
!    subroutine checkThreads()
!    implicit none
!    integer nthreads, thread_id, omp_get_thread_num, omp_get_num_threads    
! $OMP parallel private(nthreads, thread_id)
!      thread_id = omp_get_thread_num()
!      write(23,*) "Thread ", thread_id, " says: Hello World "
!
!      if (thread_id == 0) then
!         nthreads = omp_get_num_threads()
!        write(23,*), "Thread", thread_id,  "reports: the number of threads are ",  nthreads 
!      endif
! $OMP end parallel
!    return
!    end	subroutine checkThreads       	  	           
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
end
!
!
!
subroutine dsort(n,nt,x,y)
implicit none

integer, intent(in) :: n, nt
integer, intent(inout) :: x(nt)
real *8, intent(inout) :: y(nt)
integer m

integer i, j, itemp
integer jvec(1)
real *8 temp
!
m = n
if( n>=nt ) m = nt - 1
!
do i=1,m
   jvec=minloc(y(i:nt))
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
end  
!
!
!
function angle( dx, dy )
implicit none
real *8 dx, dy, dr, angle
dr = dsqrt(dx**2 + dy**2)
if( dr == 0.0d0 ) then
    angle = 0.0d0
else if( dy .ge. 0.0d0) then
    angle = dacos(dx/dr)
else
   angle =-dacos(dx/dr)
endif
return
end function angle
!
!
!
      subroutine PrintPreview()
      use dbase
      implicit none
      integer i, j
      real *8 xglobal, yglobal, zglobal
!
      open(82,file='./output/preview.dat')      
!
      write(82,"(I6)") totnode
      do i = 1 , totnode
         write(82,"(3(e16.9,3x),2x,2(i6,2x))") coord(i,1), coord(i,2), coord(i,3), plynum(i,1), ply(plynum(i,1))%imat
      enddo
      close(82)      
!
      return
      end
!
!
!          
      subroutine copy_avec_to_bvec( bvec, aval, avec, n )
      implicit none
      real *8 aval, avec(*), bvec(*), factor
      integer i, n
      do i = 1 , n
         bvec(i) = aval*avec(i)
      enddo
      return
      end
!
!
!          
      subroutine add_avec_to_bvec( bval, bvec, aval, avec, n )
      implicit none
      real *8 avec(*), bvec(*), bval, aval
      integer i, n
      do i = 1 , n
         bvec(i) = bval*bvec(i) + aval*avec(i)
      enddo
      return
      end
!
!
!          
      subroutine scaleAvec( sval, avec, n )
      implicit none
      real *8 avec(*), aval, sval
      integer i, n
      do i = 1 , n
         avec(i) = sval*avec(i)
      enddo
      return
      end
!
!
!          
      subroutine dotProduct( abval, avec, bvec, n )
      implicit none
      real *8 avec(*), bvec(*), abval
      integer i, n
      abval = 0.0d0
      do i = 1 , n
         abval = abval + avec(i)*bvec(i)
      enddo
      return
      end
!
!
!          
      subroutine getNorm( avec, anorm, n )
      implicit none
      real *8 avec(*), anorm, anprm
      integer i, n
      anorm = 0.0d0
      do i = 1 , n
         anprm = anorm + avec(i)**2
      enddo
      anorm = dsqrt(anorm)
      return
      end
!
!
!      
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
!
      subroutine inverse( n, a, y, np )
      implicit none
      integer nmax
      parameter (nmax=5000)
      integer np, indx(NMAX), n, i, j
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
      end subroutine inverse
!
!
!
      subroutine solve( n, a, y, np )      
      implicit none
      integer nmax
      parameter (nmax=5000)
      integer np, indx(NMAX), n, i, j
      real *8 a(np,np), y(np), d
      call ludcmp(a,n,np,indx,d)
      call lubksb(a,n,np,indx,y)
      return
      end subroutine solve
!
!
!
      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit none
      integer nmax
      parameter(nmax=5000)
      INTEGER n,np,indx(n)
      REAL *8 d,a(np,np),TINY
      PARAMETER (TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL *8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.) then
            print *, 'singular matrix in ludcmp'
            read(*,*)
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
!
      subroutine recondition( amat, rcvec, n )
      implicit none
      integer n, i, j
      real *8 amat(n,n), rcvec(n), cmax
      do j = 1 , n
         cmax = dabs(amat(1,j))
         do i = 1 , n
            if(dabs(amat(i,j)) > cmax) cmax = dabs(amat(i,j))
         enddo
         rcvec(j) = cmax
         do i = 1 , n
            amat(i,j) = amat(i,j)/cmax
         enddo
      enddo
!         
      return
      end
!
!
!
      subroutine recover( avec, rcvec, n )
      implicit none
      integer n, i
      real *8 avec(n), rcvec(n)
      do i = 1 , n
         avec(i) = avec(i)/rcvec(i)
      enddo
!         
      return
      end
!
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
!           