      program shktub
c
c*******************************************************
c     Shock Tube Problem: 1D Euler eqs solver
c********************************************************
      common/ grid1d/ mx,dx,x(-5:505)
      common/ cmpcnd/ cfl,dt,nlast,time,eps,g,rhol,ul,pl,rhor,ur,pr
     &              ,rags,cvgas,ecp,n
      common/ solute/ rho(-5:505),u(-5:505),p(-5:505),e(-5:505)
      common/ solver/ q(3,-5:505),qold(3,-5:505),flux(3,-5:505)
      
      integer, parameter :: nn=3
      double precision :: a(nn,nn), xx(nn,nn)
      double precision, parameter:: abserr=1.0e-09
      integer ii, jj
      
c
cc    (1)set grid
      mx   =  401
      xmin = -2.0
      xmax =  2.0
      dx   = (xmax-xmin)/float(mx-1)
      do 100 i=-2,mx+3
        x(i) = xmin+dx*float(i-1)
  100 continue
c     (2)set parameters
      atmpa= 1.013e05
      eps  = 1.0e-06
      g    = 1.4
      rgas =8.314e+03/28.96
      cvgas=rgas/(g-1.0)
c     (2-1) Left chamber
c     u(m/s),T(K),rho(kg/m3),P(atm)
      ul=0
      rhol=1.0
      platm=1.0
      pl=platm*atmpa
      tl=pl/(rgas*rhol)
c     (2-2) Right chamber
      ur=0
      rhor=0.125
      pratm=0.1
      pr=pratm*atmpa
      tr=pr/(rgas*rhor)
c     (2-3) Timestep settings
      tmsec=2.
      time= tmsec*1.0e-03
      cfl=0.1
c     (2-4) entropy
      ecp= 0.125
c
c     (3)numerical solutions
      uref=sqrt(g*rgas*amax1(tl,tr))
      dt   = cfl*dx/abs(uref)
      nlast= int(time/dt)
      time = dt*float(nlast)
c     (4) initialize properties
      do 300 i=1,mx
        if(x(i).lt.0.0) then
c     (4-1) Left chamber
          rho(i)= rhol
            u(i)=   ul
            p(i)=   pl
            e(i)= p(i)/((g-1.0)*rho(i))
        else
c     (4-2) Right chamber
          rho(i)= rhor
            u(i)=   ur
            p(i)=   pr
            e(i)= p(i)/((g-1.0)*rho(i))
        end if
c     (4-3) conserved quantitiy
        q(1,i)= rho(i)
        q(2,i)= rho(i)*u(i)
        q(3,i)= rho(i)*(e(i)+0.5*u(i)**2)
  300 continue
cc
      open(unit=60,file='sw.plt',form='formatted')
      write(60,*) 'VARIABLES="x, m","t, ms","rho, kg/m3","u, m/s",',
     &            '"p, Pa","t, K"'
      write(60,*) 'ZONE T="SW",I=',mx,',J=',nlast,
     &            ',F=POINT'
      do 410 i=1,mx
         write(60,9001) x(i),0.,rho(i),u(i),p(i),e(i)/cvgas
 410  continue
c     (5) Main loop: Time marching
      do 1000 n=1,nlast
c
        do 1010 i=1,mx
          qold(1,i)= q(1,i)
          qold(2,i)= q(2,i)
          qold(3,i)= q(3,i)
 1010   continue
c
        call calflx
        
        do 1100 i=4,mx-3
          q(1,i)= qold(1,i)-(dt/dx)*(flux(1,i)-flux(1,i-1))
          q(2,i)= qold(2,i)-(dt/dx)*(flux(2,i)-flux(2,i-1))
          q(3,i)= qold(3,i)-(dt/dx)*(flux(3,i)-flux(3,i-1))
          rho(i)= q(1,i)
            u(i)= q(2,i)/rho(i)
            p(i)= (g-1.0)*(q(3,i)-0.5*rho(i)*u(i)**2)
            e(i)= p(i)/((g-1.0)*rho(i))
 1100   continue
c
        if(mod(n,10).eq.0) write(6,*) ' .... step = ',n
        do 9010 i=1,mx
           write(60,9001) x(i),dt*n*1.e3,rho(i),u(i),p(i),e(i)/cvgas
 9010   continue
 9001   format(6f15.5)
 1000 continue
      close (unit=60)
c     waveform at t=2ms
      open(unit=70,file='waveform.txt',form='formatted')
      do 9200 i=1,mx
         write(70,9201) x(i),rho(i),u(i),p(i),e(i)/cvgas
 9200 continue
 9201 format(5e15.5)
      close (unit=70)
c
      write(6,*) 'DONE!'
c
      stop
      end
c
      subroutine calflx
c========================================================================
c     Compute Flux with Yee's Symmetric TVD Scheme
c========================================================================
      common/ grid1d/ mx,dx,x(-5:505)
      common/ cmpcnd/ cfl,dt,nlast,time,eps,g,rhol,ul,pl,rhor,ur,pr
     &              ,rags,cvgas,ecp,n
      common/ solute/ rho(-5:505),u(-5:505),p(-5:505),e(-5:505)
      common/ solver/ q(3,-5:505),qold(3,-5:505),flux(3,-5:505)
      dimension eigl(3,-5:505),rmat(3,3,-5:505),alpha(3,-5:505)
      integer, parameter :: nn=3
      double precision :: a(nn,nn), xx(nn,nn)
      double precision, parameter:: abserr=1.0e-011
      integer ii, jj
      DOUBLE PRECISION, DIMENSION(3,3) :: MAT, MATINV
      LOGICAL :: OK_FLAG      
      
c
c     Find L/R values
      do 100 i=1,mx-1
        q1l= q(1,i  )
        q1r= q(1,i+1)
        q2l= q(2,i  )
        q2r= q(2,i+1)
        q3l= q(3,i  )
        q3r= q(3,i+1)
        
c
        dq1= q(1,i+1)-q(1,i  )
        dq2= q(2,i+1)-q(2,i  )
        dq3= q(3,i+1)-q(3,i  )
        
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq1',dq1
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq2',dq2
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq3',dq3
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq4',dq4
c
        rl= q1l
        ul= q2l/rl
        pl=(g-1.0)*(q3l-.5*rl*ul**2)
        hl=(q3l+pl)/rl
        rr= q1r
        ur= q2r/rr
        pr=(g-1.0)*(q3r-.5*rr*ur**2)
        hr=(q3r+pr)/rr
c     Roe Average at i+1/2
        ubar= (sqrt(rl)*ul+sqrt(rr)*ur)/(sqrt(rl)+sqrt(rr))
        hbar= (sqrt(rl)*hl+sqrt(rr)*hr)/(sqrt(rl)+sqrt(rr))
        abar= sqrt((g-1.0)*(hbar-0.5*ubar**2))
        abar= sqrt(amax1( abar, amin1(g*pl/rl, g*pr/rr) ) )
c     Eigen values and corresponding eigen vectors:

! matrix A
      a(1,1)=0.0
      a(1,2)=1.0
      a(1,3)=0.0
      
      a(2,1)=(g-3.0)*0.5*ubar**2
      a(2,2)=(3.0-g)*ubar
      a(2,3)=g-1.0
      
      a(3,1)=(((g-1.0)*0.5*ubar**2)-hbar)*ubar
      a(3,2)=hbar-((g-1.0)*ubar**2)
      a(3,3)=g*ubar
      
! print a header and the original matrix      
      if(n.eq.1 .and. i.eq.201)  write (*,200)
      
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a11', a(1,1)
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a12', a(1,2)
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a13', a(1,3)
      
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a21', a(2,1)
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a22', a(2,2)
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a23', a(2,3)
      
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a31', a(3,1)
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a32', a(3,2)
      if(n.eq.1 .and. i.eq.201)  write (*,*)'a33', a(3,3)
      

      call Jacobi(a,xx,abserr,nn)
      
      ! print solutions
      if(n.eq.1 .and. i.eq.201)  write (*,202)
      if(n.eq.1 .and. i.eq.201)  write (*,*) (a(ii,ii),ii=1,nn)
      if(n.eq.1 .and. i.eq.201)  write (*,203)
      do ii = 1,nn
         if(n.eq.1 .and. i.eq.201)  write (*,*)  (xx(ii,jj),jj=1,nn)
      end do
      
      xxx=1
      do ii = 1,nn
         xxx= xxx+1
         if(n.eq.1 .and. i.eq.201)  write (*,*)  (xxx,jj=1,nn)
      end do
      
      
200   format (' Eigenvalues and eigenvectors (Jacobi method) ',/, 
     &      ' Matrix A')
201   format (6f12.6)
202   format (/,' Eigenvalues')
203   format (/,' Eigenvectors')

c          Matrix R&Lambda at i+1/2
          eigl(1,i)= a(3,3)
        rmat(1,1,i)= xx(1,3)
        rmat(2,1,i)= xx(2,3)
        rmat(3,1,i)= xx(3,3)
          eigl(2,i)= a(1,1)
        rmat(1,2,i)= xx(1,1)
        rmat(2,2,i)= xx(2,1)
        rmat(3,2,i)= xx(3,1)
          eigl(3,i)= a(2,2)
        rmat(1,3,i)= xx(1,2)
        rmat(2,3,i)= xx(2,2)
        rmat(3,3,i)= xx(3,2)

c          eigl(1,i)= ubar-abar
c        rmat(1,1,i)= 1.0
c        rmat(2,1,i)= ubar-abar
c        rmat(3,1,i)= hbar-ubar*abar
c          eigl(2,i)= ubar
c        rmat(1,2,i)= 1.0
c        rmat(2,2,i)= ubar
c        rmat(3,2,i)= 0.5*ubar**2
c          eigl(3,i)= ubar+abar
c        rmat(1,3,i)= 1.0
c        rmat(2,3,i)= ubar+abar
c        rmat(3,3,i)= hbar+ubar*abar
       
       
       if(n.eq.1 .and. i.eq.201)  write (*,*) 'eigl1',eigl(1,i)
       if(n.eq.1 .and. i.eq.201)  write (*,*) 'eigl2',eigl(2,i)
       if(n.eq.1 .and. i.eq.201)  write (*,*) 'eigl3',eigl(3,i)
       
       if(n.eq.1 .and. i.eq.201)  write (*,*) 'rmat(1,3,i)',rmat(1,3,i)
       if(n.eq.1 .and. i.eq.201)  write (*,*) 'rmat(2,3,i)',rmat(2,3,i)
       if(n.eq.1 .and. i.eq.201)  write (*,*) 'rmat(3,3,i)',rmat(3,3,i)
       
      MAT(1,1)=rmat(1,3,i)
      MAT(2,1)=rmat(2,3,i)
      MAT(3,1)=rmat(3,3,i)
      
      MAT(1,2)=rmat(1,1,i)
      MAT(2,2)=rmat(2,1,i)
      MAT(3,2)=rmat(3,1,i)
      
      MAT(1,3)=rmat(1,2,i)
      MAT(2,3)=rmat(2,2,i)
      MAT(3,3)=rmat(3,2,i)

!
!     Invert the input matrix.
!

      CALL M33INV (MAT, MATINV, OK_FLAG)



c     alpha = R^(-1) . dQ_(j+1/2) = R^(-1) . (Q_(j+1)-Q_j)
        bb= (g-1.0)/abar**2
        aa= bb*ubar**2/2.0
        alpha(1,i)=
     1              .5*(aa+ubar/abar)    *dq1
     2             -.5*(bb*ubar+1./abar) *dq2
     3             +.5*bb                *dq3
        alpha(2,i)=
     1              (1.-aa)  *dq1
     2              +bb*ubar *dq2
     3              -bb      *dq3
        alpha(3,i)= 
     1              .5*(aa-ubar/abar)    *dq1
     2             -.5*(bb*ubar-1./abar) *dq2
     3             +.5*bb                *dq3
     
c          alpha(1,i)=MATINV(1,1)*dq1+MATINV(1,2)*dq2+MATINV(1,3)*dq3
c          alpha(2,i)=MATINV(2,1)*dq1+MATINV(2,2)*dq2+MATINV(2,3)*dq3
c          alpha(3,i)=MATINV(3,3)*dq1+MATINV(3,2)*dq2+MATINV(3,3)*dq3

       if(n.eq.1 .and. i.eq.201) write(6,*) 'alpha(1,i)',alpha(1,i)
       if(n.eq.1 .and. i.eq.201) write(6,*) 'alpha(2,i)',alpha(2,i)
       if(n.eq.1 .and. i.eq.201) write(6,*) 'alpha(3,i)',alpha(3,i)
       if(n.eq.1 .and. i.eq.201) write(6,*) 'alpha(4,i)',0.0
       
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq1',dq1
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq2',dq2
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq3',dq3
       if(n.eq.1 .and. i.eq.201) write(6,*) 'dq4',dq4

c
  100 continue
c     compute limiter func. and flux
      do 400 i=3,mx-3
        corr=ecp*amax1(abs(eigl(1,i)), abs(eigl(2,i)), abs(eigl(3,i)) )
c     : Lambda = u-c
        abc= abs(eigl(1,i))
        if(abc.lt.corr) abc= (abc**2+corr**2)*.5/corr
        sgn = sign(1., alpha(1,i))
        qlim= sgn*amax1(0.,amin1(sgn*2.*alpha(1,i-1),sgn*2.*alpha(1,i)
     &       ,sgn*2.*alpha(1,i+1),sgn*.5*(alpha(1,i-1)+alpha(1,i+1)) ))
        ph1= -(dt/dx)*eigl(1,i)**2*qlim-abc*(alpha(1,i)-qlim)
c     : Lambda = c
        abc= abs(eigl(2,i))
        if(abc.lt.corr) abc= (abc**2+corr**2)*.5/corr
        sgn = sign(1., alpha(2,i))
        qlim= sgn*amax1(0.,amin1(sgn*2.*alpha(2,i-1),sgn*2.*alpha(2,i)
     &       ,sgn*2.*alpha(2,i+1),sgn*.5*(alpha(2,i-1)+alpha(2,i+1)) ))
        ph2= -(dt/dx)*eigl(2,i)**2*qlim-abc*(alpha(2,i)-qlim)
c     : Lambda = u+c
        abc= abs(eigl(3,i))
        if(abc.lt.corr) abc= (abc**2+corr**2)*.5/corr
        sgn = sign(1., alpha(3,i))
        qlim= sgn*amax1(0.,amin1(sgn*2.*alpha(3,i-1),sgn*2.*alpha(3,i)
     &       ,sgn*2.*alpha(3,i+1),sgn*.5*(alpha(3,i-1)+alpha(3,i+1)) ))
        ph3= -(dt/dx)*eigl(3,i)**2*qlim-abc*(alpha(3,i)-qlim)
c     Phi_(i+1/2) = R . phi
        rphi1= rmat(1,1,i)*ph1+rmat(1,2,i)*ph2+rmat(1,3,i)*ph3
        rphi2= rmat(2,1,i)*ph1+rmat(2,2,i)*ph2+rmat(2,3,i)*ph3
        rphi3= rmat(3,1,i)*ph1+rmat(3,2,i)*ph2+rmat(3,3,i)*ph3
        q1l= q(1,i  )
        q1r= q(1,i+1)
        q2l= q(2,i  )
        q2r= q(2,i+1)
        q3l= q(3,i  )
        q3r= q(3,i+1)
         rl= q1l
         ul= q2l/rl
         pl= (g-1.0)*(q3l-.5*rl*ul**2)
         hl= (q3l+pl)/rl
         rr= q1r
         ur= q2r/rr
         pr= (g-1.0)*(q3r-.5*rr*ur**2)
         hr= (q3r+pr)/rr
        e1l= rl*ul
        e2l= rl*ul**2+pl
        e3l= ul*(q3l+pl)
        e1r= rr*ur
        e2r= rr*ur**2+pr
        e3r= ur*(q3r+pr)
        flux(1,i)= .5*(e1l+e1r+rphi1)
        flux(2,i)= .5*(e2l+e2r+rphi2)
        flux(3,i)= .5*(e3l+e3r+rphi3)
        
       if(n.eq.1 .and. i.eq.201) write(6,*) 'flux(1,i)',flux(1,i)
       if(n.eq.1 .and. i.eq.201) write(6,*) 'flux(2,i)',flux(2,i)
       if(n.eq.1 .and. i.eq.201) write(6,*) 'flux(3,i)',flux(3,i)

c
  400 continue
c
      return
      end
      
      
      subroutine Jacobi(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================

      integer i, j, k, n
      double precision a(n,n),x(n,n)
      double precision abserr, b2, bar
      double precision beta, coeff, c, s, cs, sc
      
! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
      x = 0.0
      do i=1,n
        x(i,i) = 1.0
      end do
      
! find the sum of all off-diagonal elements (squared)
      b2 = 0.0
      do i=1,n
        do j=1,n
          if (i.ne.j) b2 = b2 + a(i,j)**2
        end do
      end do      
      
      if (b2 <= abserr) return
      
! average for off-diagonal elements /2
      bar = 0.5*b2/float(n*n)
      
      do while (b2.gt.abserr)
        do i=1,n-1
          do j=i+1,n
            if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
            b2 = b2 - 2.0*a(j,i)**2
            bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
            beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
            coeff = 0.5*beta/sqrt(1.0+beta**2)
            s = sqrt(max(0.5+coeff,0.0))
            c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
            do k=1,n
              cs =  c*a(i,k)+s*a(j,k)
              sc = -s*a(i,k)+c*a(j,k)
              a(i,k) = cs
              a(j,k) = sc
            end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
            do k=1,n
              cs =  c*a(k,i)+s*a(k,j)
              sc = -s*a(k,i)+c*a(k,j)
              a(k,i) = cs
              a(k,j) = sc
              cs =  c*x(k,i)+s*x(k,j)
              sc = -s*x(k,i)+c*x(k,j)
              x(k,i) = cs
              x(k,j) = sc
            end do
          end do
        end do
      end do
      return
      end      
      
!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

      SUBROUTINE M33INV (A, AINV, OK_FLAG)

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  
     &      - A(1,1)*A(2,3)*A(3,2)  
     &      - A(1,2)*A(2,1)*A(3,3)  
     &      + A(1,2)*A(2,3)*A(3,1)  
     &      + A(1,3)*A(2,1)*A(3,2)  
     &      - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV
