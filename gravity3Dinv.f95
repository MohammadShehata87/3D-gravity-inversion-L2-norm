!--------------------------------------------------------------
!      3-D INVERSION PROGRAM FOR Gravity METHOD
!
!    Gravity3Dinv_shehata by H. Mizunaga @ Kyushu Univ. 
!
!--------------------------------------------------------------
!   << CONSTANTS >>
!
!     MXP   : MAX. NUM. OF PARAMETER
!     MXD   : MAX. NUM. OF OBSERVED DATA
!     MXN   : MAX. NUM. OF NODE(X,Y AND Z DIRECTION)
!
!   << VARIABLES >>
!
!     ALOC   : LOCATION
!     MDATA  : NUM. OF DATA
!     NPARA  : NUM. OF PARAM.
!     PAR()  : PARAMETER , LOG SCALE
!     DPAR() : MODIFIED TERM FOR PAR()
!     AXIS() : X-AXIS AND Y-AXIS OF OBSERVED POINTS
!     OBS0() : OBSERVED DATA
!     WT()   : WEIGHT FACTOR
!     AJAC() : JACOBIAN MATRIX
!--------------------------------------------------------------
      Implicit None
!      Integer, Parameter :: mxp=20000, mxd=2000, mxn=200
      Integer, Parameter :: mxp=20000, mxd=30000, mxn=20000
!      Integer, Parameter :: itmax=200, mdata=500
!      Integer, Parameter :: itmax=10, mdata=2000
      Integer, Parameter :: itmax=20, mdata=28435
!      Integer, Parameter :: itmax=10, mdata=28435
      Real(kind(0d0)), Parameter :: pi=3.14159265358979D0
      Real(kind(0d0)), Parameter :: eps=1.0D-10, epl=0.5D-3
      Real(kind(0d0)) :: par(mxp), dpar(mxp), axis(mxd, 2), obs(mxd)
      Real(kind(0d0)) :: par0(mxp)
      Real(kind(0d0)) :: ajac(mxd+mxp, mxp), dx(mxd+mxp)
      Real(kind(0d0)) :: aj(mxd+mxp, mxp), obs0(mxd), wt(mxd)
!
      Real(kind(0d0)) :: xa(mxn), ya(mxn), za(mxn), wij(mxd, mxp)
      Real(kind(0d0)) :: fnew(mxd), fold(mxd), rnew(mxd), rold(mxd), ptem(mxp)
      Real(kind(0d0)) :: smooth, rmyu, snew, sold
      Real(kind(0d0)) :: rms, abf, abdf, delta, parx, calc
      Real(kind(0d0)) :: rnd1, rnd2
      Real(kind(0d0)) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,rho0,g
      Integer :: i, j, k, nxx, nyy, nzz, npara, npa, np, ierr, icond
      character*10 :: dummy
      character*20 :: fname
      character*3  :: s
!-------------------
!    for LSQR
!-------------------
      External :: aprod
      Integer :: istop, itnlim, nout, m, n, itr, ijk, ijcode, itn
      Real(kind(0d0)) :: damp, dra, rm
!
!      Integer, Parameter :: imul=30, maxm=3000*imul, maxn=3000*imul
!      Integer, Parameter :: imul=30, maxm=10000*imul, maxn=10000*imul
!      Integer, Parameter :: imul=30, maxm=20000*imul, maxn=20000*imul
      Integer, Parameter :: imul=30, maxm=30000*imul, maxn=30000*imul
!      Integer, Parameter :: imul=200, maxm=20000*imul, maxn=20000*imul
      Real(kind(0d0)) :: u(maxm), v(maxn), w(maxn), x(maxn), se(maxn)
      Real(kind(0d0)) :: atol, btol, conlim, anorm, acond, rnorm, arnorm, xnorm
!
!      Integer, Parameter :: imax=3000, jmax=3000, ijmax=imax*jmax/10
!      Integer, Parameter :: imax=10000, jmax=10000, ijmax=imax*jmax/10
!      Integer, Parameter :: imax=20000, jmax=20000, ijmax=imax*jmax/10
      Integer, Parameter :: imax=30000, jmax=30000, ijmax=imax*jmax/10
      Real(kind(0d0)) :: val(ijmax)
      Integer :: irow(imax), jcol(ijmax)
      Common /lbl1/val
      Common /lbl2/irow, jcol
!
!      call random_seed()
!------------------------------
!  Location of Observed Points
!------------------------------
!      Do i = 1, mdata
!        call random_number(rnd1)
!        call random_number(rnd2)
!        axis(i, 1)=-10.0d0+20.0d0*rnd1  ! (km)
!        axis(i, 2)=-10.0d0+20.0d0*rnd2  ! (km)
!        wt(i)=1.0d0
!      End Do
!------------------------------
      Open (1, file = 'Data.txt', Status = 'Old', Action = 'Read')
      read(1,*) dummy
      Do i = 1, mdata
        Read (1, *) axis(i,1), axis(i,2), obs0(i)
        wt(i)= 1.0d0
      End Do
      close(1)
!------------------------------
!  Location of gravity blocks
!------------------------------
      nxx=11
      nyy=21
!      nxx=14
!      nyy=25
      nzz=11
!      
      do i=1, nxx
        xa(i)=374.0d0+25.0d0*(i-1)
!        xa(i)=374.0d0+20.0d0*(i-1)
      enddo
!      
      do j=1, nyy
        ya(j)=4250.0d0+25.0d0*(j-1)
!        ya(j)=4250.0d0+20.0d0*(j-1)
      enddo
!      
      do k=1, nzz
        za(k)=1.0*(k-1)
!        za(k)=2.0*(k-1)
      enddo
!
      npara = (nxx-1)*(nyy-1)*(nzz-1) ! Caution
!
!      do i=1,npara
!        call random_number(rnd1)
!        par0(i)=1000.0d0*(rnd1-0.5d0)*2.0  ! Density of each block (kg/m^3)
!      enddo        
      par0=0.0d0
!-----------------------
!     CALC. JACOBIAN
!-----------------------
      rho0 = 1.0d0
      z0=0.0d0
!      
      npa = 0
      Do k = 1, nzz-1
        z1 = za(k)
        z2 = za(k+1)
        Do j = 1, nyy-1
          y1 = ya(j)
          y2 = ya(j+1)
          Do i = 1, nxx-1
            x1 = xa(i)
            x2 = xa(i+1)
            npa = npa + 1
            Do m = 1, mdata
              x0 = axis(m, 1)
              y0 = axis(m, 2)
!              call gbox1(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho0,g)
              call gboxnew1(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho0,g)
              wij(m, npa) = g
              ajac(m, npa) = g*sqrt(wt(m))
            End Do
!
          End Do
        End Do
      End Do
!-----------------------
!   Calc. Observed Data
!-----------------------
!      do i=1, mdata
!        obs0(i)=0.0d0
!        do j=1, npara
!          obs0(i)=obs0(i)+wij(i,j)*par0(j)
!        enddo
!      enddo
      obs=obs0            
!----------------------------
!   SMOOTHNESS FACTOR
!----------------------------
      smooth=0.05d0 ! Original
!      smooth=0.1d0
!      smooth=0.005d0
      rmyu = sqrt(smooth)
!
      Do i = 1, npara
        par(i) = 0.0d0
      End Do
!
      snew = 0.0d0
      Do i = 1, mdata
        rnew(i) = obs(i)*sqrt(wt(i))
        fnew(i) = 0.0d0
        snew = snew + rnew(i)*rnew(i)
      End Do
      
!=========================
!  LEAST SQUARE SOLUTION
!=========================
      Do itr = 1, itmax
!-------------------
!  INPUT JACOBIAN
!-------------------
        Do i = 1, mdata
          Do j = 1, npara
            aj(i, j) = ajac(i, j)
          End Do
        End Do
        Do i = 1, npara
          Do j = 1, npara
            aj(i+mdata, j) = 0.0
          End Do
        End Do
!--------------------------------------
!  INPUT SMOOTH PARAMETER TO JACOBIAN
!--------------------------------------
        ijk = 0
        rm = rmyu
        Do k = 1, nzz-1
          rm = rmyu*(1.0+(k-1)/10.0)
          Do j = 1, nyy-1
            Do i = 1, nxx-1
              ijk = ijk + 1
              If (i/=1) Then
                aj(ijk+mdata, ijk-1) = rm
                aj(ijk+mdata, ijk) = aj(ijk+mdata, ijk) - rm
              End If
              If (i/=nxx) Then
                aj(ijk+mdata, ijk+1) = rm
                aj(ijk+mdata, ijk) = aj(ijk+mdata, ijk) - rm
              End If
              If (j/=1) Then
                aj(ijk+mdata, ijk-nxx) = rm
                aj(ijk+mdata, ijk) = aj(ijk+mdata, ijk) - rm
              End If
              If (j/=nyy) Then
                aj(ijk+mdata, ijk+nxx) = rm
                aj(ijk+mdata, ijk) = aj(ijk+mdata, ijk) - rm
              End If
              If (k/=1) Then
                aj(ijk+mdata, ijk-nxx*nyy) = rm
                aj(ijk+mdata, ijk) = aj(ijk+mdata, ijk) - rm
              End If
              If (k/=nzz) Then
                aj(ijk+mdata, ijk+nxx*nyy) = rm
                aj(ijk+mdata, ijk) = aj(ijk+mdata, ijk) - rm
              End If
!
            End Do
          End Do
        End Do
!
        Do i = 1, mdata
          dx(i) = rnew(i)
        End Do
        Do i = 1, npara
          dx(i+mdata) = 0.0
        End Do
!
!        Call lls(aj, mxd+mxp, mdata+npara, npara, dx, ierr)
!
!-------  use lsqr -------------
        atol = 1.0D-10
        btol = atol
        m = mdata + npara
        n = npara
        damp = 0.0
        itnlim = m + n + 50
        nout = 6
        Do i = 1, m
          u(i) = dx(i)
        End Do
        ijcode = 0
        Do i = 1, m
          irow(i) = ijcode + 1
          Do j = 1, n
            If (aj(i,j)/=0.0) Then
              ijcode = ijcode + 1
              val(ijcode) = aj(i, j)
              jcol(ijcode) = j
            End If
          End Do
        End Do
        irow(m+1) = ijcode + 1
        Call lsqr(m, n, aprod, damp, u, v, w, x, se, atol, btol, conlim, &
          itnlim, nout, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
!-------  use lsqr -------------
!----------------------
!  CALC. NEW DATA
!----------------------
        Do i = 1, mdata
          fold(i) = fnew(i)
          rold(i) = rnew(i)
        End Do
        Do i = 1, npara
!          dpar(i) = dx(i)
          dpar(i) = x(i)
          ptem(i) = par(i) + dpar(i)
        End Do
!----------------------------
!   CALC. NEW SP VALUES
!----------------------------
        Do i = 1, mdata
          dra = 0.0
          Do j = 1, npara
            dra = dra + ptem(j)*wij(i, j)
          End Do
          fnew(i) = dra
        End Do
!----------------------------------
!   CALC. NEW RESIDUAL
!----------------------------------
        sold = snew
        snew = 0.0
        Do i = 1, mdata
          rnew(i) = (obs(i)-fnew(i))*sqrt(wt(i))
          snew = snew + rnew(i)*rnew(i)
        End Do
        rms = sqrt(snew/mdata)
!-----------------------
!   CONVERGENCE CHECK
!-----------------------
        If (snew<sold) Then
          abf = 0.0
          abdf = 0.0
          Do i = 1, mdata
            abf = abf + abs(fnew(i))
            abdf = abdf + abs(fold(i)-fnew(i))
          End Do
          delta = abdf/abf
          If (delta<=epl) Then
            Do i = 1, npara
              par(i) = par(i) + dpar(i)
            End Do
            icond = 0
            Go To 170
          End If
!-----------------------
!   SET NEW PARAMETER
!-----------------------
          Do i = 1, npara
            par(i) = par(i) + dpar(i)
          End Do
          rmyu = 0.5*rmyu
        Else
          Do i = 1, npara
            dpar(i) = 0.0
          End Do
          rmyu = 2.0*rmyu
        End If
!-------------------------------
!  PRINT INTERMEDIATE RESULT
!-------------------------------
!      WRITE(6,*) ' '
!      WRITE(6,*) ' **************************** '
!      WRITE(6,*) '    INTERMEDIATE RESULTS '
!      WRITE(6,*) ' **************************** '
!      WRITE(6,*) '  '
!      WRITE(6,*) '   ITERATION : ',ITR
!      WRITE(6,*) '   RMS       : ',RMS
        Write (6, *) '  ITERATION: ', itr, ' RMS: ', rms
!      WRITE(6,*) '  '
!
!      WRITE(6,*) ' SP POINTS  X-AXIS   Y-AXIS   Z-AXIS   STRENGTH'
      write (s,'(i3)') itr  
      fname='shehata' // s // '.txt'
      open(2,file=fname,status='new')
!      
      Write (2, *) '  ITERATION: ', itr, ' RMS: ', rms
      NP=0
      DO 981 K=1,NZZ-1
       DO 982 J=1,NYY-1
         DO 983 I=1,NXX-1
          NP=NP+1
          PARX=PAR(NP)
!          WRITE(6,112) NP,XA(I),YA(J),ZA(K),PARX
          Write (2, '(i6,6f8.2,e14.6)') np, xa(i), xa(i+1), ya(j), ya(j+1), za(k), za(k+1), parx
!          WRITE(2,'(i6,4f14.6)') NP,XA(I),YA(J),ZA(K),PARX
  983    CONTINUE
  982  CONTINUE
  981 CONTINUE
!
      close(2)
!
      End Do
!
170   Continue
!-------------------------
!    PRINT FINAL RESULT
!-------------------------
!
      fname='shehata_final.txt'
      open(3,file=fname,status='new')
!
      Write (3, *) '  '
      Write (3, *) ' **************************************************** '
      Write (3, *) '  3-D Gravity inversion by H.Mizunaga @ Kyushu Univ.  '
      Write (3, *) ' **************************************************** '
      Write (3, *) '  '
      Write (3, 180) smooth
      Write (3, *) '  '
180   Format ('    SMOOTHNESS FACTOR   : ', F10.5)
!
!      Write (6, *) ' SP POINTS  X-AXIS   Y-AXIS   Z-AXIS   STRENGTH'
!
      np = 0
      Do k = 1, nzz-1
        Do j = 1, nyy-1
          Do i = 1, nxx-1
            np = np + 1
            parx = par(np)
            Write (3, '(i6,6f8.2,e14.6)') np, xa(i), xa(i+1), ya(j), ya(j+1), za(k), za(k+1), parx
!            Write (6, '(i6,6f8.2,e14.6)') np, xa(i), xa(i+1), ya(j), ya(j+1), za(k), za(k+1), parx
!            Write (6, '(i6,6f8.2,2e14.6)') np, xa(i), xa(i+1), ya(j), ya(j+1), za(k), za(k+1), parx, par0(np)
!            Write (6, '(4i6,2e14.6)') np, i, j, k, parx, par0(np)
          End Do
        End Do
      End Do
190   Format (' ', I7, 2X, 3F8.2, 2E14.6)
!
      Write (6, *) '  '
      Write (6, *) '  NUM.        OBS.           CAL.           RES. '
      Do i = 1, mdata
        calc = fnew(i)
!        Write (6, 200) i, obs0(i), calc, obs0(i) - calc
!        Write (6, '(i6,2f10.3,3e15.7)') i, axis(i,1), axis(i,2), obs0(i), calc, obs0(i) - calc
        Write (3, '(i6,2f10.3,3e15.7)') i, axis(i,1), axis(i,2), obs0(i), calc, obs0(i) - calc
      End Do
!
      close(3)      
200   Format (' ', I6, 3E15.7)
!
    End Program
!------------------
!   Haaz(1953)
!------------------
subroutine gbox1(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho,g)
  implicit none 
  real(kind(0d0)),intent(in)  :: x0,y0,z0,x1,y1,z1,x2,y2,z2,rho
  real(kind(0d0)),intent(out) :: g
!
  real(kind(0d0)),parameter :: pi=3.14159265358979d0, twopi=2.0d0*pi
  real(kind(0d0)),parameter :: gamma=6.670d-11, si2mg=1.0d+5, km2m=1.0d+3
  real(kind(0d0)) :: x(2), y(2), z(2), sum
  real(kind(0d0)) :: rijk, arg1, arg2, arg3
!
  integer :: i, j, k, ijk 
  integer :: isign(2)
  data isign /-1, 1/
!
  if(x2-x1 <= 0.0d0) stop 'GBOX1 : x-axis error <= 0.0'
  if(y2-y1 <= 0.0d0) stop 'GBOX1 : y-axis error <= 0.0'
  if(z2-z1 <= 0.0d0) stop 'GBOX1 : z-axis error <= 0.0'
!
  x(1)=x0-x1
  y(1)=y0-y1
  z(1)=z0-z1
  x(2)=x0-x2
  y(2)=y0-y2
  z(2)=z0-z2
  sum=0.0d0
!
  do i=1,2
    do j=1,2
      do k=1,2  
        rijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
        ijk=isign(i)*isign(j)*isign(k)
        arg1=atan2(x(i)*y(j),z(k)*rijk)
        if(arg1 < 0.0d0) arg1=arg1+twopi
        arg2=rijk+y(j)
        arg3=rijk+x(i)
        sum=sum+ijk*(z(k)*arg1-x(i)*log(arg2)-y(j)*log(arg3))
      enddo
    enddo
  enddo
  g=rho*gamma*sum*si2mg*km2m        
!    
end subroutine gbox1
!
subroutine gboxnew1(x0,y0,z0,x1,y1,z1,x2,y2,z2,rho,g)
  implicit none 
  real(kind(0d0)),intent(in)  :: x0,y0,z0,x1,y1,z1,x2,y2,z2,rho
  real(kind(0d0)),intent(out) :: g
!
  real(kind(0d0)),parameter :: pi=3.14159265358979d0, twopi=2.0d0*pi
  real(kind(0d0)),parameter :: gamma=6.670d-11, si2mg=1.0d+5, km2m=1.0d+3
  real(kind(0d0)) :: x(2), y(2), z(2), sum
  real(kind(0d0)) :: rijk, arg1, arg2, arg3
!
  real(kind(0d0)) :: dx, dy, dz, xx, yy, zz, rr
  real(kind(0d0)) :: p(3), w(3)
  data p /-0.774596669241438d0,               0.0d0, 0.774596669241438d0/
  data w / 0.555555555555556d0, 0.888888888888889d0, 0.555555555555556d0/
!
  integer :: i, j, k 
!
  dx=x2-x1
  dy=y2-y1
  dz=z2-z1
  if(dx <= 0.0d0) stop 'GBOXnew : x-axis error <= 0.0'
  if(dy <= 0.0d0) stop 'GBOXnew : y-axis error <= 0.0'
  if(dz <= 0.0d0) stop 'GBOXnew : z-axis error <= 0.0'
!
  sum=0.0d0
!
  do i=1,3
    xx=0.5*(p(i)*dx+x1+x2)
    do j=1,3
      yy=0.5*(p(j)*dy+y1+y2)
      do k=1,3
        zz=0.5*(p(k)*dz+z1+z2)
        rr=sqrt((xx-x0)**2+(yy-y0)**2+(zz-z0)**2)
        sum=sum+w(i)*w(j)*w(k)*(zz-z0)/rr**3.0d0
      enddo
    enddo
  enddo
  sum=sum*dx*dy*dz/8.0d0
  g=rho*gamma*sum*si2mg*km2m        
!    
end subroutine gboxnew1
!
!     ***************
!     *     LLS     *
!     ***************
    Subroutine lls(a, max, m, n, b, irc)
      Implicit None
      Integer, Parameter :: mxp=20000
      Integer, intent(in) :: max, m, n
      Integer :: irc, i, j, k, ip
      Integer :: iw(mxp)
      Real(kind(0d0)) :: a(max, n), b(m), w(mxp)
      Real(kind(0d0)) :: sig, sum, temp, dk, bk
      Real(kind(0d0)) :: eps = 16.0d0**(-13)  !  Machine Epsilon
!
!     COMPUTE NORM FOR EACH COLUM OF ORIGINAL MATRIX,
!
!     SET UP PIVOTTING INFORMATIONS.
!
      irc = 0
      sig = -1.0d0
      Do j = 1, n
        sum = 0.0D0
        Do i = 1, m
          sum = sum + a(i, j)*a(i, j)
        End Do
        w(j) = sum
        iw(j) = j
        If (sig<sum) sig = sum
      End Do
      eps = sig*eps*eps
!
!     HOUSEHOLDER'S TRANSFORMATION
!
      Do k = 1, n
!
!     SEARCH PIVOT COLUMN WHICH HAS THE BIGGEST NORM
!
        sig = -1.0
        Do j = k, n
          If (sig>w(j)) Go To 100
          sig = w(j)
          ip = j
100     End Do
!
!     PIVOTTING (COLUMN INTERCHANGE)
!
        If (ip/=k) Then
          iw(k) = ip
          w(ip) = w(k)
          Do i = 1, m
            temp = a(i, k)
            a(i, k) = a(i, ip)
            a(i, ip) = temp
          End Do
        End If
!
        sum = 0.0D0
        Do i = k, m
          sum = sum + a(i, k)*a(i, k)
        End Do
!
!     SET UP ORTHOGONAL MATRIX OF K-TH STAGE
!
        If (sum>eps) Then
          dk = -sign(sqrt(sum), a(k,k))
          bk = 1.0D0/(sum-a(k,k)*dk)
          a(k, k) = a(k, k) - dk
          w(k) = dk
          Do j = k + 1, n
            sum = 0.0D0
            Do i = k, m
              sum = sum + a(i, j)*a(i, k)
            End Do
            w(n+j) = sum*bk
          End Do
!
          Do j = k + 1, n
            Do i = k, m
              a(i, j) = a(i, j) - a(i, k)*w(n+j)
            End Do
            w(j) = w(j) - a(k, j)*a(k, j)
          End Do
        Else
          irc = 10
          Return
        End If
      End Do
!
!     APPLY ORTHOGONAL TRANSFORM TO CONSTANT VECTOR
!
      Do j = 1, n
        sum = 0.0D0
        Do i = j, m
          sum = sum + a(i, j)*b(i)
        End Do
        sum = sum/(w(j)*a(j,j))
        Do i = j, m
          b(i) = b(i) + sum*a(i, j)
        End Do
      End Do
!
!     COMPUTE SOLUTION BY BACKWARD SUBSTITUTIONS
!
      b(n) = b(n)/w(n)
      Do i = n - 1, 1, -1
        sum = 0.0D0
        Do j = i + 1, n
          sum = sum + a(i, j)*b(j)
        End Do
        b(i) = (b(i)-sum)/w(i)
      End Do
!
!     REORDER THE SOLUTION ACCORDING TO PIVOTTING HISTORY
!
      Do j = 1, n
        i = n - j + 1
        temp = b(i)
        b(i) = b(iw(i))
        b(iw(i)) = temp
      End Do
      Return
    End Subroutine lls
!    
!-------------
!    lsqr
!-------------
    Subroutine lsqr(m, n, aprod, damp, u, v, w, x, se, atol, btol, conlim, &
        itnlim, nout, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
!
      Implicit None
      External aprod
      Integer  :: m, n, itnlim, nout, istop, itn
      Real(kind(0d0)) :: u(m), v(n), w(n), x(n), se(n), atol, btol, conlim, &
        damp, anorm, acond, rnorm, arnorm, xnorm

!     intrinsics and local variables

      Integer :: i, nconv, nstop
      Real(kind(0d0)) :: dnrm2
      Real(kind(0d0)) :: alfa, bbnorm, beta, bnorm, cs, cs1, cs2, ctol, dampsq, &
        ddnorm, delta, gamma, gambar, phi, phibar, psi, res1, res2, rho, &
        rhobar, rhbar1, rhbar2, rhs, rtol, sn, sn1, sn2, t, tau, test1, test2, &
        test3, theta, t1, t2, t3, xxnorm, z, zbar
!
      Real(kind(0d0)), Parameter :: zero=0.0D+0, one=1.0D+0
!
      Character *16 :: enter, exit
      Character *60 :: msg(0:7)
!
      Data enter/' enter lsqr.    '/, exit/' exit  lsqr.    '/

!-----------------------------------------------------------------------
      itn = 0
      istop = 0
      nstop = 0
      ctol = zero
      If (conlim>zero) ctol = one/conlim
      anorm = zero
      acond = zero
      bbnorm = zero
      dampsq = damp**2
      ddnorm = zero
      res2 = zero
      xnorm = zero
      xxnorm = zero
      cs2 = -one
      sn2 = zero
      z = zero
!
      Do i = 1, n
        v(i) = zero
        x(i) = zero
        se(i) = zero
      End Do
!     set up the first vectors u and v for the bidiagonalization.
!     these satisfy  beta*u = b,  alfa*v = a(transpose)*u.
      alfa = zero
      beta = dnrm2(m, u, 1)
      If (beta>zero) Then
        Call dscal(m, (one/beta), u, 1)
        Call aprod(2, m, n, v, u)
        alfa = dnrm2(n, v, 1)
      End If
      If (alfa>zero) Then
        Call dscal(n, (one/alfa), v, 1)
        Call dcopy(n, v, 1, w, 1)
      End If
      arnorm = alfa*beta
      If (arnorm==zero) Go To 110
      rhobar = alfa
      phibar = beta
      bnorm = beta
      rnorm = beta
      If (nout>0) Then
        test1 = one
        test2 = alfa/beta
      End If
!     ------------------------------------------------------------------
!     main iteration loop.
!     ------------------------------------------------------------------
100   itn = itn + 1
!     perform the next step of the bidiagonalization to obtain the
!     next  beta, u, alfa, v.  these satisfy the relations
!                beta*u  =  a*v  -  alfa*u,
!                alfa*v  =  a(transpose)*u  -  beta*v.
      Call dscal(m, (-alfa), u, 1)
      Call aprod(1, m, n, v, u)
      beta = dnrm2(m, u, 1)
      bbnorm = bbnorm + alfa**2 + beta**2 + dampsq
      If (beta>zero) Then
        Call dscal(m, (one/beta), u, 1)
        Call dscal(n, (-beta), v, 1)
        Call aprod(2, m, n, v, u)
        alfa = dnrm2(n, v, 1)
        If (alfa>zero) Then
          Call dscal(n, (one/alfa), v, 1)
        End If
      End If
!     use a plane rotation to eliminate the damping parameter.
!     this alters the diagonal (rhobar) of the lower-bidiagonal matrix.
      rhbar2 = rhobar**2 + dampsq
      rhbar1 = sqrt(rhbar2)
      cs1 = rhobar/rhbar1
      sn1 = damp/rhbar1
      psi = sn1*phibar
      phibar = cs1*phibar
!     use a plane rotation to eliminate the subdiagonal element (beta)
!     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
      rho = sqrt(rhbar2+beta**2)
      cs = rhbar1/rho
      sn = beta/rho
      theta = sn*alfa
      rhobar = -cs*alfa
      phi = cs*phibar
      phibar = sn*phibar
      tau = sn*phi
!     update  x, w  and the standard error estimates.
      t1 = phi/rho
      t2 = -theta/rho
      t3 = one/rho
      Do i = 1, n
        t = w(i)
        x(i) = t1*t + x(i)
        w(i) = t2*t + v(i)
        t = (t3*t)**2
        se(i) = t + se(i)
        ddnorm = t + ddnorm
      End Do
!     use a plane rotation on the right to eliminate the
!     super-diagonal element (theta) of the upper-bidiagonal matrix.
!     then use the result to estimate  norm(x).
      delta = sn2*rho
      gambar = -cs2*rho
      rhs = phi - delta*z
      zbar = rhs/gambar
      xnorm = sqrt(xxnorm+zbar**2)
      gamma = sqrt(gambar**2+theta**2)
      cs2 = gambar/gamma
      sn2 = theta/gamma
      z = rhs/gamma
      xxnorm = xxnorm + z**2
!     test for convergence.
!     first, estimate the norm and condition of the matrix  abar,
!     and the norms of  rbar  and  abar(transpose)*rbar.
      anorm = sqrt(bbnorm)
      acond = anorm*sqrt(ddnorm)
      res1 = phibar**2
      res2 = res2 + psi**2
      rnorm = sqrt(res1+res2)
      arnorm = alfa*abs(tau)
!     now use these norms to estimate certain other quantities,
!     some of which will be small near a solution.
      test1 = rnorm/bnorm
      test2 = zero
      If (rnorm>zero) test2 = arnorm/(anorm*rnorm)
      test3 = one/acond
      t1 = test1/(one+anorm*xnorm/bnorm)
      rtol = btol + atol*anorm*xnorm/bnorm
!     the following tests guard against extremely small values of
!     atol, btol  or  ctol.  (the user may have set any or all of
!     the parameters  atol, btol, conlim  to zero.)
!     the effect is equivalent to the normal tests using
!     atol = relpr,  btol = relpr,  conlim = 1/relpr.
      t3 = one + test3
      t2 = one + test2
      t1 = one + t1
      If (itn>=itnlim) istop = 7
      If (t3<=one) istop = 6
      If (t2<=one) istop = 5
      If (t1<=one) istop = 4
!     allow for tolerances set by the user.
      If (test3<=ctol) istop = 3
      If (test2<=atol) istop = 2
      If (test1<=rtol) istop = 1
!     stop if appropriate.
!     the convergence criteria are required to be met on  nconv
!     consecutive iterations, where  nconv  is set below.
!     suggested value:  nconv = 1, 2  or  3.
      If (istop==0) nstop = 0
      If (istop==0) Go To 100
      nconv = 1
      nstop = nstop + 1
      If (nstop<nconv .And. itn<itnlim) istop = 0
      If (istop==0) Go To 100
!     ------------------------------------------------------------------
!     end of iteration loop.
!     ------------------------------------------------------------------
!     finish off the standard error estimates.
      t = one
      If (m>n) t = m - n
      If (dampsq>zero) t = m
      t = rnorm/sqrt(t)
      Do i = 1, n
        se(i) = t*sqrt(se(i))
      End Do
110   Continue
      Return
!     end of lsqr
    End Subroutine lsqr
!
    Subroutine dcopy(n, x, incx, y, incy)
      Integer n, incx, incy
      Double Precision x(n), y(n)
!
      Do i = 1, n
        y(i) = x(i)
      End Do
      Return
    End Subroutine dcopy
!
    Double Precision Function dnrm2(n, x, incx)
      Integer n, incx, i
      Double Precision x(n), d
!
      d = 0.0
      Do i = 1, n
        d = d + x(i)**2
      End Do
      dnrm2 = dsqrt(d)
!
      Return
    End Function dnrm2
!
    Subroutine dscal(n, a, x, incx)
      Integer n, incx
      Double Precision a, x(n)
!
      Do i = 1, n
        x(i) = a*x(i)
      End Do
!
      Return
    End Subroutine dscal
!********************************************************
!     aprod(1,---) : y = y + a * x
!     aprod(2,---) : x = x + a(transpose) * y
!********************************************************
    Subroutine aprod(mode, m, n, x, y)
      Implicit None
      Integer :: mode, m, n
      Real(kind(0d0)) :: x(n), y(m)
!         case 1: y = y + a * x
      If (mode==1) Then
        Call aprod1(m, n, x, y)
      End If
!         case 2: x = x + a(transpose) * y
      If (mode==2) Then
        Call aprod2(m, n, x, y)
      End If
!
      Return
    End Subroutine aprod
!
    Subroutine aprod1(m, n, x, y)
!     ------------------------------------------------------------------
!     aprod1  computes  y = y + a*x  for subroutine aprod,
!     ------------------------------------------------------------------
      Implicit None
      Integer :: m, n
      Real(kind(0d0)) :: x(n), y(m)
      Integer :: i, j
!      Integer, Parameter :: imax=3000, jmax=3000, ijmax=imax*jmax/10
!      Integer, Parameter :: imax=10000, jmax=10000, ijmax=imax*jmax/10
!      Integer, Parameter :: imax=20000, jmax=20000, ijmax=imax*jmax/10
      Integer, Parameter :: imax=30000, jmax=30000, ijmax=imax*jmax/10
      Real(kind(0d0)) :: val(ijmax)
      Integer :: irow(imax), jcol(ijmax)
      Common /lbl1/val
      Common /lbl2/irow, jcol
!
      Do i = 1, m
        Do j = irow(i), irow(i+1) - 1
          y(i) = y(i) + val(j)*x(jcol(j))
        End Do
      End Do
!
    End Subroutine aprod1
!
    Subroutine aprod2(m, n, x, y)
!     ------------------------------------------------------------------
!     aprod2  computes  x = x + a(t)*y  for subroutine aprod,
!     ------------------------------------------------------------------
      Implicit None
      Integer :: m, n
      Real(kind(0d0)) :: x(n), y(m)
      Integer :: i, j
!      Integer, Parameter :: imax=3000, jmax=3000, ijmax=imax*jmax/10
!      Integer, Parameter :: imax=10000, jmax=10000, ijmax=imax*jmax/10
!      Integer, Parameter :: imax=20000, jmax=20000, ijmax=imax*jmax/10
      Integer, Parameter :: imax=30000, jmax=30000, ijmax=imax*jmax/10
      Real(kind(0d0)) :: val(ijmax)
      Integer :: irow(imax), jcol(ijmax)
      Common /lbl1/val
      Common /lbl2/irow, jcol
!
      Do j = 1, m
        Do i = irow(j), irow(j+1) - 1
          x(jcol(i)) = x(jcol(i)) + val(i)*y(j)
        End Do
      End Do
!
    End Subroutine aprod2
