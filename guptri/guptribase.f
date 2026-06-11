        subroutine upddel(total,xinc)
c       implicit none
c
c****   formal parameter declarations
        real*8 total, xinc
c
c****   accumulate root sum of squares in total with increment xinc
c       assume both arguments nonnegative
c
c****   this version dated june 16, 1987
c       authors: jim demmel and bo kagstrom
c 
        if (total.gt.xinc) then
          total = total * sqrt(1.0 + (xinc/total)**2)
        elseif (total.lt.xinc) then
          total = xinc * sqrt(1.0 + (total/xinc)**2)
        else
          total = total * sqrt(2.0)
        endif
        return
        end
c
        subroutine cident(c,ldc,n)
c       implicit none
c
c****   formal parameter declarations
        integer ldc, n
        complex*16 c(ldc,n)
c
c****   set c = n by n indentity matrix
c
c****   this version dated june 16, 1987
c
c****   internal variables
        integer i, j
c
        do 1 j=1,n
          do 2 i=1,n
            c(i,j)=0
2         continue
          c(j,j)=1
1       continue
        return
        end
c
        subroutine krnstr(m,n,kstr,kfirst,last,nisum,risum,case,
     *              nmsing,nmjord,dmjord)
c       implicit none
c
c****   debug space
        common /debug2/ idbg(20),outunit
        integer idbg, outunit
c
c****   formal parameter declarations
        integer m,n, kstr(4,*), kfirst, last, nisum, risum
        integer case, nmsing, nmjord, dmjord
c
c****   interpret null space dimensions as kronecker indices
c       there are 7 cases (for details see the code below)
c
c****   this version dated june 16, 1987
c       authors: jim demmel and bo kagstrom
c
c****   internal variables
        integer j, lastm1, ni, ri, nnew, rnew
        logical ldebug
c
c       set debug flag
        ldebug= (idbg(3).ne.0)
c****
c       in cases 2,3,6 below, adjoin column to kstr so
c       kstr(2,last)=0 in all cases
        if (kstr(1,last).eq.0) kstr(2,last)=0
        nisum=0
        risum=0
        if (last.ge.kfirst) then
          do 1 j=kfirst,last
            nisum=nisum+kstr(1,j)
            risum=risum+kstr(2,j)
1         continue
        end if
        nnew=n-nisum
        rnew=m-risum
c
        ni=kstr(1,last)
        ri=kstr(2,last)
        if (ldebug) write(outunit,100) n,m,kfirst,last,nisum,risum,
     *   nnew,rnew,ni,ri
100     format(//'entering krnstr',/
     *     'n,m,kfirst,last,nisum,risum,nnew,rnew,ni,ri=',10i3)
        if (ldebug) write(outunit,101) (j,j=1,last)
101     format('kstr='/20i4)
        if (ldebug) write(outunit,102) (kstr(1,j),j=1,last)
        if (ldebug) write(outunit,102) (kstr(2,j),j=1,last)
102     format(20i4)
c
        if (ni.eq.0 .and. nnew.gt.0 .and. rnew.gt.0) then
c****     case 1
c         rest of pencil begins at (risum+1,nisum+1)
          case=1
        else if (ri.gt.0 .and. nnew.eq.0 .and. rnew.eq.0) then
c****     case 2
c         entire pencil reduced; no indices or eigenvalues of other type
          case=2
          last=last+1
          kstr(1,last)=0
          kstr(2,last)=0
          ni=0
          ri=0
        else if (ri.gt.0 .and. nnew.gt.0 .and. rnew.eq.0) then
c****     case 3
c         entire pencil reduced; no indices or eigenvalues of other type
          case=3
          last=last+1
          kstr(1,last)=nnew
          kstr(2,last)=0
          ni=nnew
          ri=0
          nisum=nisum+nnew
          nnew=0
        else if (ni.gt.0 .and. ri.eq.0 .and. rnew.gt.0 .and. nnew.gt.0)
     *       then
c****     case 4
c         rest of pencil begins at (risum+1,nisum+1)
          case=4
        else if (ni.gt.0 .and. ri.eq.0 .and. nnew.eq.0 .and. rnew.gt.0)
     *       then
c****     case 5
c         entire pencil reduced; last rnew rows are 0
c         (i.e. there are rnew zero indices of other type)
          case=5
        else if (ni.gt.0 .and. ri.gt.0 .and. nnew.eq.0 .and. rnew.gt.0)
     *       then
c****     case 6
c         entire pencil reduced; last rnew rows are 0
c         (i.e. there are rnew zero indices of other kind)
          case=6
          last=last+1
          kstr(1,last)=0
          kstr(2,last)=0
          ni=0
          ri=0
        else
c****     cannot happen, error state, print error message
          if (ldebug) write(outunit,105)
105       format(//'error condition')
          case=7
        end if
c
        if (ldebug) then
          write(outunit,107) case,n,m,kfirst,last,nisum,risum,nnew,
     +     rnew,ni,ri
107       format(/' case,n,m,kfirst,last,nisum,risum,nnew,rnew,ni,ri=',
     +     /,11i4)
          write(outunit,101) (j,j=1,last)
          write(outunit,102) (kstr(1,j),j=1,last)
          write(outunit,102) (kstr(2,j),j=1,last)
        endif
c       compute number of singular blocks
        nmsing=nisum-risum
c
c       compute number of jordan blocks
        nmjord=-nmsing+kstr(1,kfirst)
c
c       compute dimension of jordan blocks
        dmjord=0
        if (last.gt.kfirst) then
          lastm1=last-1
          do 3 j=kfirst,lastm1
            dmjord=dmjord+(j-kfirst+1)*(kstr(2,j)-kstr(1,j+1))
3         continue
        end if
        if (ldebug) then
          write(outunit,106) case,nmsing,nmjord,dmjord
106       format(/'case,nmsing,nmjord,dmjord=',4i4)
          write(outunit,101) (j,j=1,last)
          write(outunit,102) (kstr(1,j),j=1,last)
          write(outunit,102) (kstr(2,j),j=1,last)
        endif
c
        return
        end
c
        real*8 function norme(a, ldab, m, n)
c       implicit none
c****   formal parameter declarations
        integer ldab, m, n
        complex*16 a(ldab,*)
c****   compute frobenius norm of matrix a
c
        real*8 sum
        integer i, j
c
        sum = 0.
        do 1 i = 1, m
          do 2 j = 1, n
            sum = sum + dreal(a(i,j))**2 + dimag(a(i,j))**2
2         continue
1       continue
        norme = sqrt(sum)
        return
        end

c     On this file - blas routines: double precision complex
c     zaxpy, zswap, dcabs1, dznrm2, zcopy, zdotc, zdotu, zscal,
c     zrotg, zdrot, drotg
c
      subroutine zaxpy(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c
      double complex zx(1),zy(1),za
      double precision dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end
      
      subroutine  zswap (n,zx,incx,zy,incy)
c
c     interchanges two vectors.
c     jack dongarra, 3/11/78.
c
      double complex zx(1),zy(1),ztemp
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
 
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end

 
      double precision function dznrm2( n, zx, incx)
      logical imag, scale
      integer          next
      double precision cutlo, cuthi, hitest, sum, xmax, absx, zero, one
      double complex      zx(1)
      double precision dreal,dimag
      double complex zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      data         zero, one /0.0d0, 1.0d0/
c
c     unitary norm of the complex n-vector stored in zx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson , 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dznrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      do 210 i=1,nn,incx
         absx = dabs(dreal(zx(i)))
         imag = .false.
         go to next,(30, 50, 70, 90, 110)
   30 if( absx .gt. cutlo) go to 85
      assign 50 to next
      scale = .false.
c
c                        phase 1.  sum is zero
c
   50 if( absx .eq. zero) go to 200
      if( absx .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 assign 110 to next
      sum = (sum / absx) / absx
  105 scale = .true.
      xmax = absx
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( absx .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( absx .le. xmax ) go to 115
         sum = one + sum * (xmax / absx)**2
         xmax = absx
         go to 200
c
  115 sum = sum + (absx/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
   85 assign 90 to next
      scale = .false.
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
      hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
   90 if(absx .ge. hitest) go to 100
         sum = sum + absx**2
  200 continue
c                  control selection of real and imaginary parts.
c
      if(imag) go to 210
         absx = dabs(dimag(zx(i)))
         imag = .true.
      go to next,(  50, 70, 90, 110 )
c
  210 continue
c
c              end of main loop.
c              compute square root and adjust for scaling.
c
      dznrm2 = dsqrt(sum)
      if(scale) dznrm2 = dznrm2 * xmax
  300 continue
      return
      end
 
      subroutine  zcopy(n,zx,incx,zy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 4/11/78.
c
      double complex zx(1),zy(1)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end
 
      double complex function zdotc(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c
      double complex zx(1),zy(1),ztemp
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end
 
      double complex function zdotu(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c
      double complex zx(1),zy(1),ztemp
      ztemp = (0.0d0,0.0d0)
      zdotu = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + zx(ix)*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotu = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + zx(i)*zy(i)
   30 continue
      zdotu = ztemp
      return
      end
 
      subroutine  zscal(n,za,zx,incx)
c
c    scales a vector by a constant.
c    jack dongarra, 3/11/78.
c
      double complex za,zx(1)
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increments not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increments equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
 
      subroutine zrotg(ca,cb,c,s)
      double complex ca,cb,s
      double precision c,dcabs1
      double precision norm,scale
      double complex alpha
      if (dcabs1(ca) .ne. 0.0d0) go to 10
         c = 0.0d0
         s = (1.0d0,0.0d0)
         ca = cb
         go to 20
   10 continue
         scale = dcabs1(ca) + dcabs1(cb)
         norm = scale*dsqrt((dcabs1(ca/dcmplx(scale,0.0d0)))**2 +
     *                      (dcabs1(cb/dcmplx(scale,0.0d0)))**2)
         alpha = ca /dcabs1(ca)
         c = dcabs1(ca) / norm
         s = alpha * dconjg(cb) / norm
         ca = alpha * norm
   20 continue
      return
      end
 
      subroutine  zdrot (n,zx,incx,zy,incy,c,s)
c
c     applies a plane rotation, where the cos and sin (c and s) are
c     double precision and the vectors zx and zy are double complex.
c     jack dongarra, linpack, 3/11/78.
c
      double complex zx(1),zy(1),ztemp
      double precision c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = c*zx(ix) + s*zy(iy)
        zy(iy) = c*zy(iy) - s*zx(ix)
        zx(ix) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = c*zx(i) + s*zy(i)
        zy(i) = c*zy(i) - s*zx(i)
        zx(i) = ztemp
   30 continue
      return
      end
 
      subroutine drotg(da,db,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,db,c,s,roe,scale,r,z
c
      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
      if( scale .ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         go to 20
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
   20 z = 1.0d0
      if( dabs(da) .gt. dabs(db) ) z = s
      if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
      da = r
      db = z
      return
      end
c
c*** no more on this file

c     as of june 22, 1987 this file contains
c     bound, ebdreg, gvec, pbound, blddfl, blddfu, bldrhs, prml,
c     prmlct, svdiv, evalbd, bndwsp
c
      subroutine bound(a,b,ldab,m,n,irstrt,icstrt,dimreg,
     +                 evala,evalb,edlmax,gvcond,pqnorm,ecase,
     +                 sdlmax, difl, difu, qnorm, pnorm, scase, 
     +                 work, info)
c
c     implicit none
c
c**** debug space
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
c**** formal parameter declarations
      integer ldab,m,n,irstrt,icstrt,dimreg,info,ecase,scase
      complex*16 a(ldab,*), b(ldab,*), evala(*),evalb(*)
      complex*16 work(*)
      real*8 gvcond(*), edlmax, sdlmax, qnorm, pnorm, pqnorm
      real*8 difl, difu
c
c********************************************************************
c
c     compute error bounds for selected eigenvalues of general pencil
c     and error bounds for left and right reducing subspaces
c
c     this version requires all selected eigenvalues be simple
c     input pencil a - lambda b must be in guptri form
c
c     theorems and corollaries referred to below appear in:
c     'accurate solutions of ill-posed problems in control theory'
c     proc. of the 25th ieee conference on decision and control,
c     athens, greece, december 10-12, 1986, pp 558-563
c     by j. demmel and b. kagstrom
c
c     see also:
c     j. demmel and b. kagstrom, 'computing stable eigendecompositions
c        of matrix pencils', linear algebra and its applications,
c        vol 88/89, 1987, pp 139-185
c
c     inputs
c
c       a(ldab,n), b(ldab,n) - complex*16 - input pencil in 
c                                           guptri form
c
c       lda - integer - leading dimension of a and b
c
c       m,n - integer - row, column dimensions of a and b
c
c       irstrt, icstrt - integer - starting row and column of selected 
c                        part of pencil for which eigenvalue bounds 
c                        are desired. reducing subspace bounds will be
c                        supplied for right reducing subspace spanned
c                        by leading icstrt-1 components and for left
c                        reducing subspace spanned by leading icstrt-1
c                        components.
c                        note: set icstrt=n+1 to make right reducing
c                                  subspace whole space
c                              set irstrt=m+1 to make left reducing
c                                  subspace whole space
c
c       dimreg - integer - number of selected eigenvalues;
c         if dimreg.eq.0 only subspace perturbation bounds will be
c         computed
c        (note - one can select a subset of the regular part only;
c         this gives generally different bounds for common eigenvalues
c         from a different selected subset; see paper above for 
c         discussion)
c
c     outputs
c
c       evala(dimreg), evalb(dimreg) - complex*16 - 
c          normalized selected eigenvalues;
c          evala(i)/evalb(i) is i-th eigenvalue and
c          abs(evala(i))**2 + abs(evalb(i))**2 = 1
c
c       edlmax - real*8 - maximum frobenius norm of perturbation for 
c                which eigenvalue perturbation bounds hold. 
c                if no maximum norm then edlmax=-1.
c
c       gvcond(dimreg) - real*8 - condition numbers; suppose the pencil
c         is perturbed by amount delta .le. edlmax (if edlmax=-1. then
c         delta arbitrary) such that the conditions of theorem 5 or 
c         corollary 1 hold (edlmax=-1. implies these conditions always
c         hold). then if c/s is a perturbed eigenvalue such that
c         abs(c)**2 + abs(s)**2 = 1, then for some i
c         abs(c*evalb(i)-s*evala(i)) .le. delta * gvcond(i)
c
c       pqnorm - real*8 - overall condition number; under same 
c         conditions as for gvcond, if areg - lambda breg is regular 
c         part of unperturbed pencil in guptri form, then
c         sigma-min(c*breg - s*areg) .le. delta * pqnorm
c         (sigma-min is the smallest singular value)
c
c       ecase - integer - which of 5 cases for eigenvalue bounds 
c               the pencil falls depending on input dimensions;
c               the first four cases are for dimreg.gt.0, in which
c               case the description gives:
c                  (part of KCF to above, left of selected part) and
c                  (part of KCF to below, right of selected part) 
c          ecase=1 - (right singular and/or regular part) and
c                    (left singular and/or regular part)
c          ecase=2 - (right singular and/or regular part) and (nothing)
c          ecase=3 - (nothing) and (left singular and/or regular part)
c          ecase=4 - (nothing) and (nothing)
c          ecase=5 - dimreg.eq.0 (no eigenvalue bounds)
c
c       sdlmax - real*8 - maximum frobenius norm of perturbation for 
c                which reducing subspace perturbation bounds hold
c                (if scase=4 (see below) sdlmax=-1. to indicate that
c                 this bound does not apply)
c
c       difl, difu - real*8 - difl and difu functions (used to
c                    compute sdlmax, see paper for details)
c                    (if scase=4 (see below), both set to 0)
c
c       qnorm, pnorm, - real*8 - norms of left and right projectors
c                       (used in reducing subspace bounds)
c                       (if scase=4 (see below), both set to 1)
c
c       scase - integer - which of 4 cases for reducing subspace
c               bounds the pencil falls depending on input dimensions:
c          scase=1 - both left and right subspaces nontrivial
c          scase=2 - left space trivial (0) and right space nontrivial
c          scase=3 - left space nontrivial and right space trivial
c                   (whole space)
c          scase=4 - both spaces trivial (either 0 or whole space)
c
c       the reducing subspace bounds may be calculated from 
c          scase, sdlmax, pnorm and qnorm as follows:
c          let delta be the distance in the frobenius norm from a
c          perturbed pencil with the same structure as a - lambda b
c          to a - lambda b (see the above paper by demmel and
c          kagstrom for more details). if delta.lt.sdlmax then the
c          following bounds apply, where relerr=delta/sdlmax :
c
c          upper bound on angular perturbation in left reducing subspace
c            if scase=1 (theorem 4, case 1 in paper)
c              atan(relerr/(pnorm-relerr*sqrt(pnorm**2-1)))
c            if scase=2
c              0 (since left subspace trivial)
c            if scase=3
c              atan(relerr/(1-relerr))
c            if scase=4
c              0 (since left subspace trivial)
c
c          upper bound on angular perturbation in right reducing subspace
c            if scase=1 (theorem 4, case 1 in paper)
c              atan(relerr/(qnorm-relerr*sqrt(qnorm**2-1)))
c            if scase=2
c              atan(relerr/(1-relerr))
c            if scase=3
c              0 (since right subspace trivial)
c            if scase=4
c              0 (since right subspace trivial)
c
c          lower bound on angular perturbation in left reducing subspace
c            if scase=1 (theorem 4, case 2 in paper)
c              atan(1/(sqrt(2*min(irstrt-1,m-irstrt+1))*pnorm +
c                   sqrt(pnorm**2-1)))
c            if scase=2 this bound does not apply
c            if scase=3 this bound does not apply
c            if scase=4 this bound does not apply
c
c          lower bound on angular perturbation in right reducing subspace
c            if scase=1 (theorem 4, case 2 in paper)
c              atan(1/(sqrt(2*min(icstrt-1,n-icstrt+1))*qnorm +
c                   sqrt(qnorm**2-1)))
c            if scase=2 this bound does not apply
c            if scase=3 this bound does not apply
c            if scase=4 this bound does not apply
c
c         (note: given scase, sdlmax, pnorm, qnorm, m, n, icstrt, irstrt
c          and delta (the frobenius norm of a perturbation), subroutine
c          evalbd will compute the above upper and lower subspace bounds)
c
c       info - integer - 0 if normal return
c                        1 if svd error in difu calculation in pbound
c                        2 if difu=0 in pbound
c                        3 if svd error in difl calculation in pbound
c                        4 if difl=0 in pbound
c                        5 if multiple eigenvalues
c                        6 if inconsistent input dimensions
c
c     workspace
c       work(*) - complex*16 - exact amount is complicated function of 
c                 input dimensions and depends on ecase, and computed
c                 as follows:
c
c                    irend=irstrt+dimreg-1; icend=icstrt+dimreg-1;
c       if ecase=1 - m11=irstrt-1; m21=m-m11; n11=icstrt-1; n21=n-n11;
c                    m12=irend-irstrt+1; m22=m-irend; 
c                    n12=icend-icstrt+1; n22=n-icend;
c                    workspace = max( (2*n21*m11*(n11*n21+m11*m21+
c                                     2*n21*m11+2)+n11*n21+m11*m21) ,
c                                     (2*((m21*n11+1)*(n11*n21+
c                                     m11*m21+1)-1)) ,
c                                     (2*n22*m12*(n12*n22+m12*m22+
c                                     2*n22*m12+2)+n12*n22+m12*m22) ,
c                                     (2*((m22*n12+1)*(n12*n22+
c                                     m12*m22+1)-1)) )
c       if ecase=2 or ecase=5 - 
c                    m11=irstrt-1; m21=m-m11; n11=icstrt-1; n21=n-n11;
c                    workspace = max( (2*n21*m11*(n11*n21+m11*m21+
c                                    2*n21*m11+2)+n11*n21+m11*m21) ,
c                                    (2*((m21*n11+1)*(n11*n21+
c                                    m11*m21+1)-1)) )
c       if ecase=3 - m11=irend; m21=m-m11; n11=icend; n21=n-icend;
c                    workspace = max( (2*n21*m11*(n11*n21+m11*m21+
c                                    2*n21*m11+2)+n11*n21+m11*m21) ,
c                                    (2*((m21*n11+1)*(n11*n21+
c                                    m11*m21+1)-1)) )
c       if ecase=4 - workspace = n*n
c
c       the following simple expression bounds the workspace also, but
c          may occasionally be much too large (especially if ecase=4):
c            workspace .le. 2*m*n* (n*n + m*m + 2*n + m + 2) + n*n + m*m
c*********************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c    
c     addresses:
c             jim demmel, courant institute, 251 mercer str, 
c                 new york, new york 10012, usa
c                 electronic address: demmel at nyu.edu or
c                                     na.demmel at score.stanford.edu
c              bo kagstrom, institute of information processing,
c                 university of umea, s-90187 umea, sweden
c                 electronic address: bokg at seumdc51.bitnet or
c                                     na.kagstrom at score.stanford.edu
c
c**** bound uses the following functions and subroutines
c        pbound, ebdreg, cmatpr (debug only), gvec, dznrm2 (blas),
c        blddfu, blddfl, bldrhs, prml, prmlct, svdiv, zsvdc (linpack)
c 
c**** internal variables
      integer irend,icend,idummy,i
      real*8 rdummy, difu1, difu2, difl1, difl2, pnorm1, pnorm2
      real*8 qnorm1, qnorm2, pdelta1, pdelta2, delta
c
c     test input dimensions for consistency
      info = 0
      if (irstrt.gt.icstrt .or. irstrt.le.0 .or.
     +    n-icstrt-dimreg.gt.m-irstrt-dimreg .or.
     +    n-icstrt-dimreg+1.lt.0 .or. dimreg.lt.0) then
c       inconsistent input dimensions
        info = 6
        return
      endif
      icend = icstrt+dimreg-1
      irend = irstrt+dimreg-1
      delta = 0.
c
      if (dimreg.gt.0) then
c       there are eigenvalue bounds to compute
c
c       ecase 1 - in addition to selected regular part KCF has
c       (right singular part and/or regular part) and
c       (left singular part and/or regular part)   
        if (icstrt.ne.1 .and. irend.ne.m) then
          ecase = 1
          if (irstrt.eq.1) then
            scase = 2
          else
            scase = 1
          endif
c         see corollary 1 for explanation of bounds
          call pbound(a,b,ldab,m,n,irstrt-1,icstrt-1,
     +                delta,difl1,difu1,qnorm1,pnorm1, pdelta1,
     +                rdummy,rdummy,rdummy,rdummy,idummy,work,info)
          if (info.ne.0) return
          call pbound(a(irstrt,icstrt),b(irstrt,icstrt),ldab,
     +                m-irstrt+1,n-icstrt+1,irend-irstrt+1,
     +                icend-icstrt+1,
     +                delta,difl2,difu2,qnorm2,pnorm2,pdelta2,
     +                rdummy,rdummy,rdummy,rdummy,idummy,work,info)
          if (info.ne.0) return
          edlmax = min (pdelta1, pdelta2/(sqrt(2.)*qnorm1))
          pqnorm = 2.*pnorm2*qnorm1
c
          sdlmax = pdelta1
          pnorm = pnorm1
          qnorm = qnorm1
          difl = difl1
          difu = difu1
        endif
c
c       ecase 2 - in addition to selected regular part KCF has
c                (right singular part and/or regular part) and
c                (nothing)
        if (icstrt.ne.1 .and. irend.eq.m) then
          ecase=2
          if (irstrt.eq.1) then
            scase = 2
          else
            scase = 1
          endif
c         see part 1 of theorem 5 for explanation of bounds
          call pbound(a,b,ldab,m,n,irstrt-1,icstrt-1,delta,difl1,
     +                difu1,qnorm1,pnorm1,pdelta1,rdummy,rdummy,
     +                rdummy,rdummy,idummy,work,info)
          if (info.ne.0) return
          edlmax= pdelta1
          pqnorm=1.
          if (idummy.eq.1) pqnorm=sqrt(2.)*qnorm1
c
          sdlmax = pdelta1
          pnorm = pnorm1
          qnorm = qnorm1
          difl = difl1
          difu = difu1
        endif
c
c       ecase 3 - in addition to selected regular part KCF has
c                (nothing) and
c                (left singular part  and/or regular part)
        if (icstrt.eq.1 .and. irend.ne.m) then
          ecase = 3
          scase = 4
c         see part 2 of theorem 5 for explanation of bounds
          call pbound(a,b,ldab,m,n,irend,icend,
     +                delta,difl2,difu2,qnorm2,pnorm2,pdelta2,
     +                rdummy,rdummy,rdummy,rdummy,idummy,work,info)
          if (info.ne.0) return
          edlmax = pdelta2
          pqnorm = 1.
          if (idummy.eq.1) pqnorm = sqrt(2.)*pnorm2
          difl = 0.
          difu = 0.
          pnorm = 1.
          qnorm = 1.
          sdlmax = -1.
        endif
c
c       ecase 4 - pencil regular and entire spectrum selected
        if (icstrt.eq.1 .and. irend.eq.m) then
          ecase=4
          edlmax=-1.
          pqnorm=1.
c
          scase = 4
          difl = 0.
          difu = 0.
          pnorm = 1.
          qnorm = 1.
          sdlmax = -1.
        endif
c
        call ebdreg(a,b,ldab,irstrt,icstrt,dimreg,
     +              gvcond,evala,evalb,work,info)
        if (info.ne.0) then
          info = 5
          return
        endif
        if (pqnorm.ne.1.) then
          do 1 i=1,dimreg
            gvcond(i)=gvcond(i)*pqnorm
1         continue
        endif
c
      else
c       dimreg.eq.0, so only compute subspace bounds
        ecase = 5
        call pbound(a,b,ldab,m,n,irstrt-1,icstrt-1,
     +              delta,difl,difu,qnorm,pnorm,sdlmax,
     +              rdummy,rdummy,rdummy,rdummy,scase,work,info)
      endif
c
      if (idbg(20).ne.0) then
        write(outunit,100) ldab,m,n,irstrt,icstrt,dimreg,ecase,scase
100     format(' bound - ldab,m,n,irstrt,icstrt,dimreg,ecase,scase=',
     +         /,8i5)
        if (ecase.ne.5) then
          write(outunit,101) edlmax,pqnorm
101       format(' edlmax,pqnorm=',2d15.6,/,' gvcond=')
          write(outunit,102) (gvcond(i),i=1,dimreg)
102       format(5d15.6)
          call cmatpr(work,dimreg,dimreg,dimreg,'gvec')
        endif
        if (scase.ne.4) write(outunit,103) sdlmax,pnorm,qnorm
103     format(' sdlmax,pnorm,qnorm=',3d15.6)
      endif
      return
      end        
c
c
      subroutine ebdreg(a,b,ldab,irstrt,icstrt,dimreg,
     +                  gvcond,evala,evalb,work,info)
c     implicit none
c**** formal parameter declarations
      integer ldab, dimreg, irstrt, icstrt, info
      complex*16 a(ldab,*), b(ldab,*), work(*), evala(*), evalb(*)
      real*8 gvcond(*)
c     
c*****************************************************************
c
c     compute error bounds for eigenvalues of a regular pencil
c     requires all simple eigenvalues
c
c     inputs:
c       a(ldab,*), b(ldab,*) - complex*16 - contain pencil
c       irstrt, icstrt - integer - starting row and column locations
c                        of pencil within a and b
c       dimreg - integer - dimension of regular pencil
c 
c     outputs:
c       evala(dimreg), evalb(dimreg) - complex*16 - normalized 
c                        eigenvalues:
c                        evala(i)/evalb(i) is i-th eigenvalue and
c                        abs(evala(i))**2 + abs(evalb(i))**2 =1
c       gvcond(dimreg) - real*8 - gvcond(i) is condition number of 
c                 i-th eigenvalue where if the pencil is perturbed by 
c                 frobenius norm delta and the perturbed eigenvalue 
c                 is c/s where
c                 abs(c)**2 + abs(s)**2 = 1 then for some i
c                 abs(c*evalb(i) - s*evala(i)) .le. delta * gvcond(i)
c       info - integer - returns 0 (normal) if no multiple eigenvalues, 
c                  else nonzero
c
c     workspace:
c       work(dimreg**2) - complex*16 - work space
c
c***********************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
c**** ebdreg uses the following functions and subroutines:
c      gvec
c
c**** internal variables
c
      real*8 scl
      integer i
c
c     compute eigenvectors
      call gvec(a( irstrt , icstrt ),
     +     b( irstrt , icstrt ), ldab,
     +     dimreg , work, dimreg, gvcond, info)
c
c     compute normalized eigenpairs
      do 555 i=1,dimreg
        scl=sqrt(abs(a(irstrt-1+i,icstrt-1+i))**2+
     +           abs(b(irstrt-1+i,icstrt-1+i))**2)
        evala(i) = a(irstrt-1+i,icstrt-1+i)/scl
        evalb(i) = b(irstrt-1+i,icstrt-1+i)/scl
        if (info.eq.0) gvcond(i)= dimreg * gvcond(i) / scl
555   continue
c
      return
      end
c
c
      subroutine gvec(a,b,ldab,n,vec,ldvec,gvcond,info)
c
c     implicit none
c**** debug space
      common /debug2/ idbg(20),outunit
      integer idbg,outunit
      logical ldebug
c**** formal parameter declarations
      integer ldab, n, ldvec, info
      complex*16 a(ldab,*), b(ldab,n), vec(ldvec,*)
      real*8 gvcond(*)
c
c********************************************************************
c
c     compute the left and right eigenvectors of the upper triangular
c     regular pencil a - lambda b
c     compute condition numbers of eigenvalues
c
c     inputs
c       a(ldab,n),b(ldab,n) - complex*16 - n by n matrices
c       ldab - integer - leading dimension of a, b
c       n - integer - dimension of a, b
c       ldvec - integer - leading dimension of vec
c
c       idbg(10) - if idbg(10) ne 0, print debug output
c
c     outputs
c       vec(ldvec,n) - complex*16 -  matrix containing eigenvectors
c             vec(1:i,i) contains the right eigenvector of the i-th
c               eigenvalue, normalized so vec(i,i)=1. the other
c               components of the eigenvector are zero
c             vec(i:n,i) contains the left eigenvector of the i-th
c               eigenvalue, normalized so vec(i,i)=1. the other 
c               components of the eigenvector are zero
c       gvcond(n) - real*8 - array of condition numbers of eigenvalues.
c                if right eigenvectors scaled by diagonal matrix d
c                to have unit norm, scale left eigenvectors by d**-1.
c                then condition number is norm of left eigenvector.
c       info - integer - 0 if pencil regular without multiple eigenvalues
c              nonzero index of a multiple or 0/0 eigenvalue otherwise

c***********************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
c**** gvec uses the following external function:
c     dznrm2 (blas)
      real*8 dznrm2
c**** internal variables
      integer nm1, i, im1, im2, j, jp1, k, ip1, ip2, jm1
      complex*16 alpha, beta, diag, cmul, csum
      real*8 ca, cb, dmax, dmin, d 
c
      ldebug=(idbg(10).ne.0)
      info=0
      nm1=n-1
      if (ldebug) write(outunit,99)
99    format(' entering gvec')
      do 1 i=1,n
c
        if (ldebug) write(outunit,100) i
100     format(' i=',i4)
        vec(i,i)=1.
c
c       compute alpha, beta so that zz = beta*a - alpha*b is a
c       singular matrix whose left and right null spaces are the
c       left and right eigenspaces we seek
        ca=abs(a(i,i))
        cb=abs(b(i,i))
        dmax=max(ca,cb)
        if (ldebug) write(outunit,101) a(i,i),b(i,i),ca,cb,dmax
101     format(' a(i,i)=',2d20.5,/,' b(i,i)=',2d20.5,/,' ca=',d20.5,/,
     +  ' cb=',d20.5,/,' dmax=',d20.5)
        if (dmax.eq.0.0) then
c         singular pencil
          info=i
          return
        endif 
        dmin=min(ca,cb)
        d=dmax*sqrt(1+(dmin/dmax)**2)
        alpha = a(i,i)/d
        beta = b(i,i)/d
        if (ldebug) write(outunit,102) dmin,d,alpha,beta
102     format(' dmin=',d20.5,/,' d=',d20.5,/,' alpha=',2d20.5,/,
     +  ' beta=',2d20.5)
c
c       compute right eigenvector
        if (i.ne.1) then
c
c         solve zz(1:i-1,1:i-1) * x = -zz(1:i-1,i) for
c         x = vec(1:i-1,i)
          diag=beta*a(i-1,i-1) - alpha*b(i-1,i-1)
          im1=i-1
          if (ldebug) write(outunit,103) im1,i,diag
103       format(' i,j,diag=',2i4,2d20.5)
          if (abs(diag).eq.0.0) then
c           multiple eigenvalue
            info=i-1
            return
          endif
          vec(i-1,i)=-(beta*a(i-1,i)-alpha*b(i-1,i))/diag
          if (i.ne.2) then
            im1=i-1
            im2=i-2
            do 2 j=im2,1,-1
              diag=beta*a(j,j)-alpha*b(j,j)
              if (ldebug) write(outunit,103) j,i,diag
              if (abs(diag).eq.0.0) then
c               multiple eigenvalue
                info=j
                return
              endif
              csum=-(beta*a(j,i)-alpha*b(j,i))
              jp1=j+1
              do 3 k=jp1,im1
                cmul=beta*a(j,k)-alpha*b(j,k)
                csum=csum-cmul*vec(k,i)
3             continue
              vec(j,i)=csum/diag
2           continue
          endif
        endif
c
c       compute left eigenvector
        if (i.ne.n) then
c         solve xt * zz(i+1:n,i+1:n) = -zz(i,i+1:n) for
c         x = vec(i+1:n,i)
          diag=beta*a(i+1,i+1)-alpha*b(i+1,i+1)
          ip1=i+1
          if (ldebug) write(outunit,103) i,ip1,diag
          if (abs(diag).eq.0.0) then
c           multiple eigenvalue
            info=i
            return
          endif
          vec(i+1,i)=-(beta*a(i,i+1)-alpha*b(i,i+1))/diag
          if (i.ne.nm1) then
            ip1=i+1
            ip2=i+2
            do 4 j=ip2,n
              diag=beta*a(j,j)-alpha*b(j,j)
              if (ldebug) write(outunit,103) i,j,diag
              if (abs(diag).eq.0.0) then
c               multiple eigenvalue
                info=i
                return
              endif
              csum=-(beta*a(i,j)-alpha*b(i,j))
              jm1=j-1
              do 5 k=ip1,jm1
                cmul=beta*a(k,j)-alpha*b(k,j)
                csum=csum-cmul*vec(k,i)
5             continue
              vec(j,i)=csum/diag
4           continue
          endif
        endif
1     continue
c
c     compute condition numbers
      do 6 i=1,n
        gvcond(i)=dznrm2(i,vec(1,i),1)*dznrm2(n-i+1,vec(i,i),1)
6     continue
      return
      end
c
      subroutine pbound(a,b,ldab,m,n,rowred,colred,delta,difl,difu,
     +    qnorm,pnorm,pdelta,lbndup,rbndup,lbndlw,rbndlw,scase,work,
     +    ierr)
c
c     implicit none
c
c**** formal parameter declarations
      integer ldab,m,n,rowred,colred,ierr,scase
      complex*16 a(ldab,*),b(ldab,*),work(*)
      real*8 delta,difl,difu,qnorm,pnorm,pdelta,lbndup,rbndup
      real*8 lbndlw, rbndlw
c
c*******************************************************************
c
c     compute perturbation bounds for reducing subspaces of
c     singular pencil a - lambda b
c     assume a - lambda b has been reduced to generalized upper
c     triangular form by guptri
c     need rowred .le. colred and n-colred .le. m-rowred
c       as implied by generalized upper triangular form
c
c     there are 4 cases, depending on dimension:
c
c      case 1: 0 .lt. rowred and 0 .lt. n-colred so that
c        both left and right reducing subspaces nontrivial
c
c      case 2: if rowred=0 and 0 .lt. colred .lt. n then left reducing
c        subspace 0 but right one nontrivial and bounds exist for it
c
c      case 3: if colred=n and 0 .lt. rowred .lt. m then right reducing
c        subspace is entire space but left one nontrivial with bounds
c
c      case 4: if ( (rowred=0 and colred=0) or
c                   (rowred=0 and colred=n) or
c                   (rowred=m and colred=n) ) then
c              both left and right subspaces trivial
c
c     inputs:
c
c       a(ldab,n),b(ldab,n) - complex*16 - m by n matrices
c
c       ldab - integer - leading dimension of a and b
c
c       m,n - integer - dimensions of a and b
c
c       rowred,colred - integer - number of rows and columns in 
c             (1,1) position of a,b.  dimensions of desired left 
c             and right reducing subspaces
c
c       delta - real*8 - distance of perturbed pencil from a - lambda b 
c
c       idbg(9) - integer - if idbg(9) ne 0, print debug output
c
c     outputs: (described in more detail in 
c       'accurate solutions of ill-posed problems in control theory'
c       25th conference on decision and control, 
c       j. demmel and b. kagstrom
c
c       difl - real*8 - difl function (in case 4, difl=0)
c
c       difu - real*8 - difu function (in case 4, difu=0)
c
c       qnorm - real*8 - right projector norm ( sqrt(r0**2+1) )
c                        (in case 4, qnorm=1.)
c
c       pnorm - real*8 - left projector norm ( sqrt(l0**2+1) )
c                        (in case 4, prnorm=1.)
c
c       pdelta - real*8 - radius of ball around a - lambda b within 
c                which perturbation bounds hold (in case 4, pdelta=-1.
c                to show pdelta does not apply). if delta.ge.pdelta, 
c                the following bounds are set to -1. the following
c                outputs are given in terms of relerr = delta/pdelta
c
c       lbndup - real*8 - upper bound on angular perturbation in left 
c                reducing subspace (case 1 of theorem 4 of above paper)
c                 in case 1: 
c                  lbndup=atan(relerr/(pnorm-relerr*sqrt(pnorm**2-1)))
c                 in case 2:
c                  lbndup=0
c                 in case 3:
c                  lbndup=atan(relerr/(1-relerr))
c                 in case 4: 
c                  lbndup=0
c
c       rbndup - real*8 - upper bound on angular perturbation in right 
c                reducing subspace (case 1 of theorem 4)
c                 in case 1:
c                  rbndup=atan(relerr/(qnorm-relerr*sqrt(qnorm**2-1)))
c                 in case 2:
c                  rbndup=atan(relerr/(1-relerr))
c                 in case 3:
c                  rbndup=0
c                 in case 4: 
c                  rbndup=0
c
c       lbndlw - real*8 - lower bound on angular perturbation in left 
c                reducing subspace (case 2 of theorem 4)
c                 in case 1:
c                  lbndlw=atan(1/(sqrt(2*min(rowred,m-rowred))*pnorm +
c                         sqrt(pnorm**2-1)))
c                 in case 2: lbndlw=-1 since this bound does not apply
c                 in case 3: lbndlw=-1 since this bound does not apply
c                 in case 4: lbndlw=-1 since this bound does not apply
c
c       rbndlw - real*8 - lower bound on angular perturbation in right
c                reducing subspace (case 2 of theorem 4)
c                 in case 1:
c                  rbndlw=atan(1/(sqrt(2*min(colred,n-colred))*qnorm +
c                         sqrt(qnorm**2-1)))
c                 in case 2: rbndlw=-1 since this bound does not apply
c                 in case 3: rbndlw=-1 since this bound does not apply
c                 in case 4: rbndlw=-1 since this bound does not apply
c
c       scase - integer - 1, 2, 3 or 4 as described above
c
c       ierr - integer - error flag
c              0 means no error (normal return)
c              1 means error in svd of difu
c              2 means difu = 0
c              3 means error in svd of difl
c              4 means difl = 0
c              5 means bad rowred or colred
c
c     work space
c       work - complex*16 - array of length at least
c              max ( rowdfu*coldfu+coldfu**2+2*coldfu+rowdfu ,
c                    rowdfl*coldfl+2*coldfl+rowdfl )
c            where
c              rowdfu=coldfl=colred*(n-colred)+rowred*(m-rowred)
c              coldfu=2*(n-colred)*rowred
c              rowdfl=2*(m-rowred)*colred
c
c*********************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel, courant institute, 251 mercer str, new york,  
c                 new york, 10012
c                 electronic address: demmel at nyu.edu
c              bo kagstrom, institute of information processing,
c                 university of umea, s-90187 umea, sweden
c                 electronic address: bokg at seumdc51.bitnet
c
c**** pbound uses the following subroutines and functions
c     dznrm2, blddfu, blddfl, bldrhs, prml, prmlct, svdiv, zsvdc
c
c**** internal variables
c
      complex*16 dummy
      integer rowdfu,coldfu,sstrt,wstrt,estrt,rowdfl,coldfl,vstrt
      integer isub, i, j, info, len
      real*8 r0, l0, relerr, dznrm2
c
      ierr=0
      if ((rowred.gt.colred).or.((n-colred).gt.(m-rowred))) then
c       inconsistent dimensions
        ierr = 5
      elseif ((0.lt.rowred) .and. (0.lt.n-colred)) then
c       case 1
        scase = 1
c       compute difu
c       build transposed difu matrix starting at work(1)
c       rowdfu = number of rows in difut
        rowdfu=colred*(n-colred)+rowred*(m-rowred)
c       coldfu = number of columns in difut
        coldfu=2*(n-colred)*rowred
c
        call blddfu(work,rowdfu,a,b,ldab,m,n,rowred,colred)
c
c       setup workspace for svd
c       store left singular vectors u over difu starting at work(1)
        sstrt=1+rowdfu*coldfu
c       store singular values starting at work(sstrt)
        wstrt=sstrt+coldfu
c       store work array needed for svd starting at work(wstrt)
        estrt=wstrt+rowdfu
c       store e array needed for svd starting at work(estrt)
        vstrt=estrt+coldfu
c       store right singular vectors v starting at work(vstrt)
c
c       compute svd
        call zsvdc(work(1),rowdfu,rowdfu,coldfu,work(sstrt),
     +    work(estrt),work(1),rowdfu,work(vstrt),coldfu,work(wstrt),
     +    21,info)
c
        if (info.eq.0) goto 10
          ierr=1
          return
10      continue
c
c       extract difu
        difu=dreal(work(sstrt-1+coldfu))
c
        if (difu.gt.0.) goto 20
          ierr=2
          return
20      continue
c
c       compute pnorm, qnorm
c       build rhs = (-col a12, -col b12) starting at work(wstrt)
        call bldrhs(work(wstrt),a,b,ldab,m,n,rowred,colred)
c
c       solve underdetermined least squares problem
c       premultiply rhs by v* storing result at work(estrt)
        call prmlct(work(vstrt),coldfu,coldfu,coldfu,
     +              work(wstrt),work(estrt))
c
c       premultiply by inverted singular values
        call svdiv(work(estrt),coldfu,work(sstrt))
c
c       premultiply by u storing result at work(wstrt)
        call prml(work,rowdfu,rowdfu,coldfu,work(estrt),work(wstrt))
c
        len=colred*(n-colred)
c       compute r0 = norm of leading len components
        r0=dznrm2(len,work(wstrt),1)
c
c       compute l0 = norm of remaining components
        len=rowred*(m-rowred)
        l0=dznrm2(len,work(wstrt+len),1)
c       compute pnorm, qnorm from l0, r0
        pnorm=sqrt(1+l0**2)
        qnorm=sqrt(1+r0**2)
c
c       compute difl
c       build difl matrix starting at work(1)
c       rowdfl = number of rows in difl
        rowdfl=2*colred*(m-rowred)
c       coldfl=number of columns in difl
        coldfl=rowred*(m-rowred)+colred*(n-colred)
        call blddfl(work,rowdfl,a,b,ldab,m,n,rowred,colred)
c
c       setup workspace for svd
c       do not compute any singular vectors
        sstrt=1+rowdfl*coldfl
c       store singular values starting at work(sstrt)
        wstrt=sstrt+coldfl
c       store work array needed by svd starting at work(wstrt)
        estrt=wstrt+rowdfl
c       store e array needed by svd starting at work(estrt)
c
        call zsvdc(work(1),rowdfl,rowdfl,coldfl,work(sstrt),
     +             work(estrt),dummy,1,dummy,1,work(wstrt),0,info)
c
        if (info.eq.0) goto 30
          ierr=3
          return
30      continue
c
c       extract difl
        difl=dreal(work(sstrt-1+coldfl))
        if (difl.gt.0.) goto 40
          ierr=4
          return
40      continue
c       compute perturbation bounds
        pdelta=min(difl,difu)/(sqrt(pnorm**2+qnorm**2)+
     +         2.*max(pnorm,qnorm))
        relerr=delta/pdelta
        lbndup=-1.
        rbndup=-1.
        lbndlw=-1.
        rbndlw=-1.
        if (relerr.ge.1.) goto 50
          lbndup=atan(relerr/(pnorm-relerr*sqrt(pnorm**2-1.)))
          rbndup=atan(relerr/(qnorm-relerr*sqrt(qnorm**2-1.)))
          lbndlw=atan(1./(sqrt(2.*min(rowred,m-rowred))*pnorm+
     +           sqrt(pnorm**2-1.)))
          rbndlw=atan(1./(sqrt(2.*min(colred,n-colred))*qnorm+
     +           sqrt(qnorm**2-1.)))
50      continue
      elseif (rowred.eq.0.and.colred.gt.0.and.colred.lt.n) then
c       case 2
        scase = 2
c       compute difl
c       build difl matrix ( (a**t b**t)**t ) starting at work(1)
        isub = 0
        do 100 j=colred+1, n
          do 101 i=1, m
            isub = isub +1
            work(isub) = a(i,j)
101       continue
          do 102 i=1,m
            isub = isub +1
            work(isub) = b(i,j)
102       continue
100     continue
c       compute singular values
        sstrt=1+isub
        estrt=sstrt + n-colred
        wstrt=estrt + n-colred
        call zsvdc(work,2*m,2*m,n-colred,work(sstrt),work(estrt),
     +             dummy,1,dummy,1,work(wstrt),0,info)
        if (info.ne.0) then
          ierr=3
          return
        endif
c       extract difl
        difl = abs(work(sstrt+n-colred-1))
        difu=difl
        if (difl.eq.0.) then
           ierr=4
           return
        endif
        pdelta=difl
        relerr=delta/pdelta
        pnorm = 1.
        qnorm = 1.
        lbndlw = -1.
        rbndlw = -1.
        lbndup = -1.
        rbndup = -1.
        if (relerr.lt.1.) then
          lbndup = 0.
          rbndup = atan(relerr/(1.-relerr))
        endif
      elseif (colred.eq.n.and.rowred.gt.0.and.rowred.lt.m) then
c       case 3
        scase = 3
c       compute difu
c       build difu matrix (a,b) starting at work(1)
        isub = 0
        do 104 j=1,n
          do 105 i=1,rowred
            isub = isub +1
            work(isub) = a(i,j)
105       continue
104     continue
        do 106 j=1,n
          do 107 i=1,rowred
            isub = isub +1
            work(isub) = b(i,j)
107       continue
106     continue
c       compute singular values
        sstrt=isub+1
        estrt=sstrt+rowred+1
        wstrt=estrt+2*n
        call zsvdc(work,rowred,rowred,2*n,work(sstrt),work(estrt),
     +             dummy,1,dummy,1,work(wstrt),0,info)
        if (info.ne.0) then
          ierr = 1
          return
        endif
c       extract difu
        difu=abs(work(sstrt+rowred-1))
        difl = difu
        if (difu.eq.0.0) then
          ierr = 2
          return
        endif
        pdelta = difu
        relerr = delta/pdelta
        pnorm = 1.
        qnorm = 1.
        lbndup = -1.
        rbndup = -1.
        lbndlw = -1.
        rbndlw = -1.
        if ( relerr.lt.1.0) then
          rbndup = 0.
          lbndup = atan(relerr/(1.-relerr))
        endif
      else
c       both left and right subspace trivial
        scase = 4
        lbndup = 0.
        rbndup = 0.
        lbndlw = -1.
        rbndlw = -1.        
        difl = 0.
        difu = 0.
        pdelta = -1.
        pnorm = 1.
        qnorm = 1.
      endif
      return
      end
c
c
      subroutine blddfl(work,wrow,a,b,ldab,m,n,rowred,colred)
c     implicit none
c**** formal parameter declarations
      integer ldab, m, n, rowred, colred, wrow
      complex*16 work(wrow,*),a(ldab,*),b(ldab,*)
c
c***************************************************************
c
c     build difl matrix in work
c     in matlab notation
c
c     difl matrix = < <a11' .*. eye(m-rowred) , -eye(colred) .*. a22 >;
c                     <b11' .*. eye(m-rowred) , -eye(colred) .*. b22 >>
c
c     where a11 = a(1:rowred , 1:colred) 
c           a22 = a(rowred+1 : m , colred+1 : n)
c           b11 = b(1:rowred , 1:colred)
c           b22 = b(rowred+1 : m , colred+1 : n)
c
c***************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
c**** internal variables
      integer wcol,rstrta,rstrtb,cstrt,cnt,i,j
      integer row12,col1,col2,mmrwrd,nmclrd
c
c     nmclrd = number of columns in (1,2), (2,2) blocks of a, b
      nmclrd = n-colred
c     mmrwrd = number of rows in (2,1), (2,2) blocks of a, b
      mmrwrd = m-rowred
c     row12 = numbers of rows in each subblock of difl matrix
      row12 = colred*mmrwrd
c     col1 = number of columns in (1,1), (2,1) blocks of difl
      col1 = rowred*mmrwrd
c     col2 = number of columns in (1,2), (2,2) blocks of difl
      col2 = colred*nmclrd
c     wcol = total number of columns in difl
      wcol = col1+col2
c
c     zero out difl
      do 10 j=1,wcol
        do 11 i=1,wrow
          work(i,j)=0.
11      continue
10    continue
c
c     fill in (1,1), (2,1) blocks of difl
      rstrta=0
      rstrtb=row12
      cstrt=0
      do 1 j=1,colred
        do 2 i=1,rowred
          do 3 cnt=1,mmrwrd
            work(cnt+rstrta,cnt+cstrt)=a(i,j)
            work(cnt+rstrtb,cnt+cstrt)=b(i,j)
3         continue
          cstrt=cstrt+mmrwrd
2       continue
        cstrt=0
        rstrta=rstrta+mmrwrd
        rstrtb=rstrta+row12
1     continue
c
c     fill in (1,2), (2,2) blocks of difl
      rstrta=0
      cstrt=col1
      do 4 cnt=1,colred
        rstrtb=rstrta+row12
        do 5 j=1,nmclrd
          do 6 i=1,mmrwrd
            work(rstrta+i,cstrt+j)=-a(i+rowred,j+colred)
            work(rstrtb+i,cstrt+j)=-b(i+rowred,j+colred)
6         continue
5       continue
        rstrta=rstrta+mmrwrd
        cstrt=cstrt+nmclrd
4     continue
      return
      end
c
c
      subroutine blddfu(work,wrow,a,b,ldab,m,n,rowred,colred)
c     implicit none
c**** formal parameter declarations
      integer ldab, m, n, rowred, colred, wrow
      complex*16 work(wrow,*),a(ldab,*),b(ldab,*)
c*********************************************************************
c
c     build conjugate transpose difu matrix in work
c     in matlab notation
c
c     (difu matrix)' =
c
c       < < eye(n-colred) .*. a11' , eye(n-colred) .*. b11' >;
c         < -conj(a22) .*. eye(rowred) , -conj(b22) .*. eye(rowred) >>
c
c     where a11 = a(1:rowred , 1:colred) 
c           a22 = a(rowred+1 : m , colred+1 : n)
c           b11 = b(1:rowred , 1:colred)
c           b22 = b(rowred+1 : m , colred+1 : n)
c
c*********************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
c**** internal variables
c
      integer wcol,cstrta,cstrtb,rstrt,cnt,i,j
      integer mmrwrd,nmclrd,rwrdp1,clrdp1
      integer row1, row2, col12
c
c     nmclrd = number of columns in (1,2), (2,2) entries of a, b
      nmclrd=n-colred
c     col12 = number of columns in each subblock of difuct matrix
      col12=rowred*nmclrd
c     mmrwrd = number of rows in (2,1), (2,2) entries of a, b
      mmrwrd = m-rowred
c     row1 = number of rows in (1,1), (2,1) sublocks of difu
      row1 = colred*nmclrd
c     row2 = number of rows in (1,2), (2,2) subblocks of difu
      row2 = rowred*mmrwrd
c     wcol = total number of columns in difu matrix
      wcol = 2*col12
c     initialize difu to zero
      do 1 j=1,wcol
        do 2 i=1,wrow
          work(i,j)=0.
2       continue
1     continue
c
c     fill in (1,1), (1,2) positions of difu
      cstrta=0
      rstrt=0
      do 3 cnt=1,nmclrd
        cstrtb=cstrta+col12
          do 4 j=1,colred
            do 5 i=1,rowred
              work(rstrt+j,cstrta+i)=conjg(a(i,j))
              work(rstrt+j,cstrtb+i)=conjg(b(i,j))
5           continue
4         continue
        cstrta=cstrta+rowred
        rstrt=rstrt+colred
3     continue
c
c     fill in (2,1), (2,2) positions of difuct
      rwrdp1=rowred+1
      clrdp1=colred+1
      cstrta=0
      cstrtb=col12
      rstrt=row1
      do 6 j=clrdp1,n
        do 7 i=rwrdp1,m
          do 8 cnt=1,rowred
            work(cnt+rstrt,cnt+cstrta)=-conjg(a(i,j))
            work(cnt+rstrt,cnt+cstrtb)=-conjg(b(i,j))
8         continue
          rstrt=rstrt+rowred
7       continue
        rstrt=row1
        cstrta=cstrta+rowred
        cstrtb=cstrta+col12
6     continue
      return
      end
c
c
      subroutine bldrhs(work,a,b,ldab,m,n,rowred,colred)
c     implicit none
c**** formal parameter declarations
      integer ldab, m, n, rowred, colred
      complex*16 work(*), a(ldab,*), b(ldab,*)
c
c*********************************************************************
c
c     extract a12 = (1,2) block of a and b12 = (1,2) block of b
c     and store columnwise in work=(-col a12, -col b12)
c
c*********************************************************************
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
c**** internal variables
      integer clrdp1, j, i, loc
c
      clrdp1=colred+1
      loc=0
      do 1 j=clrdp1,n
        do 2 i=1,rowred
          loc=loc+1
          work(loc)=-a(i,j)
2       continue
1     continue
      do 3 j=clrdp1,n
        do 4 i=1,rowred
          loc=loc+1
          work(loc)=-b(i,j)
4       continue
3     continue
      return
      end
c
c
      subroutine prml(u,ldu,m,n,rhs,prod)
c     implicit none
      integer ldu, m, n
      complex*16 u(ldu,n),rhs(n),prod(m)
c
c*********************************************************************
c     compute prod = u * rhs
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
      integer i, j
c
      do 1 j=1,m
        prod(j)=rhs(1)*u(j,1)
1     continue
      if (n.eq.1) return
      do 2 i=2,n
        call zaxpy(m,rhs(i),u(1,i),1,prod,1)
2     continue
      return
      end
c
c
      subroutine prmlct(u,ldu,m,n,rhs,prod)
c     implicit none
      integer ldu, m, n
      complex*16 u(ldu,n),rhs(m),prod(n),zdotc
c
c*********************************************************************
c     compute prod = (conjugate transpose u) * rhs
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
      integer j
c
      do 1 j=1,n
        prod(j)=zdotc(m,u(1,j),1,rhs,1)
1     continue
      return
      end
c
c
      subroutine svdiv(z,n,s)
c     implicit none
      integer n
      complex*16 z(n),s(n)
c
c*********************************************************************
c     divide one array by another
c
c**** this version dated 16 june 1987
c     authors: jim demmel and bo kagstrom
c
      integer j
c
      do 1 j=1,n
        z(j)=z(j)/s(j)
1     continue
      return
      end
c
      subroutine evalbd(delta, sdlmax, qnorm, pnorm, scase,
     +                  m, n, irstrt, icstrt, 
     +                  lbndup, rbndup, lbndlw, rbndlw)
c
c     implicit none
c**** formal parameter declarations
c
      real*8 delta, sdlmax, qnorm, pnorm
      real*8 lbndup, rbndup, lbndlw, rbndlw
      integer scase, m, n, irstrt, icstrt
c
c******************************************************************
c
c     evaluate reducing subspace angular perturbation bounds computed
c     by subroutine bound for a perturbation of frobenius
c     norm delta. see documentation to subroutine bound for more details.
c
c     inputs:
c
c       sdlmax, qnorm, pnorm and scase are computed by bound. 
c       m, n, irstrt and icstrt are dimensions also input to bound
c       in order to compute sdlmax, qnorm, pnorm and scase.
c
c     outputs:
c
c       lbndup - real*8 - upper bound on angular perturbation in 
c                         left reducing subspace 
c                         (0 if space trivial and -1 if inapplicable)
c
c       rbndup - real*8 - upper bound on angular perturbation in
c                         right reducing subspace 
c                         (0 if space trivial and -1 if inapplicable)
c
c       lbndlw - real*8 - lower bound on angular perturbation in
c                         left reducing subspace (-1 if inapplicable)
c
c       rbndlw - real*8 - lower bound on angular perturbation in
c                         right reducing subspace (-1 if inapplicable)
c
c************************************************************************
c
c**** this version dated 16 june 87
c     authors: jim demmel and bo kagstrom
c
c**** internal variables
      real*8 relerr
c
      if (scase .ne. 4) relerr = delta/sdlmax
      if (scase.eq.1) then
        lbndup = atan(relerr/(pnorm-relerr*sqrt(pnorm**2-1.)))
        rbndup = atan(relerr/(qnorm-relerr*sqrt(qnorm**2-1.)))
        lbndlw = atan(1./(sqrt(2.*min(irstrt-1,m-irstrt+1))*pnorm +
     +           sqrt(pnorm**2-1.)))
        rbndlw = atan(1./(sqrt(2.*min(icstrt-1,n-icstrt+1))*qnorm +
     +           sqrt(qnorm**2-1.)))
      elseif (scase.eq.2) then
        lbndup = 0.
        rbndup = atan(relerr/(1.-relerr))
        lbndlw = -1.
        rbndlw = -1.
      elseif (scase.eq.3) then
        lbndup = atan(relerr/(1.-relerr))
        rbndup = 0.
        lbndlw = -1.
        rbndlw = -1.
      elseif (scase.eq.4) then
        lbndup = 0.
        rbndup = 0.
        lbndlw = -1.
        rbndlw = -1.
      endif
      return
      end
c
      subroutine bndwsp(m,n,irstrt,icstrt,dimreg,ecase,space,info)
c
c     implicit none
c
c**** debug space
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
c**** formal parameter declarations
      integer m,n,irstrt,icstrt,dimreg,info,ecase,space
c
c********************************************************************
c
c     compute work space needed by subroutine bound
c
c     inputs
c
c       m,n - integer - row, column dimensions of a and b
c
c       irstrt, icstrt - integer - starting row and column of selected 
c                        part of pencil for which eigenvalue bounds 
c                        are desired. reducing subspace bounds will be
c                        supplied for right reducing subspace spanned
c                        by leading icstrt-1 components and for left
c                        reducing subspace spanned by leading icstrt-1
c                        components.
c                        note: set icstrt=n+1 to make right reducing
c                                  subspace whole space
c                              set irstrt=m+1 to make left reducing
c                                  subspace whole space
c
c       dimreg - integer - number of selected eigenvalues;
c         if dimreg.eq.0 only subspace perturbation bounds will be
c         computed
c        (note - one can select a subset of the regular part only;
c         this gives generally different bounds for common eigenvalues
c         from a different selected subset; see paper above for 
c         discussion)
c
c     outputs
c
c       ecase - integer - which of 5 cases for eigenvalue bounds 
c               the pencil falls depending on input dimensions;
c               the first four cases are for dimreg.gt.0, in which
c               case the description gives:
c                  (part of KCF to above, left of selected part) and
c                  (part of KCF to below, right of selected part) 
c          ecase=1 - (right singular and/or regular part) and
c                    (left singular and/or regular part)
c          ecase=2 - (right singular and/or regular part) and (nothing)
c          ecase=3 - (nothing) and (left singular and/or regular part)
c          ecase=4 - (nothing) and (nothing)
c          ecase=5 - dimreg.eq.0 (no eigenvalue bounds)
c
c       space - integer - amount of workspace (double precision complex
c                         words) needed by subroutine bound
c       (the following simple expression bounds the workspace also, but
c          may occasionally be much too large (especially if ecase=4):
c            workspace .le. 2*m*n* (n*n + m*m + 2*n + m + 2) + n*n + m*m)
c
c       info - integer - 0 if normal return
c                        1 if inconsistent input dimensions
c
c*************************************************************************
c
c**** this version dated 22 june 1987
c     authors: jim demmel, courant institute, 251 mercer str, 
c                 new york, new york, 10012
c                 electronic address: demmel at nyu.edu
c              bo kagstrom, institute of information processing,
c                 university of umea, s-90187 umea, sweden
c                 electronic address: bokg at seumdc51.bitnet 
c
c**** internal variables
      integer irend,icend,m11,m21,m12,m22,n11,n12,n21,n22
c
c     test input dimensions for consistency
      info = 0
      icend = icstrt+dimreg-1
      irend = irstrt+dimreg-1
      if (irstrt.gt.icstrt .or. irstrt.le.0 .or.
     +    n-icstrt-dimreg.gt.m-irstrt-dimreg .or.
     +    n-icstrt-dimreg+1.lt.0 .or. dimreg.lt.0) then
c       inconsistent input dimensions
        info = 1
      else
        if (dimreg.gt.0) then
c         there are eigenvalue bounds to compute
c
c         ecase 1 - in addition to selected regular part KCF has
c         (right singular part and/or regular part) and
c         (left singular part and/or regular part)   
          if (icstrt.ne.1 .and. irend.ne.m) then
            ecase = 1
          endif
c
c         ecase 2 - in addition to selected regular part KCF has
c                  (right singular part and/or regular part) and
c                  (nothing)
          if (icstrt.ne.1 .and. irend.eq.m) then
            ecase=2
          endif
c
c         ecase 3 - in addition to selected regular part KCF has
c                  (nothing) and
c                  (left singular part  and/or regular part)
          if (icstrt.eq.1 .and. irend.ne.m) then
            ecase = 3
          endif
c
c         ecase 4 - pencil regular and entire spectrum selected
          if (icstrt.eq.1 .and. irend.eq.m) then
            ecase=4
          endif
c
        else
c         dimreg.eq.0, so only compute subspace bounds
          ecase = 5
        endif
c
        if (ecase .eq. 1) then
          m11=irstrt-1
          m21=m-m11
          n11=icstrt-1
          n21=n-n11
          m12=irend-irstrt+1
          m22=m-irend
          n12=icend-icstrt+1 
          n22=n-icend
          space = max( (2*n21*m11*(n11*n21+m11*m21+
     +                  2*n21*m11+2)+n11*n21+m11*m21) ,
     +                  (2*((m21*n11+1)*(n11*n21+
     +                  m11*m21+1)-1)) ,
     +                  (2*n22*m12*(n12*n22+m12*m22+
     +                  2*n22*m12+2)+n12*n22+m12*m22) ,
     +                  (2*((m22*n12+1)*(n12*n22+
     +                  m12*m22+1)-1)) )
        elseif (ecase .eq. 2 .or. ecase .eq. 5) then
          m11=irstrt-1
          m21=m-m11
          n11=icstrt-1
          n21=n-n11
          space = max( (2*n21*m11*(n11*n21+m11*m21+
     +                 2*n21*m11+2)+n11*n21+m11*m21) ,
     +                 (2*((m21*n11+1)*(n11*n21+
     +                 m11*m21+1)-1)) )
        elseif (ecase .eq. 3) then
          m11=irend
          m21=m-m11
          n11=icend
          n21=n-icend
          space = max( (2*n21*m11*(n11*n21+m11*m21+
     +                 2*n21*m11+2)+n11*n21+m11*m21) ,
     +                 (2*((m21*n11+1)*(n11*n21+
     +                 m11*m21+1)-1)) )
        elseif (ecase .eq. 4) then
          space = n*n
        endif
      endif
c
      if (idbg(19).ne.0) then
        write(outunit,100) m,n,irstrt,icstrt,dimreg,ecase,
     +  space,info
100     format(' bndwsp - m,n,irstrt,icstrt,dimreg'
     +         ',ecase,space,info=',/,8i5)
      endif
      return
      end        

c   on this file june 7, 1987: cmatml, cmatmr
c
      subroutine cmatml(a,lda,rowa,cola,b,ldb,rowb,c,ldc,work,job)
c
c     implicit none
      integer       lda,rowa,cola,ldb,rowb,ldc,job
      complex*16    a(lda,lda),b(ldb,ldb),c(ldc,ldc),work(*)
      complex*16    zdotu,zdotc
c
c***********************************************************************
c
c     cmatml performs  a complex (left) matrix multiplication  b * a,
c     or b' * a (' = transpose ,conjugate) where a is rowa * cola,
c     b is rowb * rowa. the result is stored in c or overwritten in a.
c     note the extra restrictions on dimensions of b when job = 3 or 4.
c
c     on entry
c         
c         a         complex(lda,cola), where lda>=rowa.
c
c         lda       integer
c                   lda is the leading dimension of the array a.
c                    
c         rowa      integer
c                   rowa is the number of rows of a, which is also
c                   the number of columns of b.
c         cola      integer
c                   cola is the number of columns of a, which is also
c                   the number columns of the resulting matrix.
c
c         b         complex(ldb,rowa), ldb>=rowb.
c                   
c         ldb       integer
c                   ldb is the leading dimension of the array b.
c                                           
c         rowb      integer
c                   rowb is the number of rows of the array, which
c                   is also the number of rows of the resulting matrix.
c
c         ldc       integer
c                   ldc is the leading dimension of the array c
c
c         work      complex(rowa)
c                   work is a scratch array.
c
c         job       integer
c                   job controls the matrix multiplication, and has
c                   the following meaning
c                   job=1       a = b * a
c                   job=2       c = b * a
c                   job=3       a = b' * a
c                   job=4       c = b' * a
c
c    on return
c
c         c         complex(ldc,cola), where ldc>=rowb.
c                   c is the matrix product of a and b. if rowa (=colb)
c                   = rowb then it is possible to call cmatml with c 
c                   equals to a, and the result is overwritten in a.
c
c*********************************************************************
c
c         this version dated june 7, 1987
c         authors: jim demmel and bo kagstrom
c
c*****    internal variables
c
      integer       i,j
c
c*****    cmatml uses the following functions and subroutines
c
c         blas      zcopy, zdotc, zdotu
c
c*****    determine what is to be computed via nested if-then -else's
c
      do 20 j = 1, cola
          do 10 i = 1, rowb
            if     (job .eq. 1) then
               work(i) = zdotu(rowa,b(i,1),ldb,a(1,j),1)
            elseif (job .eq. 2) then
               c(i,j) = zdotu(rowa,b(i,1),ldb,a(1,j),1)
            elseif (job .eq. 3) then
               work(i) = zdotc(rowa,b(1,i),1,a(1,j),1)
            else
c                  (job .eq. 4)
               c(i,j) = zdotc(rowa,b(1,i),1,a(1,j),1)
            endif
   10     continue
          if (job .eq. 1 .or. job .eq. 3) then
             call zcopy(rowa,work,1,a(1,j),1)
          endif
   20 continue
      return
      end


      subroutine cmatmr(a,lda,rowa,cola,b,ldb,colb,c,ldc,work,job)
c
c     implicit none
      integer       lda,rowa,cola,ldb,colb,ldc,job
      complex*16    a(lda,lda),b(ldb,ldb),c(ldc,ldc),work(*)
      complex*16    zdotu,zdotc
c
c***********************************************************************
c
c     cmatmr performs  a complex (right) matrix multiplication  a * b,
c     or a * b' ,(' = transpose ,conjugate), where a is rowa * cola,
c     b is cola * colb. the result is stored in c or overwritten in a.
c     note the extra restrictions in dimension of b when job = 3 or 4.
c
c     on entry
c         
c         a         complex(lda,cola), where lda>=rowa.
c
c         lda       integer
c                   lda is the leading dimension of the array a.
c                    
c         rowa      integer
c                   rowa is the number of rows of a, which is also
c                   the number of rows in the resulting matrix.
c         cola      integer
c                   cola is the number of columns of a, which is also
c                   the number of rows of b.
c
c         b         complex(ldb,colb), ldb>=cola.
c                   
c         ldb       integer
c                   ldb is the leading dimension of the array b.
c                                           
c         colb      integer
c                   colb is the number of columns  of b, which is
c                   also the number of columns of the resulting matrix
c
c         ldc       integer
c                   ldc is the leading dimension of the array c
c
c         work      complex(cola)
c                   work is a scratch array.
c
c         job       integer
c                   job controls the matrix multiplication, and has
c                   the following meaning
c                   job=1       a = a * b
c                   job=2       c = a * b 
c                   job=3       a = a * b'
c                   job=4       c = a * b'
c
c    on return
c
c         c         complex(ldc,colb), where ldc>=rowa.
c                   c is the matrix product of a and b. if cola(=rowb)
c                   = colb then it is possible to call cmatmr with c 
c                   equals to a, and the result is overwritten in a.
c
c*********************************************************************
c
c         this version dated june 7, 1987
c         authors: jim demmel and bo kagstrom 
c
c*****    internal variables
c
      integer       i,j
c
c*****    cmatmr uses the following functions and subroutines
c
c         blas      zcopy, zdotc, zdotu
c
c*****    determine what is to be computed via nested if-then -else's
c
      do 20 i = 1, rowa
          do 10 j = 1, colb
            if     (job .eq. 1) then
               work(j) = zdotu(cola,a(i,1),lda,b(1,j),1)
            else if (job .eq. 2) then
               c(i,j) = zdotu(cola,a(i,1),lda,b(1,j),1)
            else if (job .eq. 3) then
               work(j) = zdotc(cola,b(j,1),ldb,a(i,1),lda)
            else
c                  (job .eq. 4)
               c(i,j) = zdotc(cola,b(j,1),ldb,a(i,1),lda)
            end if
   10     continue
          if (job .eq. 1 .or. job .eq. 3) then
             call zcopy(cola,work,1,a(i,1),lda)
          end if
   20 continue
      return
      end

      subroutine zqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      integer ldx,n,p,job
      integer jpvt(1)
      complex*16 x(ldx,1),qraux(1),work(1)
c
c     zqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       complex*16(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    complex*16(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the unitary part of the decomposition

c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   complex*16(p).
c                qraux contains further information required to recover
c                the unitary part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     zqrdc uses the following functions and subprograms.
c
c     blas zaxpy,zdotc,zscal,zswap,dznrm2
c     fortran dabs,dmax1,cdabs,dcmplx,cdsqrt,min0
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dznrm2,tt
      complex*16 zdotc,nrmxl,t
      logical negj,swapj
c
      complex*16 csign,zdum,zdum1,zdum2
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2))
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call zswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call zswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dcmplx(dznrm2(n,x(1,j),1),0.0d0)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (dreal(qraux(j)) .le. maxnrm) go to 90
                  maxnrm = dreal(qraux(j))
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call zswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = (0.0d0,0.0d0)
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dcmplx(dznrm2(n-l+1,x(l,l),1),0.0d0)
            if (cabs1(nrmxl) .eq. 0.0d0) go to 180
               if (cabs1(x(l,l)) .ne. 0.0d0)
     *            nrmxl = csign(nrmxl,x(l,l))
               call zscal(n-l+1,(1.0d0,0.0d0)/nrmxl,x(l,l),1)
               x(l,l) = (1.0d0,0.0d0) + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -zdotc(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call zaxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (cabs1(qraux(j)) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (cdabs(x(l,j))/dreal(qraux(j)))**2
                     tt = dmax1(tt,0.0d0)
                     t = dcmplx(tt,0.0d0)
                     tt = 1.0d0
     *                    + 0.05d0*tt
     *                      *(dreal(qraux(j))/dreal(work(j)))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*cdsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dcmplx(dznrm2(n-l,x(l+1,j),1),0.0d0)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end


      subroutine zqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      complex*16 x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),xb(1)
c
c     zqrsl applies the output of zqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to zqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  zqrdc produces a factored unitary matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      complex*16(ldx,p).
c               x contains the output of zqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in zqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to zqrdc.
c
c        qraux  complex*16(p).
c               qraux contains the auxiliary output from zqrdc.
c
c        y      complex*16(n)
c               y contains an n-vector that is to be manipulated
c               by zqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     complex*16(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    complex*16(n).
c               qty contains ctrans(q)*y, if its computation has
c               been requested.  here ctrans(q) is the conjugate
c               transpose of the matrix q.
c
c        b      complex*16(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in zqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into zqrdc.)
c
c        rsd    complex*16(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     complex*16(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call zqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)

c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     zqrsl uses the following functions and subprograms.
c
c     blas zaxpy,zcopy,zdotc
c     fortran dabs,min0,mod
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      complex*16 zdotc,t,temp
      logical cb,cqy,cqty,cr,cxb
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (cabs1(x(1,1)) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = (0.0d0,0.0d0)
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call zcopy(n,y,1,qy,1)
         if (cqty) call zcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (cabs1(qraux(j)) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -zdotc(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call zaxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute ctrans(q)*y.
c
            do 90 j = 1, ju
               if (cabs1(qraux(j)) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -zdotc(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call zaxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call zcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call zcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call zcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = (0.0d0,0.0d0)
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = (0.0d0,0.0d0)
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (cabs1(x(j,j)) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call zaxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (cabs1(qraux(j)) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -zdotc(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call zaxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -zdotc(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call zaxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end


c    On this file june 13, 1987:
c    listr, ppcj
      subroutine listr (opt, a, b, ldab, m, n, rowb, rowe,
     *                  colb, cole, first, zero, epsua, epsub, gap,
     *                  pp, ldpp, qq, ldqq, kstr, kfirst, step,
     *                  adlsvd, bdlsvd,
     *                  work, x, sx, ex, q, arow, brow, w, qraux, y,
     *                  qty, info)
c
c     implicit none
c**** debug space
c     the common-block declarations assume that the dimension of the
c     input matrix pencil a - lambda b is not larger than abdim.
c     the debug space is used for producing debug outputs (optional,
c     see below)
c
      integer abdim
      parameter (abdim = 30)
      common /debug1/ acopy(abdim,abdim),bcopy(abdim,abdim),
     *                atest(abdim,abdim),btest(abdim,abdim),swap
      common /debug2/ idbg(20), outunit
      complex*16 acopy, bcopy, atest, btest
      logical swap 
      integer idbg, outunit
c
c**** formal parameter declarations
      character*(*) opt
      integer ldab, m, n, rowb, rowe, colb, cole, ldpp, ldqq,
     *        kstr(4,*), step, kfirst, info
      logical first, zero
      real*8 adlsvd, bdlsvd, epsua, epsub, gap
      complex*16 a(ldab,*), b(ldab,*),pp(ldpp,*), qq(ldqq,*),
     *           work(*)
c
c**** work space
c
        complex*16   x(m,n), sx(*), ex(*), q(m,m), 
     *               arow(*), brow(*), w(m,m), qraux(*), y(*),
     *               qty(*)
c
c*******************************************************************
c
c     listr computes the kronecker left (row) structure and
c     the jordan structure of the infinite-eigenvalue of a singular
c     pencil a-lambda*b. for details concerning the listr-kernel see 
c     the following papers:
c     
c        b.kagstrom, rgsvd - an algorithm for computing the kronecker
c             structure and reducing subspaces of singular a - lambda b
c             pencils, siam j.sci.stat.comput., vol. 7, 1986, pp 185-211
c
c        j.demmel and b.kagstrom, stably computing the kronecker 
c             structure and reducing subspaces of singular pencils
c             a - lambda b for uncertain data, in large scale eigenvalue
c             problems (cullum, willoughby eds), north holland, 1986,
c             pp 283-323.
c
c
c     formal parameters
c
c     on entry
c
c        opt*(*) character, if opt = 'cind' listr computes indices
c                           if opt = 'rind' already computed indices
c                           are reused in the reduction
c
c        a(ldab,*) complex*16, input matrix a of order m by n
c
c        b(ldab,*) complex*16, input matrix b of order m by n
c
c        ldab      integer, leading dimension of a and b
c
c        m         integer, current row dimension of a and b
c
c        n         integer, current column dimension of and b
c
c        rowb      integer, first row of the subpencil
c
c        rowe      integer, last row of the subpencil
c
c        colb      integer, first column of the subpencil
c
c        cole      integer, last column of the subpencil
c
c        first     logical, first should be 'true' if first call to 
c                  listr, else 'false'
c
c        zero      logical, if 'true', zero out small singular values
c                  so returned pencil really has structure described
c                  in kstr (see below), else returned pencil is a
c                  true equivalence transformation of input pencil
c                  (no singular values are deleted)
c
c        epsua     real*8, threshold for deleting singular values of a
c                  (used when compressing rows of a)
c
c        epsub     real*8, threshold for deleting singular values of b
c                  (used when compressing rows of b)
c
c        gap       real*8, should be at least 1 and nominally 1000.
c                  used by subroutine rcsvdc to make rank decisions
c                  by searching for adjacent singular values whose
c                  ratio exceeds gap.
c
c        ldpp      integer, leading dimension of pp
c
c        ldqq      integer, leading dimension of qq
c
c        kfirst    integer, index to the first location in kstr
c                  where structure-index information is stored
c                  from this reduction (see below)
c
c     on exit
c
c        pp(ldpp,*)complex*16, left unitary transformation matrix 
c                  pp of order m by m
c
c        qq(ldqq,*)complex*16, right unitary transformation matrix
c                  qq of order n by n
c
c        a(ldab,*) transformed matrix a (pp**h * a * qq)
c
c        b(ldab,*) transformed matrix b (pp**h * b * qq)
c
c        kstr(4,*) integer, stores information concerning left 
c                  kronecker indices and the jordan structure of
c                  the infinite eigenvalue.
c                  kstr(1,kfirst-1+j) - kstr(2,kfirst-1+j) =
c                  number of l(j-1)**t blocks (left indices of
c                  degree j-1).
c                  kstr(2,kfirst-1+j) - kstr(1,kfirst+j) = 
c                  number of jordan blocks of the infinite
c                  eigenvalue of dimension j.
c                  index j goes from 1 to step (see below)
c                  note: rows 3 and 4 of kstr are not used inside 
c                  listr.
c
c        step      integer,  the number of deflation-steps in this
c                  reduction
c
c        adlsvd    real*8, root sum of squares of deleted singular
c                  values of a (independent of the input zero)
c
c        bdlsvd    real*8, root sum of squares of deleted singular
c                  values of b (independent of the input zero)
c
c        info      integer, zero if normal return,
c                           1 if svd does not converge
c
c        on exit form listr a and b will be in block upper triangular form:
c
c
c                 a = ( a11   *  )        b = ( b11    *  )
c                     (  0   ali )            (  0    bli )
c
c       the block structure of the pencil ali - lambda*bli describes 
c       the kronecker row (left) structure and the jordan structure 
c       of the infinite eigenvalue. if ni and ri denote the dimension of
c       the diagonal blocks in ali and bli (see example below),
c       then they have the following interpretation:
c
c         ni - ri = the number of l(i-1)**t -blocks of order i by i-1
c         ri - ni+1 = the number of j(inf)-blocks of order i by i
c
c       note that if a - lambda*b is a regular pencil then ni=ri.
c       the listr reduction stops when an ni.eq.0 or ni.ne.0 but ri.eq.0. 
c       then b11 will have full row rank. a11 - lambda*b11 might
c       still be a singular pencil (can have right (column) indices). 
c       an example illustrates the two cases (see papers for details):
c       case 1 - n4.eq.0:
c
c               ( a11 a12 a13 ) n3           (  0  b12 b13 ) n3
c         ali = (  0  a22 a23 ) n2     bli = (  0   0  b23 ) n2
c               (  0   0  a33 ) n1           (  0   0   0  ) n1
c                  r3  r2  r1                   r3  r2  r1
c
c        case 2 - n4.ne.0 and r4.eq.0:
c
c                ( a11 a12 a13 ) n4          ( b11  b12 b13 ) n4
c          ali = ( a21 a22 a23 ) n3    bli = (  0   b22 b23 ) n3
c                (  0  a32 a33 ) n2          (  0   0   b33 ) n2
c                (  0   0  a43 ) n1          (  0   0    0  ) n1
c                   r3  r2  r1                   r3  r2  r1
c
c        the ni by ri subdiagonal blocks ai+1i of ali are in the form
c        (rii)
c        ( 0 ), where rii is ri by ri, nonsingular and upper triangular.
c
c        if kfirst = 1 on input then case 2 above cause the following
c        output for step and kstr:
c          step = 4
c          kstr(1,1) = n1    kstr(2,1) = r1
c          kstr(1,2) = n2    kstr(2,2) = r2
c          kstr(1,3) = n3    kstr(2,3) = r3
c          kstr(1,4) = n4    kstr(2,4) = 0
c
c        note that on output (ali,bli) or (a11,b11) can be nonexistent
c        in the block upper triangular form (a,b). (ali,bli) does not 
c        exist if n1=r1=0. (a11,b11) does not exist if the input pencil
c        a -lambda*b has no right (column) singular structure and no
c        finite eigenvalues.
c
c
c***     work space including size (all variables complex*16)
c        work(*)           max(m,n) (if linpack svd is used)
c                          or 2*min(m,n)+6*max(m,n) (lapack svd, 920708)
c        x(m,*)            m by n
c        sx(*)             min(m,n) + 1
c        ex(*)             n
c        q(m,*)            m by m
c        arow(*)           max(m,n)
c        brow(*)           max(m,n)
c        w(m,*)            m by m
c        qraux(*)          max(m,n)
c        y(*)              max(m,n)
c        qty(*)            max(m,n)
c
c*****************************************************************
c
c****    this version dated june 16, 1987
c        authors: jim demmel, bo kagstrom
c
c****    listr uses the following functions and subroutines
c        kcfpack - cmatml, cmatmr, cmatpr, cmcopy, ppcj,
c                  rcsvdc, upddel
c        linpack - zqrdc, zqrsl
c
c****    internal variables
c
        logical ldebug
        integer mrow, ncol, i, j, sn1, sr1, rep, rowsn1, colsr1, xrow
     *          , xcol, job, ldx, ldq, n1, rnull, ldw, cnull, r1
     *          , rowbm1, colbm1, idummy, ikstr, mxrc, k, iii, jjj
c
        real*8 del, difa, difb
c
        complex*16 dummy
c
c**** set leading dimensions of x, q, w
        ldx = m
        ldq = m
        ldw = m
c       set debug switch
        ldebug = idbg(5) .ne. 0
c****   compute the order of the pencil in action (mrow * ncol)
        mrow = rowe - rowb + 1
        ncol = cole - colb + 1
c
c*+*+*+ accumulate deleted singular values in adlsvd and bdlsvd
        adlsvd = 0.0
        bdlsvd = 0.0
c
      if (ldebug) then
         write (outunit,1000) 'epsua=', epsua
         write (outunit,1000) 'epsub=', epsub
      endif
c
c
c**** set rep depending on what option
c
      if ( opt .eq. 'cind' ) then
c         perhaps not enough !!
          rep = rowe * cole
      else
c         the number of deflation steps
c          rep = step - kfirst + 1
c*****      Changes made 1986-06-17
         rep = step
      endif
c***  6/18/87
      if (ldebug) write(outunit,2000) 'kfirst=',kfirst,
     +            'step=',step,'rep=',rep
      sn1 = 0
      sr1 = 0
      step = 0
c**** while rep > 0 do
   30 continue
      if (ldebug) write(outunit,2000) 'rep at top of loop=',rep
      if (rep .eq. 0) go to 500
c     jump when while - loop satisfied
c
c     while - clause
        step = step + 1
        if (ldebug) then
           write( outunit, 2000 ) 'Results from step = ', step
 2000      format( t5, a, i3/)
           write(outunit,2005) opt
 2005      format(t5,a)
        endif
c
c**** set n1 and r1 if we are reusing kstr
c
      if ( opt .eq. 'rind' ) then
         ikstr = kfirst + step - 1
         n1 = kstr(1, ikstr)
         r1 = kstr(2, ikstr)
         rnull = n1 -r1
      endif
c
c**** step 1 - compress rows of b (gives n1 = dimension of the
c              row nullspace)
c* 1.1
c      rows, rowb:rowe-sn1
c      cols, colb:cole-sr1
c------------------------------
        rowsn1 = rowe - sn1
        colsr1 = cole - sr1
c
        rowbm1 = rowb - 1
        colbm1 = colb - 1
        xrow = mrow - sn1
        xcol = ncol - sr1
c****   6/18/87 fix 
        if (opt .eq. 'rind') cnull = n1 - (xrow - xcol)
c****
        do 40 i = 1, xrow
           do 35 j = 1, xcol
              x(i, j) = b(rowbm1 + i, colbm1 + j)
   35      continue
   40   continue
        job = 1000
      if (ldebug) then
        write(outunit,5000) 'rowsn1=',rowsn1,'colsr1=',colsr1,
     *                      'xrow=',xrow
        write(outunit,5000) 'xcol=',xcol,'rowb=',rowb,'rowe=',rowe
        write(outunit,5000) 'colb=',colb,'cole=',cole,'sr1=',sr1,
     *                      'sn1=',sn1
        if ( opt .eq. 'rind') then
           write(outunit,5000) 'cnull=',cnull,'rnull=',rnull
           write(outunit,5000) 'n1=',n1,'r1=',r1
        endif
      endif
c       put m*n in info before calling (why ? 870608)
        info = m*n
        call rcsvdc (x, ldx, xrow, xcol, sx, ex, q, ldq, dummy, 1, opt,
     *               epsub, gap, cnull, n1, del, work, job, info )
c
c
        call upddel(bdlsvd, del)
c
      if (ldebug) then
        write(outunit,1000) 'bdlsvd=', bdlsvd, 'del=', del
        mxrc = min0( xrow, xcol)
        call cmatpr( sx, 1, 1, mxrc,
     *               'singular values - row compress b')
        call cmatpr( ex, 1, 1, mxrc,'sub diagonal - should be zero')
        write (outunit, 1005) 'info=', info, '(rownullity) n1=', n1
 1005   format(t5, a, i3/ )
      endif
c
        if (info .ne. 0 ) then
c***    6/18/87
          if (ldebug) write(outunit,2007) info
 2007     format('listr - after first call to rcsvdc, info =',i4)
          info = 1
          return
        endif
c     
c       if n1=0, we are done
        if (n1 .eq. 0) then
           r1=0
           goto 450
        end if
c
c* 1.2 - apply left transformation q to a and b (the full matrices)
c        cols in a and b: colb:n
c        rows in a:  rowb:rowe-sn1 ( xrow row's)
c        rows in b:  rowb:rowe-sn1
c-----------------------------------
        do 70 i = colb, n
           do 50 j = 1, xrow
              arow(j) = 0.d0
              brow(j) = 0.d0
              do 45 k = 1, xrow
                 arow(j) = arow(j) + a(rowbm1 + k, i) * conjg(q(k,j))
                 brow(j) = brow(j) + b(rowbm1 + k, i) * conjg(q(k,j))
   45         continue
   50      continue
           do 60 j = 1, xrow
                 a(rowbm1 + j, i) = arow(j)
                 b(rowbm1 + j, i) = brow(j)
   60      continue
   70   continue
c
c*         zero part of b
c          rows, rowe-sn1-n1+1:rowe-sn1
c          cols, colb:cole-sr1
c----------------------------------------
        if (zero) then
          do 80 i = rowe - sn1 - n1 + 1, rowe - sn1
             do 75 j = colb, cole - sr1
                b(i, j) = 0.d0
   75        continue
   80     continue
        endif
c
c**** Step 2 - row compress part of A ( gives n1 - r1 =
c              dimension of the common nullspace)
c
c* 2.1
c       rows, rowe-sn1-n1+1:rowe-sn1
c       cols, colb:cole-sr1
c-----------------------------------
        xrow = n1
        xcol = ncol - sr1
        do 90 i = 1, xrow
           do 85 j = 1, xcol
              x(i, j) = a( rowsn1 - n1 + i, colbm1 + j)
   85      continue
   90   continue
c
        job = 1000
        info = m*n
c****   6/18/87 fix
        if (opt .eq. 'rind') cnull = xcol - r1
c
        call rcsvdc ( x, ldx, xrow, xcol, sx, ex, w, ldw, dummy, 1, opt,
     *                epsua, gap, cnull, rnull, del, work, job, info )
c
c
        if ( opt .eq. 'cind' ) r1 = n1 - rnull
c
c       if r1 = 0 then we are done ! zero part in a and update qq
c
      if (ldebug) then
        write (outunit, 1005) 'info=', info, 'rnull=', rnull,
     *                  'n1=', n1,'r1=', r1
c
        mxrc = min0( xrow, xcol)
        call cmatpr( sx, 1, 1, mxrc,
     *               'singular values - row compress a')
        call cmatpr( ex, 1, 1, mxrc,'sub diagonal - should be zero')
      endif
c
        if (info .ne. 0) then
c****   6/18/87
          if (ldebug) write(outunit,2008) info
 2008     format('listr - after second call to rcsvdc, info= ',i4)
          info = 1
          return
        endif
c
       call upddel(adlsvd, del)
c
c
       if (r1 .eq. 0) goto 3500
c
c* 2.2
c      update left transformation q
c                rows, 1:mrow-sn1 (xrow)
c                cols, xrow-n1+1:xrow
c___________________________________________________________________
c
       xrow = mrow - sn1
       do 110 i = 1, xrow
          do 100 j = 1, n1
             arow(j) = 0.d0
             do 95 k = 1, n1
                arow(j) = arow(j) + q(i, xrow - n1 + k) * w(k, j)
   95        continue
  100     continue
c
          do 105 j = 1, n1
             q(i, xrow - n1 + j) = arow(j)
  105     continue
  110   continue
c
 1000      format(t5, a, d13.5)
        if (idbg(5) .gt. 1) then
          call cmatpr(q,ldq,mrow-sn1,mrow-sn1,'q after step 2.2')
          call cmatpr(a,ldab,m,n,'a before step 2.2')
          call cmatpr(b,ldab,m,n,'b before step 2.2')
        endif
c       
c****   now a and b ....with w too
c           rows, rowe-sn1-n1+1:rowe-sn1
c           cols, colb:n
c
c       note that we do not make use of that some of the elements
c       in b are zero
c
        do 120 i = colb,n
           do 114 j = 1, n1
              arow(j) = 0.d0
              brow(j) = 0.d0
              do 112 k = 1, n1
                 arow(j) = arow(j) + a(rowsn1-n1+k,i) * conjg(w(k,j))
                 brow(j) = brow(j) + b(rowsn1-n1+k,i) * conjg(w(k,j))
  112         continue
  114      continue
           do 116 j = 1, n1
              a(rowsn1-n1+j,i) = arow(j)
              b(rowsn1-n1+j,i) = brow(j)
  116      continue
  120   continue
        if (idbg(5) .gt. 1) then
          call cmatpr(a,ldab,m,n,'a after step 2.2')
          call cmatpr(b,ldab,m,n,'b after step 2.2')
        endif
c      
c*        zero part of a
c         rows, rowe-sn1-(n1-r1)+1:rowe-sn1
c         cols, colb:cole-sr1
c--------------------------------------------
c
 3500  continue
       if (zero) then
         if (ldebug) then
             write(outunit, 4005) 'loop indices in 130',
     *                       rowsn1 - (n1 - r1) + 1,rowsn1
             write(outunit, 4005) 'loop indices in 125', colb, colsr1
 4005       format(t5, a, 2i5)
         endif 
         do 130 i = rowsn1 - (n1 - r1) + 1,rowsn1
            do 125 j = colb, colsr1
               a(i,j) = 0.d0
  125       continue
  130    continue
       endif
       if (r1 .eq. 0) go to 350
c
c**** Step 3 - Triangularize A by a rq-decomposition ( using qr)
c
c* 3.1
c         rows, rowb:rowe-sn1-(n1-r1)
c         cols, colb:cole-sr1
c---------------------------------------------
c
          xrow = mrow - sn1 - (n1 - r1)
          xcol = ncol - sr1
c         move a(trans,conjg) with permuted columns (n,n-1,...1)
c
          do 140 i = 1, xcol
             do 135 j = 1, xrow
                 x(i, j) = conjg( a(rowsn1 - (n1-r1)+1-j, colbm1 + i))
  135        continue
  140     continue
          job = 0
         if (idbg(5) .gt. 1) then
          call cmatpr(x,ldx,xcol,xrow,'part of a before qr-decomp')
         endif
          call zqrdc( x, ldx, xcol, xrow, qraux, idummy, dummy, job)
c
c****     move the upper triangular part to a
c
          if (ldebug) then
             write(outunit, 5000) 'xrow=', xrow, 'xcol=', xcol
             write(outunit, 5000) 'rowb=',rowb, 'colb=', colb
c             call cmatpr(x,ldx,xcol,xrow,'x after call to zqrdc')
             write(outunit, 1010) 'a(rowb,colb)', a(rowb,colb)
1010         format(t5, a, 2d15.5)
          endif
          call ppcj( x, ldx, 1, xcol, 1, xrow, a(rowb, colb), ldab)
c
c         zero elements in a to make it upper triangular!
c             rows, rowe-sn1-(n1-r1)-(xcol-2):rowe-sn1-(n1-r1)
c             cols, colb:cole-sr1-1

          do 150 i = colb, cole - sr1 -1
             do 148 j = i-colb+( rowsn1-(n1-r1)-xcol+2),rowsn1-(n1-r1)
                a(j, i) = 0.d0
  148        continue
  150    continue
         if (idbg(5) .gt. 1) then
           call cmatpr(a,ldab,m,n,'A after triangularization')
         endif
c
c* 3.2
c         apply v (xcol*xcol)to remaining rows of a
c              rows, 1:rowb-1
c              cols, colb:cole-sr1
c
c-----------------------------------------------------
c
         do 170 j = 1, rowb - 1

            do 160 i = 1, xcol
               y(i) = conjg(a( j, colbm1 + i))
  160       continue
            job = 01000
            call zqrsl(x, ldx, xcol, xrow, qraux, y, dummy, qty,
     *                 dummy, dummy, dummy, job, info)
            do 165 i = 1, xcol
               a(j, colsr1-i+1) = conjg( qty(i))
  165       continue
  170    continue
         if (idbg(5) .gt. 1) then
            call cmatpr(a, ldab, m, n,
     *              ' A after triangularization - step 3.1')
         endif
c
c        apply v to b from right (xcol*xcol)
c            rows, 1:rowe-sn1-n1
c            cols, colb:cole-sr1
c----------------------------------------------------------
c
         do 185 j = 1, rowe - sn1 - n1
            do 180 i = 1, xcol
               y(i) = conjg(b(j, colbm1 + i))
  180       continue
            job = 01000
            call zqrsl(x, ldx, xcol, xrow, qraux, y, dummy, qty,
     *                 dummy, dummy, dummy, job, info)
            do 175 i = 1, xcol
               b(j, colsr1 - i + 1) = conjg(qty(i))
  175       continue
  185    continue
         if (idbg(5) .gt. 1) then
            call cmatpr(b, ldab, m, n,
     *              ' B after triangularization - step 3.1')
         endif
c
c****    update right transformation matrix qq ( n*n )
c        rows, 1:n
c        cols, colb:cole-sr1
c-----------------------------------------------------
c
         do 200 j = 1, n
            do 195 i = 1, xcol
               y(i) = conjg(qq(j ,colbm1 + i))
  195       continue
            job = 01000
            call zqrsl(x, ldx, xcol, xrow, qraux, y, dummy, qty,
     *                 dummy, dummy, dummy, job, info)
            do 198 i = 1, xcol
               qq(j, colsr1 - i + 1) = conjg(qty(i))
  198       continue
  200    continue
         if ( idbg(5) .gt. 1) then
           call cmatpr(qq,ldqq,n,n,'qq after updating')
         endif
  350    continue
c
c****    update left transformation matrix pp ( m*m )
c        rows, 1:m
c        cols, rowb:rowe-sn1
c-----------------------------------------------------
c
         xrow = mrow - sn1
         if (first) then
            do 210 i = 1, m
               do 205 j = 1, m
                  pp(i, j) = q(i, j)
  205          continue
  210       continue
        else
            do 240 i = 1, m
               do 230 j = 1, xrow
                  arow(j) = 0.d0
                  do 220 k = 1, xrow
                     arow(j) = arow(j) + pp(i, rowbm1 + k) * q(k, j)
  220             continue
  230          continue
               do 235 j = 1, xrow
                  pp(i, rowbm1 + j) = arow(j)
  235          continue
  240       continue
         endif
c       
c****    update indices
c
         sn1 = sn1 + n1
         sr1 = sr1 + r1
       if (ldebug) then
        write(outunit,5000) 'rowsn1=',rowsn1,'colsr1=',colsr1,
     *                      'xrow=',xrow
        write(outunit,5000) 'xrow=',xrow,'rowb=',rowb,'rowe=',rowe
        write(outunit,5000) 'colb=',colb,'cole=',cole,'sr1=',sr1,
     *                      'sn1=',sn1
       endif
c*       monitoring of the r1 and n1 in kstr
c
  450    continue
       if (ldebug) then
         if (swap) then
           call cmcopy(bcopy,20,m,n,atest)
           call cmcopy(acopy,20,m,n,btest)
         else
           call cmcopy(acopy,20,m,n,atest)
           call cmcopy(bcopy,20,m,n,btest)
         end if
         call cmatml(atest,20,m,n,pp,ldpp,m,atest,20,work,3)
         call cmatmr(atest,20,m,n,qq,ldqq,n,atest,20,work,1)
         call cmatml(btest,20,m,n,pp,ldpp,m,btest,20,work,3)
         call cmatmr(btest,20,m,n,qq,ldqq,n,btest,20,work,1)
         difa=0
         difb=0
         do 1234 iii=1,m
           do 5678 jjj=1,n
             difa=difa+abs(atest(iii,jjj)-a(iii,jjj))
             difb=difb+abs(btest(iii,jjj)-b(iii,jjj))
 5678      continue
 1234    continue
         write(outunit,201) 'difa=',difa
 201     format(t5,a,d13.6/)
         call cmatpr(atest,20,m,n,'atest')
         write(outunit,201) 'difb=',difb
         call cmatpr(btest,20,m,n,'btest')
        endif
c
c****    compute rep depending on what option is used
c
         if ( opt .eq. 'cind') then
            kstr(1, kfirst - 1 + step) = n1
            kstr(2, kfirst - 1 + step) = r1
c*****      changed 1986-06-17
c            rep = n1 * r1 * (mrow - sr1) * (ncol - sn1)
            rep = n1 * r1 * (ncol - sr1) * (mrow - sn1)
         else
            rep =rep - 1
         endif
         if (ldebug) then
           write(outunit,5000) 'sn1=',sn1,'sr1=',sr1,'rep=',rep
 5000      format(t5,a,i4/)
         endif
         first = .false.
         go to 30
c
c**** end of while clause
  500 continue
c
      return
      end
c
        subroutine ppcj(from,ldfrom,rowb,rowe,colb,cole,to,ldto)
c
c       take from(rowb:rowe, colb:cole), reverse the columns, reverse
c       the rows, take its conjugate transpose, and store in
c       to(1:cole-colb+1, 1:rowe-rowb+1)
        complex*16 from(ldfrom,*), to(ldto,*)
        integer rowb,rowe,colb,cole,rsum,csum
        rsum=rowe+1
        csum=cole+1
        nrow=rowe-rowb+1
        ncol=cole-colb+1
        do 1 i=1,ncol
          do 2 j=1,nrow
            to(i,j)=conjg(from(rsum-j,csum-i))
2         continue
1       continue
        return
        end

c   in this file june 12 1987:
c   cmatg1, pertb1, cmcopy, cdife, cnorm, cond, cmatpr, matblk1
c
c   routines here are only for debug and test, but called by guptri
c
      subroutine    cmatg1(a,lda,b,ldb,m,n,acopy,bcopy,work,job,
     *                    epsper, trpose, binfile)
c     implicit none
c**
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
      integer       lda,ldb,m,n,job, mtemp
      real*8        epsper 
      complex*16    a(lda,lda),b(ldb,ldb),acopy(lda,lda),
     *              work(*),bcopy(ldb,ldb)
      character*80 binfile
c
c******************************************************************
c
c     this routine reads input pencil from a binary file
c     created by bpenc1.for or bpenc2.for
c     revised:    870612, 920702
c
c******************************************************************
      integer i, j
      logical wanta, wantb, pertur, prints, trpose
      real*8  nea, neb, cnorm, drandu
c
c note: if lda or ldb>= 20, then the dimensions of p and qinv have
c       to be changed
c
c
c*****    determine what is to be computed
c         only prints, pertur and trpose are used  ****
      wanta = job / 1000 .ne. 0
      wantb = mod(job,1000) / 100 .ne. 0
      pertur = mod(job,100) / 10 .ne. 0
      prints = mod(job,10 ) .ne. 0
      write( outunit, 400) 'trpose=', trpose
      write (outunit, 400) 'pertur=', pertur
      write (outunit, 400) 'prints=', prints
  400 format(t5, a, l1)
c
c     read acopy and bcopy from binary file
c
      open(15, file = binfile, form='unformatted', status='old')
      read(15) m, n
      read(15) ((acopy(i,j), j = 1, n), i = 1, m)
      read(15) ((bcopy(i,j), j = 1, n), i = 1, m)
      close(15, status ='keep')
c
c
      if (trpose) then
        do 750 i=1,m
          do 751 j=1,n
            a(j,i)=acopy(i,j)
            b(j,i)=bcopy(i,j)
751       continue
750     continue
        mtemp = m
        m = n
        n = mtemp
      else
        do 752 i=1,m
          do 753 j=1,n
            a(i,j)=acopy(i,j)
            b(i,j)=bcopy(i,j)
753       continue
752     continue
      endif
c
      if (pertur) then
c
      nea = cnorm(a,lda,m,n,0,work) * epsper
      neb = cnorm(b,ldb,m,n,0,work) * epsper
c     add perturbations to a and b
          do 50 i = 1, m
              do 50 j = 1, n
                  a(i,j) = a(i,j) + drandu(1.d0) * nea
                  b(i,j) = b(i,j) + drandu(1.d0) * neb
   50     continue
	  
	  endif
c
c     compute norm(a,e) and norm(b,e)
c
      nea = cnorm(a,lda,m,n,0,work)
      neb = cnorm(b,ldb,m,n,0,work)
      write(outunit,350) 'epsper=', epsper
      write(outunit,350) 'norm(a,e)=', nea, 'norm(b,e)=', neb
  350 format(t5,a,d12.5,tr5,a,d12.5,tr5,a,d12.5)
c
c     copy a and b to acopy and bcopy, respectively
c
      call cmcopy(a,lda,m,n,acopy)
      call cmcopy(b,ldb,m,n,bcopy)
c
      if (prints) then
        call cmatpr(a,lda,m,n,'final version of a input')
        call cmatpr(b,lda,m,n,'final version of b input')
      endif
      return
c
      end
c
      subroutine pertb1(aorig, borig, a, b, ldab,m ,n , epsbnd,
     +                  work, job, nostat)
c     implicit none
c***  debug space
      integer abdim
      parameter (abdim = 30)
      common /debug1/ acopy(abdim,abdim),bcopy(abdim,abdim), 
     +                 atest(abdim, abdim), btest(abdim,abdim), swap
      common /debug2/ idbg(20), outunit
      complex*16 acopy, bcopy, atest, btest
      integer idbg, outunit
      logical swap
c**** formal parameter declarations
      integer ldab, m, n, job
      complex*16   aorig(ldab,*), borig(ldab,*), a(ldab,*), 
     *             b(ldab,*), work(*)
      real*8 epsbnd
      logical nostat
c
c**** add random noise of relative size epsbnd to aorig and borig
c     and store in a and b 
c     job cotrols the structure of the perturbations
c     job = 1 add random perturbations to a and b
c     job = 2 add random perturbations to a
c             add random perturbations to the last n-m columns of b
c     job = 3 add general random perturbations to a only
c
c     if (idbg(11) .ne. 0 ) print out perturbed a and b
c
      real*8 nea, neb, cnorm, drandu
      integer i, j, colb
c
      nea = cnorm(aorig, ldab, m, n, 0, work) *epsbnd
      neb = cnorm(borig, ldab, m, n, 0, work) *epsbnd
c     add perturbations to a
          do 50 i = 1, m
              do 50 j = 1, n
                  a(i,j) = aorig(i,j) + drandu(1.d0) * nea
   50     continue
           
c      add perturbations to b in columns colb to n
       if (job .eq. 1) then
          colb = 1
       elseif (job .eq. 2) then
          colb = n - m + 1
       else
c         job .eq. 3
          colb = n + 1
       endif
       if (job .eq. 1 .or. job .eq. 2) then
          if ( colb .ge. 1) then
             do 70 i = 1, m
                do 60 j = colb, n
                   b(i,j) = borig(i,j) + drandu(1.d0) * neb
   60           continue
                if (job .eq. 2) then
                   do 65 j = 1, colb - 1
                      b(i,j) = borig(i,j)
   65             continue
                endif
   70        continue
          else
            write(outunit,300) 'wrong dimensions! m,n=', m,n
          endif
       else
           call cmcopy(borig, ldab, m, n, b) 
       endif
c     compute  and norm(a,e),norm(b,e)
c
      if (nostat) then
        nea = cnorm(a,ldab,m,n,0,work)
        neb = cnorm(b,ldab,m,n,0,work)
        write(outunit,350) 'epsbnd=', epsbnd
        write(outunit,350) 'norm(aper,e)=', nea, 'norm(bper,e)=', neb
      endif
c
c     copy a and b to acopy and bcopy, respectively
c
      call cmcopy(a,ldab,m,n,acopy)
      call cmcopy(b,ldab,m,n,bcopy)
c
      if (idbg(11) .gt. 0) then
           call cmatpr(a,ldab,m,n,' perturbed a for input to guptri')
           call cmatpr(b,ldab,m,n,' perturbed b for input to guptri')
      endif
      return
c
  100 format(2i4)
  300 format(t5,a,2i5)
  350 format(t5,a,d12.5,tr5,a,d12.5,tr5,a,d12.5)
      end
c
      subroutine    cmcopy(a,lda,m,n,acopy)
c     implicit none
      integer       lda,m,n
      complex*16    a(lda,1),acopy(lda,1)
c
c***  the routine cmcopy copies matrix a to acopy
c
      integer       i,j
c
      do 10 i = 1, m
          do 10 j =1, n
              acopy(i,j) = a(i,j)
   10 continue
c
      return
      end
c
      real*8 function cdife(a,b,ldab,m,n)
c     implicit none
      integer  ldab, m, n
      complex*16 a(ldab,*), b(ldab,*), z
c
c**** the routine computes the frobenius norm of the
c     difference between a and b
c
      integer i,j
      real*8 sum
c
      sum=0.0
      do 10 i = 1, m
         do 5 j = 1, n
            z=a(i,j)-b(i,j)
            sum=sum+dreal(z)**2 + dimag(z)**2
    5    continue
   10 continue
      cdife = sqrt(sum)
      return
      end
c
      real*8 function cnorm(a,lda,m,n,job,work)
c
c     implicit none
c
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
      integer       lda,m,n,job
      complex*16    a(lda,1),work(1)
c
      integer       joba, info, ss, se, sx, sw, i, j
      complex*16    u,v
c
      if (job .eq. 2) then
c         compute the 2-norm
c         allocate space for s(min(m+1,n)),e(n) and x(lda,n) 
          ss = 1
          se = ss + min(m+1,n)
          sx = se + n
          sw = sx + lda*n
          call cmcopy(a,lda,m,n,work(sx))
          joba = 00
          call zsvdc(work(sx),lda,m,n,work(ss),work(se),
     *               u,1,v,1,work(sw),joba,info)
          if (info .ne. 0) then
             write(outunit,100) 
     +            'csvdc did not converge, called from cnorm'
  100        format(t5,a/)
             call cmatpr(work(ss),1,1,n,
     *                     'singular values - main diagonal')
             call cmatpr(work(se),1,1,n,
     *                     'sub-diagonal - should be zero')
          else
c                  = s(1)
             cnorm = work(ss)
          endif
c              ( info .eq. 2)
      else
c         (job .eq. 0), compute the frobenius norm
          cnorm = 0.0
          do 20 i = 1,m
             do 20 j = 1,n
                  cnorm = cnorm + conjg(a(i,j)) * a(i,j)
   20     continue
          cnorm = sqrt(cnorm)
      endif
c         ( job .eq. 2)
      return
      end
c

      real*8 function cond(a,lda,m,n,work)
c  
c     implicit none
c
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
      integer       lda,m,n
      complex*16    a(lda,1),work(1)
c
      integer       joba, info, nn, ss, se, sx, sw
      complex*16    u, v
c
c         allocate space for s(min(m+1,n)),e(n) and x(lda,n) 
          ss = 1
          se = ss + min(m+1,n)
          sx = se + n
          sw = sx + lda*n
          call cmcopy(a,lda,m,n,work(sx))
          joba = 00
          call zsvdc(work(sx),lda,m,n,work(ss),work(se),
     *               u,1,v,1,work(sw),joba,info)
c
          if (info .ne. 0) then
             write(outunit,100) 
     +             'csvdc did not converge, called from cond'
  100        format(t5,a/)
             call cmatpr(work(ss),1,1,n,
     *                    'singular values - main diagonal')
             call cmatpr(work(se),1,1,n,
     *                    'sub-diagonal - should be zero')
          else
             nn = min(m,n)
             if (work(nn) .eq. (0.,0.)) then
                   cond = 0.0
                   write(outunit,100) 
     +                   'cond = the matrix is singular'
             else
                   cond = dreal(work(ss)/work(nn))
c                                       s(1)/s(nn))
             endif
          endif
          return
          end


      subroutine cmatpr(a,lda,m,n,text)
c     implicit none
c
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
      integer       lda,m,n, k
      complex*16    a(lda,*)
      character*(*) text
c
      write(outunit, 300) 'lda=',lda, 'm=', m, 'n=', n
  300 format(3(5x,a,i3))
      write(outunit,100)text
  100 format(t5,a)
      write(outunit,200) ('-',k=1,70)
  200 format(t5,70a)
c
      if (lda .eq. 1) then
          call matblk1(a,lda,1,1,1,n)
      else
          call matblk1(a,lda,1,m,1,n)
      endif
      return
      end
c
      subroutine matblk1(a,lda,rf,rs,kf,ks)
c     implicit none
c
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
      integer lda, rf, rs, kf, ks
      real*8 a(2,lda,*)
c
      integer tpr,blk,ifirst,bl,ilast,i,j,k,l
      real*8 aim
c
c     tpr is the number of elements per output-row
c
c     is a real or complex ? yes if aim = 0.d0
      aim = 0.d0
      do 20 i = rf, rs
       do 10 j = kf, ks
          if (a(2,i,j) .ne. 0.0d0) aim = 1.
   10  continue
   20 continue 
      tpr = 3
      blk = (ks - kf) / tpr + 1
      ifirst = kf
      do 40 bl = 1, blk
         if (bl .ne. blk) then
            ilast = ifirst + tpr - 1
         else
            ilast = ks
         endif
         do 30 k = rf, rs
             if ( aim .eq. 0.d0) then
c               a is real
                write(outunit,50) (a(1,k,i), i=ifirst, ilast)
             else
c                a is complex 
               do 25 l = 1, 2
                  write(outunit,50) (a(l,k,i), i=ifirst, ilast)
   25          continue
               write(outunit,60)
             endif
   30    continue
         ifirst = ifirst +tpr
      write(outunit,60)
   40 continue
      return
c
   50 format(t3,3d24.17)
   60 format(/)
      end
c     File zqz.for contains: zqz, rcopy, icopy, ricopy,
c                            zqzhs1, zqzvl1
c     Date: 12 june, 1987
c
      subroutine zqz(a, b, ldab, dimreg, rowb, colb, q, ldq,
     *               ph, ldp, ierr, work)
c
c     implicit none
c
c**** debug space
      common/debug2/ idbg(20), outunit
      integer idbg, outunit
      logical ldebug
c
c**** formal parameter declarations
      integer ldab, dimreg, rowb, colb, ldq, ldp, ierr
      complex*16 a(ldab,*), b(ldab,*), q(ldq,*), ph(ldp,*)
      complex*16 work(*)
c
c*********************************************************************
c
c     this routine reduces the remaining regular part (corresponding to
c     the nonzero and finite eigenvalues) to upper
c     triangular form by using the qz algorithm.
c     this routine is necessary since there is no complex*16
c     version of the qz-routine in eispack or linpack.
c
c     on entry
c
c     a       complex*16(ldab,*), where ldab >= dimreg
c
c     b       complex*16(ldab,*)
c
c     ldab    integer
c             leading dimension of the arrays a and b
c
c     dimreg  integer
c             dimension of the remaining regular part
c
c     rowb    integer
c             first row in a and b of remaining regular part
c
c     colb    integer
c             first column in a and b of remaining regular part
c
c     ldq     integer
c             leading dimension of the array q
c
c     ldp     integer
c             leading dimension of the array p
c
c     work    complex*16(2*dimreg*dimreg+4 + 3*dimreg)
c             scratch array
c
c     idbg(12) integer
c             if nonzero, turn on debug output
c
c     on exit
c
c     a       changed, contains the upper triangular a-part of the
c             qz-decomposition
c
c     b       changed, contains the upper triangular b-part of the
c             qz-decomposition
c
c     q       complex*16(ldq,*), where ldq >= dimreg
c             right transformation matrix
c
c     ph      complex*16(ldp,*), where ldp >= dimreg
c             conjugate transpose of the left transformation matrix
c
c     ierr    integer
c             error messages from qz-algorithm
c             zero for normal return
c             nonzero if eigenvalue has not converged in 50 iterations
c             (for more details see routine zqzvl1)
c
c*********************************************************************
c
c     this version dated june 10, 1987
c     authors: jim demmel and bo kagstrom
c
c**** zqz uses the following functions and subroutines
c
c     cmatpr, icopy, rcopy, ricopy, zqzhs1, zqzvl1
c
c**** internal variables
c
      integer dimsqr, strtar, strtai, strtbr, strtbi
      integer j
      integer alfarb, alfaib, betab
c
      ldebug = idbg(12) .ne. 0
c
c**** reduce the pencil a(rowb:rowb+dimreg-1,colb:colb+dimreg-1)
c       - lambda b(rowb:rowb+dimreg-1,colb:colb+dimreg-1)
c     to upper triangular form with the qz algorithm
c     copy a, b to separate arrays for real and imaginary parts in
c     preparation for using zqzhs1, zqzvl1
c
      dimsqr = dimreg * dimreg / 2 + 1
      strtar = 1
      strtai = strtar + dimsqr
      strtbr = strtai + dimsqr
      strtbi = strtbr + dimsqr
      alfarb = strtbi + dimsqr
      alfaib = alfarb + dimreg
      betab = alfaib + dimreg
      call rcopy(a(rowb,colb), ldab, dimreg, work(strtar))
      call icopy(a(rowb,colb), ldab, dimreg, work(strtai))
      call rcopy(b(rowb,colb), ldab, dimreg, work(strtbr))
      call icopy(b(rowb,colb), ldab, dimreg, work(strtbi))
c
      if (ldebug) then
        write(outunit,100) 'entering qz'
        write(outunit,100) 'ldab=',ldab,'dimreg=',dimreg,
     +                     'rowb=',rowb,'colb=',colb,
     +                     'ldq=',ldq,'ldp=',ldp,
     +                     'dimsrq=',dimsqr,'strtar=',strtar,
     +                     'strtai=',strtai,'strtbr=',strtbr,
     +                     'strtbi=',strtbi,'alfarb=',alfarb,
     +                     'alfaib=',alfaib,'betab=',betab
100     format(3x,a,1x,i4)
        write(outunit,100) 'areal'
        write(outunit,101) (work(strtar-1+j),j=1,dimsqr)
        write(outunit,100) 'aimag'
        write(outunit,101) (work(strtai-1+j),j=1,dimsqr)
        write(outunit,100) 'breal'
        write(outunit,101) (work(strtbr-1+j),j=1,dimsqr)
        write(outunit,100) 'bimag'
        write(outunit,101) (work(strtbi-1+j),j=1,dimsqr)
101     format(3d23.16)
      endif
      call zqzhs1(dimreg, dimreg, work(strtar),work(strtai),
     *            work(strtbr),work(strtbi),
     *              .true., q, ldq, .true., ph, ldp )
c
      if (ldebug) then
        write(outunit,100) 'after zqzhs1'
        write(outunit,100) 'areal'
        write(outunit,101) (work(strtar-1+j),j=1,dimsqr)
        write(outunit,100) 'aimag'
        write(outunit,101) (work(strtai-1+j),j=1,dimsqr)
        write(outunit,100) 'breal'
        write(outunit,101) (work(strtbr-1+j),j=1,dimsqr)
        write(outunit,100) 'bimag'
        write(outunit,101) (work(strtbi-1+j),j=1,dimsqr)
        call cmatpr(q,ldq,dimreg,dimreg,'q after zqzhs1')
        call cmatpr(ph,ldp,dimreg,dimreg,'ph after zqzhs1')
      endif     
c
      call zqzvl1(dimreg, dimreg, work(strtar),work(strtai),
     *            work(strtbr),work(strtbi),
     *              0.0d0, work(alfarb), work(alfaib), work(betab),
     *              .true., q, ldq, .true., ph, ldp, ierr)
c
      if (ldebug) then
        write(outunit,100) 'after zqzvl1, ierr=',ierr
        write(outunit,100) 'areal'
        write(outunit,101) (work(strtar-1+j),j=1,dimsqr)
        write(outunit,100) 'aimag'
        write(outunit,101) (work(strtai-1+j),j=1,dimsqr)
        write(outunit,100) 'breal'
        write(outunit,101) (work(strtbr-1+j),j=1,dimsqr)
        write(outunit,100) 'bimag'
        write(outunit,101) (work(strtbi-1+j),j=1,dimsqr)
        write(outunit,100) 'alfarb'
        write(outunit,101) (work(alfarb-1+j),j=1,dimreg)
        write(outunit,100) 'alfaib'
        write(outunit,101) (work(alfaib-1+j),j=1,dimreg)
        write(outunit,100) 'betab'
        write(outunit,101) (work(betab-1+j),j=1,dimreg)
        call cmatpr(q,ldq,dimreg,dimreg,'q after zqzvl1')
        call cmatpr(ph,ldp,dimreg,dimreg,'ph after zqzvl1')
      endif     
c
c        if (idbg(2) .gt. 1) then
c            call cmatpr(q,ldq,dimreg,dimreg,'q from qz')
c            call cmatpr(ph,ldp,dimreg,dimreg,'ph from qz')
c        endif
        if (ierr.ne.0) return
c
c     copy the real and imaginary parts of the qz-decomposition
c     to a and b, respectively
c                                                          
      call ricopy(a(rowb,colb), ldab, dimreg, work(strtar),
     *            work(strtai))
      call ricopy(b(rowb,colb), ldab, dimreg, work(strtbr),
     *            work(strtbi))
c
      if (ldebug) then
        call cmatpr(a(rowb,colb),ldab,dimreg,dimreg,'a after qz')
        call cmatpr(b(rowb,colb),ldab,dimreg,dimreg,'b after qz')
      endif
      return
      end

      subroutine rcopy(a, lda, dimreg, acopy)
c     implicit none
c
c**** formal parameter declarations
      integer lda, dimreg
      complex*16 a(lda,*)
      real*8 acopy(*)
c
c***  copy the real parts of a to the real vector acopy
c
c     this version dated june 10, 1987
c     authors: jim demmel and bo kagstrom
c     
c***  internal variables
      integer i, j
c
      do 20 i = 1, dimreg
         do 10 j = 1, dimreg
              acopy(i + (j - 1) * dimreg) = dreal(a(i,j))
   10    continue
   20 continue
      return
      end 
 
      subroutine icopy(a, lda, dimreg, acopy)
c     implicit none
c**** formal parameter declarations
      integer lda, dimreg
      complex*16 a(lda,*)
      real*8 acopy(*)
c
c***  copy the imaginary parts of a to the real vector acopy
c
c     this version dated june 10, 1987
c     authors: jim demmel and bo kagstrom
c     
c***  internal variables
c     
      integer i, j
      do 20 i = 1, dimreg
         do 10 j = 1, dimreg
              acopy(i + (j - 1) * dimreg) = dimag(a(i,j))
   10    continue
   20 continue
      return
      end 

      subroutine ricopy(a, lda, dimreg, arcopy, aicopy)
c     implicit none
c
c**** formal parameter declarations
      integer lda, dimreg
      complex*16 a(lda,*)
      real*8 arcopy(*), aicopy(*)
c
c***  copy arcopy and aicopy to the real and imaginary parts of a,
c     respectively
c
c     this version dated june 10, 1987
c     authors: jim demmel and bo kagstrom
c     
c***  internal variables
c     
      integer i, j
      do 20 i = 1, dimreg
         do 10 j = 1, dimreg
              a(i,j) = dcmplx(arcopy(i + (j - 1) * dimreg),
     *                        aicopy(i + (j - 1) * dimreg))
   10    continue
   20 continue
      return
      end 

c
c     ------------------------------------------------------------------
c
      subroutine zqzhs1(nm,n,ar,ai,br,bi,matz,z,ldz,matzl,zl,ldzl)
c
c
c     modified by demmel,6/23/86 to compute left, right transformations
c     in complex arithmetic
c  
      integer i,j,k,l,n,k1,lb,l1,nm,nk1,nm1
      real*8 ar(nm,n),ai(nm,n),br(nm,n),bi(nm,n)
      complex*16 z(ldz,*),zl(ldzl,*)
      complex*16 zu,tz,zll,zll1
      real*8 r,s,t,ti,u1,u2,xi,xr,yi,yr,rho,u1i
      logical matz,matzl
c
c     this subroutine is a complex analogue of the first step of the
c     qz algorithm for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of complex general matrices and
c     reduces one of them to upper hessenberg form with real (and non-
c     negative) subdiagonal elements and the other to upper triangular
c     form using unitary transformations.  it is usually followed by
c     cqzval  and possibly  cqzvec.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrices,
c
c        a=(ar,ai) contains a complex general matrix,
c
c        b=(br,bi) contains a complex general matrix,
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        matzl same as matz for left hand transformations
c
c     on output-
c
c        a has been reduced to upper hessenberg form.  the elements
c          below the first subdiagonal have been set to zero, and the
c          subdiagonal elements have been made real (and non-negative),
c
c        b has been reduced to upper triangular form.  the elements
c          below the main diagonal have been set to zero,
c
c        z contains the product of the right hand
c          transformations if matz has been set to .true.
c          otherwise, z is not referenced.
c
c        zl same as z for left transformations
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** initialize z **********
      if (.not. matz) go to 10
c
      do 3 i = 1, n
c
         do 2 j = 1, n
            z(i,j) = dcmplx(0.0d0,0.0d0)
    2    continue
c
         z(i,i) = dcmplx(1.0d0,0.0d0)
    3 continue
c
c     ********** initialize zl **********
      if (matzl) then
        do 300 i=1,n
          do 200 j=1,n
            zl(i,j)=0.
200       continue
          zl(i,i)=1.
300     continue
      endif
c
c     ********** reduce b to upper triangular form with
c                temporarily real diagonal elements **********
   10 if (n .le. 1) go to 170
      nm1 = n - 1
c
      do 100 l = 1, nm1
         l1 = l + 1
         s = 0.0
c
         do 20 i = l, n
            s = s + abs(br(i,l)) + abs(bi(i,l))
   20    continue
c
         if (s .eq. 0.0) go to 100
         rho = 0.0
c
         do 25 i = l, n
            br(i,l) = br(i,l) / s
            bi(i,l) = bi(i,l) / s
            rho = rho + br(i,l)**2 + bi(i,l)**2
   25    continue
c
         r = sqrt(rho)
         xr = abs(dcmplx(br(l,l),bi(l,l)))
         if (xr .eq. 0.0) go to 27
         rho = rho + xr * r
         u1 = -br(l,l) / xr
         u1i = -bi(l,l) / xr
         yr = r / xr + 1.0
         br(l,l) = yr * br(l,l)
         bi(l,l) = yr * bi(l,l)
         go to 28
c
   27    br(l,l) = r
         u1 = -1.0
         u1i = 0.0
c
   28    do 50 j = l1, n
            t = 0.0
            ti = 0.0
c
            do 30 i = l, n
               t = t + br(i,l) * br(i,j) + bi(i,l) * bi(i,j)
               ti = ti + br(i,l) * bi(i,j) - bi(i,l) * br(i,j)
   30       continue
c
            t = t / rho
            ti = ti / rho
c
            do 40 i = l, n
               br(i,j) = br(i,j) - t * br(i,l) + ti * bi(i,l)
               bi(i,j) = bi(i,j) - t * bi(i,l) - ti * br(i,l)
   40       continue
c
            xi = u1 * bi(l,j) - u1i * br(l,j)
            br(l,j) = u1 * br(l,j) + u1i * bi(l,j)
            bi(l,j) = xi
   50    continue
c
         do 80 j = 1, n
            t = 0.0
            ti = 0.0
c
            do 60 i = l, n
               t = t + br(i,l) * ar(i,j) + bi(i,l) * ai(i,j)
               ti = ti + br(i,l) * ai(i,j) - bi(i,l) * ar(i,j)
   60       continue
c
            t = t / rho
            ti = ti / rho
c
            do 70 i = l, n
               ar(i,j) = ar(i,j) - t * br(i,l) + ti * bi(i,l)
               ai(i,j) = ai(i,j) - t * bi(i,l) - ti * br(i,l)
   70       continue
c
            xi = u1 * ai(l,j) - u1i * ar(l,j)
            ar(l,j) = u1 * ar(l,j) + u1i * ai(l,j)
            ai(l,j) = xi
c
c        update zl
         if (matzl) then
           t=0.
           ti=0.
           do 600 i=l,n
             t= t + br(i,l)*dreal(zl(i,j)) + bi(i,l)*dimag(zl(i,j))
             ti=ti+ br(i,l)*dimag(zl(i,j))- bi(i,l)*dreal(zl(i,j))
600        continue
           tz=dcmplx(t/rho,ti/rho)
           do 700 i=l,n
             zl(i,j)=zl(i,j)-tz*dcmplx(br(i,l),bi(i,l))
700        continue
           zl(l,j)=zl(l,j)*dcmplx(u1,-u1i)
         endif
80       continue
c
         br(l,l) = r * s
         bi(l,l) = 0.0
c
         do 90 i = l1, n
            br(i,l) = 0.0
            bi(i,l) = 0.0
   90    continue
c
  100 continue
c     ********** reduce a to upper hessenberg form with real subdiagonal
c                elements, while keeping b triangular **********
      do 160 k = 1, nm1
         k1 = k + 1
c     ********** set bottom element in k-th column of a real **********
         if (ai(n,k) .eq. 0.0) go to 105
         r = abs(dcmplx(ar(n,k),ai(n,k)))
         u1 = ar(n,k) / r
         u1i = ai(n,k) / r
         ar(n,k) = r
         ai(n,k) = 0.0
c
         do 103 j = k1, n
            xi = u1 * ai(n,j) - u1i * ar(n,j)
            ar(n,j) = u1 * ar(n,j) + u1i * ai(n,j)
            ai(n,j) = xi
  103    continue
c
c        update zl
         if (matzl) then
           do 1030 j=1,n
             zl(n,j)=zl(n,j)*dcmplx(u1,-u1i)
1030       continue
         endif
c
         xi = u1 * bi(n,n) - u1i * br(n,n)
         br(n,n) = u1 * br(n,n) + u1i * bi(n,n)
         bi(n,n) = xi
  105    if (k .eq. nm1) go to 170
         nk1 = nm1 - k
c     ********** for l=n-1 step -1 until k+1 do -- **********
         do 150 lb = 1, nk1
            l = n - lb
            l1 = l + 1
c     ********** zero a(l+1,k) **********
            s = abs(ar(l,k)) + abs(ai(l,k)) + ar(l1,k)
            if (s .eq. 0.0) go to 150
            u1 = ar(l,k) / s
            u1i = ai(l,k) / s
            u2 = ar(l1,k) / s
            r = sqrt(u1*u1+u1i*u1i+u2*u2)
            u1 = u1 / r
            u1i = u1i / r
            u2 = u2 / r
            ar(l,k) = r * s
            ai(l,k) = 0.0
            ar(l1,k) = 0.0
c
            do 110 j = k1, n
               xr = ar(l,j)
               xi = ai(l,j)
               yr = ar(l1,j)
               yi = ai(l1,j)
               ar(l,j) = u1 * xr + u1i * xi + u2 * yr
               ai(l,j) = u1 * xi - u1i * xr + u2 * yi
               ar(l1,j) = u1 * yr - u1i * yi - u2 * xr
               ai(l1,j) = u1 * yi + u1i * yr - u2 * xi
  110       continue
c
c           update zl
            if (matzl) then
              zu=dcmplx(u1,-u1i)
              do 1100 j=1,n
                zll=zl(l,j)
                zll1=zl(l1,j)
                zl(l,j)= zu*zll+u2*zll1
                zl(l1,j)=conjg(zu)*zll1-u2*zll
1100          continue
            endif
c
            xr = br(l,l)
            br(l,l) = u1 * xr
            bi(l,l) = -u1i * xr
            br(l1,l) = -u2 * xr
c
            do 120 j = l1, n
               xr = br(l,j)
               xi = bi(l,j)
               yr = br(l1,j)
               yi = bi(l1,j)
               br(l,j) = u1 * xr + u1i * xi + u2 * yr
               bi(l,j) = u1 * xi - u1i * xr + u2 * yi
               br(l1,j) = u1 * yr - u1i * yi - u2 * xr
               bi(l1,j) = u1 * yi + u1i * yr - u2 * xi
  120       continue
c     ********** zero b(l+1,l) **********
            s = abs(br(l1,l1)) + abs(bi(l1,l1)) + abs(br(l1,l))
            if (s .eq. 0.0) go to 150
            u1 = br(l1,l1) / s
            u1i = bi(l1,l1) / s
            u2 = br(l1,l) / s
            r = sqrt(u1*u1+u1i*u1i+u2*u2)
            u1 = u1 / r
            u1i = u1i / r
            u2 = u2 / r
            br(l1,l1) = r * s
            bi(l1,l1) = 0.0
            br(l1,l) = 0.0
c
            do 130 i = 1, l
               xr = br(i,l1)
               xi = bi(i,l1)
               yr = br(i,l)
               yi = bi(i,l)
               br(i,l1) = u1 * xr + u1i * xi + u2 * yr
               bi(i,l1) = u1 * xi - u1i * xr + u2 * yi
               br(i,l) = u1 * yr - u1i * yi - u2 * xr
               bi(i,l) = u1 * yi + u1i * yr - u2 * xi
  130       continue
c
            do 140 i = 1, n
               xr = ar(i,l1)
               xi = ai(i,l1)
               yr = ar(i,l)
               yi = ai(i,l)
               ar(i,l1) = u1 * xr + u1i * xi + u2 * yr
               ai(i,l1) = u1 * xi - u1i * xr + u2 * yi
               ar(i,l) = u1 * yr - u1i * yi - u2 * xr
               ai(i,l) = u1 * yi + u1i * yr - u2 * xi
  140       continue
c
            if (.not. matz) go to 150
c
            zu=dcmplx(u1,-u1i)
            do 145 i = 1, n
              zll1=z(i,l1)
              zll=z(i,l)
              z(i,l1)=zu*zll1+u2*zll
              z(i,l)= conjg(zu)*zll-u2*zll1
  145       continue
c
  150    continue
c
  160 continue
c
  170 return
c     ********** last card of zqzhes **********
      end
c
c     ------------------------------------------------------------------
c
      subroutine zqzvl1(nm,n,ar,ai,br,bi,eps1,alfr,alfi,beta,
     x              matz,z,ldz,matzl,zl,ldzl,ierr)
c
c     modified by demmel, 6/23/86 to compute left and right 
c     transformations using complex arithmetic
c
      integer i,j,k,l,n,en,k1,k2,ll,l1,na,nm,its,km1,lm1,
     x        enm2,ierr,lor1,enorn
      real*8 ar(nm,n),ai(nm,n),br(nm,n),bi(nm,n),alfr(n),alfi(n),
     x       beta(n)
      complex*16 z(ldz,*),zl(ldzl,*)
      complex*16 zu,zll,zll1
      real*8 r,s,a1,a2,ep,sh,u1,u2,xi,xr,yi,yr,ani,a1i,a33,a34,a43,a44,
     x       bni,b11,b33,b44,shi,u1i,a33i,a34i,a43i,a44i,b33i,b44i,
     x       epsa,epsb,eps1,anorm,bnorm,b3344,b3344i
      integer max0
      logical matz,matzl
      complex*16 z3
c
c
c
c
c
c     this subroutine is a complex analogue of steps 2 and 3 of the
c     qz algorithm for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
c     as modified in technical note nasa tn e-7305(1973) by ward.
c
c     this subroutine accepts a pair of complex matrices, one of them
c     in upper hessenberg form and the other in upper triangular form,
c     the hessenberg matrix must further have real subdiagonal elements.
c     it reduces the hessenberg matrix to triangular form using
c     unitary transformations while maintaining the triangular form
c     of the other matrix and further making its diagonal elements
c     real and non-negative.  it then returns quantities whose ratios
c     give the generalized eigenvalues.  it is usually preceded by
c     cqzhes  and possibly followed by  cqzvec.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrices,
c
c        a=(ar,ai) contains a complex upper hessenberg matrix
c          with real subdiagonal elements,
c
c        b=(br,bi) contains a complex upper triangular matrix,
c
c        eps1 is a tolerance used to determine negligible elements.
c          eps1 = 0.0 (or negative) may be input, in which case an
c          element will be neglected only if it is less than roundoff
c          error times the norm of its matrix.  if the input eps1 is
c          positive, then an element will be considered negligible
c          if it is less than eps1 times the norm of its matrix.  a
c          positive value of eps1 may result in faster execution,
c          but less accurate results,
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise,
c
c        z=(zr,zi) contains, if matz has been set to .true., the
c          transformation matrix produced in the reduction
c          by  cqzhes, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output-
c
c        a has been reduced to upper triangular form.  the elements
c          below the main diagonal have been set to zero,
c
c        b is still in upper triangular form, although its elements
c          have been altered.  in particular, its diagonal has been set
c          real and non-negative.  the location br(n,1) is used to
c          store eps1 times the norm of b for later use by  cqzvec,
c
c        alfr and alfi contain the real and imaginary parts of the
c          diagonal elements of the triangularized a matrix,
c
c        beta contains the real non-negative diagonal elements of the
c          corresponding b.  the generalized eigenvalues are then
c          the ratios ((alfr+i*alfi)/beta),
c
c        z contains the product of the right hand transformations
c          (for both steps) if matz has been set to .true.,
c
c        ierr is set to
c          zero       for normal return,
c          j          if ar(j,j-1) has not become
c                     zero after 50 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     ********** compute epsa,epsb **********
      anorm = 0.0
      bnorm = 0.0
c
      do 30 i = 1, n
         ani = 0.0
         if (i .ne. 1) ani = abs(ar(i,i-1))
         bni = 0.0
c
         do 20 j = i, n
            ani = ani + abs(ar(i,j)) + abs(ai(i,j))
            bni = bni + abs(br(i,j)) + abs(bi(i,j))
   20    continue
c
         if (ani .gt. anorm) anorm = ani
         if (bni .gt. bnorm) bnorm = bni
   30 continue
c
      if (anorm .eq. 0.0) anorm = 1.0
      if (bnorm .eq. 0.0) bnorm = 1.0
      ep = eps1
      if (ep .gt. 0.0) go to 50
c     ********** compute roundoff level if eps1 is zero **********
      ep = 1.0d0
   40 ep = ep / 2.0d0
      if (1.0d0 + ep .gt. 1.0d0) go to 40
   50 epsa = ep * anorm
      epsb = ep * bnorm
c     ********** reduce a to triangular form, while
c                keeping b triangular **********
      lor1 = 1
      enorn = n
      en = n
c     ********** begin qz step **********
   60 if (en .eq. 0) go to 1001
      if (.not. matz) enorn = en
      its = 0
      na = en - 1
      enm2 = na - 1
c     ********** check for convergence or reducibility.
c                for l=en step -1 until 1 do -- **********
   70 do 80 ll = 1, en
         lm1 = en - ll
         l = lm1 + 1
         if (l .eq. 1) go to 95
         if (abs(ar(l,lm1)) .le. epsa) go to 90
   80 continue
c
   90 ar(l,lm1) = 0.0
c     ********** set diagonal element at top of b real **********
   95 b11 = abs(dcmplx(br(l,l),bi(l,l)))
      if (b11     .eq. 0.0) go to 98
      u1 = br(l,l) / b11
      u1i = bi(l,l) / b11
c
      do 97 j = l, enorn
         xi = u1 * ai(l,j) - u1i * ar(l,j)
         ar(l,j) = u1 * ar(l,j) + u1i * ai(l,j)
         ai(l,j) = xi
         xi = u1 * bi(l,j) - u1i * br(l,j)
         br(l,j) = u1 * br(l,j) + u1i * bi(l,j)
         bi(l,j) = xi
   97 continue
c
c     update zl
      if (matzl) then
        do 970 j=1,n
          zl(l,j)=zl(l,j)*dcmplx(u1,-u1i)
970     continue
      endif
c
      bi(l,l) = 0.0
   98 if (l .ne. en) go to 100
c     ********** 1-by-1 block isolated **********
      alfr(en) = ar(en,en)
      alfi(en) = ai(en,en)
      beta(en) = b11
      en = na
      go to 60
c     ********** check for small top of b **********
  100 l1 = l + 1
      if (b11 .gt. epsb) go to 120
      br(l,l) = 0.0
      s = abs(ar(l,l)) + abs(ai(l,l)) + abs(ar(l1,l))
      u1 = ar(l,l) / s
      u1i = ai(l,l) / s
      u2 = ar(l1,l) / s
      r = sqrt(u1*u1+u1i*u1i+u2*u2)
      u1 = u1 / r
      u1i = u1i / r
      u2 = u2 / r
      ar(l,l) = r * s
      ai(l,l) = 0.0
c
      do 110 j = l1, enorn
         xr = ar(l,j)
         xi = ai(l,j)
         yr = ar(l1,j)
         yi = ai(l1,j)
         ar(l,j) = u1 * xr + u1i * xi + u2 * yr
         ai(l,j) = u1 * xi - u1i * xr + u2 * yi
         ar(l1,j) = u1 * yr - u1i * yi - u2 * xr
         ai(l1,j) = u1 * yi + u1i * yr - u2 * xi
         xr = br(l,j)
         xi = bi(l,j)
         yr = br(l1,j)
         yi = bi(l1,j)
         br(l1,j) = u1 * yr - u1i * yi - u2 * xr
         br(l,j) = u1 * xr + u1i * xi + u2 * yr
         bi(l,j) = u1 * xi - u1i * xr + u2 * yi
         bi(l1,j) = u1 * yi + u1i * yr - u2 * xi
  110 continue
c
c     update zl
      if (matzl) then
        zu=dcmplx(u1,-u1i)
        do 1110 j=1,n
          zll=zl(l,j)
          zll1=zl(l1,j)
          zl(l,j)=zll*zu+zll1*u2
          zl(l1,j)=zll1*conjg(zu)-zll*u2
1110    continue
      endif
c
      lm1 = l
      l = l1
      go to 90
c     ********** iteration strategy **********
  120 if (its .eq. 50) go to 1000
      if (its .eq. 10) go to 135
c     ********** determine shift **********
      b33 = br(na,na)
      b33i = bi(na,na)
      if (abs(dcmplx(b33,b33i)) .ge. epsb) go to 122
      b33 = epsb
      b33i = 0.0
  122 b44 = br(en,en)
      b44i = bi(en,en)
      if (abs(dcmplx(b44,b44i)) .ge. epsb) go to 124
      b44 = epsb
      b44i = 0.0
  124 b3344 = b33 * b44 - b33i * b44i
      b3344i = b33 * b44i + b33i * b44
      a33 = ar(na,na) * b44 - ai(na,na) * b44i
      a33i = ar(na,na) * b44i + ai(na,na) * b44
      a34 = ar(na,en) * b33 - ai(na,en) * b33i
     x    - ar(na,na) * br(na,en) + ai(na,na) * bi(na,en)
      a34i = ar(na,en) * b33i + ai(na,en) * b33
     x     - ar(na,na) * bi(na,en) - ai(na,na) * br(na,en)
      a43 = ar(en,na) * b44
      a43i = ar(en,na) * b44i
      a44 = ar(en,en) * b33 - ai(en,en) * b33i - ar(en,na) * br(na,en)
      a44i = ar(en,en) * b33i + ai(en,en) * b33 - ar(en,na) * bi(na,en)
      sh = a44
      shi = a44i
      xr = a34 * a43 - a34i * a43i
      xi = a34 * a43i + a34i * a43
      if (xr .eq. 0.0 .and. xi .eq. 0.0) go to 140
      yr = (a33 - sh) / 2.0
      yi = (a33i - shi) / 2.0
      z3 = sqrt(dcmplx(yr**2-yi**2+xr,2.0*yr*yi+xi))
      u1 = dreal(z3)
      u1i = dimag(z3)
      if (yr * u1 + yi * u1i .ge. 0.0) go to 125
      u1 = -u1
      u1i = -u1i
  125 z3 = (dcmplx(sh,shi) - dcmplx(xr,xi) / dcmplx(yr+u1,yi+u1i))
     x   / dcmplx(b3344,b3344i)
      sh = dreal(z3)
      shi = dimag(z3)
      go to 140
c     ********** ad hoc shift **********
  135 sh = ar(en,na) + ar(na,enm2)
      shi = 0.0
c     ********** determine zeroth column of a **********
  140 a1 = ar(l,l) / b11 - sh
      a1i = ai(l,l) / b11 - shi
      a2 = ar(l1,l) / b11
      its = its + 1
      if (.not. matz) lor1 = l
c     ********** main loop **********
      do 260 k = l, na
         k1 = k + 1
         k2 = k + 2
         km1 = max0(k-1,l)
c     ********** zero a(k+1,k-1) **********
         if (k .eq. l) go to 170
         a1 = ar(k,km1)
         a1i = ai(k,km1)
         a2 = ar(k1,km1)
  170    s = abs(a1) + abs(a1i) + abs(a2)
         u1 = a1 / s
         u1i = a1i / s
         u2 = a2 / s
         r = sqrt(u1*u1+u1i*u1i+u2*u2)
         u1 = u1 / r
         u1i = u1i / r
         u2 = u2 / r
c
         do 180 j = km1, enorn
            xr = ar(k,j)
            xi = ai(k,j)
            yr = ar(k1,j)
            yi = ai(k1,j)
            ar(k,j) = u1 * xr + u1i * xi + u2 * yr
            ai(k,j) = u1 * xi - u1i * xr + u2 * yi
            ar(k1,j) = u1 * yr - u1i * yi - u2 * xr
            ai(k1,j) = u1 * yi + u1i * yr - u2 * xi
            xr = br(k,j)
            xi = bi(k,j)
            yr = br(k1,j)
            yi = bi(k1,j)
            br(k,j) = u1 * xr + u1i * xi + u2 * yr
            bi(k,j) = u1 * xi - u1i * xr + u2 * yi
            br(k1,j) = u1 * yr - u1i * yi - u2 * xr
            bi(k1,j) = u1 * yi + u1i * yr - u2 * xi
  180    continue
c
c        update zl
         if (matzl) then
           zu=dcmplx(u1,-u1i)
           do 1800 j=1,n
             zll=zl(k,j)
             zll1=zl(k1,j)
             zl(k,j)=zu*zll+u2*zll1
             zl(k1,j)=conjg(zu)*zll1-u2*zll
1800       continue
         endif
c
         if (k .eq. l) go to 240
         ai(k,km1) = 0.0
         ar(k1,km1) = 0.0
         ai(k1,km1) = 0.0
c     ********** zero b(k+1,k) **********
  240    s = abs(br(k1,k1)) + abs(bi(k1,k1)) + abs(br(k1,k))
         u1 = br(k1,k1) / s
         u1i = bi(k1,k1) / s
         u2 = br(k1,k) / s
         r = sqrt(u1*u1+u1i*u1i+u2*u2)
         u1 = u1 / r
         u1i = u1i / r
         u2 = u2 / r
         if (k .eq. na) go to 245
         xr = ar(k2,k1)
         ar(k2,k1) = u1 * xr
         ai(k2,k1) = -u1i * xr
         ar(k2,k) = -u2 * xr
c
  245    do 250 i = lor1, k1
            xr = ar(i,k1)
            xi = ai(i,k1)
            yr = ar(i,k)
            yi = ai(i,k)
            ar(i,k1) = u1 * xr + u1i * xi + u2 * yr
            ai(i,k1) = u1 * xi - u1i * xr + u2 * yi
            ar(i,k) = u1 * yr - u1i * yi - u2 * xr
            ai(i,k) = u1 * yi + u1i * yr - u2 * xi
            xr = br(i,k1)
            xi = bi(i,k1)
            yr = br(i,k)
            yi = bi(i,k)
            br(i,k1) = u1 * xr + u1i * xi + u2 * yr
            bi(i,k1) = u1 * xi - u1i * xr + u2 * yi
            br(i,k) = u1 * yr - u1i * yi - u2 * xr
            bi(i,k) = u1 * yi + u1i * yr - u2 * xi
  250    continue
c
         bi(k1,k1) = 0.0
         br(k1,k) = 0.0
         bi(k1,k) = 0.0
         if (.not. matz) go to 260
c
         zu=dcmplx(u1,-u1i)
         do 255 i = 1, n
           zll=z(i,k)
           zll1=z(i,k1)
           z(i,k)=conjg(zu)*zll-u2*zll1
           z(i,k1)=zu*zll1+u2*zll
  255    continue
c
  260 continue
c     ********** set last a subdiagonal real and end qz step **********
      if (ai(en,na) .eq. 0.0) go to 70
      r = abs(dcmplx(ar(en,na),ai(en,na)))
      u1 = ar(en,na) / r
      u1i = ai(en,na) / r
      ar(en,na) = r
      ai(en,na) = 0.0
c
      do 270 j = en, enorn
         xi = u1 * ai(en,j) - u1i * ar(en,j)
         ar(en,j) = u1 * ar(en,j) + u1i * ai(en,j)
         ai(en,j) = xi
         xi = u1 * bi(en,j) - u1i * br(en,j)
         br(en,j) = u1 * br(en,j) + u1i * bi(en,j)
         bi(en,j) = xi
  270 continue
c
c     update zl
      if (matzl) then
        zu=dcmplx(u1,-u1i)
        do 2700 j=1,n
          zl(en,j)=zu*zl(en,j)
2700    continue
      endif
c
      go to 70
c     ********** set error -- bottom subdiagonal element has not
c                become negligible after 50 iterations **********
 1000 ierr = en
c     ********** save epsb for use by cqzvec **********
 1001 if (n .gt. 1) br(n,1) = 0.
c     if (n .gt. 1) br(n,1) = epsb
      return
c     ********** last card of zqzval **********
      end



      subroutine rcsvdc(x, ldx, m, n, s, e, u, ldu, v, ldv,
     *                  opt, epsu, gap, cnull, rnull, del,
     *                  work, job, info)
c
c     implicit none
c**** debug space
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
c**** formal parameter declarations
      integer ldx,m,n,ldu,ldv,cnull,rnull,job,info
      real*8 epsu, gap, del
      complex*16 x(ldx,n),s(m),e(m),u(ldu,1),v(ldv,n),work(n)
      character*(*) opt
c
c*********************************************************************
c
c     rcsvdc computes the singular value decomposition (svd)
c     of a m by n matrix x, and its numerical column and row
c     nullities, respectively. the diagonal elements s(i) are
c     the singular values. the user controls the ordering and 
c     the placing of the singular values. 
c     the columns of the unitary matrices u and v correspond
c     to the left and right singular vectors, respectively.
c
c     on entry
c
c         x         complex(ldx,n), where ldx>=m.
c
c         ldx       integer
c                   ldx is the leading dimension of the array x.
c
c         m         integer
c                   m is the number of rows of x.
c
c         n         integer
c                   n is the number of columns of x.            
c
c         ldu       integer
c                   ldu is the leading dimension of the array u.
c
c         ldv       integer
c                   ldv is the leading dimension of the array v.
c
c         work      complex(n)
c                   work is a scratch array.
c
c         job       integer
c                   job controls the computations to be done. it has
c                   the decimal expansion abcd with the following
c                   meaning                   
c                      a=0    do not compute the left singular vectors.
c                      a=1    return the m left singular vextors in u.
c                      a=2    return the first min(m,n) left singular
c                             vectors in u.
c                      b=0    do not compute the right singular vectors.
c                      b=1    return the right singular vectors in v.
c                      c=0    singular values are ordered in decreasing
c                             order.
c                      c=1    singular values are ordered in increasing
c                             order.
c                      d=0    diagonal of singular values starts in
c                             position (1,1).
c                      d=1    diagonal of singular values ends in
c                             position (m,n).
c
c     on return
c
c         s         complex(mm), where mm=min(m+1,n).     ??????
c                   the first min(m,n) entries of s contain the 
c                   singular values of x.      
c                                
c         e         complex(m)
c                   e contains the subdiagonal from computing
c                   the svd. should ordinarily be zeros.
c
c         u         complex(ldu,k), where ldu>=m.
c                   if joba=1 then k=m, if joba=2 then k=min(m,n).
c                   u contains the matrix of left singular
c                   vectors of x.
c                   u is not referenced if joba=0. if m<=n or if
c                   joba=2, then u may be identified with x in the
c                   subroutine call.
c
c         v         complex(ldv,n), where ldv>=n.
c                   v contains the matrix of right singular
c                   vectors of x.
c                   v is not referenced if jobb=0. if n<=m,
c                   then v may be identified with x in the 
c                   subroutine call.
c
c         cnull     integer
c                   cnull contains the numerical column nullity of x.
c
c         rnull     integer
c                   rnull contains the numerical row nullity of x.
c
c         del       real*8
c                   del contains the squareroot of the sum of the
c                   squares of the singular values interpreted as zeros.
c
c         info      integer
c                   info tells the user what has been done.
c                   info=0, all the singular values and
c                   vectors are correct. if info .ne.o, then
c                   cnull, rnull and del contain no meaningful
c                   information. for more details see the
c                   linpack routine csvdc.
c
c********************************************************************
c
c         this version dated june 13, 1987
c         authors: jim demmel and bo kagstrom 
c
c*****    rcsvdc uses the following functions and subroutines
c
c         linpack    zsvdc
c         blas       zswap
c
c*****    internal variables
c
c
c***  if idbg(6) .eq. 0 then debug output is switched off
c     on input info contains the product of the row and
c     column dimensions of the original a and b
c
      integer       jobu, ncu, jobx, nsvd, i, j, n1, mn, mpn, k
      logical       wantu, wantv, incr, posmn, ldebug
      real*8        t1, t2
      complex*16    cell
c
c     save m*n (=info) in mpn
      mpn = info
c
c*****    determine what is to be computed
      ldebug = idbg(6) .ne. 0
      jobu = job/1000
      wantu = jobu .ne. 0
c
c     ncu is the number of columns in u
      ncu = m
      if (jobu .eq. 2) ncu = min0(m,n)
      wantv = mod(job,1000)/100 .ne. 0
      incr = mod(job,100)/10 .ne. 0
      posmn = mod(job,10) .ne. 0
c
c*****    compute the svd of x
c     singular values in decraesing order
c
      jobx = job/100
      call zsvdc(x,ldx,m,n,s,e,u,ldu,v,ldv,work,jobx,info)
c**** 6/18/87
c      if( info .ne. 0)return
       if (info .ne. 0) then
         if (ldebug) write(outunit,101) info
 101     format('rcsvdc - after zsvd, info= ',i4)
         return
       endif
c
c*****    compute the column and row nullities of x
c         n1 = number of singular values interpreted as zeros
c
c         we seek n1 so that
c             s(nsvd-n1) >= t2 > t1 >= s(nsvd-n1+1 )
c         if this relation does not hold n1 is decreased by one
c         until we have a gap t2/t1 (=gap) between the singular
c         values we interpret as zeros and the others.
c*****    works only if singular values in increasing order
c
      t1 = epsu
      t2 = gap * t1
      if (ldebug) then
           write(outunit,100) 't1= ', t1, 't2= ', t2
  100      format(t5,a,d12.5,tr5,a,d12.5)
      endif
c

      nsvd = min0(m,n)
c
c**** shall we compute cnull and rnull or not?
c
      if (opt .eq. 'cind') then
c
c**** note that if only one singular value then we interpret it
c     as zero if it is less than t2
         if (nsvd .eq. 1) then
            n1 = 0
            if ( abs(s(1)) .le. t2 ) n1 = 1
         else
c
            do 20 i = nsvd, 1 , -1
               if (abs(s(i)) .ge. t1) go to 25
   20       continue
            n1 = nsvd
            go to 35
   25       continue
            if ( i .ge. 1) then
               if (abs(s(i)) .gt. t2) go to 30
               i = i - 1
               go to 25
            endif
   30       continue
            n1 = nsvd - i
         endif
   35    continue
c
         if ( m .ge. n ) then
            cnull = n1
            rnull = (m - n) + n1
         else
            cnull = (n - m) + n1
            rnull = n1
         endif
       else
c
c        cnull and rnull are alreday known from earlier computations
         if ( m.ge. n) then
             n1 = cnull
         else
             n1 = cnull - (n -m)
         endif
       endif
       del = 0.
       do 40 i = nsvd, (nsvd - n1 + 1), -1
c*     accumulate square root of sum of squares
           call upddel(del, abs(s(i)))
   40  continue
c
      if (incr) then
c
c      reorder the singular values (and  the corresponding vectors)
c      into increasing order
c
        do 50 i = 1, nsvd/2
           j = nsvd - i + 1
           if (wantu)
     *       call zswap(m,u(1,i),1,u(1,j),1)
           if (wantv)
     *       call zswap(n,v(1,i),1,v(1,j),1)
           cell = s(i)
           s(i) = s(j)
           s(j) = cell
   50   continue
      endif
c**       incr
      if (posmn) then
c
c      move the columns of u and v, such that the diagonal of
c      singular values (of u'*x*v where '=transpose conjugate)
c      ends at position (m,n)
c
        if ( (jobu .eq. 1) .and.  (m .gt. n)) then
c
c         move the last m-n columns of u to the first positions in u,
c         and adjust the remaining col's accordingly.
c         (remember the case when a=2, ncu= number of col's of u)
c
          mn = m - n
          do 70 k = 1, mn
              do 70 i = 1, m
                 cell = u(i,ncu)
                  do 60 j = ncu, 2, -1
                     u(i,j) = u(i,j-1)
   60         continue
              u(i,1) = cell
   70     continue
        endif
c**         (jobu = 1)
        if (wantv .and. (m .lt. n)) then
c
c         move the last n-m columns of v to the first positions in v
c         and adjust the remaining columns accordingly.
c         (n= the of col's of v)
c
          mn = n - m
          do 90 k = 1, mn
             do 90 i = 1, n
                cell = v(i,n)
                do 80 j = n, 2, -1
                   v(i,j) = v(i,j-1)
   80        continue
             v(i,1) = cell
   90     continue
        endif
c**         wantv
      endif
c**       posmn
c
      return
      end




c     In this file June 13, 1987: reordr, exchng, cgiv, zcsrot
c     
c
      subroutine reordr (a, b, ldab, m, n, rowb, colb, rowe, cole,
     *                   ftest, ndim, ind, pp, ldpp, qq, ldqq)
c
c     implicit none
c***  debug space
      common /debug2/ idbg(20), outunit
      external ftest
      integer idbg, outunit
c 
      integer ldab, m, n, rowb, colb, rowe, cole, ftest
      integer ndim, ind(*), ldpp, ldqq
      complex*16 a(ldab,*), b(ldab,*), pp(ldpp,*), qq(ldqq,*)
c
c***********************************************************************
c     given that the specified regular part of a - lambda*b is in
c     upper triangular form reordr reorders the 1 by 1 diagonal blocks
c     (the generalized eigenvalues) by constructing
c     equivalence transformations (pairs of left and right givens
c     transformations). the givens transformations that perform the
c     reordering are accumulated in the left and right transformation
c     matrices pp and qq, respectively. normally pp and qq result
c     from previous reductions or are initialized to the identity 
c     matrix before the call.
c
c     after the reordering the eigenvalues specified by the function
c     ftest (provided by the user) appear at the top north-west corner
c     of the specified regular part of a - lambda*b.
c     if ndim is the number of eigenvalues in the spectrum specified
c     by ftest then the rowb+ndim-1 first columns of pp, and the 
c     colb+ndim-1 first columns of qq, respectively, 
c     span a pair of reducing subspaces corresponding to this 
c     part of the spectrum of a - lambda*b. for algorithmic details of
c     the reordering of eigenvalues see p. van dooren: algorithm 590:
c     dsubsp and exchng, fortran routines for computing deflating 
c     subspaces with specified spectrum, acm toms, vol.4, 1982,
c     pp 376-382
c 
c     if idb(7) .eq. 0 then debug output is switched off
c
c**** formal parameters
c 
c    on entry
c
c     a(ldab,*) complex*16, input matrix a in upper triangular form
c
c     b(ldab,*) complex*16, input matrix b in upper triangular form
c
c     ldab    integer, leading dimension of a and b
c
c     m       integer, current row dimension of a and b
c
c     n       integer, current column dimension of a and b
c
c     rowb    integer, first row of the regular part of a-lambda*b
c
c     colb    integer, first column of the regular part of a-lambda*b
c
c     rowe    integer, last row of the regular part of a-lambda*b
c
c     cole    integer, last column of the regular part of a-lambda*b
c
c     ftest(alpha, beta)  integer function describing the spectrum
c             of the deflating subspace to be computed. if alpha/beta
c             is in that spectrum then ftest = 1, otherwise ftest = -1.
c
c     ldpp    integer, leading dimension of pp
c
c     ldqq    integer, leading dimension of qq
c
c    on exit
c
c     ndim    integer, the dimension of the computed pair of
c             deflating subspace
c
c     ind(*)  integer array, working array of dimension at least
c             min(rowe-rowb+1)
c
c     pp(ldpp,*) complex*16, array, the unitary right hand transformation
c             matrix of order m by m.
c             accumulates all right hand givens transformations.
c
c     qq(ldqq,*) complex*16, the unitary left hand transformation
c             matrix of order n by n.
c             accumulates all left hand givens transformations.
c
c     a(ldab,*) in upper tringular form with reordered diagonal
c               elements
c
c     b(ldab,*) in upper triangular form with reordered diagonal
c               elements
c     note: the reordered eigenvalues are a(i,i)/b(i,i) (see also above)
c
cc************************************************************************
c
c**** this version dated 14 june, 1987
c     authors: jim demmel and bo kagstrom
c
c**** reordr uses the following functions and subroutines
c     cmatpr, exchng
c     ftest (user written)
c
c**** internal variables
      integer dimr, i, k, j, inside, rfirst, kfirst, jj, nswap
      integer indk, ii
      logical ldebug
c
c     set debug flag
      ldebug = idbg(7) .ne. 0
c
      if (ldebug) then
         write(outunit, 2005) 'eigenvalues before reordering'
         do 770 i = rowb, rowe
           j = colb + i - rowb
           if (abs(b(i ,j)) .eq. 0. ) then
               write(outunit, 2005) 'infinite eigenvalue',a(i,j), b(i,j)
 2005          format(t5,a,4d15.5)
           else
               write(outunit, 2005) 'eigenvalue=', a(i,j)/b(i,j)
           endif
  770    continue
      endif
c***
c
      dimr = rowe - rowb +1
      if (ldebug) then
           write(outunit, 500) 'dimr=', dimr
  500      format(t5,a,i3)
      endif
c
      if( dimr .ge. 2) then
c**** search through the eigenvalues and note down in ind(*) which
c     eigenvalues are in the spectrum determined by ftest
c
      ndim = 0
      dimr = 0
      do 10 i = rowb, rowe
         dimr = dimr + 1
         inside = ftest( a(i, colb + dimr - 1), b(i, colb + dimr - 1))
         if ( inside .eq. 1 ) ndim = ndim + 1
         ind(dimr) = inside
   10 continue
      if (ldebug) then
          write(outunit, 700) 'ind(*) before reordering',
     *                  (ind(i), i = 1, dimr)
  700     format(t5,a,20i3)
      endif
c
c**** reorder the blocks (eigenvalues) such that those that belong
c     to the specified spectrum appear first at the top north-west corner
c     of the specified regular part of a-lambda*b
c
      do 100 i = 1, dimr
         if ( ind(i) .lt. 0) then
c
c           search for the first block to be moved ( first ind(k)
c           that is positive)
            do 60 k = i + 1, dimr
                 if ( ind(k) .gt. 0) go to 70
   60       continue
c           no more blocks to test or to move, go to exit
            go to 110
         else
c           continue the search
            go to 100
         endif
c
c        make k-i interchanges so that block k appear before block i
   70    continue
         nswap = k - i
         if (ldebug) write(outunit, 500) 'nswap=',nswap
         indk = ind(k)
         do 80 j =1, nswap
              jj = k - j
              rfirst = rowb + i - 1 + nswap - j
              kfirst = colb + i - 1 + nswap - j
              call exchng(a, b, ldab, m, n, rfirst, kfirst,
     *                    pp, ldpp, qq, ldqq)
              ind(jj + 1) = ind(jj)
   80    continue
         ind(i) = indk
         if (ldebug) then
              write(outunit, 700) 'ind(*) after do 80',
     *                      (ind(ii),ii=1,dimr)
         endif
c
c     continue to search for eigenvalues that should be reordered
  100 continue
c
c     exit
  110 continue
      if (ldebug) then
           write(outunit, 700) 'final ind(*) from reorder',
     *                   (ind(ii),ii=1,dimr)
      endif
c
c     end of if ( dimr .ge. 2)
      endif
c
      if (idbg(2) .gt. 1) then
            call cmatpr(qq,ldqq,n,n,'qq after reordr')
            call cmatpr(pp,ldpp,m,m,'pp after reordr')
      endif
c
      if (ldebug) then
         write(outunit, 2005) ' eigenvalues at exit from reordr'
         do 75 i = rowb, rowe
           j = colb + i - rowb
           if (abs(b(i ,j)) .eq. 0. ) then
               write(outunit, 2005) 'infinite eigenvalue',a(i,j), b(i,j)
           else
               write(outunit, 2005) 'eigenvalue=', a(i,j)/b(i,j)
           endif
   75    continue
       endif
c
      return
      end
 
      subroutine exchng(a, b, ldab, m, n, rowb, colb,
     *                  pp, ldpp, qq, ldqq)
c
c     implicit none
c***  debug space
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c
      integer ldab, m, n, rowb, colb, ldpp, ldqq
      complex*16 a(ldab,*), b(ldab,*), pp(ldpp,*), qq(ldqq,*)
c
c***********************************************************************
c     given that the regular part of a - lambda*b is on upper
c     triangular form exchng computes a unitary equivalence
c     transformation that exchanges the 1 by 1 diagonal blocks
c     at positions (rowb, colb) and (rowb+1, colb+1), respectively,
c     along with their generalized eigenvalues.
c     the givens rotations that perform the exchange are
c     accumulated in the left and right transformation matrices
c     pp and qq, respectively.
c
c     if idbg(8) .eq. 0 then debug output is switched off
c
c     formal parameters
c     
c    on entry
c
c     a(ldab,*) complex*16, input matrix a in upper triangular form
c
c     b(ldab,*) complex*16, input matrix b in upper triangular form
c
c     ldab    integer, leading dimension of a and b 
c
c     m       integer, current row dimension of a and b
c
c     n       integer, current column dimension of a and b
c
c     rowb    integer, first row of the regular part of a-lambda*b
c
c     colb    integer, forst column of the regular part of a-lambda*b
c
c     ldpp    integer, leading dimension of pp
c
c     ldqq    integer, leading dimension of qq
c
c    on exit
c
c     pp(ldpp,*) complex*16, the unitary right hand transformation
c             of order m by m
c
c     qq(ldqq,*) complex*16, the unitary left hand transformation
c             of order n by n
c
c     a(ldab,*)  in upper triangular form with two diagonal elements
c             exchanged
c
c     b(ldab,*)  in upper trinagular form with two diagonal elements
c             exchanged
cc
c************************************************************************
c
c**** this version dated june 14, 1986
c     authors: jim demmel and bo kagstrom
c
c**** exchng uses the following functions and subroutines
c     cgiv, cmatpr, zcsrot
c
c**** internal variables
      logical altb, ldebug
      integer rbp1, cbp1
      real*8    maxab1
      complex*16 sa1, sb1, f, g, s, c, ctemp
c
      ldebug = idbg(8) .ne. 0
c
      rbp1 = rowb + 1
      cbp1 = colb + 1
      if (ldebug) then
        write (outunit, 2000) 'results from exchng: rowb, colb', 
     *                        rowb, colb
 2000   format( t5, a, 2i3)
        write(outunit, 2000) 'eigenvalues before exchange'
        write(outunit, 2000) 'rbp1,cbp1=', rbp1, cbp1
        if (abs(b(rowb,colb)) .gt. 0.) then
              ctemp = a(rowb,colb)/b(rowb,colb)
              write(outunit, 3000) ctemp
        else
              write(outunit, 3000) a(rowb,colb), b(rowb,colb)
        endif
        if (abs(b(rbp1,cbp1)) .gt. 0.) then
             write(outunit, 3000) a(rbp1,cbp1), b(rbp1,cbp1)
             ctemp = a(rbp1,cbp1)/b(rbp1,cbp1)
             write(outunit, 3000) ctemp
        else
             write(outunit, 3000) a(rbp1,cbp1), b(rbp1,cbp1)
        endif
 3000   format( t5,d15.5)
c     end of output for debugging
      endif
      maxab1 = max(abs(a(rbp1, cbp1)), abs(b(rbp1, cbp1)))
      altb = .true.
      if (abs(a(rbp1, cbp1)) .ge. maxab1) altb = .false.
      if (ldebug) then
           write(outunit, 310)  'maxab1=', maxab1
  310      format(t5, a, d15.5)
           write(outunit,305) 'altb=', altb
  305      format(t5,a,l1)
      endif
      sa1 = a(rbp1, cbp1) / maxab1
      sb1 = b(rbp1, cbp1) / maxab1
      f = sa1 * b(rowb, colb) - sb1 * a(rowb, colb)
      g = sa1 * b(rowb, cbp1) - sb1 * a(rowb, cbp1)
c
c**** construct the right hand transformation (affects the columns
c     colb and colb + 1 of a, b and qq)
      call cgiv(f, g, c, s)
      call zcsrot(rbp1, a(1, colb), 1, a(1, cbp1),1, conjg(s), -c)
      call zcsrot(rbp1, b(1, colb), 1, b(1, cbp1),1, conjg(s), -c)
      call zcsrot(n, qq(1, colb), 1, qq(1, cbp1), 1, conjg(s), -c)
      if (ldebug) then
           call cmatpr( a,ldab,m,n, ' A after right transf.')
           call cmatpr( b,ldab,m,n, ' B after right transf.')
      endif
c
c**** construct the left hand transformation (affects the rows
c     rowb and rowb + 1 of a, b, and pp(conjg,trans))
      if (altb) then
         call cgiv(b(rowb, colb), b(rbp1, colb), c, s)
      else
         call cgiv(a(rowb, colb), a(rbp1, colb), c, s)
      endif
      call zcsrot(n-colb+1, a(rowb,colb), ldab, a(rbp1, colb),
     *           ldab, c, s)
      call zcsrot(n-colb+1, b(rowb,colb), ldab, b(rbp1, colb),
     *           ldab, c, s)
      call zcsrot(m, pp(1, rowb), 1, pp(1, rbp1), 1, c, conjg(s))
      if (ldebug) then
           call cmatpr( a,ldab,m,n, ' A after left transf.')
           call cmatpr( b,ldab,m,n, ' B after left transf.')
      endif
c
      a(rbp1, colb) = (0.d0, 0.d0)
      b(rbp1, colb) = (0.d0, 0.d0)
      if (ldebug) then
           write (outunit, 2000) 'eigenvalues after exchange'
           if (abs(b(rowb,colb)) .gt. 0.) then
              write(outunit, 3000) a(rowb,colb)/b(rowb,colb)
           else
              write(outunit, 3000) a(rowb,colb), b(rowb,colb)
           endif
           if (abs(b(rbp1,cbp1)) .gt. 0.) then
              write(outunit, 3000) a(rbp1,cbp1)/b(rbp1,cbp1)
           else
              write(outunit, 3000) a(rbp1,cbp1), b(rbp1,cbp1)
           endif
           call cmatpr( a,ldab,m,n, 'Final A after one exchange')
           call cmatpr( b,ldab,m,n, 'Final B after one exchange')
c     end of outputs for debugging
      endif
      return
      end

      subroutine cgiv( a, b, c, s)
c
c     implicit none
c***  debug space
      common /debug2/ idbg(20), outunit
      integer idbg, outunit
c***  formal parameter declarations
      complex*16 a, b, s , c
c
c**** cgiv constructs a complex givens transformation
c
c                c      s
c        g =                  c*c + s*conjg(s) = 1
c            -conjg(s)  c
c
c     which zeros the second entry of the 2-vector (a,b)**t:
c               a   aprim
c           g * b =   0   
c
c     cgiv leaves the arguments a and b unchanged,
c     (aprim is computed but no returned in this version).
c     note that the resulting c could have been chosen real
c     (but not for our application since we interchange c and s
c      when applying the the transformation in an equivalence
c      transformation)
c
c     if idbg(8) .eq. 0 then debug output is  withed off
c
c**** this version dated june, 1986
c
c**** internal variables
c
      real*8 sigma, delta, absa 
      complex*16  aprim, alfa
      logical ldebug
      ldebug = idbg(8) .ne. 0
c
      absa = abs(a)
      if ( absa .eq. 0) then
         c = 0.d0
         s = (1.d0, 0.d0)
         aprim = b
      else
         sigma = absa + abs(b)
         delta = sigma*sqrt(abs(a/sigma)**2 + abs(b/sigma)**2)
         alfa = a / absa
         c = absa /delta
         s = alfa * conjg(b) / delta
         aprim = alfa * delta
      endif
      if (ldebug) then
          write(outunit, 100) 'cos=', c, 'sin=', s
  100     format (t5, a, 2d12.5)
          write(outunit, 100) 'cos-sin-identity', 
     +                        c*conjg(c)+s*conjg(s)
      endif
      return
      end

      subroutine  zcsrot (n,cx,incx,cy,incy,c,s)
c
c     implicit none
      complex*16 cx(*), cy(*), c, s
      integer incx,incy,n
c
c**** zcsrot
c     applies a givens transformation where cos (c) and sin (s)
c     are complex as well as the vectors cx and cy.
c     the transformation is computed by cgiv.
c     note that c can be chosen real. however since we
c     will be able to interchange the values of c and s when
c     calling zcsrot we have to declare c as complex too.
c 
c     zcsrot is a modification of csrot
c     deal with complex sin (s) and cos (c)
c
c**** this version dated june, 1986
c
      integer i, ix, iy
      complex*16 ctemp
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = c*cx(ix) + s*cy(iy)
        cy(iy) = conjg(c)*cy(iy) - conjg(s)*cx(ix)
        cx(ix) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 continue
      do 30 i = 1,n
        ctemp = c*cx(i) + s*cy(i)
        cy(i) = conjg(c)*cy(i) - conjg(s)*cx(i)
        cx(i) = ctemp
   30 continue
      return
      end



      DOUBLE PRECISION FUNCTION DRANDU(DLIM)
C*************************************
C     REVISED: 920702                *
C*************************************
      IMPLICIT NONE
      DOUBLE PRECISION DLIM
C     DRANDU COMPUTES UNIFORM PSEUDO RANDOM NUMBERS IN THE
C     INTERVAL (-DLIM, DLIM)
C
      INTEGER*4 SEED(4)
      LOGICAL START
	  DOUBLE PRECISION DLARAN
	  EXTERNAL DLARAN
      DATA START/.TRUE./
      SAVE START,SEED
      IF (START) THEN
         SEED(1) = 1
         SEED(2) = 2 
         SEED(3) = 3
         SEED(4) = 4011
      ENDIF
      DRANDU = DLIM*DBLE(2.0D0*DLARAN(SEED)-1.0D0)
      START = .FALSE.
      RETURN
      END

      INTEGER FUNCTION IRANDU(ILIM)
C**************************************************
C     The routine computes a integer random number
C     in the interval (-LIM,LIM)
C     REVISED:   860731 13:00, 920702
C***************************************************
c     IMPLICIT NONE
      INTEGER ILIM
      DOUBLE PRECISION DRANDU
      IRANDU = NINT(DRANDU(DBLE(ILIM)))
      RETURN
      END


      DOUBLE PRECISION FUNCTION DLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 1.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
* *  This routine uses a multiplicative congruential method with modulus *  2**48 and multiplier 33952834046453 (see G.S.Fishman, *  'Multiplicative congruential random number generators with modulus *  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      DLARAN = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $         ( DBLE( IT4 ) ) ) ) )
      RETURN
*
*     End of DLARAN
*
      END
      subroutine codim(kstr,dimreg)
c     subroutine codim(a,ldab,b,kstr,colb,cole,rowb,rowe)
c     Erik Elmroth 921119
      
      implicit none
      
C     complex*16 a(ldab,*), b(ldab,*)
C     integer ldab,colb,cole,rowb,rowe,kstr(4,*)
      integer dimreg,kstr(4,*)
      
C     Computes the codimension of the matrix pencil a - lambda b,
C     where the matrix pencil on entry is on guptri form and the
C     kstr contains information abpout the Kronecker structure.
C     
C     This is a preliminary version where it is assumed that
C     the KCF for the regular nozero finite part only consists
C     simple eigenvalues.
      
      integer maxblk, i, p, k, Cjor, Cright, Cleft, Cjorsing
      integer Csing, cod
      integer regsiz, nosing, ii, jj, nn, Czero, Cinf
      parameter (maxblk = 20)
      integer L(0:maxblk), J(maxblk), R(maxblk), N(maxblk), LT(0:maxblk)
      
      L(0) = 0
      LT(0) = 0
      do 5 i = 1,maxblk
         L(i) = 0
         J(i) = 0
         R(i) = 0
         N(i) = 0
         LT(i) = 0
 5    continue
      
      
C     Compute L(i) = number of Li blocks (i = 0, 1, 2,... )
C     Compute J(i) = number of Ji blocks (i = 1, 2, 3,... )
      i = 1
      if (kstr(1,i) .ne. -1) then
         L(0) = kstr(1,i) - kstr(2,i)
 10      if (kstr(1,i+1) .ne. -1) then
            L(i) = kstr(1,i+1) - kstr(2,i+1)
            J(i) = kstr(2,i) - kstr(1,i+1)
            i = i + 1
            goto 10
         end if
         i = i + 1
      end if
      i = i + 1
      
C     Compute LT(i) = number of LTi blocks (i = 0, 1, 2,... )
C     Compute N(i) = number of Ni blocks (i = 1, 2, 3,... )
      ii = 1
      if (kstr(1,i) .ne. -1) then
         LT(0) = kstr(1,i) - kstr(2,i)
 20      if (kstr(1,i+1) .ne. -1) then
            LT(ii) = kstr(1,i+1) - kstr(2,i+1)
            N(ii) = kstr(2,i) - kstr(1,i+1)
            i = i + 1
            ii = ii + 1
            goto 20
         end if
         i = i + 1
      end if
      i = i + 1
      
C     Determine the number of singular blocks and
C     the size of the complete regular part
      regsiz = dimreg
      nosing = L(0) + LT(0)
      do 50 i = 1, maxblk
         nosing = nosing + L(i)
         regsiz = regsiz + J(i)*i
 50   continue
      do 51 i = 1, maxblk
         nosing = nosing + LT(i)
         regsiz = regsiz + N(i)*i
 51   continue

      Cjor = 0
      Czero = 0
      jj = 1
      do 70 i = maxblk, 1, -1
         if (J(i) .ne. 0) then
            do 80 k = 1, J(i)
               Czero = Czero + jj*i
               jj = jj + 2
 80         continue
         end if
 70   continue

      Cinf = 0
      nn = 1
      do 71 i = maxblk, 1, -1
         if (N(i) .ne. 0) then
            do 90 k = 1, N(i)
               Cinf = Cinf + nn*i
               nn = nn + 2
 90         continue
         end if
 71   continue
      
      Cright = 0
      do 30 k = 0, maxblk
         do 30 p = k+1, maxblk
            if (L(p) .ne. 0 .and. L(k) .ne. 0)
     $           Cright = Cright + (p - k - 1) * (L(p) * L(k))
 30      continue
         
         Cleft = 0
         do 40 k = 0, maxblk
            do 40 p = k+1, maxblk 
               if (LT(p) .ne. 0 .and. LT(k) .ne. 0)
     $              Cleft = Cleft + (p - k - 1) * (LT(p) * LT(k))
 40         continue
            
            Cjorsing = nosing*regsiz
            
            Csing = 0
            do 60 p = 0, maxblk
               do 60 k = 0, maxblk
                  if (L(p) .ne. 0 .and. LT(k) .ne. 0)
     $                 Csing = Csing + (p + k + 2) * (L(p) * LT(k))
 60            continue
               
               cod = Cjor + Czero + Cinf + Cright + Cleft + Cjorsing + 
     $              Csing
               
               return
               end
      
c   In this file June 7, 1987:Linpack routines - zsvdc, zqrdc, zqrsl
c
      subroutine zsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
      integer ldx,n,p,ldu,ldv,job,info
      complex*16 x(ldx,p),s(p),e(p),u(ldu,n),v(ldv,p),work(n)
c
c
c     zsvdc is a subroutine to reduce a complex*16 nxp matrix x by
c     unitary transformations u and v to diagonal form.  the
c     diagonal elements s(i) are the singular values of x.  the
c     columns of u are the corresponding left singular vectors,
c     and the columns of v the right singular vectors.
c
c     on entry
c
c         x         complex*16(ldx,p), where ldx.ge.n.
c                   x contains the matrix whose singular value
c                   decomposition is to be computed.  x is
c                   destroyed by zsvdc.
c
c         ldx       integer.
c                   ldx is the leading dimension of the array x.
c
c         n         integer.
c                   n is the number of columns of the matrix x.
c
c         p         integer.
c                   p is the number of rows of the matrix x.
c
c         ldu       integer.
c                   ldu is the leading dimension of the array u
c                   (see below).
c
c         ldv       integer.
c                   ldv is the leading dimension of the array v
c                   (see below).
c
c         work      complex*16(n).
c                   work is a scratch array.
c
c         job       integer.
c                   job controls the computation of the singular
c                   vectors.  it has the decimal expansion ab
c                   with the following meaning
c
c                        a.eq.0    do not compute the left singular
c                                  vectors.
c                        a.eq.1    return the n left singular vectors
c                                  in u.
c                        a.ge.2    returns the first min(n,p)
c                                  left singular vectors in u.
c                        b.eq.0    do not compute the right singular
c                                  vectors.
c                        b.eq.1    return the right singular vectors
c                                  in v.
c
c     on return
c
c         s         complex*16(mm), where mm=min(n+1,p).
c                   the first min(n,p) entries of s contain the
c                   singular values of x arranged in descending
c                   order of magnitude.
c
c         e         complex*16(p).
c                   e ordinarily contains zeros.  however see the
c                   discussion of info for exceptions.
c
c         u         complex*16(ldu,k), where ldu.ge.n.  if joba.eq.1
c                                   then k.eq.n, if joba.ge.2 then
c
c                                   k.eq.min(n,p).
c                   u contains the matrix of right singular vectors.
c                   u is not referenced if joba.eq.0.  if n.le.p
c                   or if joba.gt.2, then u may be identified with x
c                   in the subroutine call.
c
c         v         complex*16(ldv,p), where ldv.ge.p.
c                   v contains the matrix of right singular vectors.
c                   v is not referenced if jobb.eq.0.  if p.le.n,
c                   then v may be identified whth x in the
c                   subroutine call.
c
c         info      integer.
c                   the singular values (and their corresponding
c                   singular vectors) s(info+1),s(info+2),...,s(m)
c                   are correct (here m=min(n,p)).  thus if
c                   info.eq.0, all the singular values and their
c                   vectors are correct.  in any event, the matrix
c                   b = ctrans(u)*x*v is the bidiagonal matrix
c                   with the elements of s on its diagonal and the
c                   elements of e on its super-diagonal (ctrans(u)
c                   is the conjugate-transpose of u).  thus the
c                   singular values of x and b are the same.
c
c     linpack. this version dated 03/19/79 .
c              correction to shift calculation made 2/85.
c     g.w. stewart, university of maryland, argonne national lab.
c
c     zsvdc uses the following functions and subprograms.
c
c     external zdrot
c     blas zaxpy,zdotc,zscal,zswap,dznrm2,drotg
c     fortran dabs,dmax1,cdabs,dcmplx
c     fortran dconjg,max0,min0,mod,dsqrt
c
c     internal variables
c
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,
     *        mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      complex*16 zdotc,t,r
      double precision b,c,cs,el,emm1,f,g,dznrm2,scale,shift,sl,sm,sn,
     *                 smm1,t1,test,ztest
      logical wantu,wantv
c
      complex*16 csign,zdum,zdum1,zdum2
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
      csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2))
c
c     set the maximum number of iterations.
c
c**** 6/21/87
c      maxit = 30
       maxit = 100
c
c     determine what is to be computed.
c
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
c
c     reduce x to bidiagonal form, storing the diagonal elements
c     in s and the super-diagonal elements in e.
c
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
c
c           compute the transformation for the l-th column and
c           place the l-th diagonal in s(l).
c
            s(l) = dcmplx(dznrm2(n-l+1,x(l,l),1),0.0d0)
            if (cabs1(s(l)) .eq. 0.0d0) go to 10
               if (cabs1(x(l,l)) .ne. 0.0d0) s(l) = csign(s(l),x(l,l))
               call zscal(n-l+1,1.0d0/s(l),x(l,l),1)
               x(l,l) = (1.0d0,0.0d0) + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (cabs1(s(l)) .eq. 0.0d0) go to 30
c
c              apply the transformation.
c
               t = -zdotc(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
               call zaxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
c
c           place the l-th row of x into  e for the
c           subsequent calculation of the row transformation.
c
            e(j) = dconjg(x(l,j))
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
c
c           place the transformation in u for subsequent back
c           multiplication.
c
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
c
c           compute the l-th row transformation and place the
c           l-th super-diagonal in e(l).
c
            e(l) = dcmplx(dznrm2(p-l,e(lp1),1),0.0d0)
            if (cabs1(e(l)) .eq. 0.0d0) go to 80
               if (cabs1(e(lp1)) .ne. 0.0d0) e(l) = csign(e(l),e(lp1))
               call zscal(p-l,1.0d0/e(l),e(lp1),1)
               e(lp1) = (1.0d0,0.0d0) + e(lp1)
   80       continue
            e(l) = -dconjg(e(l))
            if (lp1 .gt. n .or. cabs1(e(l)) .eq. 0.0d0) go to 120
c
c              apply the transformation.
c
               do 90 i = lp1, n
                  work(i) = (0.0d0,0.0d0)
   90          continue
               do 100 j = lp1, p
                  call zaxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call zaxpy(n-l,dconjg(-e(j)/e(lp1)),work(lp1),1,
     *                       x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
c
c              place the transformation in v for subsequent
c              back multiplication.
c
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
c
c     set up the final bidiagonal matrix or order m.
c
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = (0.0d0,0.0d0)
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = (0.0d0,0.0d0)
c
c     if required, generate u.
c
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = (0.0d0,0.0d0)
  180       continue
            u(j,j) = (1.0d0,0.0d0)
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (cabs1(s(l)) .eq. 0.0d0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -zdotc(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call zaxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call zscal(n-l+1,(-1.0d0,0.0d0),u(l,l),1)
               u(l,l) = (1.0d0,0.0d0) + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = (0.0d0,0.0d0)
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = (0.0d0,0.0d0)
  260          continue
               u(l,l) = (1.0d0,0.0d0)
  270       continue
  280    continue
  290    continue
  300 continue
c
c     if it is required, generate v.
c
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (cabs1(e(l)) .eq. 0.0d0) go to 320
               do 310 j = lp1, p
                  t = -zdotc(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call zaxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = (0.0d0,0.0d0)
  330       continue
            v(l,l) = (1.0d0,0.0d0)
  340    continue
  350 continue
c
c     transform s and e so that they are double precision.
c
      do 380 i = 1, m
         if (cabs1(s(i)) .eq. 0.0d0) go to 360
            t = dcmplx(cdabs(s(i)),0.0d0)
            r = s(i)/t
            s(i) = t
            if (i .lt. m) e(i) = e(i)/r
            if (wantu) call zscal(n,r,u(1,i),1)
  360    continue
c     ...exit
         if (i .eq. m) go to 390
         if (cabs1(e(i)) .eq. 0.0d0) go to 370
            t = dcmplx(cdabs(e(i)),0.0d0)
            r = t/e(i)
            e(i) = t
            s(i+1) = s(i+1)*r
            if (wantv) call zscal(p,r,v(1,i+1),1)
  370    continue
  380 continue
  390 continue
c
c     main iteration loop for the singular values.
c
      mm = m
      iter = 0
c**** 6/23/87 added code to ensure convergence
c     compute norm of matrix
      test = abs(s(m))
      do 975 i=1,m-1
        test = test + abs(s(i)) + abs(e(i))
975   continue
      test = test * m * 100.
c****
  400 continue
c
c        quit if all the singular values have been found.
c
c     ...exit
         if (m .eq. 0) go to 660
c
c        if too many iterations have been performed, set
c        flag and return.
c
         if (iter .lt. maxit) go to 410
            info = m
c     ......exit
            go to 660
  410    continue
c
c        this section of the program inspects for
c        negligible elements in the s and e arrays.  on
c        completion the variables kase and l are set as follows.
c
c           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
c           kase = 2     if s(l) is negligible and l.lt.m
c           kase = 3     if e(l-1) is negligible, l.lt.m, and
c                        s(l), ..., s(m) are not negligible (qr step).
c           kase = 4     if e(m-1) is negligible (convergence).
c
         do 430 ll = 1, m
            l = m - ll
c        ...exit
            if (l .eq. 0) go to 440
c****       6/24/87, nonconvergence fix
c            test = cdabs(s(l)) + cdabs(s(l+1))
c****
            ztest = test + cdabs(e(l))
            if (ztest .ne. test) go to 420
               e(l) = (0.0d0,0.0d0)
c        ......exit
               go to 440
  420       continue
  430    continue
  440    continue
         if (l .ne. m - 1) go to 450
            kase = 4
         go to 520
  450    continue
            lp1 = l + 1
            mp1 = m + 1
            do 470 lls = lp1, mp1
               ls = m - lls + lp1
c           ...exit
               if (ls .eq. l) go to 480
c****          6/24/87, nonconvergence fix
c               test = 0.0d0
c               if (ls .ne. m) test = test + cdabs(e(ls))
c               if (ls .ne. l + 1) test = test + cdabs(e(ls-1))
c****
               ztest = test + cdabs(s(ls))
               if (ztest .ne. test) go to 460
                  s(ls) = (0.0d0,0.0d0)
c           ......exit
                  go to 480
  460          continue
  470       continue
  480       continue
            if (ls .ne. l) go to 490
               kase = 3
            go to 510
  490       continue
            if (ls .ne. m) go to 500
               kase = 1
            go to 510
  500       continue
               kase = 2
               l = ls
  510       continue
  520    continue
         l = l + 1
c
c        perform the task indicated by kase.
c
         go to (530, 560, 580, 610), kase
c
c        deflate negligible s(m).
c
  530    continue
            mm1 = m - 1
            f = dreal(e(m-1))
            e(m-1) = (0.0d0,0.0d0)
            do 550 kk = l, mm1
               k = mm1 - kk + l
               t1 = dreal(s(k))
               call drotg(t1,f,cs,sn)
               s(k) = dcmplx(t1,0.0d0)
               if (k .eq. l) go to 540
                  f = -sn*dreal(e(k-1))
                  e(k-1) = cs*e(k-1)
  540          continue
               if (wantv) call zdrot(p,v(1,k),1,v(1,m),1,cs,sn)
  550       continue
         go to 650
c
c        split at negligible s(l).
c
  560    continue
            f = dreal(e(l-1))
            e(l-1) = (0.0d0,0.0d0)
            do 570 k = l, m
               t1 = dreal(s(k))
               call drotg(t1,f,cs,sn)
               s(k) = dcmplx(t1,0.0d0)
               f = -sn*dreal(e(k))
               e(k) = cs*e(k)
               if (wantu) call zdrot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  570       continue
         go to 650
c
c        perform one qr step.
c
  580    continue
c
c           calculate the shift.
c
            scale = dmax1(cdabs(s(m)),cdabs(s(m-1)),cdabs(e(m-1)),
     *                    cdabs(s(l)),cdabs(e(l)))
            sm = dreal(s(m))/scale
            smm1 = dreal(s(m-1))/scale
            emm1 = dreal(e(m-1))/scale
            sl = dreal(s(l))/scale
            el = dreal(e(l))/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
            c = (sm*emm1)**2
            shift = 0.0d0
            if (b .eq. 0.0d0 .and. c .eq. 0.0d0) go to 590
               shift = dsqrt(b**2+c)
               if (b .lt. 0.0d0) shift = -shift
               shift = c/(b + shift)
  590       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
c
c           chase zeros.
c
            mm1 = m - 1
            do 600 k = l, mm1
               call drotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = dcmplx(f,0.0d0)
               f = cs*dreal(s(k)) + sn*dreal(e(k))
               e(k) = cs*e(k) - sn*s(k)
               g = sn*dreal(s(k+1))
               s(k+1) = cs*s(k+1)
               if (wantv) call zdrot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call drotg(f,g,cs,sn)
               s(k) = dcmplx(f,0.0d0)
               f = cs*dreal(e(k)) + sn*dreal(s(k+1))
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*dreal(e(k+1))
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n)
     *            call zdrot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  600       continue
            e(m-1) = dcmplx(f,0.0d0)
            iter = iter + 1
         go to 650
c
c        convergence.
c
  610    continue
c
c           make the singular value  positive
c
            if (dreal(s(l)) .ge. 0.0d0) go to 620
               s(l) = -s(l)
               if (wantv) call zscal(p,(-1.0d0,0.0d0),v(1,l),1)
  620       continue
c
c           order the singular value.
c
  630       if (l .eq. mm) go to 640
c           ...exit
               if (dreal(s(l)) .ge. dreal(s(l+1))) go to 640
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p)
     *            call zswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n)
     *            call zswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 630
  640       continue
            iter = 0
            m = m - 1
  650    continue
      go to 400
  660 continue
      return
      end
      subroutine rzstr (opt, a, b, ldab, m, n, rowb, rowe,
     *                  colb, cole, first, zero, epsua, epsub, gap,
     *                  pp, ldpp, qq, ldqq, kstr, kfirst, step,
     *                  adlsvd, bdlsvd,
     *                  work, x, sx, ex, q, arow, brow, w, qraux, y,
     *                  qty, info)
c
c     implicit none
c**** debug space
c     the common-block declarations assume that the dimension of the
c     input matrix pencil a - lambda b is not larger than abdim.
c     the debug space is used for producing debug outputs (optional,
c     see below)
c
      integer abdim
      parameter (abdim = 30)
      common /debug1/ acopy(abdim,abdim),bcopy(abdim,abdim),
     *              atest(abdim,abdim), btest(abdim,abdim), swap
      common /debug2/ idbg(20), outunit
      complex*16 acopy,bcopy,atest,btest
      logical swap
      integer idbg, outunit
c
c**** formal parameter declarations
      character*(*) opt
      integer ldab, m, n, rowb, rowe, colb, cole, ldpp, ldqq,
     *        kstr(4,*), step, kfirst, info
      logical first, zero
      real*8 adlsvd, bdlsvd, epsua, epsub, gap
      complex*16 a(ldab,*), b(ldab,*), pp(ldpp,*), qq(ldqq,*),
     *        work(*)
c
c****    workspace 
c
        complex*16   x(m,*), sx(*), ex(*), q(n,*), 
     *           arow(*), brow(*), w(n,*), qraux(*),
     *           y(*), qty(*)
c
c*******************************************************************
c
c     rzstr computes the kronecker right (column) structure and
c     the jordan structure of the zero-eigenvalue of a singular
c     pencil a-lambda*b.for details concerning the listr-kernel see 
c     the following papers:
c     
c        b.kagstrom, rgsvd - an algorithm for computing the kronecker
c             structure and reducing subspaces of singular a - lambda b
c             pencils, siam j.sci.stat.comput., vol. 7, 1986, pp 185-211
c
c        j.demmel and b.kagstrom, stably computing the kronecker 
c             structure and reducing subspaces of singular pencils
c             a - lambda b for uncertain data, in large scale eigenvalue
c             problems (cullum, willoughby eds), north holland, 1986,
c             pp 283-323.
c
c
c     formal parameters
c
c     on entry
c
c        opt*(*) character, if opt = 'cind' rzstr computes indices
c                           if opt = 'rind' already computed indices
c                           are reused in the reduction
c
c        a(ldab,*) complex*16, input matrix a of order m by n
c
c        b(ldab,*) complex*16, input matrix b of order m by n
c
c        ldab      integer, leading dimension of a and b
c
c        m         integer, current row dimension of a and b
c
c        n         integer, current column dimension of and b
c
c        rowb      integer, first row of the subpencil
c
c        rowe      integer, last row of the subpencil
c
c        colb      integer, first column of the subpencil
c
c        cole      integer, last column of the subpencil
c
c        first     logical, first should be 'true' if first call to 
c                  rzstr, else 'false'
c
c        zero      logical, if 'true', zero out small singular values
c                  so returned pencil really has structure described
c                  in kstr (see below), else returned pencil is a
c                  true equivalence transformation of input pencil
c                  (no singular values are deleted)
c
c        epsua     real*8, threshold for deleting singular values of a
c                  (used when compressing columns of a)
c
c        epsub     real*8, threshold for deleting singular values of b
c                  (used when compressing columns of b)
c
c        gap       real*8, should be at least 1 and nominally 1000.
c                  used by subroutine rcsvdc to make rank decisions
c                  by searching for adjacent singular values whose
c                  ratio exceeds gap.
c
c        ldpp      integer, leading dimension of pp
c
c        ldqq      integer, leading dimension of qq
c
c        kfirst    integer, index to the first location in kstr
c                  where structure-index information is stored
c                  from this reduction (see below)
c
c     on exit
c
c        pp(ldpp,*)complex*16, left unitary transformation matrix 
c                  pp of order m by m
c
c        qq(ldqq,*)complex*16, right unitary transformation matrix
c                  qq of order n by n
c
c        a(ldab,*) transformed matrix a (pp**h * a * pp)
c
c        b(ldab,*) transformed matrix b (pp**h * b * pp)
c
c        kstr(4,*) integer, stores information concerning right 
c                  kronecker indices and the jordan structure of
c                  the zero eigenvalue.
c                  kstr(1,kfirst-1+j) - kstr(2,kfirst-1+j) =
c                  number of l(j-1) blocks (right indices of
c                  degree j-1).
c                  kstr(2,kfirst-1+j) - kstr(1,kfirst+j) = 
c                  number of jordan blocks of the zero
c                  eigenvalue of dimension j.
c                  index j goes from 1 to step (see below)
c                  note: rows 3 and 4 of kstr are not used inside 
c                  rzstr.
c
c        step      integer,  the number of deflation-steps in this
c                  reduction
c
c        adlsvd    real*8, root sum of squares of deleted singular
c                  values of a (independent of the input zero)
c
c        bdlsvd    real*8, root sum of squares of deleted singular
c                  values of b (independent of the input zero)
c
c        info      integer, zero if normal return,
c                           1 if svd does not converge
c
c        on exit from rzstr a and b will be in block upper triangular form:
c
c
c              a = ( arz   *  )        b = ( brz    *  )
c                  (  0   a22 )            (  0    b22 )
c
c        the block structure of arz - lambda*brz describes the 
c        kronecker column (right) structure and the jordan structure 
c        of the zero eigenvalue. if ni and ri denote the dimension of
c        the diagonal blocks in arz and brz (see example below), then 
c        they have the following interpretation:
c
c          ni - ri = the number of l(i-1) -blocks of order (i-1) by i
c          ri - ni+1 = the number of j(0)-blocks of order i by i
c
c        note that if a - lambda*b is a regular pencil then ni=ri.
c        the rzstr reduction stops when an ni.eq.0 or ni.ne.0 but ri.eq.0. 
c        then a22 will have full column rank. a22 - lambda*b22 might
c        still be a singular pencil (can have row (left) indices). 
c        an example illustrates the two cases (see papers for details):
c        case 1 - n4.eq.0:
c
c                ( 0  a12 a13 ) r1           ( b11 b12 b13 ) r1
c          arz = ( 0   0  a23 ) r2     brz = (  0  b22 b23 ) r2
c                ( 0   0   0  ) r3           (  0   0  b33 ) r3
c                  n1  n2  n3                   n1  n2  n3
c
c        case 2 - n4.ne.0 and r4.eq.0:
c
c                ( 0  a12 a13 a14 ) r1        ( b11 b12 b13 b14 ) r1
c          arz = ( 0   0  a23 a24 ) r2  brz = (  0  b22 b23 b24 ) r2
c                ( 0   0   0  a34 ) r3        (  0   0  b33 b34 ) r3
c                  n1  n2  n3  n4                n1  n2  n3  n4
c
c       the ri by ni diagonal blocks bii of brz are in the form
c       ( 0 rii), where rii is ri by ri, nonsingular and upper
c       triangular.
c
c       if kfirst = 1 on input then case 2 above cause the following
c       output for step and kstr:
c         step = 4
c         kstr(1,1) = n1   kstr(2,1) = r1
c         kstr(1,2) = n2   kstr(2,2) = r2
c         kstr(1,3) = n3   kstr(2,3) = r3
c         kstr(1,4) = n4   kstr(2,4) = 0
c
c       note that on output (arz,brz) or (a22,b22) can be nonexistent
c       in the block upper triangular form (a,b). (arz,brz) does not 
c       exist if n1=r1=0. (a22,b22) does not exist if the input pencil
c       a -lambda*b has no left (row) singular structure, no
c       infinite eigenvalue and no nonzero eigenvalues.
c
c***     work space including size (all variables complex*16)
c        work(*)           max(m,n) (if linpack svd is used)
c                          or 2*min(m,n)+6*max(m,n) (lapack svd, 920708)
c        x(m,*)            m by n
c        sx(*)             min(m,n) + 1
c        ex(*)             n
c        q(n,*)            n by n
c        arow(*)           max(m,n)
c        brow(*)           max(m,n)
c        w(n,*)            n by n
c        qraux(*)          max(m,n)
c        y(*)              max(m,n)
c        qty(*)            max(m,n)
c
c*****************************************************************
c
c****    this version dated june 16, 1987
c        authors: jim demmel and bo kagstrom
c
c****    rzstr uses the following functions and subroutines
c        kcfpack  - cmatml, cmatmr, cmatpr, cmcopy, rcsvdc, upddel 
c        linpack  - zqrdc, zqrsl
c
c****    internal variables
c
        logical ldebug
        integer mrow, ncol, i, j, sn1, sr1, rep, rowsr1, colsn1, xrow
     *          , xcol, job, ldx, ldq, n1, rnull, ldw, cnull, r1,
     *          colsnb, jend, idummy, ikstr, mxrc, k, iii, jjj
c
        real*8 del, difa, difb
c
        complex*16 dummy

c
c****   set leading dimensions of x, q, and w
c
        ldx = m
        ldq = n
        ldw = n
c       set debug switch
        ldebug= (idbg(4).ne.0)
c****   compute the order of the pencil in action (mrow * ncol)
c
        mrow = rowe - rowb + 1
        ncol = cole - colb + 1
c
c*+*+*  accumulate deleted singular values in adlsvd, bdlsvd 
        adlsvd = 0.0
        bdlsvd = 0.0
c
      if (ldebug) write (outunit,1001) 'epsua=', epsua
      if (ldebug) write (outunit,1001) 'epsub=', epsub
1001  format(t5,a,d13.6)
c
c
c**** set rep depending on what option
c
      if ( opt .eq. 'cind' ) then
c         perhaps not enough !!
          rep = rowe * cole
      else
          rep = step - kfirst + 1
      endif
c***  6/18/87
      if (ldebug) write(outunit,2000) 'kfirst=',kfirst,
     +            'step=',step,'rep=',rep
c
      sn1 = 0
      sr1 = 0
      step = 0
c**** while rep > 0 do
   30 continue
      if (ldebug) write(outunit,2000) 'rep at top of loop=',rep
      if (rep .eq. 0) go to 500
c     jump when while - loop satisfied
c
c     while - clause
        step = step + 1
        if (ldebug) write(outunit,2000) 'Results from step = ', step
 2000   format( t5, a, i3/)
        if (ldebug) write(outunit,2005) opt
 2005   format(t5,a)
c
c**** set n1 and r1 if we are reusing kstr
c
      if ( opt .eq. 'rind' ) then
         ikstr = kfirst + step - 1
         n1 = kstr(1, ikstr)
         r1 = kstr(2, ikstr)
         cnull = n1 -r1
      endif
c
c**** step 1 - compress columns of a (gives n1 = dimension of the
c              column nullspace)
c* 1.1
c      rows, rowb+sr1:rowe
c      cols, colb+sn1:cole
c
        rowsr1 = rowb + sr1 - 1
        colsn1 = colb + sn1 - 1
        xrow = mrow - sr1
        xcol = ncol - sn1
        do 40 i = 1, xrow
           do 35 j = 1, xcol
              x(i, j) = a(rowsr1 + i, colsn1 + j)
   35      continue
   40   continue
        if ( xrow .ge. xcol ) then
           job = 0110
        else
           job = 0111
        endif
        if (ldebug) then
          write(outunit,5000) 'rowsr1=',rowsr1,'colsn1=',colsn1,
     +                        'xrow=',xrow
          write(outunit,5000) 'xcol=',xcol,'rowb=',rowb,'rowe=',rowe
          write(outunit,5000) 'colb=',colb,'cole=',cole,'sr1=',sr1,
     +                        'sn1=',sn1
        endif
c
c       put m*n in info before calling
        if (idbg(4) .gt. 2) then
          call cmatpr(x ,ldx, xrow, xcol,'a-input rcsvdc')
        endif
        info = m*n

        call rcsvdc (x, ldx, xrow, xcol, sx, ex, dummy, 1, q, ldq, opt,
     *               epsua, gap, n1, rnull, del, work, job, info )
c
c
         call upddel(adlsvd, del)
c
        mxrc = min0( xrow, xcol)
        if (ldebug) call cmatpr( sx, 1, 1, mxrc,
     *               'singular values - column compress a')
       if (idbg(4) .gt. 1 .or. (info .ne. 0 .and. ldebug) ) then
         call cmatpr( ex, 1, 1, mxrc,'sub diagonal - should be zero')
         call cmatpr(q, ldq, xcol, xcol,
     *               'step1.1: right singular vectors of A')
       endif
        if (ldebug) write (outunit,1005) 'info=', info, 'n1=', n1
 1005   format(t5, a, i3/ )
c
         if (info .ne. 0) then
c****      6/18/87
           if (ldebug) write(outunit,2007) info
 2007      format('rzstr - after first call to rcsvdc, info= ',i4)
           info = 1
           return
         endif
c
c       if n1=0, we are done
        if (n1 .eq. 0) then
           r1=0
           goto 450
        end if
c
c* 1.2 - apply right transformation q to a and b (the full matrices)
c        rows in a and b: 1:rowe
c        columns in a: colb+sn1:cole ( xcol col's)
c        columns in b: colb+sn1:cole
c
        do 70 i = 1, rowe
           do 50 j = 1, xcol
              arow(j) = 0.d0
              brow(j) = 0.d0
              do 45 k = 1, xcol
                 arow(j) = arow(j) + a(i, colsn1 + k) * q(k, j)
                 brow(j) = brow(j) + b(i, colsn1 + k) * q(k, j)
   45         continue
   50      continue
           do 60 j = 1, xcol
                 a(i, colsn1 + j) = arow(j)
                 b(i, colsn1 + j) = brow(j)
   60      continue
   70   continue
c
c*         zero part of a
c          rows, rowb+sr1:rowe
c          cols, colb+sn1:colb+sn1+n1-1
c
        if (zero) then
          do 80 i = rowb + sr1, rowe
             do 75 j = colb + sn1, colsn1 + n1
                a(i, j) = 0.d0
   75        continue
   80     continue
        endif
c
c**** Step 2 - column compress part of B ( gives n1 - r1 =
c              dimension of the common nullspace)
c
c* 2.1
c       rows, rowb+sr1:rowe
c       cols, colb+sn1:colb+sn1+n1-1
c
        xrow = mrow - sr1
        xcol = n1
        do 90 i = 1, xrow
           do 85 j = 1, xcol
              x(i, j) = b( rowsr1 + i, colsn1 + j)
   85      continue
   90   continue

        if (xrow .ge. xcol) then
           job = 0110
        else
           job = 0111
        endif
        if (idbg(4) .gt. 2) then
          call cmatpr(x ,ldx, xrow, xcol,'b-input rcsvdc')
        endif
        info = m*n
        call rcsvdc ( x, ldx, xrow, xcol, sx, ex, dummy, 1, w, ldw,
     *           opt, epsub, gap, cnull, rnull, del, work, job, info )
c
        if ( opt .eq. 'cind' ) r1 = n1 - cnull
c
c       if r1 = 0 then we are done ! Zero part in b and then update qq
c
c
        if (ldebug) write(outunit,1005) 'info=', info, 'cnull=', cnull,
     *                  'n1=', n1,'r1=', r1
c
        mxrc = min0( xrow, xcol)
        if (ldebug) call cmatpr( sx, 1, 1, mxrc,
     *               'singular values - column compress b')
       if (idbg(4) .gt. 1 .or. (info .ne. 0 .and. ldebug) ) then
         call cmatpr( ex, 1, 1, mxrc,'sub diagonal - should be zero')
         call cmatpr ( w, ldw, xcol, xcol,
     *                'step 2.1: right singular vectors of b')
       endif
c     
       call upddel(bdlsvd, del)
c
        if (info .ne. 0) then
c****     6/18/87
          if (ldebug) write(outunit,2008) info
 2008     format('rzstr - after second call to rcsvdc, info= ',i4)
          info = 1
          return
        endif
c
       if (r1 .eq. 0) goto 3500
c
c* 2.2
c      update q  rows, 1:ncol-sn1
c                cols, 1:n1
c         a, b   rows, 1:rowe
c                cols, colb+sn1:colb+sn1+n1-1
c
c      note that we do not make use of that some of the elements in a
c      are zero
c      first q
       xcol = ncol - sn1
       do 110 i = 1, xcol
          do 100 j = 1, n1
             arow(j) = 0.d0
             do 95 k = 1, n1
                arow(j) = arow(j) + q(i, k) * w(k, j)
   95        continue
  100     continue
c
          do 105 j = 1, n1
             q(i, j) = arow(j)
  105     continue
  110   continue
c
        if (idbg(4) .gt. 2) then
           call cmatpr(q, ldq, xcol, xcol,
     *                'updated q after second column compress')
        endif
c        
c       now a and b ....
        do 120 i = 1, rowe
           do 114 j = 1, n1
              arow(j) = 0.d0
              brow(j) = 0.d0
              do 112 k = 1, n1
                 arow(j) = arow(j) + a(i, colsn1 + k) * w(k, j)
                 brow(j) = brow(j) + b(i, colsn1 + k) * w(k, j)
  112         continue
  114      continue
           do 116 j = 1, n1
              a(i, colsn1 + j) = arow(j)
              b(i, colsn1 + j) = brow(j)
  116      continue
  120   continue
c
c*        zero part of b
c         rows, rowb+sr1:rowe
c         cols, colb+sn1:colb+sn1+(n1-r1)-1
c
 3500  continue
       if (zero) then
         do 130 i = rowb + sr1, rowe
            do 125 j = 1, n1 - r1
               b(i, colsn1 + j) = 0.d0
  125       continue
  130    continue
       endif
c
       if (r1 .eq. 0 )go to 350 
c
c**** Step 3 - Triangularize b ( using qr)
c
c* 3.1
c         rows, rowb+sr1:rowe
c         cols, colb+sn1+(n1-r1):cole
c
          xrow = mrow - sr1
          xcol = ncol - sn1 - (n1 - r1)
          colsnb = colsn1 + (n1-r1)
          do 140 i = 1, xrow
             do 135 j = 1, xcol
                x(i, j) = b( rowsr1 + i, colsnb + j)
  135        continue
  140     continue
          job = 0
          call zqrdc( x, ldx, xrow, xcol, qraux, idummy, dummy, job)
c
c         move the upper triangular part to b
c
          do 150 i = 1, xrow
             do 145 j = i, xcol
                b(rowsr1 + i, colsnb + j) = x(i, j)
  145        continue
             jend = min0(xcol, i - 1)
             do 148 j = 1, jend
                b(rowsr1 + i, colsnb + j) = 0.d0
  148        continue
  150    continue
c
c* 3.2
c        apply v(conj,trans) to remaining cols of b
c        from the left (xrow*xrow)
c        rows, rowb+sr1:rowe
c        cols, cole+1:n
c
         do 170 j = cole+1, n
            do 160 i = 1, xrow
               y(i) = b(rowsr1 + i, j)
  160       continue
            job = 01000
            call zqrsl(x, ldx, xrow, xcol, qraux, y, dummy, qty,
     *                 dummy, dummy, dummy, job, info)
            do 165 i = 1, xrow
               b(rowsr1 + i, j) = qty(i)
  165       continue
  170    continue
c        if (ldebug) call cmatpr(b, ldab, m, n,
c    *              'B after triangularization - step 3.1')
c
c        apply v(conj,trans) to a from the left (xrow*xrow)
c        rows, rowb+sr1:rowe
c        cols, colb+sn1+n1:n
c
         do 185 j = colb + sn1 + n1, n
            do 180 i = 1, xrow
               y(i) = a(rowsr1 + i, j)
  180       continue
            job = 01000
            call zqrsl(x, ldx, xrow, xcol, qraux, y, dummy, qty,
     *                 dummy, dummy, dummy, job, info)
            do 175 i = 1, xrow
               a(rowsr1 + i, j) = qty(i)
  175       continue
  185    continue
c        if (ldebug) call cmatpr(a,ldab,m,n,'A after step 3.2')
c
c****    update left transformation matrix pp ( m*m )
c        rows, 1:m
c        cols, rowb+sr1:rowe
c
         do 200 i = 1, m
            do 190 j = 1, xrow
               y(j) = conjg( pp(i, rowsr1+j) )
  190       continue
            job = 01000
            call zqrsl( x, ldx, xrow, xcol, qraux, y, dummy, qty,
     *                 dummy, dummy, dummy, job, info)
            do 195 j = 1, xrow
               pp(i, rowsr1 + j) = conjg( qty(j) )
  195       continue
  200    continue
c        
         if (idbg(4) .gt. 1) then
            call cmatpr(pp, ldpp, m, m,
     *                  'step 3.2: pp after updating with w from qr')
         endif
  350    continue
c
c****    update right transformation matrix qq (n*n)
c        rows, 1:n
c        cols, colb+sn1:cole
c
         xcol = ncol - sn1
         if (first) then
            do 210 i = 1, n
               do 205 j = 1, n
                  qq(i, j) = q(i, j)
  205          continue
  210       continue
        else
            do 240 i = 1, n
               do 230 j = 1, xcol
                  arow(j) = 0.d0
                  do 220 k = 1, xcol
                     arow(j) = arow(j) + qq(i, colsn1 + k) * q(k, j)
  220             continue
  230          continue
               do 235 j = 1, xcol
                  qq(i, colsn1 + j) = arow(j)
  235          continue
  240       continue
         endif
         if (idbg(4) .gt. 2) then
            call cmatpr(qq,ldqq, m, m, 'updated qq')
         endif
c        
c****    update indices
c
         sn1 = sn1 + n1
         sr1 = sr1 + r1
      if (ldebug) then
        write(outunit,5000) 'rowsr1=',rowsr1,'colsn1=',colsn1,
     +                      'xrow=',xrow
        write(outunit,5000) 'xcol=',xcol,'rowb=',rowb,'rowe=',rowe
        write(outunit,5000) 'colb=',colb,'cole=',cole,'sr1=',sr1,
     +                      'sn1=',sn1
      endif
c*       monitoring of the r1 and n1 in kstr
c
  450    continue
c****    added 060787 to match zlistr
       if (ldebug) then
         if (swap) then
           call cmcopy(bcopy,20,m,n,atest)
           call cmcopy(acopy,20,m,n,btest)
         else
           call cmcopy(acopy,20,m,n,atest)
           call cmcopy(bcopy,20,m,n,btest)
         end if
         call cmatml(atest,20,m,n,pp,ldpp,m,atest,20,work,3)
         call cmatmr(atest,20,m,n,qq,ldqq,n,atest,20,work,1)
         call cmatml(btest,20,m,n,pp,ldpp,m,btest,20,work,3)
         call cmatmr(btest,20,m,n,qq,ldqq,n,btest,20,work,1)
         difa=0
         difb=0
         do 1234 iii=1,m
           do 5678 jjj=1,n
             difa=difa+abs(atest(iii,jjj)-a(iii,jjj))
             difb=difb+abs(btest(iii,jjj)-b(iii,jjj))
 5678      continue
 1234    continue
         write(outunit,201) 'difa=',difa
 201     format(t5,a,d13.6/)
c        call cmatpr(atest,20,m,n,'atest')
         write(outunit,201) 'difb=',difb
c        call cmatpr(btest,20,m,n,'btest')
       endif
c
c****    compute rep depending on what option is used
c
         if ( opt .eq. 'cind') then
            kstr(1, step) = n1
            kstr(2, step) = r1
            rep = n1 * r1 * (mrow - sr1) * (ncol - sn1)
         else
            rep =rep - 1
         endif
         if (ldebug) write(outunit,5000) 'sn1=',sn1,'sr1=',sr1,
     +                                   'rep=',rep
 5000    format(t5,a,i4/)
         first = .false.
         go to 30
c
c**** end of while clause
  500 continue
c
      return
      end
c     last line of zrstr


