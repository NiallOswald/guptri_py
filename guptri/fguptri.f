C  This file is derived from historical GUPTRI source code.
C  Original authorship is attributed to James Demmel and Bo Kagstrom.
C  This copy has been modified for inclusion in this package.

C
C     Changed 1995-02-14 torkel
C     line 17/18

      subroutine fguptri(a,b,m,n,epsu,gap,zero, work, lwork,
     $                     pp, qq, kstr, lkstr, info)
      implicit none
      integer m, n, info, lwork, lkstr 
      integer kstr(4,lkstr)
      complex*16 a(m,n), b(m,n), pp(m,m), qq(n,n), work(lwork)
      double precision  epsu, gap, adelta, bdelta
      
      integer ldab, ldpp, ldqq, pstruc(4)
      integer rtre, rtce, zrre, zrce, fnre, fnce, inre, ince
      logical zero

!f2py intent(inplace) :: a, b
!f2py intent(in) :: epsu, gap, zero
!f2py intent(out) :: pp, qq, info, kstr
!f2py intent(hide) :: work
!f2py integer intent(hide), depend(a) :: m, n
!f2py integer intent(hide), depend(m,n) :: lkstr = max(m,n)+6
!f2py integer intent(hide), depend(m,n) :: lwork = 2*(max(m,n)*max(m,n))
!f2py+ +m*n+(min(m,n)*min(m,n))+12*max(m,n)+3*min(m,n)+1

       integer idbg(20), outunit
       common /debug2/ idbg, outunit
       data idbg/20*0/
      
      ldab = m
      ldpp = m
      ldqq = n
      
      call guptri(a, b, ldab, m, n,  epsu, gap, zero,
     $     pp, ldpp, qq, ldqq, adelta, bdelta, rtre, rtce,
     $     zrre, zrce, fnre, fnce, inre, ince, pstruc,
     $     work, work(max(m,n)+1), kstr, info)

      return
      end
      

      subroutine convertkstr(ikstr, dkstr, work, kstrcols)

      implicit none
      integer ikstr(4,*), work(4,*), kstrcols, i, j
      double precision dkstr(4,*)
C
C     Copies the integer array ikstr to the double precision
C     array dkstr. The operation is performed in two steps
C     whith first a copy to work, then another copy to dkstr,
C     which enables call where ikstr and dkstr refers to the
C     the same memory locations.
C

      do i = 1, 4
         do j = 1, kstrcols
            work(i,j) = ikstr(i,j)
         end do
      end do

      do i = 1, 4
         do j = 1, kstrcols
            dkstr(i,j) = dble(work(i,j))
         end do
      end do

      end
