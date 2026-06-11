C
C     Changed 1995-02-14 torkel
C     line 17/18

      subroutine fguptri(a,b,m,n,epsu,gap,zero, work, lwork,
     $                     pp, qq, kstr, info)
      implicit none
      integer m, n, info, lwork, lkstr 
      integer kstr(4,*)
      complex*16 a(*), b(*), pp(*), qq(*), work(*)
      double precision  epsu, gap, adelta, bdelta
      
      integer ldab, ldpp, ldqq, pstruc(4)
      integer rtre, rtce, zrre, zrce, fnre, fnce, inre, ince
      logical zero

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
