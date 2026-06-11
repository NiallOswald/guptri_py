c    On this file June 13, 1987:
c    guptri, upddel, cident, krnstr, norme  
c
        subroutine guptri(a, b, ldab, m, n, epsu, gap, zero,
     *                    pp, ldpp, qq, ldqq,
     *                    adelta, bdelta, rtre, rtce, zrre, zrce,
     *                    fnre, fnce, inre, ince, pstruc, struc,
     *                    work, kstr, info)
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
     *                  atest(abdim,abdim),btest(abdim,abdim),swap
        common /debug2/ idbg(20),outunit
        complex*16 acopy,bcopy,atest,btest
        logical swap
        integer idbg, outunit
c
c**** formal parameter declarations
        integer ldab, ldpp, ldqq, m, n
        complex*16 a(ldab,*), b(ldab,*), pp(ldpp,*), qq(ldqq,*)
        real*8 epsu, gap, adelta, bdelta
        integer rtre, rtce, zrre, zrce, fnre, fnce, inre, ince
        integer pstruc(4), struc(*), info
        logical zero
c
c****   work space
        integer kstr(4,*)
        complex*16 work(*)
c
c***********************************************************************
c
c       guptri reduces the pencil a - lambda b to generalized upper 
c       triangular (guptri) form via unitary equivalence transformations.
c       the guptri reduction is based on an improved version of the
c       rgqzd algortihm (a unitary version of the rgsvd algorithm).
c       for details see the papers:
c        b.kagstrom, rgsvd - an algorithm for computing the kronecker
c             structure and reducing subspaces of singular a - lambda b
c             pencils, siam j.sci.stat.comput., vol. 7, 1986, pp 185-211
c
c        j.demmel and b.kagstrom, stably computing the kronecker 
c             structure and reducing subspaces of singular pencils
c             a - lambda b for uncertain data, in large scale eigenvalue
c             problems (cullum, willoughby eds), north holland, 1986
c             pp 283-323.
c
c      debug switch for guptri is idbg(2)
c       - if idbg(2) ne 0, print debug output, else no output
c
c       on entry
c
c        a(ldab,*)      complex*16, input matrix a of order m by n
c
c        b(ldab,*)      complex*16, input matrix b of order m by n
c
c        ldab           integer, leading dimension of a and b
c
c        m              integer, current row dimension of a and b
c
c        n              integer, current column dimension of a and b
c
c        epsu           real*8, relative uncertainty in data 
c                       (should be at least about macheps). used by
c                       subroutine rcsvdc to make rank decisions
c
c        gap            real*8, should be at least 1 and nominally 
c                       1000. used by subroutine rcsvdc to make rank
c                       decisions by searching for adjacent singular
c                       values whose ratio exceeds gap
c
c        zero           logical, if true, zero out small singular values
c                       so returned pencil really has structure described
c                       in pstruc and struc (see below), else returned
c                       pencil is a true equivalence transformation of
c                       input pencil (no singular values are deleted)
c
c        ldpp           integer, leading dimension of pp
c
c        ldqq           integer, leading dimension of qq
c
c       on exit
c
c        pp(ldpp,*)     complex*16, left unitary transformation matrix pp
c                       of order m by m such that
c                       pp**h * (a - lambda b) * qq is in guptri form
c                       (described below)
c
c        qq(ldqq,*)     complex*16, right unitary transformation matrix qq
c                       of order m by m such that
c                       pp**h * (a - lambda b) * qq is in guptri form
c
c        a(ldab,*)      transformed matrix a (pp**H * a * qq) in 
c                       guptri form
c
c        b(ldab,*)      transformed matrix b (pp**H * b * qq) in
c                       guptri form
c
c        guptri (generalized upper triangular) form is described as 
c           follows: on output
c
c               ( art  *   *   *   *  )      ( brt  *   *   *   *  )
c               (  0  azr  *   *   *  )      (  0  bzr  *   *   *  )
c           a = (  0   0  afn  *   *  ), b = (  0   0  bfn  *   *  )
c               (  0   0   0  ain  *  )      (  0   0   0  bin  *  )
c               (  0   0   0   0  alt )      (  0   0   0   0  blt )
c
c           the diagonal blocks describe the kronecker canonical form 
c           (kcf) of the pencil a - lambda b as follows:
c
c             art - lambda brt has all right singular structure
c             azr - lambda bzr has all jordan structure for 0 eigenvalue
c             afn - lambda bfn has all jordan structure for finite
c                              nonzero eigenvalues
c             ain - lambda bin has all jordan structure for infinite
c                              eigenvalue
c             alt - lambda blt has all left singular structure
c
c           any subset of these blocks may not appear in a - lambda b.
c           the dimensions of these blocks are given by the following
c           integer output parameters:
c
c             rtre, rtce - last row and column of art, brt blocks
c                          (if both are zero, no right singular 
c                           structure)
c                          (if rtre.eq.0 and rtce.gt.0 then only l(0) 
c                           blocks in kcf)
c
c             zrre, zrce - last row and column of azr, bzr blocks
c                          (if zrre.eq.rtre and zrce.eq.rtce then no 
c                           0 eigenvalue)
c
c             fnre, fnce - last row and column of afn, bfn blocks
c                          (if fnre.eq.zrre and fnce.eq.zrce then no
c                           finite nonzero eigenvalues)
c
c             inre, ince - last row and column of ain, bin blocks
c                          (if inre.eq.fnre and ince.eq.fnce then no
c                           infinite eigenvalues)
c
c          notes: m, n are last row and column of alt, blt blocks.
c                 if inre.eq.m and ince.eq.n then no left singular
c                      structure.
c                 if inre.lt.m and ince.eq.n then only l(0)**t blocks
c                      in kcf.
c                 ince-rtce = inre-rtre = dimension of regular part.
c                 ince-fnce = inre-fnre = multiplicity of infinite
c                                         eigenvalue.
c                 fnce-zrce = fnre-zrce = total multiplicity of finite
c                                         nonzero eigenvalues
c                 zrce-rtce = zrre-rtre = multiplicity of 0 eigenvalue
c
c           the block structure of all the blocks (except afn and bfn)
c           are described by the integer output parameters:
c
c           pstruc(4)   integer, see below
c           struc(*)    integer, see below
c
c           (for more details about the block structure of (art,brt)
c            and (azr,bzr) see the output from routine rzstr. for 
c            more details about the block structure of (ain,bin)
c            and  (alt,blt) see the output from routine listr.)
c
c             struc(1 : pstruc(1)) describes the structure of art, brt
c                 (if pstruc(1).eq.0 then art and brt are not present).
c                 art and brt are both block upper triangular.
c                 the number of column blocks are pstruc(1) of
c                 dimensions struc(1) ... struc(pstruc(1)).
c                 the number of row blocks are pstruc(1)-1 of
c                 dimensions struc(2) ... struc(pstruc(1)).
c                 if pstruc(1).eq.1 then art and brt are 
c                 '0 by struc(1)' representing struc(1) zero columns
c                 (l(0) blocks in the kcf). 
c                 the number of l(j) blocks in the kcf is given by
c                 struc(j+1) - struc(j+2) for
c                 j.le.pstruc(1)-2 and struc(pstruc(1)) for 
c                 j.eq.pstruc(1)-1.
c
c             struc(pstruc(1)+1 : pstruc(2)) describes the structure 
c                 of azr, bzr (if pstruc(2).eq.pstruc(1) then azr and
c                 bzr are not present). 
c                 azr and bzr are both block upper triangular with
c                 pstruc(2)-pstruc(1) column and row blocks of
c                 dimensions struc(pstruc(1)+1) ... struc(pstruc(2)).
c                 the number of j by j jordan blocks for the zero
c                 eigenvalue in the kcf is given by
c                 struc(pstruc(1)+j)-struc(pstruc(1)+j+1) for
c                 j.le.pstruc(2)-pstruc(1)-1, and struc(pstruc(2))
c                 for j.eq.pstruc(2)-pstruc(1)
c
c             afn and bfn are both upper triangular. the finite nonzero
c                 eigenvalues of a - lambda b are given by the ratios
c                 afn(i,i)/bfn(i,i) of the diagonal entries of afn 
c                 and bfn.
c
c             struc(pstruc(2)+1 : pstruc(3)) describes the structure 
c                 of ain, bin (if pstruc(3).eq.pstruc(2) then ain and
c                 bin are not present).
c                 ain and bin are both block upper triangular with
c                 pstruc(3)-pstruc(2) column and row blocks of
c                 dimensions struc(pstruc(2)+1) ... struc(pstruc(3)).
c                 the number of j by j jordan blocks for the infinite 
c                 eigenvalue in the kcf is given by
c                 struc(pstruc(2)+j)-struc(pstruc(2)+j+1) for
c                 j.le.pstruc(3)-pstruc(2)-1, and struc(pstruc(3))
c                 for j.eq.pstruc(3)-pstruc(2)
c           
c             struc(pstruc(3)+1 : pstruc(1)) describes the structure 
c                 of alt, blt (if pstruc(3).eq.pstruc(4)
c                 then alt and blt are not present).
c                 alt and blt are both block upper triangular.
c                 the number of row blocks are pstruc(4)-pstruc(3)
c                 of dimensions struc(pstruc(3)+1) ... struc(pstruc(4)).
c                 the number of column blocks are pstruc(4)-pstruc(3)-1 
c                 of dimensions struc(pstruc(3)+2) ... struc(pstruc(4)).
c                 if pstruc(4).eq.pstruc(3)+1 then alt and blt are 
c                 'struc(pstruc(4)) by 0' representing struc(pstruc(4)) 
c                 zero rows (l(0)**t blocks in the kcf).
c                 the number of l(j)**t blocks in the kcf is given by
c                 struc(pstruc(3)+j+1) - struc(pstruc(3)+j+2) for
c                 j.le.pstruc(4)-pstruc(3)-2, and struc(pstruc(4)) for 
c                 j.eq.pstruc(4)-pstruc(3)-1.
c                      
c
c        adelta         real*8, relative distance from input matrix a
c                 to output a (if zero true).
c                 should be no larger than about epsu
c                 (otherwise pencil has ill-conditioned structure)
c
c        bdelta         real*8, relative distance from input matrix b 
c                 to output b (if zero true). 
c                 should be no larger than about epsu
c                 (otherwise pencil has ill-conditioned structure)
c
c        info - 0 if normal return
c                1 if svd failed to converge somewhere 
c                2 if qz failed to converge
c                3 if failed index error
c                  (should never occur. if it does contact either
c                   author below)
c           (if more detailed debug info needed, turn on appropriate 
c            idbg flags)
c
c
c*****  work space 
c
c    size needed for use of linpack svd  
c***       work(*)    complex*16 - 2*(max(m,n)*max(m,n)) + m*n +
c***                               min(m,n)*min(m,n) + 6*max(m,n) +
c***                               min(m,n) + 1 locations
c
c    size needed for use of lapack svd, 920708
c       work(*)    complex*16 - 2*(max(m,n)*max(m,n)) + m*n +
c                               min(m,n)*min(m,n) + 11*max(m,n) +
c                               3*min(m,n) + 1 locations
c
c       kstr(4,*) - integer - 4*max(m,n) + 24 locations
c
c***********************************************************************
c
c****    this version dated june 16, 1987
c        authors: jim demmel and bo kagstrom
c
c        addresses:
c             jim demmel, courant institute, new york university,
c             215 mercer str., new york, ny 10012, usa
c             ( phone int: country code 01 -(212)998 3391) 
c             ( email: demmel at nyu.edu or
c                      na.demmel at score.stanford.edu )
c 
c             bo kagstrom, institute of information processing,
c             university of umea, s-901 87 umea, sweden
c             (phone int - country code 46 - 90165419)
c             (email: bokg at seumdc51.bitnet or
c                     na.kagstrom at score.stanford.edu )
c
c****    guptri uses the following functions and subroutines
c
c        kcfpack -  cident, cmatml, cmatmr, cmatpr, krnstr, listr
c                   norme, rzstr, updel, zqz
c
        real*8 norme
c       
c***** internal variables
c
        logical ldebug, first
        integer mnmin, mnmax, stwork, stx, stsx, stex, stq
        integer starow, stbrow, stw, stqrax, sty, stqty
        integer rzcase, rowb, colb, rowe, cole 
        integer i, j, nsingr, lastm1, kfirst, last
        integer nstep, licase, nsingl 
        integer ierr, nlast, nsqrd, mtimn, msqrd, strtph, strtq
        integer stck
        integer nsumrz, rsumrz, nsumli, rsumli, djordz, djordi, dimreg
        integer njordz, njordi
        real*8 addlta, bddlta, epsua, epsub, anorme, bnorme
        complex*16 dummy
c
c       set debug flag
        ldebug= (idbg(2).ne.0)

c
c****   initialize pp and qq to identity matrices
        call cident(pp,ldpp,m)
        call cident(qq,ldqq,n)
c
c**     accumulate total perturbation in adelta, bdelta
        adelta = 0.
        bdelta = 0.
c**     compute norms and thresholds
        anorme = norme(a, ldab, m, n)    
        bnorme = norme(b, ldab, m, n)
        epsua = anorme * epsu
        epsub = bnorme * epsu
c*****  allocate workspace
         mnmin = min0(m,n)
         mnmax = max0(m,n)
         nsqrd = n * n
         mtimn = m * n
         stwork = 1
c
c****    modified for use of lapack svd, 920708
c        stx = stwork + mnmax
         stx = stwork + 2*mnmin + 6*mnmax
         stsx = stx + mtimn
         stex = stsx + mnmin + 1
c****    6/18/87 fix
c         stq = stex + mnmax
c         starow = stq + nsqrd
         starow = stex + mnmax
c
         stbrow = starow + mnmax
c****    6/18/87 fix
c         stw = stbrow + mnmax
c         stqrax = stw + nsqrd
         stqrax = stbrow + mnmax
c
         sty = stqrax + mnmax
         stqty = sty + mnmax
c****    6/18/87
         stq = stqty + mnmax
         stw = stq + nsqrd
c
         if (ldebug) then
            write(outunit,1642) m,n,mnmin,mnmax,stwork,
     *      stx,stsx,stex,stq,starow,stbrow,stw,stqrax,sty,stqty
1642        format(' guptri - workspace for rzstr -', 5i5,/,1x,10i5)
         endif
c
c****   reduction 1: 
c       find and put the Jordan structure of the zero eigenvalue
c       and the right singular structure in upper left corner
c       of (a,b)
c
        if (ldebug) write(outunit,100) m,n,epsu
100     format(//'guptri - m,n,epsu=',2i3,d13.6,//,'reduction 1')
        first = .true.
        swap = .false.
        call rzstr('cind', a, b, ldab, m, n, 1, m, 1, n,
     *             first, zero, epsua, epsub, gap,
     *             pp, ldpp, qq, ldqq, kstr, 1, last, addlta,
     *             bddlta,
     *             work(stwork), work(stx), work(stsx), work(stex),
     *             work(stq), work(starow), work(stbrow), work(stw),
     *             work(stqrax), work(sty), work(stqty), info)
c        if (info.ne.0) return
c****    6/18/87
         if (info.ne.0) then
           if (ldebug) write(outunit,1030) 'after reduction 1, info=',
     +                                     info
           return
         endif
        if (ldebug) then
          write(outunit,102) last
102       format(/'kstr, last=',i3)
          write(outunit,103) (j,j=1,last)
          write(outunit,103) (kstr(1,j),j=1,last)
          write(outunit,103) (kstr(2,j),j=1,last)
103       format(20i4)
        endif
c
c**     update total perturbation
        call upddel(adelta, addlta)
        call upddel(bdelta, bddlta)
        if (ldebug) write(outunit,101) adelta,bdelta
101     format('accumulated perturbations in a,b = ',2d15.6)
c
c**     convert computed null space dimensions into kronecker indices
        call krnstr(m, n, kstr, 1, last, nsumrz, rsumrz, rzcase,
     *              nsingr, njordz, djordz)
c
c       check for error condition
        if (rzcase .eq. 7) then
c***      6/18/87
          if (ldebug) write(outunit,1030) 
     +                'after first krnstr, rzcase=',rzcase
          info = 3
          return
        endif
c
c****   reductions 2 and 3:
c       if there are both right singular blocks and jordan blocks
c       corresponding to the zero eigenvalue,  reduce again to
c       separate them
c
c*****  6/15/87
        if (nsingr.eq.0 .and. djordz.eq.0) then
c         no right singular or zero structure
          pstruc(1) = 0
          pstruc(2) = 0
        elseif (nsingr.gt.0 .and. djordz.eq.0) then
c         right structure but no zero structure
          do 7352 j = 1, last
            struc(j) = kstr(1,j)
 7352     continue
          pstruc(1) = last
          pstruc(2) = last
        elseif (nsingr.eq.0 .and. djordz.gt.0) then
c         no right structure but zero structure
          do 7353 j = 1, last-1
            struc(j) = kstr(1,j)
 7353     continue
          pstruc(1) = 0
          pstruc(2) = last-1
c****     6/15/87
c        elseif (nsingr.gt.0 .and. njordz.gt.0) then
        elseif (nsingr.gt.0 .and. djordz.gt.0) then
c
c****   reduction 2:
c         separate the right and zero structures
c         reduce first rsumrz rows, nsumrz columns, swapping roles
c         of a,b. insist on computing same right singular structure
c         as in reduction 1
c
          lastm1=last-1
          nlast=last
          kstr(3,last)=kstr(1,last)
          kstr(4,last)=kstr(2,last)
          if (kstr(3,last) .eq.0) nlast=nlast-1
          if (last.gt.1) then
            do 2 j=lastm1,1,-1
              kstr(4,j)=kstr(3,j+1)
              kstr(3,j)=kstr(4,j)+kstr(1,j)-kstr(2,j)
              if (kstr(3,j) .eq. 0) nlast=nlast-1
2           continue
          end if
c
c****     6/15/87
          pstruc(1) = nlast
          do 7354 j = 1, nlast
            struc(j) = kstr(3,j)
 7354     continue
c
          if (ldebug) then
            write(outunit,104) rsumrz,nsumrz 
104         format(/'reduction 2, rsumrz,nsumrz=',2i4/'newkst')
            write(outunit,103) (j,j=1,nlast)
            write(outunit,103) (kstr(3,j),j=1,nlast)
            write(outunit,103) (kstr(4,j),j=1,nlast)
          endif
          first = .false.
          swap = .true.
          call rzstr('rind', b, a, ldab, m, n, 1, rsumrz, 1, nsumrz,
     *              first, zero, epsub, epsua, gap,
     *              pp, ldpp, qq, ldqq, kstr(3,1), 1, nlast, bddlta,
     *              addlta,
     *              work(stwork), work(stx), work(stsx), work(stex),
     *              work(stq), work(starow), work(stbrow), work(stw),
     *              work(stqrax), work(sty), work(stqty), info)
c          if (info.ne.0) return
c****      6/18/87
           if (info .ne. 0) then
             if (ldebug) write(outunit,1030) 
     +                   'after reduction 2, info=',info
             return
           endif
c
c**       update total perturbation
          call upddel(adelta,addlta)
          call upddel(bdelta,bddlta)
          if (ldebug) write(outunit,101) adelta,bdelta
c
c****   reduction 3:
c         recompute the block structure of the zero eigenvalue.
c         insist on computing the same jordan structure as in
c         reduction 1
c
          if (djordz.gt.1) then
            kstr(3,last)=0
            kstr(4,last)=0
            nlast=last-1
            if (last.gt.1) then
              do 4 j=lastm1,1,-1
                kstr(4,j)=kstr(3,j+1)+kstr(2,j)-kstr(1,j+1)
                kstr(3,j)=kstr(4,j)
                if (kstr(3,j) .eq. 0) nlast=nlast-1
4             continue
            end if
c
c*****      6/15/87
            pstruc(2) = pstruc(1) + nlast
            do 7355 j = 1, nlast
              struc(pstruc(1)+j) = kstr(3,j)
 7355       continue
c
            rowb=rsumrz-djordz+1
            colb=nsumrz-djordz+1
            if (ldebug) then
              write(outunit,105) rowb,colb
105           format(/'reduction 3, rowb,colb=',2i4/'newkst')
              write(outunit,103) (j,j=1,nlast)
              write(outunit,103) (kstr(3,j),j=1,nlast)
              write(outunit,103) (kstr(4,j),j=1,nlast)
            endif
            first = .false.
            swap = .false.
c****       6/18/87 bug fix, 'nlast' used to be 'last'
            call rzstr('rind', a, b, ldab, m, n, rowb, rsumrz, colb,
     *               nsumrz, first, zero, epsua, epsub, gap,
     *               pp, ldpp, qq, ldqq, kstr(3,1), 1, nlast, addlta,
     *               bddlta,
     *               work(stwork), work(stx), work(stsx), work(stex),
     *               work(stq), work(starow), work(stbrow), work(stw),
     *               work(stqrax), work(sty), work(stqty), info)
c            if (info.ne.0) return
c****        6/18/87
             if (info .ne. 0 ) then
               if (ldebug) write(outunit,1030) 
     +                     'after reduction 3, info=', info
               return
             endif     
c
c**         update total perturbation
            call upddel(adelta,addlta)
            call upddel(bdelta,bddlta)
            if (ldebug) write(outunit,101) adelta,bdelta
          else
c**       only a single zero eigenvalue, zero out the a-part
c
            if (zero) a(rsumrz,nsumrz) = 0.
c
c*****      6/15/87
            pstruc(2) = pstruc(1) + 1
            struc(pstruc(2)) = 1
c
          end if
c
c***    end of reductions 2 and 3
        end if
c
c**     if reduction complete, clean up kstr
c
        if (rzcase.ne.1 .and. rzcase.ne.4) then
          last=last+1
          kstr(1,last)=-1
          kstr(2,last)=-1
c
          nsumli=0
          rsumli=0
c*+
          djordi = 0
          dimreg = 0
c
c*****    6/15/87
          pstruc(3) = pstruc(2)
          pstruc(4) = pstruc(3)
c
c+*
c         if there is a common row nullspace, update kstr
          if (rzcase.eq.5 .or. rzcase.eq.6) then
            last=last+1
c
            kstr(1,last)=m-rsumrz
            kstr(2,last)=0
c*+
            nsumli=m-rsumrz
            rsumli=0
c
c****       6/15/87
            pstruc(4) = pstruc(3) + 1
            struc(pstruc(4)) = nsumli
c
c+*
          end if
          last=last+1
          kstr(1,last)=-1
          kstr(2,last)=-1
          last=last+1
          kstr(1,last)=-1
          kstr(2,last)=-1
        else
c
c         if no right or zero structure, fix kstr
          if (last.eq.1 .and. kstr(1,1).eq.0) last=0
c         put -1s at end of right, zero part of kstr
          last=last+1
          kstr(1,last)=-1
          kstr(2,last)=-1
c
c****     reduce the rest of the pencil
c
c**       allocate workspace for listr
c
          msqrd = m*m
c****     6/18/87
c          starow = stq + msqrd
c          stqrax = stw + msqrd
          stw = stq + msqrd
c
c
c****     reduction 4: 
c         find and put the jordan structure of the infinite
c         eigenvalue and the left singular structure in
c         lower right corner of (a,b)
c
          kfirst=last+1
          if (ldebug) write(outunit,107) kfirst
107       format(/'reduction 4, kfirst=',i4)
          rowb = rsumrz + 1
          colb = nsumrz + 1
          first = .false.
          swap = .false.
          call listr('cind', a, b, ldab, m, n, rowb, m, colb, n,
     *             first, zero, epsua, epsub, gap,
     *             pp, ldpp, qq, ldqq, kstr, kfirst, nstep, addlta,
     *             bddlta,
     *             work(stwork), work(stx), work(stsx), work(stex),
     *             work(stq), work(starow), work(stbrow), work(stw),
     *             work(stqrax), work(sty), work(stqty), info)
c          if (info.ne.0) return
c****      6/18/87
           if (info .ne. 0) then
             if (ldebug) write(outunit,1030)
     +                   'after reduction 4, info=',info
             return
           endif
c
          last=nstep+kfirst-1
          if (ldebug) then
            write(outunit,103) (j,j=1,last)
            write(outunit,103) (kstr(1,j),j=1,last)
            write(outunit,103) (kstr(2,j),j=1,last)
          endif
c
c**       update total perturbation
          call upddel(adelta, addlta)
          call upddel(bdelta, bddlta)
          if (ldebug) write(outunit,101) adelta,bdelta
c
c**       convert computed null space dimensions into kronecker indices
          call krnstr(n-nsumrz, m-rsumrz, kstr, kfirst, last,
     *               nsumli, rsumli, licase, nsingl, njordi, djordi)
c
          if (licase.eq.5 .or. licase.eq.6 .or. licase.eq.7) then
c           error condition - this should not happen because it would
c           mean there was right singular structure in this part
            if (ldebug) write(outunit,108) licase
108         format(//'error condition, licase=',i4)
            info = 3
            return
          end if
c
c****     reductions 5 and 6:
c         if there are both left singular blocks and jordan blocks
c         corresponding to the infinite eigenvalue,  reduce again 
c         to separate them
c
c*****    6/15/87
          if (nsingl.eq.0 .and. djordi.eq.0) then
c           no left or infinity structure
            pstruc(3) = pstruc(2)
            pstruc(4) = pstruc(3)
          elseif (nsingl.gt.0 .and. djordi.eq.0) then
c           left but no infinity structure
            pstruc(3) = pstruc(2)
            do 7356 j = kfirst, last
              struc(pstruc(3)+j-kfirst+1) = kstr(1,j)
 7356       continue
            pstruc(4) = pstruc(3) + last-kfirst+1
          elseif (nsingl.eq.0 .and. djordi.gt.0) then
c           no left but infinity structure
            do 7357 j = kfirst, last-1
              struc(pstruc(2)+j-kfirst+1) = kstr(1,j)
 7357       continue
            pstruc(3) = pstruc(2)+last-kfirst
            pstruc(4) = pstruc(3)
          elseif (nsingl.gt.0 .and. djordi.gt.0) then
c
c****     reduction 5:
c           separate the left and infinite structures.
c           reduce last rsumli columns and nsumli rows, swapping 
c           roles of a,b. insist on computing same left singular 
c           structure as in reduction 4
c
            lastm1=last-1
            kstr(3,last)=kstr(1,last)
            kstr(4,last)=kstr(2,last)
            nlast=last
            if (kstr(3,last) .eq. 0) nlast=nlast-1
            if (last.gt.kfirst) then
              do 6 j=lastm1,kfirst,-1
                kstr(4,j)=kstr(3,j+1)
                kstr(3,j)=kstr(4,j)+kstr(1,j)-kstr(2,j)
                if (kstr(3,j) .eq. 0) nlast=nlast-1
6             continue
            end if
c
c*****      6/15/87
c           temporarily put left structure in struc before infinity
            pstruc(3) = pstruc(2) + nlast-kfirst+1
            do 7358 j = kfirst, nlast
              struc(pstruc(2)+j-kfirst+1) = kstr(3,j)
 7358       continue
c
            rowb = m-nsumli+1
            colb = n-rsumli+1
            if (ldebug) then
              write(outunit,109) rowb,colb
109           format(/'reduction 5, rowb,colb=',2i4/'newkst')
              write(outunit,103) (j,j=1,nlast)
              write(outunit,103) (kstr(3,j),j=1,nlast)
              write(outunit,103) (kstr(4,j),j=1,nlast)
            endif
            nstep = nlast-kfirst+1
            first = .false.
            swap = .true.
            call listr('rind', b, a, ldab, m, n, rowb, m, colb, n,
     *               first, zero, epsub, epsua, gap,
     *               pp, ldpp, qq, ldqq, kstr(3,1), kfirst,
     *               nstep, bddlta, addlta,
     *               work(stwork), work(stx), work(stsx), work(stex),
     *               work(stq), work(starow), work(stbrow), work(stw),
     *               work(stqrax), work(sty), work(stqty), info)
c            if (info .ne. 0) return
c****        6/18/87
             if (info .ne. 0 ) then
               if (ldebug) write(outunit,1030)
     +                     'after reduction 5, info=',info
               return
             endif
c
c**         update total perturbation
            call upddel(adelta, addlta)
            call upddel(bdelta, bddlta)
            if (ldebug) write(outunit,101) adelta,bdelta
c
c****     reduction 6:
c           recompute the block structure of the infinite eigenvalue.
c           insist on computing the same jordan structure as
c           in reduction 4.
c
            if (djordi.gt.1) then
              kstr(3,last)=0
              kstr(4,last)=0
              nlast=last-1
              if (last.gt.kfirst) then
                do 8 j=lastm1,kfirst,-1
c*+
c                 kstr(4,j)=kstr(3,j+1)+kstr(2,j)-kstr(1,j)
                  kstr(4,j)=kstr(3,j+1)+kstr(2,j)-kstr(1,j+1)
c
                  kstr(3,j)=kstr(4,j)
                  if (kstr(3,j) .eq. 0) nlast=nlast-1
8               continue
              end if
c
c*****        6/15/87
c             move left structure right nlast-kfirst+1 places
              do 7359 j = pstruc(3),pstruc(2)+1,-1
                struc(j+nlast-kfirst+1) = struc(j)
 7359         continue
              pstruc(4) = pstruc(3) + nlast - kfirst +1
              pstruc(3) = pstruc(2) + nlast - kfirst +1
              do 7360 j = kfirst, nlast
                struc(j+pstruc(2)-kfirst+1) = kstr(3,j)
 7360         continue
c
              rowb = m-nsumli+1
              rowe = rowb+djordi-1
              colb = n-rsumli+1
              cole = colb+djordi-1
              if (ldebug) then
                write(outunit,111) rowb,colb,rowe,cole
111             format(/'reduction 6, rowb,colb,rowe,cole=',
     +                 4i4/'newkst')
                write(outunit,103) (j,j=1,nlast)
                write(outunit,103) (kstr(3,j),j=1,nlast)
                write(outunit,103) (kstr(4,j),j=1,nlast)
              endif
              nstep = nlast-kfirst+1
              first = .false.
              swap = .false.
              call listr('rind', a, b, ldab, m, n, rowb, rowe, colb,
     *               cole, first, zero, epsua, epsub, gap,
     *               pp, ldpp, qq, ldqq, kstr(3,1), kfirst,
     *               nstep, addlta, bddlta,
     *               work(stwork), work(stx), work(stsx), work(stex),
     *               work(stq), work(starow), work(stbrow), work(stw),
     *               work(stqrax), work(sty), work(stqty), info)
c              if (info .ne. 0) return
c****          6/18/87
               if (info .ne. 0) then
                 if (ldebug) write (outunit,1030)
     +                       'after reduction 6, info=',info
                 return
               endif
c
c**           update total perturbation
              call upddel(adelta,addlta)
              call upddel(bdelta,bddlta)
              if (ldebug) write(outunit,101) adelta,bdelta
c
            else
c**           only single infinite eigenvalue, zero out the b-part
              if (zero) b(m-nsumli+1, n-rsumli+1) = 0.
c
c*****        6/15/87
c             move struc left one place
              do 7361 j = pstruc(3),pstruc(2)+1,-1
                struc(j+1) = struc(j)
 7361         continue
              pstruc(4) = pstruc(3) +1
              pstruc(3) = pstruc(2) +1
              struc(pstruc(3)) = 1
c
            end if
c
c***      end of reductions 5 and 6
          end if
c
c***      change adelta and bdelta to relative perturbations
          if (anorme .ne. 0.) adelta = adelta / anorme
c         otherwise both anorme and adelta are 0.
          if (bnorme .ne. 0.) bdelta = bdelta / bnorme
c         otherwise both bnorme and bdelta are 0.
c
c***      clean up kstr
c         if there are no left or infinite indices, shorten kstr
          if (kfirst.eq.last .and. kstr(1,last).eq.0) last=last-1
c
          last=last+1
          kstr(1,last)=-1
          kstr(2,last)=-1
c         if there is a regular part with nonzero, noninfinite entries,
c         update kstr
          dimreg=0
          if (licase.eq.1 .or. licase.eq.4) then
c
c*+
            last=last+1
c
            kstr(1,last)=m-rsumrz-nsumli
            kstr(2,last)=n-nsumrz-rsumli
            dimreg=kstr(1,last)
          end if
          last=last+1
          kstr(1,last)=-1
          kstr(2,last)=-1
c
        end if
c
c****  6/15/87
c      compute output indices
       rtre = rsumrz - djordz
       rtce = nsumrz - djordz
       zrre = rsumrz
       zrce = nsumrz
       fnre = zrre + dimreg
       fnce = zrce + dimreg
       inre = fnre + djordi
       ince = fnce + djordi
c
        if (ldebug) then
          write(outunit,112)
112       format(//'final kstr=')
          write(outunit,103) (j,j=1,last)
          write(outunit,103) (kstr(1,j),j=1,last)
          write(outunit,103) (kstr(2,j),j=1,last)
          write(outunit, 1030) 'nsumrz= ',nsumrz, 'rsumrz=',rsumrz,
     *     'djordz=', djordz,'nsumli= ', nsumli, 'rsumli=',rsumli,
     *     'djordi=', djordi, 'dimreg=', dimreg
 1030     format(t5,a,i5)
c
c****     6/15/87
          write(outunit,1031) (pstruc(j),j=1,4)
 1031     format(//'final pstruc= ',4i4,/,'final struc =')
          if (pstruc(4).gt.0) write(outunit,1032) 
     +                        (struc(j),j=1,pstruc(4))
 1032     format(15i4)
          write(outunit,1033) 'rtce=',rtce,'zrce=',zrce,'fnce=',fnce,
     +                        'ince=',ince,'rtre=',rtre,'zrre=',zrre,
     +                        'fnre=',fnre,'inre=',inre
 1033     format(4(3x,a,i4),/,4(3x,a,i4))
c
        endif
c
c**** reduction 7:
c     reduce remaining regular part (corresponding to
c     the nonzero and finite eigenvalues) to upper
c     triangular form by using the qz algorithm
c
      if (dimreg .gt. 1) then
c       rowb = first row of remaining regular part
c       colb = first column of remaining regular part
c
        rowb = rsumrz  + 1
        colb = nsumrz +1
c
c****   reduce the pencil a(rowb:rowb+dimreg-1,colb:colb+dimreg-1)
c       - lambda b(rowb:rowb+dimreg-1,colb:colb+dimreg-1)
c       to upper triangular form with the qz algorithm
c
c**     allocate workspace for transformation matrices phtemp and qtemp
        strtq = 1
        strtph = strtq + dimreg * dimreg
        stck = strtph + dimreg * dimreg
        call zqz(a, b, ldab, dimreg, rowb, colb, work(strtq),
     *           dimreg, work(strtph), dimreg, ierr, work(stck))
        if (ierr .ne. 0) then
c****     6/18/87
          if (ldebug) write (outunit,1030) 'after qz, ierr=',ierr
          info = 2
          return
        endif
c**     update rows 1 to rowb-1 above the remaining regular part
c       in columns colb to colb+dimreg-1 
c       by postmultiplying with qtemp (dimreg*dimreg)
c
        call cmatmr( a( 1, colb), ldab, rowb-1, dimreg,
     *              work(strtq), dimreg, dimreg, dummy, 1,
     *              work(stck), 1)
        call cmatmr( b( 1, colb), ldab, rowb-1, dimreg,
     *              work(strtq), dimreg, dimreg, dummy, 1, 
     *              work(stck), 1)
c
c**     update (rows 1 to n in) columns colb to colb+dimreg-1
c       of qq by postmultiplying with qtemp
c
CEE 930610        call cmatmr( qq( 1, colb), ldab, n, dimreg,
        call cmatmr( qq( 1, colb), ldqq, n, dimreg,
     *              work(strtq), dimreg, dimreg, dummy, 1,
     *              work(stck), 1)
c
c**     update columns colb+dimreg to n to the right of the remaining
c       regular part in rows rowb to rowb+dimreg-1
c       by premultiplying by phtemp (dimreg*dimreg)
c
        call cmatml( a( rowb, colb+dimreg), ldab, dimreg,
     *              n-colb-dimreg+1, work(strtph),
     *              dimreg, dimreg, dummy, 1, work(stck), 1)
c
        call cmatml( b( rowb, colb+dimreg), ldab, dimreg,
     *              n-colb-dimreg+1, work(strtph),
     *              dimreg, dimreg, dummy, 1, work(stck), 1)
c
c**     update (rows 1 to m in) columns rowb to rowb+dimreg-1
c       of pp by postmultiplying with phtemp**h
c
        call cmatmr( pp( 1, rowb), ldpp, m, dimreg,
     *              work(strtph), dimreg, dimreg, dummy, 1,
     *              work(stck), 3)
      endif
	  call codim(kstr,dimreg)
c
        if (idbg(2) .gt. 1) then
            call cmatpr(qq,ldqq,n,n,'qq at exit from guptri')
            call cmatpr(pp,ldpp,m,m,'pp at exit from guptri')
        endif
c
c
      if (ldebug) then
         write(outunit, 2005) 'computed eigenvalues'
 2005    format( t5, a, 4d15.5)
c****    6/19/87
         rowb = rtre+1
         rowe = inre
         colb = rtce+1
         cole = ince
c
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
c***   end of reduction 7

       return
	   end
