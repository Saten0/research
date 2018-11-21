c---------------------------------------------------------------
c  GA法を用いた直接結合型λ/4共振器の設計
c  入出力結合：キャパシタ      
c  未知数： 長さ al11, 特性インピーダンス z0
c  入力側キャパシタ：aCin1, 出力側キャパシタ：aCout1
c                            written by T.Ohno(Kiyomi)
c---------------------------------------------------------------

      parameter(ne=10, ncl=800, mpop=800, mfr=90, nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)!,fre(mfr)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
c !!!      
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max

      open(unit=1,file='parameter_open',status='old')
      open(unit=4,file='result.dat')

      call ingeom               !Read geometry
      call strtpop              !Start population of designs
      call optim                !
      call outputf              !Output file - save results

      stop
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine ingeom         !ファイルの読み込みとビット化

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
c      common/frequency/freq
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max

      read(1,*)jran             !random number seed
      read(1,*)pcross,pmut      !Probabiliy of cross-over and mutation
      read(1,*)temp             !Temperature
      read(1,*)npop,naCin1,naCout1,nZa,nZb,nal1,nal2 !# of designs, Chromo. length
c !!!
      read(1,*)aCin1_min,aCin1_max
      read(1,*)aCout1_min,aCout1_max
      read(1,*)Za_min,Za_max
      read(1,*)Zb_min,Zb_max
      read(1,*)al1_min,al1_max
      read(1,*)al2_min,al2_max

c      read(1,*)epsi_min,epsi_max

c      npop  - number of elements in the population
c      nal11   - n. of bits for thickness of leyer #1.
c      nz0 -  

c ビット化
c !!!
      al11_step   = (al11_max - al11_min)/(2**nal11-1)
      z0_step = (z0_max - z0_min)/(2**nz0-1)
      aCin1_step = (aCin1_max - aCin1_min)/(2**naCin1-1)
      aCout1_step = (aCout1_max - aCout1_min)/(2**naCout1-1)
      Za_step = (Za_max - Za_min)/(2**nZa-1)
      Zb_step = (Zb_max - Zb_min)/(2**nZb-1)
      al1_step = (al1_max - al1_min)/(2**nal1-1)
      al2_step = (al2_max - al2_min)/(2**nal2-1)
c      epsi_step = (epsi_max - epsi_min)/(2**nepsi-1)
c !!!
      write(*,*)'--------------------------------------------'
      write(*,*)'nal1,al1_step    ',nal1,al1_step
      write(*,*)'nal2,al2_step      ',nal2,al2_step
      write(*,*)'naCin1,aCin1_step    ',naCin1,aCin1_step
      write(*,*)'naCout1,aCout1_step    ',naCout1,aCout1_step
      write(*,*)'nZa,Za_step              ',nZa,Za_step
      write(*,*)'nZb,Zb_step                ',nZb,Zb_step
c      write(*,*)'nepsi,epsi_step',nepsi,epsi_step
      write(*,*)'--------------------------------------------'
      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine strtpop        !初期集団の生成(0 or 1の付与)

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max


      jlo=0
      jhi=1
c      nbit=nal11+nz0+nepsi      !一つの個体の総bit数
c !!!
      nbit=nal1+nal2+nZa+nZb+naCin1+naCout1       !一つの個体の総bit数

      do i=1,npop               !All the designs
         do j=1,nbit            !All the chromosome length
            i1=nranu(jran)
            if (i1.EQ.1) then
               bpop1(i,j)=.TRUE.
            else
               bpop1(i,j)=.FALSE. !Most "false"
            endif
         enddo
      enddo

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine optim          !個集団の最適化

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
      real r0(1000)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max
      common /obj_fun/d_0,r_0,d_glob,r_glob
      common /iter/iter
      rmaov=-1000000.

c   反復回数：n1回      
      call object               !Performance of each design

      write(*,*)'  rmaobj    rtoobj/npop   rmaov'
      write(*,*)rmaobj,rtoobj/npop,rmaov
      write(*,*)'Give the Number of Initial Iterations'
      read(*,*)n1

 101  do iter=1,n1
         call reprod
         call cross
         call mutat
         call object
         if (rmaobj.GT.rmaov) then
            rmaov=rmaobj        !最大値を残す
            do i=1,nbit
               bopt(i)=bpop1(jopt,i)
            enddo
            call rtf(jopt,r0,1)
         endif

c         write(*,*)'iter  rmaobj  rtoobj/npop   rmaov'
c         write(*,*)iter,rmaobj,rtoobj/npop,rmaov

      enddo

      write(*,*)'Give the number of additional iterations'
      read(*,*)n1
!      write(*,*)'Give mutation probability'
!      read(*,*)pmut
c      pmut=0.01
!      write(*,*)'Give temperature'
!      read(*,*)temp
c      temp=0.02
      write(*,*)'New seed?'
      read(*,*)nseed

      jlo=1
      jhi=npop

      if (nseed.GT.0) then
         do iseed=1,nseed
            nsch=nranu(jran)
            do it=1,nbit
               bpop1(nsch,it)=bopt(it)
            enddo
         enddo
      endif

      if (n1.GT.0) goto 101

      return
      end
      
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine object

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
      real r0(1000)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max
      common /obj_fun/d_0,r_0,d_glob,r_glob

      rmiobj=-1000000.
      rmaobj=0.
      rtoobj=0                  !Total performance
      jopt=0

      do i=1,npop
         call rtf(i,r0,0)
         do k=1,490
            if (r0(k).GT.10) then
               goto 12
            endif
         enddo
         do k=510,700
            if (r0(k).GT.10) then
               goto 12
            endif
         enddo
         robj(i)=r0(240)
         if (r0(240).GT.rmaobj) then
            rmaobj=r0(240)           !Search for the best design
            jopt=i
         endif
         if (r0(240).LT.rmiobj) rmiobj=r0(240)
 12      rtoobj=rtoobj+r0(240)
      enddo

      return
      end
      
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine reprod         !増殖

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max


      jlo=1
      jhi=npop
      do i=1,npop
c         n1=nranu(jran)
         n1=i
         n2=nranu(jran)
         r1=robj(n1)
         r2=robj(n2)
         pch=ranu(jran)
         pcc=0.5*exp(-ABS(r2-r1)/temp)
         if (r1.GT.r2) pcc=1-pcc
         n3=n1
         if (pch.GT.pcc) n3=n2
         do j=1,nbit
            bpop2(i,j)=bpop1(n3,j)
         enddo
      enddo

c      write(*,*)'Reprod. finished'

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              交差(1点交差,2点交差,複数点交差)
c          (個集団NPOPに対して交差する確率crate=0.20)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cross

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max


      jlo=1
      jhi=nbit-1

      do i=1,npop/2
         i1=2*i-1
         i2=i1+1
         r1=ranu(jran)
         if (r1.LT.pcross) then
            n1=nranu(jran)
            do j=1,n1
               bpop1(i1,j)=bpop2(i1,j)
               bpop1(i2,j)=bpop2(i2,j)
            enddo
            do j=n1+1,nbit
               bpop1(i1,j)=bpop2(i2,j)
               bpop1(i2,j)=bpop2(i1,j)
            enddo
         else
            do j=1,nbit
               bpop1(i1,j)=bpop2(i1,j)
               bpop1(i2,j)=bpop2(i2,j)
            enddo
         endif
      enddo

c      write(*,*)'Cross finished'

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine mutat          !突然変異 確率：0<pmult<1

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max


      ntot=nbit*npop
      rmut=ntot*pmut

      if (rmut.LT.4) then
         r1=ranu(jran)
         r0=exp(-rmut)
         nmut=0
         if (r0.GT.R1) return
         r0t=r0
 8       nmut=nmut+1
         r0=r0*rmut/nmut
         r0t=r0t+r0
         if (r0t.LT.r1) goto 8
      else
         nmut=int(rmut)
      endif

      do i=1,nmut
         jlo=1
         jhi=npop
         n1=nranu(jran)
         jlo=1
         jhi=nbit
         n2=nranu(jran)
         if (bpop1(n1,n2)) then
            bpop1(n1,n2)=.FALSE.
         else
            bpop1(n1,n2)=.TRUE.
         endif
      enddo

c      write(*,*)'Mutat. finished'

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine outputf

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
      real r0(1000)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max

      do it=1,nbit
         bpop1(1,it)=bopt(it)
      enddo

      iond=1
      iout=1
      
      call rtf(iond,r0,iout)

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine rtf(iond,r0,iout) !求めた群の判定

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl),bopt(ncl)
      dimension robj(mpop)
      real r0(1000)
c !!!
c      common /mat/npop,nal11,nz0,nbit
      common /mat/npop,nal11,nz0,naCin1,naCout1,nbit,nZa,nZb,nal1,nal2
      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt
      common /range/jhi,jlo
      common /randm/jran,pmult
      common /prob/pcross,pmut,temp
      common /iter/iter
      common /sopt/bopt,rmaov
      common /comp/al11_min,al11_max,al11_step,z0_min,z0_max,z0_step
      common /comq/aCin1_min,aCin1_max,aCin1_step
      common /comr/aCout1_min,aCout1_max,aCout1_step
      common /coms/Za_min,Zb_min,al1_min,al2_min
      common /comt/Za_step,Zb_step,al1_step,al2_step
      common /comu/Za_max,Zb_max,al1_max,al2_max

      i=1
      nn=0


      call trans(iond,nn,naCin1,nv)
      aaCin1= aCin1_min + aCin1_step*nv

      nn=nn+naCin1

      call trans(iond,nn,naCout1,nv)
      aaCout1= aCout1_min + aCout1_step*nv

      nn=nn+naCout1

      call trans(iond,nn,nZa,nv)
      Za= Za_min + Za_step*nv

      nn=nn+nZa

      call trans(iond,nn,nZb,nv)
      Zb= Zb_min + Zb_step*nv

      nn=nn+nZb

      call trans(iond,nn,nal1,nv)
      al1= al1_min + al1_step*nv

      nn=nn+nal1


      call trans(iond,nn,nal2,nv)
      al2= al2_min + al2_step*nv

      nn=nn+nal2


c      call trans(iond,nn,nepsi,nv)
c      ei=epsi_min + epsi_step*nv

c      write(*,*)'z0,al11,aaCin1,aaCout1',az0,al11,aaCin1,aaCout1

      call analyze(aaCin1,aaCout1,r0,Za,Zb,al1,al2)

c ===============================
c            Output
c ===============================
      
      if (iout.EQ.1) then
         write(*,*)
         write(*,*)'-------------------------------------------'
         write(*,*)'Za[ohm],Zb[ohm],al1[mm],al2[mm]'
         write(*,*)'aaCin[pF],aaCout1[pF],r0[-],iter[-]'
         write(*,*)Za,Zb,al1,al2
         write(*,*)aaCin1,aaCout1,r0(240),iter
         write(*,*)'-------------------------------------------'
         write(4,*)'Za[ohm],Zb[ohm],al1[mm],al2[mm]'
         write(4,*)'aaCin[pF],aaCout1[pF],r0[-],iter[-]'
         write(4,*)Za,Zb,al1,al2
         write(4,*)aaCin1,aaCout1,r0(240),iter
         write(4,*)'-------------------------------------------'
      endif

      return
      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine analyze(aaCin1,aaCout1,r0,Za,Zb,al1,al2) !解析

      implicit real (a-z)
      complex j
      complex s11a,s21a,s12a,s22a
      complex s11b,s21b,s12b,s22b
      complex s11c,s21c,s12c,s22c
      complex Zin1,Zin2,y0
      real r0(1000)
      data v,pi/2.99792458E+08, 3.141592653589/

      l= al11*1E-3
      l1= al1*1E-3
      l2= al2*1E-3
      Cin= aaCin1*1E-12
      Cout= aaCout1*1E-12
      j=(0.0,1.0)
      
c      do fff=0.01,2.41,0.01
c      fff=2.40
c 520  ff= fff*1E9
      do fff=1,700,1
         ff=fff*1E7
         ww= 2*pi*ff
         beta= ww/v

ccc Input-side capacitor
         s11a= 1/(1+j*ww*Cin*100)
         s21a= 1-s11a
         s22a= s11a
         s12a= s21a

ccc SIR
         Zin2= j*Zb*TAN(beta*l2)
         Zin1= Za*(Zin2+j*Za*TAN(beta*l1))/(Za+j*Zin2*TAN(beta*l1))
         y0= 1/Zin1

         s11b= -y0/(0.04+y0)
         s21b= 1+s11b
         s22b= s11b
         s12b= s21b

ccc Cal. for S-parameters
         call series(s11a,s21a,s22a,s12a,s11b,s21b,s22b,s12b,
     &        s11c,s21c,s22c,s12c)
         s11a= s11c
         s21a= s21c
         s12a= s12c
         s22a= s22c
         
ccc Output-side capacitor
         s11b= 1/(100*j*ww*Cout+1)
         s21b= 1-s11b
         s22b= s11b
         s12b= s21b
         
ccc Results of S-parameters
         call series(s11a,s21a,s22a,s12a,s11b,s21b,s22b,s12b,
     &        s11c,s21c,s22c,s12c)
         s11a= s11c
         s21a= s21c
         s12a= s12c
         s22a= s22c

         as11= ABS(s11a)
         as21= ABS(s21a)

         as11= (20*log10(as11))
         as21= (20*log10(as21))

         r0(fff)= -as11           !評価値

      enddo
      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine series(s11a,s21a,s22a,s12a,s11b,s21b,s22b,s12b,
     &     s11c,s21c,s22c,s12c)

      complex s11a,s21a,s12a,s22a
      complex s11b,s21b,s12b,s22b
      complex s11c,s21c,s12c,s22c
      
      s11c= s11a+(s21a*s11b*s12a)/(1-s11b*s22a)
      s21c= (s12a*s12b)/(1-s11b*s22a)
      s12c= (s21a*s21b)/(1-s11b*s22a)
      s22c= s22b+(s12b*s22a*s21b)/(1-s11b*s22a)

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine trans(iond,nn,ng,ns) !bit列からの再変換

      parameter(ne=10,ncl=800,mpop=800,mfr=90,nuk=150)
      implicit real(a-b,d-h,o-z),complex(c),integer(i-n)
      logical bpop1(mpop,ncl),bpop2(mpop,ncl)
      dimension robj(mpop)

      common /pop/bpop1,bpop2,robj,rmaobj,rmiobj,rtoobj,jopt

      ns=0                      !Conversion from logical
      
      do jt=1,ng                !bit-strings to integers
         if (bpop1(iond,nn+jt)) ns=ns+2**(jt-1)
      enddo

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function ranu(jran)

      data im/53125/
      data ia/171/
      data ic/11213/

      jran=MOD(jran*ia+ic,im)
      ranu=FLOAT(jran)/FLOAT(im)

      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function nranu(jran)
      common /range/jhi,jlo

      data im/53125/
      data ia/171/
      data ic/11213/

      jran=MOD(jran*ia+ic,im)
      nranu=jlo+((jhi-jlo+1)*jran)/im

c      write(*,*)nranu
      return
      end

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
