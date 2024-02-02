c---------------------------------------------------------------------
c author:  Jave G. Pascual, math&stats - Washington State University
c
c These codes were written for Fortran77.  
c They have not been tested on other Fortran versions.
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c Demo of computations of RFL (version 2) pdf, cdf
c
c The calculations are based on a standardized version of the 
c   random-fatigue limit (RFL) model described in the paper:
c     Pascual, F. G. (2003). A Standardized Form of the Random 
c     Fatigue-Limit Model. Communications in Statistics - Simulation 
c     and Computation, 32, 1209-1228.
c
c---------------------------------------------------------------------
c      implicit real*8(a-h,o-z)
c      data ndist1/1/
c      data ndist2/1/
c      data beta1/-2.5d00/
c      data sdgamma/.044d00/
c
c compute pdf/cdf
c
c      xiincre=17.d00/39
c      do 200 j=5,5
c         xi=-1.d00+(j-1.d00)*xiincre
c         zincre=40.d00/9
c         print *, xi
c         do 100 i=1,10
c            z=-10.d00+(i-1.d00)*zincre
c            call rfl0pdf(ndist1,ndist2,z,xi,beta1,sdgamma,
c     +           apdf,ierip,ierp1,ierp2)
c           isum=ierip+ierp1+ierp2
c            if (isum.gt.0) then
c               print *,"PDF:", apdf,ierip,ierp1,ierp2
c            end if
c            call rfl0cdf(ndist1,ndist2,z,xi,beta1,sdgamma,
c     +           acdf,ieric,iersc)
c            isum=ieric+iersc
c            if (isum.gt.0) then
c               print *,"CDF:", acdf,ieric,iersc
c            end if
c 100     continue
c 200  continue
c      stop
c      end

c----------------------------------------------------------------
      subroutine srfl0pdf(ndist1,ndist2,z,xi,beta1,
     +     sdgamma,num,answer,ieri,ier1,ier2)
c----------------------------------------------------------------
c
c R interface with Fortran to compute RFL pdf
c
c-------------------------------------
      implicit real*8(a-h,o-z)
      dimension ndist1(num),ndist2(num),z(num),xi(num),
     +     beta1(num),sdgamma(num),answer(num),
     +     ieri(num),ier1(num),ier2(num)
      do 10 i=1,num
         call rfl0pdf(ndist1(i),ndist2(i),z(i),xi(i),beta1(i),
     +        sdgamma(i),answer(i),ieri(i),ier1(i),ier2(i))
 10   continue
      return
      end

c----------------------------------------------------------------
      subroutine srfl0cdf(ndist1,ndist2,z,xi,beta1,sdgamma,
     +     num,answer,ieri,iers)
c----------------------------------------------------------------
c
c R interface with Fortran to compute RFL cdf
c
c-------------------------------------
      implicit real*8(a-h,o-z)
      dimension ndist1(num),ndist2(num),z(num),xi(num),
     +     beta1(num),sdgamma(num),answer(num),
     +     ieri(num),iers(num)
      do 10 i=1,num
         call rfl0cdf(ndist1(i),ndist2(i),z(i),xi(i),beta1(i),
     +        sdgamma(i),answer(i),ieri(i),iers(i))
 10   continue
      return
      end
	  
c----------------------------------------------------------------
      subroutine srfl0quan(ndist1,ndist2,xi,alpha,beta1,
     +     sdgamma,bd1,bd2,num,answer)
c----------------------------------------------------------------
c
c R interface with Fortran to compute (standard) RFL quantile
c
c-------------------------------------
      implicit real*8(a-h,o-z)
      dimension ndist1(num),ndist2(num),xi(num),alpha(num),
     +     beta1(num),sdgamma(num),answer(num)
      do 10 i=1,num
         xbd1=bd1
         xbd2=bd2
         call rfl0quan(ndist1(i),ndist2(i),xi(i),alpha(i),
     +        beta1(i),sdgamma(i),xbd1,xbd2,answer(i))
 10   continue
      return
      end
	  
c----------------------------------------------------------------
      subroutine sq2b(ndist1,ndist2,stress,alpha,quan,
     +     beta1,sigma,ugamma,sdgamma,bd1,bd2,beta0)
c----------------------------------------------------------------
c
c    reparameterize by switching a quantile with beta0
c
      implicit real*8(a-h,o-z)
      call quan2beta0(ndist1,ndist2,stress,alpha,quan,
     +     beta1,sigma,ugamma,sdgamma,bd1,bd2,beta0)
      return
      end

c----------------------------------------------------------------
      subroutine quan2beta0(ndist1,ndist2,stress,alpha,quan,
     +     beta1,sigma,ugamma,sdgamma,bd1,bd2,beta0)
c----------------------------------------------------------------
c
c       this is a subroutine to compute beta0 given the values
c       of beta1, sigma, ugamma, sdgamma and a quantile
c
c inputs:
c
c       ndist1, ndist2  integer codes for distributions
c                       1=sev, 2=normal, 3=logistic
c
c       beta1       slope of mean log(time) model eqn
c
c       xstr        stress level
c
c       sigma       std devn of log(time)
c
c       ugamma      mean of the fatigue limits
c
c       sdgamma     std devn of the fatigue limits
c 
c       stress      stress level
c 
c       alpha       proportion corresponding to quantile
c
c       quan        quantile value
c
c       b1, b2      bounds for beta0
c
c outputs:
c
c       beta0       
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer6/beta1p,sigmap,ugammap,sdgammap,ndist1p,
     +     ndist2p,stressp,alphap,quanp
      external zeroin, fcndiff6, stdcdf
c
c     constants that one might want to change to
c     achieve a higher degree of accuracy from the algorithm
c
      data tol/1.0d-8/
      data zero/0.d00/
      data top/500.d00/
	data bottom/-500.d00/
      xlog=dlog(stress)
c
c     common stuff
c
      beta1p=beta1
      sigmap=sigma
      ugammap=ugamma
      sdgammap=sdgamma
      stressp=stress
      alphap=alpha
      quanp=quan
      ndist1p=ndist1
      ndist2p=ndist2
c
      xi=(xlog-ugamma)/sdgamma
	beta1t=beta1/sigma
c
      check=stdcdf(xi,ndist2)
      if (check.gt.alpha) go to 200
c   quantile is not possible
         beta0=-1000000.0d00
         return

c
c     check if bounds are okay
c
 200	bd1=bd1-5.d00
	if (bd1.le.bottom) then
         beta0=bottom
         return
      end if
      bd1t=bd1+beta1*ugamma
	z=(log(quan)-bd1t)/sigma
	call rfl0cdf(ndist1,ndist2,z,xi,beta1t,sdgamma,ans,ieri,iers)
      if (ans.lt.alpha) go to 200
 201  bd2=bd2+5.d00
      if (bd2.ge.top) then
         beta0=top
         return
      end if
	bd2t=bd2+beta1*ugamma
      z=(log(quan)-bd2t)/sigma
      call rfl0cdf(ndist1,ndist2,z,xi,beta1t,sdgamma,ans,ieri,iers)
      if (ans.gt.alpha) go to 201
c
c     zero in on beta0 value
c

      beta0=zeroin(bd1,bd2,fcndiff6,tol)
      return
      end
	  
c----------------------------------------------------------------
      double precision function fcndiff6(x)
c----------------------------------------------------------------
c     
c     function to compute difference between cdf and alpha 
c
      implicit real*8(a-h,o-z)
      common/passer6/beta1p,sigmap,ugammap,sdgammap,ndist1p,
     +     ndist2p,stressp,alphap,quanp
	xi=(dlog(stressp)-ugammap)/sdgammap
	beta1t=beta1p/sigmap
      xt=x+beta1p*ugammap
 	z=(log(quanp)-xt)/sigmap
      call rfl0cdf(ndist1p,ndist2p,z,xi,beta1t,sdgammap,ans,ieri,iers)
      fcndiff6=ans-alphap		 
      return
      end


c----------------------------------------------------------------
      subroutine rfl0pdf(ndist1,ndist2,z,xi,beta1,sdgamma,
     +     answer,ieri,ier1,ier2)
c----------------------------------------------------------------
c
c          subroutine to compute the Random Fatigue-Limit
c          model (standardized) pdf where W is the loglife
c          to failure under stress xstr, W|V is has a 
c          location-scale distribution with
c            location=beta0+beta1*log(exp(x)-exp(v))
c               scale=sigma
c          where V=log(fatigue limit) and V is location-scale
c          with
c            location=ugamma
c               scale=sdgamma
c
c         The computations below use some transformations on
c         the model parameters.  Essentially, the algorithm
c         uses 2 parameters instead of 5.
c
c inputs:
c
c       ndist1, ndist2  integer codes for distributions
c                       1=sev, 2=normal
c
c       ndist1 = distribution for W|V 
c       nidst2 = distribution for V 
c
c       beta1     standardized slope = beta1/sigma
c
c       sdgamma   scale parameter of log(fatigue limit)
c
c       xi        standardized logstress = (x-ugamma)/sdgamma
c                 
c       z         standardized loglife = (w-beta0*)/sigma
c
c outputs:
c
c       answer      density value
c
c       ier1, ier2, ieri   return condition indicator
c                          0 if no errors were
c                          detected in dqags, dqagi, resp.
c                          see dqags, dqagi documentation for 
c                          meaning of values of ier > 0
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer/zp,xip,beta1p,sdgammap,ndist1p,
     +     ndist2p
c      common/vpass/valtolp
      external rfl0int
c      external pdfdiff
      dimension bds(3),iwork1(100),work1(400),iwork2(100),
     +     work2(400),iwork3(100),work3(400)
c
c     constants that one might want to change to
c     achieve a higher degree of accuracy from the algorithm
c
c      data eps/1.0d-30/
      data eps/1.0d-15/
      data difftol/1.0d-2/
c
      data limit/100/
      data zero/0.d00/
      ier1=0
      ier2=0
      ieri=0
c
c     variables for passing
c
      zp=z
      xip=xi
      beta1p=beta1
      sdgammap=sdgamma
      ndist1p=ndist1
      ndist2p=ndist2
c
c     do the integration (break into 3 intervals:
c     -infty - bds(1), bds(1) - bds(2), bds(2) - bds(3)
c     
      answer1=zero
      answer2=zero
      answer3=zero
c
c     sort possible bounds for integration
c
      valtol=1.d-25    
      if (abs(xi).gt.(6.d00)) then
         rlb=-3.d00
      else
         rlb=-abs(xi)/2.d00
      end if
 30   rlb=rlb-.01d00
      rflintval=rfl0int(rlb)
      if (rflintval.gt.valtol) goto 30
c      eps2=1.d-30
c      rlb=zeroin(rlb,zero,pdfdiff,eps2)
      bds(1)=rlb
      xtemp=1.d00-dexp(z/beta1)
      if (xtemp.le.zero) then
         bds(2)=rlb/2.d00
      else
         bds(2)=dlog(xtemp)/sdgamma
      end if
      bds(3)=zero
c
c sort the bounds
c
      do 10 i=1,2
         do 20 j=(i+1),3
            if(bds(i).gt.bds(j)) then
               vdummy=bds(i)
               bds(i)=bds(j)
               bds(j)=vdummy
            end if
 20      continue
 10   continue
      diff1=bds(2)-bds(1)
      difftol=sdgamma/2.d00
      if(diff1.lt.difftol) then
         if(bds(2).eq.zero) then
            bds(1)=zero-difftol
         end if
      end if
      diff2=bds(3)-bds(2)
      if(diff2.lt.difftol) then
         if(bds(3).eq.zero) then
            bds(2)=zero
         end if
      endif
      call dqagi(rfl0int,bds(1),-1,eps,eps,answer1,abserr,
     +     neval,ieri,limit,4*limit,last,iwork1,work1)
      if(zero.gt.bds(2)) then
         call dqags(rfl0int,bds(1),bds(2),eps,eps,answer2,
     +        abserr,neval,ier1,limit,4*limit,last,iwork2,work2)
      end if
      if(zero.gt.bds(2)) then
         call dqags(rfl0int,bds(2),bds(3),eps,eps,answer3,
     +        abserr,neval,ier2,limit,4*limit,last,iwork3,work3)
      end if
      answer=answer1+answer2+answer3
      return
      end

c----------------------------------------------------------------
      double precision function rfl0int(x)
c----------------------------------------------------------------
c
c       function to compute integrand for computing the RFL 
c       (standardized) pdf
c
c       the following variables are communicated through common
c
c          beta1p     slope of the mean loglife
c          zp         standardized loglife
c          xip        standardized log(stress)
c          ndist1p    distn of W|V
c          ndist2p    distn of V
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer/zp,xip,beta1p,sdgammap,ndist1p,
     +     ndist2p
      external stdpdf
      data zero/0.d00/
      if (x.eq.zero) then
         rfl0int=zero
      else
         z=zp-beta1p*(xip*sdgammap+dlog(1.d00-dexp(sdgammap*x)))
         zx=x+xip
         rfl0int=stdpdf(z,ndist1p)*stdpdf(zx,ndist2p)
      end if
      return
      end
	  
c----------------------------------------------------------------
      double precision function stdpdf(x,ndist)
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      data cval/.39894228040143d00/
      data half/.5d00/
      data one/1.d00/
      if (ndist.eq.1) then
c  sev distribution
         stdpdf=dexp(x-dexp(x)) 
      end if
      if (ndist.eq.2) then
c  normal distribution
         stdpdf=cval*dexp(-half*x*x)
      end if
      if (ndist.eq.3) then
c  logistic distribution
         stdpdf=dexp(x)/((one+dexp(x))**2.d00)
      end if
      return
      end
	  	  
c----------------------------------------------------------------
      subroutine rfl0cdf(ndist1,ndist2,z,xi,beta1,
     +     sdgamma,answer,ieri,iers)
c----------------------------------------------------------------
c
c          subroutine to compute the Random Fatigue-Limit
c          model (standardized) cdf where W is the loglife
c          to failure under stress xstr, W|V is has a 
c          location-scale distribution with
c            location=beta0+beta1*log(exp(x)-exp(v))
c               scale=sigma
c          where V=log(fatigue limit) and V is location-scale
c          with
c            location=ugamma
c               scale=sdgamma
c
c         The computations below use some transformations on
c         the model parameters.  Essentially, the algorithm
c         uses 2 parameter instead of 5.
c
c inputs:
c
c       ndist1, ndist2  integer codes for distributions
c                       1=sev, 2=normal
c
c       ndist1 = distribution for loglife given fatigue limit
c       nidst2 = distribution for log(fatigue limit)
c
c       beta1     standardized slope = beta1/sigma
c
c       sdgamma   scale parameter of log(fatigue limit)
c
c       xi        standardized logstress = (x-ugamma)/sdgamma
c                 
c       z         standardized loglife
c                    = (w-beta0*)/sigma
c
c outputs:
c
c       answer             cdf value
c
c       iers, ieri         return condition indicator
c                          0 if no errors were
c                          detected in dqags, dqagi, resp.
c                          see dqags, dqagi documentation for 
c                          meaning of values of ier > 0
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer1/zp,xip,beta1p,sdgammap,ndist1p,
     +     ndist2p
      external rfl0int1
      dimension iwork(100),work(400)
c
c     constants that one might want to change to
c     achieve a higher degree of accuracy from the algorithm
c
      data eps/1.0d-30/
c
      data limit/100/
      data zero/0.0d00/
      iers=0
      ieri=0
c
c     variables for passing
c
      zp=z
      xip=xi
      beta1p=beta1
      sdgammap=sdgamma
      ndist1p=ndist1
      ndist2p=ndist2
c
c     do the integration (break into 2 intervals:
c     -infty - bd1, bd1 - zero)
c     
      answer1=zero
      answer2=zero
c
c     bounds for integration
c
c      valtol=1d-25
      bd1=-dabs(xi)
c      rlb=-3.d00
c 30   rlb=rlb-.1d00
c      rintval=rfl0int1(rlb)
c      if (rintval.gt.valtol) goto 30
      call dqagi(rfl0int1,bd1,-1,eps,eps,answer1,abserr,
     +     neval,ieri,limit,4*limit,last,iwork,work)
      call dqags(rfl0int1,bd1,zero,eps,eps,answer2,
     +     abserr,neval,iers,limit,4*limit,last,iwork,work)
      answer=answer1+answer2
      return
      end
	  
c----------------------------------------------------------------
      double precision function rfl0int1(x)
c----------------------------------------------------------------
c
c       function to compute integrand for computing the RFL 
c       (standardized) cdf
c
c       the following variables are communicated through common
c
c          beta1p     slope of the mean loglife
c          xip        standardized log stress
c          zp         standardized loglife
c          ndist1p    distn of W|V
c          ndist2p    distn of log(fatigue limit)
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer1/zp,xip,beta1p,sdgammap,ndist1p,
     +     ndist2p
      external stdcdf, stdpdf
      data zero/0.d00/
      if (x.eq.zero) then
         rfl0int1=zero
      else
         z=zp-beta1p*(xip*sdgammap+dlog(1.d00-dexp(sdgammap*x)))
         zx=x+xip
         rfl0int1=stdcdf(z,ndist1p)*stdpdf(zx,ndist2p)
      end if
      return
      end
	  
c----------------------------------------------------------------
      double precision function stdcdf(x,ndist)
c----------------------------------------------------------------
c 
c  compute standard SEV, normal, or logistic cdf
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      external fcdfn
      if (ndist.eq.1) then
c  sev distribution
         stdcdf=1.d00-dexp(-dexp(x))
      end if
      if (ndist.eq.2) then
c  normal distribution
         stdcdf=fcdfn(x)
      end if
      if (ndist.eq.3) then
c  logistic distribution
         stdcdf=dexp(x)/(1.d00+dexp(x))
      end if
      return
      end
	  
c----------------------------------------------------------------
      double precision function fcdfn(z)
c----------------------------------------------------------------
c
c        standard normal cdf
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      external derfc

      data half,root/.5d00,.7071067811865475d00/
      zroot=z*root
      fcdfn=half*derfc(-zroot)
      return
      end
	  
	  
c----------------------------------------------------------------
      subroutine rfl0quan(ndist1,ndist2,xi,alpha,beta1,sdgamma,
     +     bd1,bd2,answer)
c----------------------------------------------------------------
c
c       this is a subroutine to compute the standardized RFL
c       quantile 
c
c inputs:
c
c       ndist1, ndist2  integer codes for distributions
c                       1=sev, 2=normal
c
c       ndist1 = distribution for loglife given fatigue limit
c       nidst2 = distribution for log(fatigue limit)
c
c       beta1     standardized slope = beta1/sigma
c
c       sdgamma   log(fatigue limit) scale
c
c       xi        standardized stress = (x-ugamma)/sdgamma
c                 
c       alpha     proportion failing corresponding to quantile
c
c       b1, b2    bounds for quantile
c
c outputs:
c
c       answer    quantile value
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer2/ndist1p,ndist2p,xip,alphap,beta1p,
     +     sdgammap
      external zeroin, alphadiff, stdcdf
c
c     constants that one might want to change to
c     achieve a higher degree of accuracy from the algorithm
c
      data tol/1.0d-12/
      data top/2000.d00/
      data bottom/-2000.d00/
c
c     common stuff
c
      ndist1p=ndist1
      ndist2p=ndist2
      beta1p=beta1
      sdgammap=sdgamma
      xip=xi
      alphap=alpha
      bd1=bd1+.5d00
      bd2=bd2-.5d00
c
c     check if quantile is possible
c
      check=stdcdf(xi,ndist2)
      if (check.gt.alpha) goto 200
      print *, 'quantile is not possible'
      answer=bottom
      return
c
c     check if bounds are okay
c
 200  bd1=bd1-.5d00
      if (bd1.le.bottom) then
         answer=bottom
         return
      end if
      call rfl0cdf(ndist1,ndist2,bd1,xi,beta1,sdgamma,
     +     check1,ieri,iers)
      if (check1.gt.alpha) go to 200
 201  bd2=bd2+.5d00
      if (bd2.ge.top) then
         answer=top
         return      
      end if      
      call rfl0cdf(ndist1,ndist2,bd2,xi,beta1,sdgamma,
     +     check2,ieri,iers)
      if (check2.lt.alpha) go to 201
c
c     zero in on quantile value
c
      answer=zeroin(bd1,bd2,alphadiff,tol)
      return
      end
	  
c----------------------------------------------------------------
      double precision function alphadiff(x)
c----------------------------------------------------------------
c     
c     function to compute difference between cdf and alpha 
c
      implicit real*8(a-h,o-z)
      common/passer2/ndist1p,ndist2p,xip,alphap,beta1p,
     +     sdgammap
      call rfl0cdf(ndist1p,ndist2p,x,xip,beta1p,sdgammap,
     +     answer,ieri,iers)
      alphadiff=answer-alphap
      return
      end

c----------------------------------------------------------------
      subroutine srfl0stress(ndist1,ndist2,zalpha,alpha,beta1,
     +     sdgamma,bd1,bd2,num,answer)
c----------------------------------------------------------------
!MS$ATTRIBUTES DLLEXPORT,ALIAS:'srfl0stress'::srfl0stress
      implicit real*8(a-h,o-z)
      dimension ndist1(num),ndist2(num),zalpha(num),alpha(num),
     +     beta1(num),sdgamma(num),answer(num)
      do 10 i=1,num
         xbd1=bd1
         xbd2=bd2
         call rfl0stress(ndist1(i),ndist2(i),zalpha(i),alpha(i),
     +        beta1(i),sdgamma(i),xbd1,xbd2,answer(i))
 10   continue
      return
      end
	  
c----------------------------------------------------------------
      subroutine rfl0stress(ndist1,ndist2,zalpha,alpha,beta1,
     +     sdgamma,bd1,bd2,answer)
c----------------------------------------------------------------
c
c       this is a subroutine to compute RFL (version 2) 
c       stress to yield desired quantile
c
c inputs:
c
c       ndist1, ndist2  integer codes for distributions
c                       1=sev, 2=normal
c
c       ndist1 = distribution for log(lifetime) given fatigue
c                limit
c       nidst2 = distribution for log(fatigue limit)
c
c       beta1     standardized slope = beta1/sigma
c
c       sdgamma   log(fatigue limit) scale
c
c       zalpha    standardized quantile = (w-beta0*)/sigma
c                 
c       alpha     proportion corresponding to quantile
c
c       b1, b2    stress bounds
c
c outputs:
c
c       answer    centered stress value = (x-mugamma)
c
c----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/passer3/ndist1p,ndist2p,zalphap,alphap,beta1p,
     +     sdgammap
      external zeroin, xalphadiff, stdcdf
c
c     constants that one might want to change to
c     achieve a higher degree of accuracy from the algorithm
c
      data tol/1.0d-10/
c
c     common stuff
c
      top=10.d00
      bottom=-10.d00
      ndist1p=ndist1
      ndist2p=ndist2
      beta1p=beta1
      sdgammap=sdgamma
      zalphap=zalpha
      alphap=alpha
      bd1=bd1+.5d00
      bd2=bd2-.5d00
c
c     check if bounds are okay
c
 200  bd1=bd1-.5d00
      if (bd1.le.bottom) then
         answer=bottom
         return
      end if
      call rfl0cdf(ndist1,ndist2,zalpha,bd1,beta1,sdgamma,
     +     check1,ieri,iers)
      if (check1.gt.alpha) go to 200
 201  bd2=bd2+.5d00
      if (bd2.ge.top) then
         answer=top
         return      
      end if      
      call rfl0cdf(ndist1,ndist2,zalpha,bd2,beta1,sdgamma,
     +     check2,ieri,iers)
      if (check2.lt.alpha) go to 201
c
c     zero in on centered stress value
c
      answer=zeroin(bd1,bd2,xalphadiff,tol)
      return
      end
	  
c----------------------------------------------------------------
      double precision function xalphadiff(x)
c----------------------------------------------------------------
c     
c     function to compute difference between cdf and alpha 
c
      implicit real*8(a-h,o-z)
      common/passer3/ndist1p,ndist2p,zalphap,alphap,beta1p,
     +     sdgammap
      call rfl0cdf(ndist1p,ndist2p,zalphap,x,beta1p,sdgammap,
     +     answer,ieri,iers)
      xalphadiff=answer-alphap
      return
      end
	  