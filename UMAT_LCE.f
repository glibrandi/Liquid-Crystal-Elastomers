ccccccccccccccccccccccccccccccccccccccccccccccccccc
c Author Gabriele Librandi, PhD 
c Materials Science & Mechanical Engineering 
c Harvard University
c written in Oct 2018
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 MATERL


      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),DFGRD0(3,3),DFGRD1(3,3),
     3 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3)

C

C    LOCAL ARRAYS      

	DIMENSION RCauchy(3,3),RtauT1(3,3),RtauT(3,3),RTau(3,3),
     *RDummy(3,3),RF(3,3),RDF(3,3),RDF1(3,3),RDF2(3,3),RQref(3,3)

C
C ----------------------------------------------------------------
C    CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------

c     matrix properties
      RYoung=PROPS(1)
      RPoisson=PROPS(2)
      Rn1ref=PROPS(3)
      Rn2ref=PROPS(4)
      Rn3ref=PROPS(5)
      
	Rmu=RYoung/(2.d00*(RPoisson+1.d00))
	RBulk=(RYoung*RPoisson)/((1.d00+RPoisson)*(1.d00-2.d00*RPoisson))
c    Please note that RBulk is the Lame' constant lambda in this case; this because I was getting errors
	
	REps=1.d-8
	Rss=1.d00-TIME(1)

c	Rss = 1.d00-TIME(1)
C	Rss = (1./10.)-TIME(1)
C    Rss = 1./3.
c	Building matrix Q

	  RQref=0.d0
	RQref(1,1)=(3./2.)*(Rn1ref*Rn1ref-(1./3.))
	RQref(1,2)=(3./2.)*(Rn1ref*Rn2ref)
	RQref(1,3)=(3./2.)*(Rn1ref*Rn3ref)
	RQref(2,1)=(3./2.)*(Rn2ref*Rn1ref)
	RQref(2,2)=(3./2.)*(Rn2ref*Rn2ref-(1./3.))
	RQref(2,3)=(3./2.)*(Rn2ref*Rn3ref)
	RQref(3,1)=(3./2.)*(Rn3ref*Rn1ref)
	RQref(3,2)=(3./2.)*(Rn3ref*Rn2ref)
	RQref(3,3)=(3./2.)*(Rn3ref*Rn3ref-(1./3.))
	
c	Rn1=DFGRD1(1,1)*Rn1ref+DFGRD1(1,2)*Rn2ref+DFGRD1(1,3)*Rn3ref
c	Rn2=DFGRD1(2,1)*Rn1ref+DFGRD1(2,2)*Rn2ref+DFGRD1(2,3)*Rn3ref
c	Rn3=DFGRD1(3,1)*Rn1ref+DFGRD1(3,2)*Rn2ref+DFGRD1(3,3)*Rn3ref
	
c	  RQ=0.d0
c	RQ(1,1)=(3./2.)*(Rn1*Rn1-(1./3.))
c	RQ(1,2)=(3./2.)*(Rn1*Rn2)
c	RQ(1,3)=(3./2.)*(Rn1*Rn3)
c	RQ(2,1)=(3./2.)*(Rn2*Rn1)
c	RQ(2,2)=(3./2.)*(Rn2*Rn2-(1./3.))
c	RQ(2,3)=(3./2.)*(Rn2*Rn3)
c	RQ(3,1)=(3./2.)*(Rn3*Rn1)
c	RQ(3,2)=(3./2.)*(Rn3*Rn2)
c	RQ(3,3)=(3./2.)*(Rn3*Rn3-(1./3.))

c	WRITE(6,*)'TIME_1',TIME(1)
c	WRITE(6,*)'alpha',Ral

c	CALL  XIT

      CALL KStressNH(DFGRD1,Rmu,RBulk,RQref,Rss,RCauchy)

c	WRITE(6,*)'RCauchy',RCauchy

	STRESS(1)=RCauchy(1,1)
	STRESS(2)=RCauchy(2,2)
	STRESS(3)=RCauchy(3,3)
	STRESS(4)=RCauchy(1,2)
	STRESS(5)=RCauchy(1,3)
	STRESS(6)=RCauchy(2,3)

c ------------------------------------------------------------------
c     Calculation of the Material Jacobian also called Tangent Modulus Tensor
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c-- JACOBIAN --
c-------------------------------------------------------------------
	CALL KDETERM(DFGRD1,RJ) 

	RTau=RJ*(RCauchy)

C We compute from Cauchy Stress the Kirchhoff one
C RJ is the determinant of F, if the material is incompressible it should be equal to 1	

	DO Ii2=1,6

	   IF (Ii2.EQ.1)THEN
	     Ii=1
	     Ij=1
	   ELSEIF (Ii2.EQ.2)THEN
	     Ii=2
	     Ij=2
	   ELSEIF (Ii2.EQ.3)THEN
	     Ii=3
	     Ij=3
	   ELSEIF (Ii2.EQ.4)THEN
	     Ii=1
	     Ij=2
	   ELSEIF (Ii2.EQ.5)THEN
	     Ii=1
	     Ij=3
	   ELSEIF (Ii2.EQ.6)THEN
	     Ii=2
	     Ij=3
	   ENDIF


         RDummy=0.d00
	   RDummy(Ii,Ij)=1.d00
	   CALL KMULT(RDummy,DFGRD1,RDF1) 
	   RDummy=0.d00
	   RDummy(Ij,Ii)=1.d00
	   CALL KMULT(RDummy,DFGRD1,RDF2)     
         RDF=REps/2.d00*(RDF1+RDF2)
	   RF=DFGRD1+RDF

	   CALL KStressNH(RF,Rmu,RBulk,RQref,Rss,RtauT1)
	   RtauT=RJ*(RtauT1)

C here we produce Kirchhoff stress


         DDSDDE(Ii2,1)=1.d0/(REps*RJ)*
     *                  (RTauT(1,1)-RTau(1,1))	
	   DDSDDE(Ii2,2)=1.d0/(REps*RJ)*
     *                  (RTauT(2,2)-RTau(2,2))
	   DDSDDE(Ii2,3)=1.d0/(REps*RJ)*
     *                  (RTauT(3,3)-RTau(3,3))
	   DDSDDE(Ii2,4)=1.d0/(REps*RJ)*
     *                  (RTauT(1,2)-RTau(1,2))
	   DDSDDE(Ii2,5)=1.d0/(REps*RJ)*
     *                  (RTauT(1,3)-RTau(1,3))
	   DDSDDE(Ii2,6)=1.d0/(REps*RJ)*
     *                  (RTauT(2,3)-RTau(2,3))

	  


	ENDDO


c  	WRITE(6,*)'alpha',Ralpha/RJ



c	CALL XIT


  
	RETURN
	END

c=====/================================================================/    




!*********************************************************
!-----------------------------------------------------------------------
      SUBROUTINE condenseETens(cc4, NDI, NSHR, NTENS, DDSDDE)
!     reorganize 4th order elasticity tensor into matrix

      INCLUDE 'ABA_PARAM.INC'                        
      DIMENSION cc4(3,3,3,3), DDSDDE(NTENS,NTENS), 
     & idx1ik(6), idx1jl(6), idx2ik(6), idx2jl(6)
      
      DATA idx1ik / 1,2,3,1,1,2 /, 
     &     idx1jl / 1,2,3,2,3,3 /,
     &     idx2ik / 1,2,3,2,3,3 /,
     &     idx2jl / 1,2,3,1,1,2 /
      
      
      do I = 1, NTENS
        do J = 1, NTENS
           DDSDDE(I,J) = ( cc4(idx1ik(I),idx1jl(I),idx1ik(J),idx1jl(J))
     &               + cc4(idx2ik(I),idx2jl(I),idx2ik(J),idx2jl(J)) )
     &             * 0.5d0
        end do
      end do
         
      RETURN
      END 

c----------------------------------------------------------------------
c-- STRESSES AS A.H. Gelebart et al. Nature 2017  --
c---------------------------------------------------------------------

      SUBROUTINE KStressNH(RF,Rmu,RBulk,RQref,Rss,RC)

	INCLUDE 'ABA_PARAM.INC'  
      
	DIMENSION RS(3,3),RF(3,3),RC(3,3),RQ(3,3)
	DIMENSION RFTr(3,3),RB(3,3),RQref(3,3)

c     J=det(F)
	CALL KDETERM(RF,RJ)

c     B=F(F)^T  
      CALL KTRANS(RF,RFTr) 
      CALL KMULT(RF,RFTr,RB)
      

c-------------------------------------------------------------------
c-- STRESSES AS A.H. Gelebart et al. Nature 2017  --
c-------------------------------------------------------------------
    
c      Ral=(1.44/43.)
      Ral=(200./10000.)
	  Rs0=1.d00 
	  RS=0.d0
  
		RS(1,1)=RF(1,1)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*RBulk+RF(1,1)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*RBulk
		RS(1,1)=2.d00*RF(1,2)*(RF(1,1)*RF(1,2)+RF(2,1)*RF(2,2)+RF(3,1)*RF(3,2))*Rmu+RS(1,1)  
		RS(1,1)=2.d00*RF(1,3)*(RF(1,1)*RF(1,3)+RF(2,1)*RF(2,3)+RF(3,1)*RF(3,3))*Rmu+RS(1,1)
		RS(1,1)=RF(1,1)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*(RBulk+2.d00*Rmu)+RS(1,1)
c	
		RS(1,1)=Ral*(-9.d00*RF(1,1)*RQref(1,1)*(Rss-Rs0))+(1./2.)*RS(1,1)
		RS(1,1)=Ral*((-9./2.)*RF(1,2)*RQref(1,2)*(Rss-Rs0))+RS(1,1)
		RS(1,1)=Ral*((-9./2.)*RF(1,3)*RQref(1,3)*(Rss-Rs0))+RS(1,1)
		RS(1,1)=Ral*((-9./2.)*RF(1,2)*RQref(2,1)*(Rss-Rs0))+RS(1,1)
		RS(1,1)=Ral*((-9./2.)*RF(1,3)*RQref(3,1)*(Rss-Rs0))+RS(1,1)
		RS(1,1)=RS(1,1)/RJ
		
			 
		RS(1,2)=RF(1,2)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*RBulk+RF(1,2)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*RBulk  
		RS(1,2)=2.d00*RF(1,1)*(RF(1,1)*RF(1,2)+RF(2,1)*RF(2,2)+RF(3,1)*RF(3,2))*Rmu+RS(1,2)  
		RS(1,2)=2.d00*RF(1,3)*(RF(1,2)*RF(1,3)+RF(2,2)*RF(2,3)+RF(3,2)*RF(3,3))*Rmu+RS(1,2)
		RS(1,2)=RF(1,2)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*(RBulk+2.d00*Rmu)+RS(1,2)
c		
		RS(1,2)=Ral*(-9.d00*RF(1,2)*RQref(2,2)*(Rss-Rs0))+(1./2.)*RS(1,2)
		RS(1,2)=Ral*((-9./2.)*RF(1,1)*RQref(1,2)*(Rss-Rs0))+RS(1,2)
		RS(1,2)=Ral*((-9./2.)*RF(1,1)*RQref(2,1)*(Rss-Rs0))+RS(1,2)
		RS(1,2)=Ral*((-9./2.)*RF(1,3)*RQref(2,3)*(Rss-Rs0))+RS(1,2)
		RS(1,2)=Ral*((-9./2.)*RF(1,3)*RQref(3,2)*(Rss-Rs0))+RS(1,2)
		RS(1,2)=RS(1,2)/RJ
		
		
		RS(1,3)=RF(1,3)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*RBulk+RF(1,3)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*RBulk
		RS(1,3)=2.d00*RF(1,1)*(RF(1,1)*RF(1,3)+RF(2,1)*RF(2,3)+RF(3,1)*RF(3,3))*Rmu+RS(1,3)  
		RS(1,3)=2.d00*RF(1,2)*(RF(1,2)*RF(1,3)+RF(2,2)*RF(2,3)+RF(3,2)*RF(3,3))*Rmu+RS(1,3)
		RS(1,3)=RF(1,3)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*(RBulk+2.d00*Rmu)+RS(1,3)
c
		RS(1,3)=Ral*(-9.d00*RF(1,3)*RQref(3,3)*(Rss-Rs0))+(1./2.)*RS(1,3)
		RS(1,3)=Ral*((-9./2.)*RF(1,1)*RQref(1,3)*(Rss-Rs0))+RS(1,3)
		RS(1,3)=Ral*((-9./2.)*RF(1,2)*RQref(2,3)*(Rss-Rs0))+RS(1,3)
		RS(1,3)=Ral*((-9./2.)*RF(1,1)*RQref(3,1)*(Rss-Rs0))+RS(1,3)
		RS(1,3)=Ral*((-9./2.)*RF(1,2)*RQref(3,2)*(Rss-Rs0))+RS(1,3)
		RS(1,3)=RS(1,3)/RJ
		
		
		RS(2,1)=RF(2,1)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*RBulk+RF(2,1)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*RBulk
		RS(2,1)=2.d00*RF(2,2)*(RF(1,1)*RF(1,2)+RF(2,1)*RF(2,2)+RF(3,1)*RF(3,2))*Rmu+RS(2,1)  
		RS(2,1)=2.d00*RF(2,3)*(RF(1,1)*RF(1,3)+RF(2,1)*RF(2,3)+RF(3,1)*RF(3,3))*Rmu+RS(2,1)
		RS(2,1)=RF(2,1)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*(RBulk+2.d00*Rmu)+RS(2,1)
c
		RS(2,1)=Ral*(-9.d00*RF(2,1)*RQref(1,1)*(Rss-Rs0))+(1./2.)*RS(2,1)
		RS(2,1)=Ral*((-9./2.)*RF(2,2)*RQref(1,2)*(Rss-Rs0))+RS(2,1)
		RS(2,1)=Ral*((-9./2.)*RF(2,3)*RQref(1,3)*(Rss-Rs0))+RS(2,1)
		RS(2,1)=Ral*((-9./2.)*RF(2,2)*RQref(2,1)*(Rss-Rs0))+RS(2,1)
		RS(2,1)=Ral*((-9./2.)*RF(2,3)*RQref(3,1)*(Rss-Rs0))+RS(2,1)
		RS(2,1)=RS(2,1)/RJ
		   
		   
		RS(2,2)=RF(2,2)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*RBulk+RF(2,2)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*RBulk
		RS(2,2)=2.d00*RF(2,1)*(RF(1,1)*RF(1,2)+RF(2,1)*RF(2,2)+RF(3,1)*RF(3,2))*Rmu+RS(2,2)  
		RS(2,2)=2.d00*RF(2,3)*(RF(1,2)*RF(1,3)+RF(2,2)*RF(2,3)+RF(3,2)*RF(3,3))*Rmu+RS(2,2)
		RS(2,2)=RF(2,2)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*(RBulk+2.d00*Rmu)+RS(2,2)
c
		RS(2,2)=Ral*(-9.d00*RF(2,2)*RQref(2,2)*(Rss-Rs0))+(1./2.)*RS(2,2)
		RS(2,2)=Ral*((-9./2.)*RF(2,1)*RQref(1,2)*(Rss-Rs0))+RS(2,2)
		RS(2,2)=Ral*((-9./2.)*RF(2,1)*RQref(2,1)*(Rss-Rs0))+RS(2,2)
		RS(2,2)=Ral*((-9./2.)*RF(2,3)*RQref(2,3)*(Rss-Rs0))+RS(2,2)
		RS(2,2)=Ral*((-9./2.)*RF(2,3)*RQref(3,2)*(Rss-Rs0))+RS(2,2)
		RS(2,2)=RS(2,2)/RJ
		
		
		RS(2,3)=RF(2,3)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*RBulk+RF(2,3)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*RBulk 
		RS(2,3)=2.d00*RF(2,1)*(RF(1,1)*RF(1,3)+RF(2,1)*RF(2,3)+RF(3,1)*RF(3,3))*Rmu+RS(2,3)  
		RS(2,3)=2.d00*RF(2,2)*(RF(1,2)*RF(1,3)+RF(2,2)*RF(2,3)+RF(3,2)*RF(3,3))*Rmu+RS(2,3)
		RS(2,3)=RF(2,3)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*(RBulk+2.d00*Rmu)+RS(2,3)
c
		RS(2,3)=Ral*(-9.d00*RF(2,3)*RQref(3,3)*(Rss-Rs0))+(1./2.)*RS(2,3)
		RS(2,3)=Ral*((-9./2.)*RF(2,1)*RQref(1,3)*(Rss-Rs0))+RS(2,3)
		RS(2,3)=Ral*((-9./2.)*RF(2,2)*RQref(2,3)*(Rss-Rs0))+RS(2,3)
		RS(2,3)=Ral*((-9./2.)*RF(2,1)*RQref(3,1)*(Rss-Rs0))+RS(2,3)
		RS(2,3)=Ral*((-9./2.)*RF(2,2)*RQref(3,2)*(Rss-Rs0))+RS(2,3)
		RS(2,3)=RS(2,3)/RJ
		
		
		RS(3,1)=RF(3,1)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*RBulk+RF(3,1)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*RBulk 
		RS(3,1)=2.d00*RF(3,2)*(RF(1,1)*RF(1,2)+RF(2,1)*RF(2,2)+RF(3,1)*RF(3,2))*Rmu+RS(3,1)  
		RS(3,1)=2.d00*RF(3,3)*(RF(1,1)*RF(1,3)+RF(2,1)*RF(2,3)+RF(3,1)*RF(3,3))*Rmu+RS(3,1)
		RS(3,1)=RF(3,1)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*(RBulk+2.d00*Rmu)+RS(3,1)
c
		RS(3,1)=Ral*(-9.d00*RF(3,1)*RQref(1,1)*(Rss-Rs0))+(1./2.)*RS(3,1)
		RS(3,1)=Ral*((-9./2.)*RF(3,2)*RQref(1,2)*(Rss-Rs0))+RS(3,1)
		RS(3,1)=Ral*((-9./2.)*RF(3,3)*RQref(1,3)*(Rss-Rs0))+RS(3,1)
		RS(3,1)=Ral*((-9./2.)*RF(3,2)*RQref(2,1)*(Rss-Rs0))+RS(3,1)
		RS(3,1)=Ral*((-9./2.)*RF(3,3)*RQref(3,1)*(Rss-Rs0))+RS(3,1)
		RS(3,1)=RS(3,1)/RJ
		  
		  
		RS(3,2)=RF(3,2)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*RBulk+RF(3,2)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*RBulk
		RS(3,2)=2.d00*RF(3,1)*(RF(1,1)*RF(1,2)+RF(2,1)*RF(2,2)+RF(3,1)*RF(3,2))*Rmu+RS(3,2)  
		RS(3,2)=2.d00*RF(3,3)*(RF(1,2)*RF(1,3)+RF(2,2)*RF(2,3)+RF(3,2)*RF(3,3))*Rmu+RS(3,2)
		RS(3,2)=RF(3,2)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*(RBulk+2.d00*Rmu)+RS(3,2)
c
		RS(3,2)=Ral*(-9.d00*RF(3,2)*RQref(2,2)*(Rss-Rs0))+(1./2.)*RS(3,2)
		RS(3,2)=Ral*((-9./2.)*RF(3,1)*RQref(1,2)*(Rss-Rs0))+RS(3,2)
		RS(3,2)=Ral*((-9./2.)*RF(3,1)*RQref(2,1)*(Rss-Rs0))+RS(3,2)
		RS(3,2)=Ral*((-9./2.)*RF(3,3)*RQref(2,3)*(Rss-Rs0))+RS(3,2)
		RS(3,2)=Ral*((-9./2.)*RF(3,3)*RQref(3,2)*(Rss-Rs0))+RS(3,2)
		RS(3,2)=RS(3,2)/RJ
	
	
		RS(3,3)=RF(3,3)*(-1.d00+RF(1,1)**2+RF(2,1)**2+RF(3,1)**2)*RBulk+RF(3,3)*(-1.d00+RF(1,2)**2+RF(2,2)**2+RF(3,2)**2)*RBulk 
		RS(3,3)=2.d00*RF(3,1)*(RF(1,1)*RF(1,3)+RF(2,1)*RF(2,3)+RF(3,1)*RF(3,3))*Rmu+RS(3,3)  
		RS(3,3)=2.d00*RF(3,2)*(RF(1,2)*RF(1,3)+RF(2,2)*RF(2,3)+RF(3,2)*RF(3,3))*Rmu+RS(3,3)
		RS(3,3)=RF(3,3)*(-1.d00+RF(1,3)**2+RF(2,3)**2+RF(3,3)**2)*(RBulk+2.d00*Rmu)+RS(3,3)
c
		RS(3,3)=Ral*(-9.d00*RF(3,3)*RQref(3,3)*(Rss-Rs0))+(1./2.)*RS(3,3)
		RS(3,3)=Ral*((-9./2.)*RF(3,1)*RQref(1,3)*(Rss-Rs0))+RS(3,3)
		RS(3,3)=Ral*((-9./2.)*RF(3,2)*RQref(2,3)*(Rss-Rs0))+RS(3,3)
		RS(3,3)=Ral*((-9./2.)*RF(3,1)*RQref(3,1)*(Rss-Rs0))+RS(3,3)
		RS(3,3)=Ral*((-9./2.)*RF(3,2)*RQref(3,2)*(Rss-Rs0))+RS(3,3)
		RS(3,3)=RS(3,3)/RJ

      CALL KMULT(RS,RFTr,RC)
      
	RETURN
	END



c=====/================================================================/    
C                   
      SUBROUTINE KDETERM(RA,RdetA)                                               
C                   
C     KDETERM(A,d): CALCULATES THE DETERMINANT, d, OF MATRIX A(3,3)         
C                   
C                   
C     variables used:                                                       
C     A = INPUT tensor                                                
C     d = scalar value of determinant                                 
C                   
C     SUBROUTINE CALLS: - none                                              
C                   
C-----/----------------------------------------------------------------/    
      INCLUDE 'ABA_PARAM.INC'                                               
      DIMENSION RA(3,3)                                                      
                    
      RdetA=  RA(1,1)*(RA(2,2)*RA(3,3) - RA(2,3)*RA(3,2))                           
     1    -RA(1,2)*(RA(2,1)*RA(3,3) - RA(2,3)*RA(3,1))                           
     2    +RA(1,3)*(RA(2,1)*RA(3,2) - RA(2,2)*RA(3,1))  
                              
      RETURN        
      END 

c=====/================================================================/    
      SUBROUTINE KMULT(RA,RB,RC)                                                                   
C                   
C     variables used:                                                       
C     RC=RA*RB                                               
C-----/----------------------------------------------------------------/    
      INCLUDE 'ABA_PARAM.INC'                                               
      DIMENSION RA(3,3),RB(3,3),RC(3,3)


	                                                                                                    
        RC(1,1)=RA(1,1)*RB(1,1)+RA(1,2)*RB(2,1)+RA(1,3)*RB(3,1)                        
        RC(1,2)=RA(1,1)*RB(1,2)+RA(1,2)*RB(2,2)+RA(1,3)*RB(3,2)                         
        RC(1,3)=RA(1,1)*RB(1,3)+RA(1,2)*RB(2,3)+RA(1,3)*RB(3,3)
	                              
        RC(2,1)=RA(2,1)*RB(1,1)+RA(2,2)*RB(2,1)+RA(2,3)*RB(3,1)                          
        RC(2,2)=RA(2,1)*RB(1,2)+RA(2,2)*RB(2,2)+RA(2,3)*RB(3,2)                            
        RC(2,3)=RA(2,1)*RB(1,3)+RA(2,2)*RB(2,3)+RA(2,3)*RB(3,3)
	                            
        RC(3,1)=RA(3,1)*RB(1,1)+RA(3,2)*RB(2,1)+RA(3,3)*RB(3,1)                           
        RC(3,2)=RA(3,1)*RB(1,2)+RA(3,2)*RB(2,2)+RA(3,3)*RB(3,2)                          
        RC(3,3)=RA(3,1)*RB(1,3)+RA(3,2)*RB(2,3)+RA(3,3)*RB(3,3) 
                        
C                   
      RETURN        
      END 
c=====/================================================================/    
C                   
      SUBROUTINE KTRANS(RA,RAT)                                               
C                   
C     RAT=transpose(RA)                                                                      
C                   
C-----/----------------------------------------------------------------/    
      INCLUDE 'ABA_PARAM.INC'                                               
      DIMENSION RA(3,3),RAT(3,3)                                                      
                    
      RAT(1,1)=RA(1,1)
	RAT(1,2)=RA(2,1)
	RAT(1,3)=RA(3,1)
	RAT(2,1)=RA(1,2)
	RAT(2,2)=RA(2,2)
	RAT(2,3)=RA(3,2)
	RAT(3,1)=RA(1,3)
	RAT(3,2)=RA(2,3)
	RAT(3,3)=RA(3,3)
                          
      RETURN        
      END  

c=====/================================================================/
c=====/================================================================/

c=====/================================================================/    
      SUBROUTINE KINVERT(RA,RB)                                          
C                   
C     KINVERT(A,detA,B): CALCULATES THE INVERSE OF MATRIX A(3,3)            
C                        WHICH HAS DETERMINANT detA (INPUT PARAMETER)       
C 			           AND STORES THE RESULT IN MATRIX B = A^(-1)                            
C                   
C     variables used:                                                       
C     A    = INPUT tensor                                                
C     detA = INPUT scalar, determinant (A)                               
C     B    = OUTPUT :inverse of input tensor                             
C                   
C     SUBROUTINE CALLS: - none                                              
C                   
C-----/----------------------------------------------------------------/    
      INCLUDE 'ABA_PARAM.INC'                                               
      DIMENSION RA(3,3),RB(3,3)
	
	  CALL KDETERM(RA,RdetA)                                               
                    
        Rdi=1.0/RdetA                                                         
        RB(1,1)=(RA(2,2)*RA(3,3)-RA(3,2)*RA(2,3))*Rdi                             
        RB(2,2)=(RA(1,1)*RA(3,3)-RA(3,1)*RA(1,3))*Rdi                             
        RB(3,3)=(RA(1,1)*RA(2,2)-RA(2,1)*RA(1,2))*Rdi                             
        RB(1,2)=(RA(3,2)*RA(1,3)-RA(1,2)*RA(3,3))*Rdi                             
        RB(1,3)=(RA(1,2)*RA(2,3)-RA(2,2)*RA(1,3))*Rdi                             
        RB(2,3)=(RA(2,1)*RA(1,3)-RA(1,1)*RA(2,3))*Rdi                             
        RB(2,1)=(RA(3,1)*RA(2,3)-RA(2,1)*RA(3,3))*Rdi                             
        RB(3,1)=(RA(2,1)*RA(3,2)-RA(3,1)*RA(2,2))*Rdi                             
        RB(3,2)=(RA(3,1)*RA(1,2)-RA(1,1)*RA(3,2))*Rdi                             
C                   
      RETURN        
      END 
