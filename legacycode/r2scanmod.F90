!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module constains the active subroutines
!     r2scaninit    r2scanfun
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module r2scanmod
  Use globalmath      
  IMPLICIT NONE

  REAL(8), PRIVATE, PARAMETER :: cx0 = 1.d0
  REAL(8), PRIVATE, PARAMETER :: cx1 = -0.667d0
  REAL(8), PRIVATE, PARAMETER :: cx2 = -0.4445555d0
  REAL(8), PRIVATE, PARAMETER :: cx3 = -0.663086601049d0
  REAL(8), PRIVATE, PARAMETER :: cx4 = 1.451297044490d0
  REAL(8), PRIVATE, PARAMETER :: cx5 = -0.887998041597d0
  REAL(8), PRIVATE, PARAMETER :: cx6 = 0.234528941479d0
  REAL(8), PRIVATE, PARAMETER :: cx7 = -0.023185843322d0
  REAL(8), PRIVATE, PARAMETER :: SCANc1x = 0.667d0
  REAL(8), PRIVATE, PARAMETER :: SCANc2x = 0.8d0
  REAL(8), PRIVATE, PARAMETER :: SCANdx = 1.24d0
  REAL(8), PRIVATE, PARAMETER :: k0 = 0.174d0
  REAL(8), PRIVATE, PARAMETER :: k1 = 0.065d0
  REAL(8), PRIVATE, PARAMETER :: mu = 10.d0/81.d0
  REAL(8), PRIVATE, PARAMETER :: SCANa1 = 4.9479d0
!  REAL(8), PRIVATE, PARAMETER :: eta = 0.001d0
  REAL(8) :: eta
  REAL(8), PRIVATE, PARAMETER :: dp2 = 0.361d0
!  REAL(8), PRIVATE, PARAMETER :: C2Ceta = (20.d0/27.d0 + (5*eta)/3.d0)*(-0.162742d0)
  REAL(8) :: C2Ceta
  REAL(8), PRIVATE, PARAMETER :: cc0 = 1.d0
  REAL(8), PRIVATE, PARAMETER :: cc1 = -0.64d0
  REAL(8), PRIVATE, PARAMETER :: cc2 = -0.4352d0
  REAL(8), PRIVATE, PARAMETER :: cc3 = -1.535685604549d0
  REAL(8), PRIVATE, PARAMETER :: cc4 = 3.061560252175d0
  REAL(8), PRIVATE, PARAMETER :: cc5 = -1.915710236206d0
  REAL(8), PRIVATE, PARAMETER :: cc6 = 0.516884468372d0
  REAL(8), PRIVATE, PARAMETER :: cc7 = -0.051848879792d0
  REAL(8), PRIVATE, PARAMETER :: SCANc1c = 0.64d0
  REAL(8), PRIVATE, PARAMETER :: SCANc2c = 1.5d0
  REAL(8), PRIVATE, PARAMETER :: SCANdc = 0.7d0
  REAL(8), PRIVATE, PARAMETER :: b1c = 0.0285764d0
  REAL(8), PRIVATE, PARAMETER :: b2c = 0.0889d0
  REAL(8), PRIVATE, PARAMETER :: b3c = 0.125541d0
  REAL(8), PRIVATE, PARAMETER :: betaMB = 0.066725d0
!  REAL(8), PRIVATE, PARAMETER :: chiinfinity = 0.128025d0 !from Tulane code
  REAL(8), PRIVATE, PARAMETER :: chiinfinity = 0.12802585262625815d0
!  REAL(8), PRIVATE, PARAMETER :: Sgam = 0.031091d0   !from Tulane code
  REAL(8), PRIVATE, PARAMETER :: Sgam = 0.031090690869655d0
!  REAL(8), PRIVATE, PARAMETER :: Dfc2 = -0.711402d0   !from Tulane code
  REAL(8), PRIVATE, PARAMETER :: Dfc2=cc1+2*cc2+3*cc3+4*cc4+5*cc5+6*cc6+7*cc7
! Parameters for the Perdew-Wang (PRB 45,13244 (1992)) LDA
! correlation
  REAL(8), PRIVATE, PARAMETER :: LDAA = 0.03109070d0
  REAL(8), PRIVATE, PARAMETER :: LDAa1 = 0.21370d0
  REAL(8), PRIVATE, PARAMETER :: LDAb1 = 7.59570d0
  REAL(8), PRIVATE, PARAMETER :: LDAb2 = 3.58760d0
  REAL(8), PRIVATE, PARAMETER :: LDAb3 = 1.63820d0
  REAL(8), PRIVATE, PARAMETER :: LDAb4 = 0.49294d0


Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  r2scaninit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE r2scaninit(etain)
   implicit none
   
   REAL(8), INTENT(IN) :: etain
   REAL(8) :: o

   write(std_out,*) 'r2scan calculation with eta ', etain
   eta=etain   
   o=cx1+2*cx2+3*cx3+4*cx4+5*cx5+6*cx6+7*cx7
   write(std_out,*) 'r2scan check C2' , o*k0
   C2Ceta = (20.d0/27.d0 + (5*eta)/3.d0)*o*k0
!!!!!  C2Ceta = (20.d0/27.d0 + (5*eta)/3.d0)*(-0.162742d0)
 END SUBROUTINE r2scaninit
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! r2scan functional
!  James W. Furness,  Aaron D. Kaplan, Jinliang Ning, John P.
!  Perdew, and Jianwei Sun J. Phys. Chem. Lett. 2020, 11, 8208−8215
!
!  length units == Bohr
!  energy units -- Hartree
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE r2scanfun(rho,grad,tau,exc,vtau,vxcn,vxcs)
IMPLICIT NONE
  REAL(8), INTENT(IN) :: rho,grad, tau
  REAL(8), INTENT(INOUT) :: exc,vtau,vxcn,vxcs

  REAL(8) :: rr,ss,tt,sigma
  REAL(8) :: kF,ks,rs,U,s,p,t,W,baralpha
  REAL(8) :: gx,x,h0x,h1x,fx,ffx,exarg,ddfx,vxtau,ex,vx
  REAL(8) :: ELDA,ELDA0,betc,fc,ddfc,w1,w0,ginfinity,H0c,H1c,ec0,ec1,y,yy
  REAL(8) :: beta,srs,ddLDA,difddLDA,gg,corarg,ddy,vctau
  REAL(8) :: alphadenom,dalphadn,dalphads
  REAL(8) :: gxp,dpdn,dgxdn,dpds,dgxds,drsdn,dt2ds,dt2dn,t2
  REAL(8) :: h1xx,xp,dh1xdn,dh1xds,vxn,vxs
  REAL(8) :: vcn,vcs,dLDA0dn,dH0cdn,dw0dn,dgindn,dw1dn
  REAL(8) :: dH0cds,dginds,dbetadn,dif2ddLDA,dLDAdn
  REAL(8) :: pp,dppds,dppdn,dyds,dydn,dUdn,dUds,dWds,dWdn
  REAL(8) :: term1,term4,dddyds,dddydnp,dddydnw,dddydnr,dH1cds
  REAL(8) :: stf1,stf2,stf3,stf4,stf5,stf6,stf7
  REAL(8) :: dddydnn,dddydn,dH1cdwn,dH1cdn,dec1dn,dec1ds
  REAL(8) :: dec0ds,dec0dn,bigdenom
  REAL(8) :: eex,eec,vtx,vtc,vvxn,vvcn,vvxs,vvcs

  eex=0.d0;eec=0.d0;vtx=0.d0;vtc=0.d0;vvxn=0.d0;vvcn=0.d0;vvxs=0.d0;vvcs=0.d0
  if(rho<(machine_zero**0.333333333333333333333333333333333d0)) return
  sigma=grad*grad
  rr=rho;ss=sigma;tt=tau
  rr=max(machine_zero,rr)
  rr=min(machine_infinity,rr)
  ss=max(machine_zero,ss)
  ss=min(machine_infinity,ss)
  tt=max(machine_zero,tt)
  tt=min(machine_infinity,tt)

!  Note that all derivatives with respect to sigma are also multiplied
!      by grad

  kF=3.0936677262801359310d0*(rr**0.33333333333333333333333333333d0)
  ks=1.9846863952198559283d0*(rr**0.16666666666666666666666666667d0)
  rs=0.62035049089940001665d0/(rr**0.33333333333333333333333333333d0)
  drsdn=-rs/(3*rr)
  srs=sqrt(rs)
  U=2.8712340001881918160d0*(rr**1.66666666666666666666666666667d0)
  dUdn=(1.666666666666666666666667d0)*2.8712340001881918160d0*(rr**0.66666666666666666666666666667d0)
  s=sqrt(ss)*0.16162045967399548133d0/(rr**1.33333333333333333333333333333d0)
  p=s*s
  dpdn=-2.666666666666666666666666667d0*p/rr
  dpds=0.026121172985233599568d0*grad/(rr**2.666666666666666666666666666666667d0)
  t=sqrt(ss)*0.25192897034224488769d0/(rr**1.166666666666666666666666666667d0)
  t2=t*t
  dt2ds=0.063468206097703704205d0*grad/(rr**2.3333333333333333333333333333d0)
  dt2dn=-t2*(2.333333333333333333333333333333333333d0)/rr
  W=0.125d0*ss/rr
  dWds=0.125d0*grad/rr
  dWdn=-W/rr
  baralpha=(tt-W)/(U+eta*W)
!!!!!!exchange part

  gx=1.d0-ddexp(-SCANa1/(p**0.25d0))
  x=(c2ceta*ddexp(-(p**2)/dp2**4)+mu)*p
  xp=c2ceta*ddexp(-(p**2)/dp2**4)*(-2*(p**2)/dp2**4+1.d0)+mu
  h0x=1+k0
  h1x=1.d0+k1-k1/(1+x/k1)
  if(baralpha<1.d-13) then
     fx=ddexp(-SCANc1x*baralpha/(1.d0-baralpha))
  elseif (baralpha.lt.2.5d0) then
     fx=cx0+baralpha*(cx1+baralpha*(cx2+baralpha*(cx3+baralpha*( &
&        cx4+baralpha*(cx5+baralpha*(cx6+baralpha*cx7))))))
  else if (baralpha.ge.2.5d0) then
     fx=-SCANdx*ddexp(SCANc2x/(1.d0-baralpha))
  endif
  if(baralpha<1.d-13) then
     ddfx=-(SCANc1x/((1.d0-baralpha)**2))*ddexp(-SCANc1x*baralpha/(1.d0-baralpha))
  else if (baralpha.lt.2.5d0) then
     ddfx=(cx1+baralpha*(2*cx2+baralpha*(3*cx3+baralpha*( &
&        4*cx4+baralpha*(5*cx5+baralpha*(6*cx6+baralpha*7*cx7))))))
  else if (baralpha.ge.2.5d0) then
     ddfx=-(SCANdx*SCANc2x/((1.d0-baralpha)**2))*ddexp(SCANc2x/(1.d0-baralpha))
  endif
  ffx=(h1x+fx*(h0x-h1x))*gx  
  ex=-0.73855876638202240587d0*(rr**1.3333333333333333333333333333333333d0)
  vx=-0.98474502184269654116d0*(rr**0.3333333333333333333333333333333333d0)
  exarg=ex*ffx
  vxtau=ex*ddfx*((h0x-h1x)*gx)/(U+eta*W)          
! density and sigma derivative terms
  alphadenom=(U+eta*W)**2
  dalphads=-(U+eta*tau)*dWds/alphadenom
  dalphadn=-((U+eta*tau)*dWdn+(tt-W)*dUdn)/alphadenom
  gxp=-0.25d0*SCANa1*ddexp(-SCANa1/(p**0.25d0))/(p**1.25d0)
  dgxdn=gxp*dpdn
  dgxds=gxp*dpds
  h1xx=1.d0/((1.d0+x/k1)**2) 
  dh1xdn=h1xx*xp*dpdn
  dh1xds=h1xx*xp*dpds
  vxn=vx*ffx+ex*(ddfx*dalphadn*gx*(h0x-h1x)+dgxdn*(h1x+fx*(h0x-h1x)) &
&      +gx*dh1xdn*(1.d0-fx))          
  vxs=ex*(ddfx*dalphads*gx*(h0x-h1x)+dgxds*(h1x+fx*(h0x-h1x)) &
&      +gx*dh1xds*(1.d0-fx))          

!
!!!! correlation part

  ELDA=-2*LDAA*(1.d0 + LDAa1*rs)*ddlog(1.d0 + 0.5d0 &
&   /(LDAA*(sqrt(rs)*(LDAb1 + LDAb3*rs) + rs*(LDAb2 + LDAb4*rs))))
  bigdenom=srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs + LDAb2)
  dLDAdn=((-2*LDAA*LDAa1*ddlog(1.d0 + &
&    1.d0/(2*LDAA*bigdenom))) &
&    + (LDAa1*rs + 1.d0)*((LDAb3*rs + LDAb1)/(2*srs) + srs*LDAb3 +  &
&    2*LDAb4*rs + LDAb2)/((bigdenom**2)+bigdenom/(2*LDAA)))*drsdn
  ELDA0= -b1c/(1.d0 + b2c*sqrt(rs) + b3c*rs)
  dLDA0dn=b1c*(0.5d0*b2c/sqrt(rs)+b3c)/((1.d0 + b2c*sqrt(rs) + b3c*rs)**2) &
&    *drsdn
  ddLDA=ELDA0-ELDA
  beta= betaMB*(1.d0 + 0.1*rs)/(1.d0 + 0.1778d0*rs)
  dbetadn=-(0.0778d0*betaMB/(1.d0 + 0.1778d0*rs)**2)*drsdn
  if(baralpha<1.d-13) then
     fc=ddexp(-SCANc1c*baralpha/(1.d0-baralpha))
  else if (baralpha.lt.2.5d0) then
     fc=cc0+baralpha*(cc1+baralpha*(cc2+baralpha*(cc3+baralpha*( &
&        cc4+baralpha*(cc5+baralpha*(cc6+baralpha*cc7))))))
  else if (baralpha.ge.2.5d0) then
     fc=-SCANdc*ddexp(SCANc2c/(1.d0-baralpha))
  endif
  if(baralpha<1.d-13) then
     ddfc=-(SCANc1c/((1.d0-baralpha)**2))*ddexp(-SCANc1c*baralpha/(1.d0-baralpha))
  else if (baralpha.lt.2.5d0) then
     ddfc=(cc1+baralpha*(2*cc2+baralpha*(3*cc3+baralpha*( &
&        4*cc4+baralpha*(5*cc5+baralpha*(6*cc6+baralpha*7*cc7))))))
  else if (baralpha.ge.2.5d0) then
     ddfc=-(SCANdc*SCANc2c/((1.d0-baralpha)**2))*ddexp(SCANc2c/(1.d0-baralpha))
  endif
  w1= ddexp(-ELDA/Sgam) - 1.d0
  dw1dn=-(ddexp(-ELDA/Sgam)/Sgam)*dLDAdn
  w0= ddexp(-ELDA0/b1c) - 1.d0
  dw0dn=-(ddexp(-ELDA0/b1c)/b1c)*dLDA0dn
  ginfinity= 1.d0/(1.d0 + 4*chiinfinity*(p))**0.25d0
  dgindn=-(chiinfinity*ginfinity/(1.d0 + 4*chiinfinity*(p)))*dpdn
  dginds=-(chiinfinity*ginfinity/(1.d0 + 4*chiinfinity*(p)))*dpds
  H0c=  b1c*ddlog(1.d0 + w0*(1.d0 - ginfinity))
  dH0cdn=(b1c/(1.d0+w0*(1.d0-ginfinity)))*((1.d0-ginfinity)*dw0dn-w0*dgindn)
  dH0cds=(-b1c*w0*dginds/(1.d0+w0*(1.d0-ginfinity)))
  ec0=  ELDA0 + H0c
  dec0dn=dLDA0dn+dH0cdn
  dec0ds=dH0cds
  difddLDA=b1c*(b2c/(2*srs) + b3c)/(1.d0 + b2c*srs + b3c*rs)**2 + &
&  2*LDAA*LDAa1*ddlog(1.d0 + &
&    1.d0/(2*LDAA*(srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs +  LDAb2)))) &
&    - (LDAa1*rs + 1.d0)*((LDAb3*rs + LDAb1)/(2*srs) + srs*LDAb3 +  &
&    2*LDAb4*rs + LDAb2)/(((srs*(LDAb3*rs + LDAb1) + &  
&    rs*(LDAb4*rs + LDAb2))**2)*(1.d0       +    1.d0 &
&  /(2*LDAA*(srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs + LDAb2)))))
  y=beta*(t2)/(Sgam*w1)
  dyds=beta*dt2ds/(Sgam*w1)
  dydn=(1.d0/(Sgam*w1))*(dbetadn*t2+beta*dt2dn-beta*t2*dw1dn/w1)
  ddy=(Dfc2/(27*Sgam*w1))*(20*rs*difddLDA-45*eta*ddLDA)*p*ddexp(-(p**2)/(dp2**4))
  gg=1.d0/(1.d0+4*(y-ddy))**0.25d0
  H1c=Sgam*ddlog(1.d0+w1*(1.d0-gg))
  ec1=ELDA+H1c
  pp=p*ddexp(-(p**2)/(dp2**4))
  dppds=((ddexp(-p**2/dp2**4))*(dp2**4 - 2*p**2)/dp2**4)*dpds
  dppdn=((ddexp(-p**2/dp2**4))*(dp2**4 - 2*p**2)/dp2**4)*dpdn
  term1=1.d0+4*(y-ddy)
  term4=term1**0.25d0
  dddyds=(Dfc2/(27*Sgam*w1))*(20*rs*difddLDA-45*eta*ddLDA)*dppds
  dddydnp=(Dfc2/(27*Sgam*w1))*(20*rs*difddLDA-45*eta*ddLDA)*dppdn
  dddydnw=-(Dfc2/(27*Sgam*w1))*(20*rs*difddLDA-45*eta*ddLDA)*pp*dw1dn/w1
  dddydnr=(Dfc2/(27*Sgam*w1))*(20*difddLDA)*pp*drsdn
  dH1cds=(Sgam*w1/(term1*(term4*(1.d0+w1)-w1)))*(dyds-dddyds)
  stf1=-b1c*(8*srs*rs*(b3c**2) + 9*b2c*b3c*rs + 3*srs*(b2c**2) + b2c)&
&          /(4*rs*srs*((1.d0 + b2c*srs + b3c*rs)**3))
  stf2=srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs + LDAb2)
  stf3=1.d0+1.d0/(2*LDAA*stf2)
  stf4=(2*LDAa1*((LDAb3*rs+LDAb1)/(2*srs)+srs*LDAb3+2*LDAb4*rs+LDAb2)) &
&      /((stf2**2)*stf3)          
  stf5=(2*(LDAa1*rs+1.d0)*(((LDAb3*rs+LDAb1)/(2*srs)+srs*LDAb3+  &
&     2*LDAb4*rs + LDAb2)**2))/((stf2**3)*stf3)
  stf6=(LDAa1*rs+1.d0)*(-(LDAb3*rs+LDAb1)/(4*srs*rs)+LDAb3/srs + 2*LDAb4) &
&   /((stf2**2)*stf3)          
  stf7=(LDAa1*rs+1.d0)*& 
&  (((LDAb3*rs+LDAb1)/(2*srs)+srs*LDAb3+2*LDAb4*rs+LDAb2)**2) &
&   /(2*LDAA*(stf2**4)*(stf3**2))
  dif2ddLDA=stf1-stf4+stf5-stf6-stf7
  dddydnn=(Dfc2/(27*Sgam*w1))*(20*rs*dif2ddLDA-45*eta*difddLDA)*pp*drsdn
  dddydn=dddydnw+dddydnp+dddydnr+dddydnn
  dH1cdwn=Sgam*((term4-1.d0)/(term4*(w1+1.d0)-w1))*dw1dn
  dH1cdn=dH1cdwn+(Sgam*w1/(term1*(term4*(1.d0+w1)-w1)))*(dydn-dddydn)
  dec1dn=dLDAdn+dH1cdn
  dec1ds=dH1cds 

  corarg=rr*(ec1 + fc*(ec0 - ec1))
  vctau=rr*ddfc*(ec0 - ec1)/(U + eta*W)
! density and sigma derivative terms

  vcn=ec1+fc*(ec0-ec1)+ddfc*dalphadn*rr*(ec0-ec1)+dec0dn*rr*fc+dec1dn*rr*(1.d0-fc)
  vcs=ddfc*dalphads*rr*(ec0-ec1)+dec0ds*rr*fc+dec1ds*rr*(1.d0-fc)
  eex=exarg;eec=corarg;               exc=exarg+corarg
  vtx=vxtau;vtc=vctau;                vtau=vxtau+vctau
  vvxn=vxn; vvcn=vcn;                 vxcn=vxn+vcn
  vvxs=vxs; vvcs=vcs;                 vxcs=vxs+vcs
  

END SUBROUTINE  

END Module
