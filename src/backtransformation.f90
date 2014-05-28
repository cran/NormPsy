


!====================================================================

! PROGRAMME PRINCIPAL

!-----------------------------------------------------------
!           CALCUL DE H-1(Y)
!------------------------------------------------------------



!===================================================================


subroutine backtransformation(mu,VC0,VC1,maxmes,spl &
     ,nbzitr,zitr0,nsim,methInteg,Ymarg)


  IMPLICIT NONE

  ! in input
  double precision,dimension(maxmes),intent(in) ::mu
  integer,intent(in)::maxmes,nbzitr,nsim,methInteg
  double precision,dimension(nbzitr),intent(in)::zitr0
  double precision,dimension(nbzitr+2),intent(in)::spl
  double precision, dimension(maxmes),intent(in)::VC0
  double precision, dimension(maxmes*(maxmes+1)/2),intent(in)::VC1
  
  ! for computation
  integer ::j,k,ntrtot,ier,l,niter,m
  double precision,dimension(:),allocatable::ysim,usim
  double precision,dimension(:,:),allocatable::V2
  double precision,dimension(maxmes*(maxmes+1)/2)::V1
  double precision :: eps,ytemp2,x22
  double precision ::ytemp,diff,SX,bb
  double precision,dimension(:),allocatable::zitr,splaa
  double precision::INV_ISPLINES
  double precision,dimension(2,51)::gauss
  
  ! for output
  double precision,dimension(maxmes),intent(out) ::Ymarg



  !============= recup des places de parametres

  allocate(ysim(maxmes),usim(maxmes))!,V1(maxmes*(maxmes+1)/2))
  allocate(V2(maxmes,maxmes))

  ntrtot=nbzitr+2
  allocate(zitr(-1:(ntrtot)),splaa(-1:(ntrtot-3)))
  zitr(1:nbzitr)=zitr0(1:nbzitr)
  zitr(-1)=zitr(1)
  zitr(0)=zitr(1)
  zitr(ntrtot-1)=zitr(ntrtot-2)
  zitr(ntrtot)=zitr(ntrtot-1)

  ysim=0.d0
  usim=0.d0
  !V1=0.d0
  V2=0.d0
  splaa=0.d0

  ymarg=0.d0
  eps=1.d-20
  ier=0
  x22=0.d0
  
  bb=spl(1)
  splaa(-1:(ntrtot-3))=spl(2:ntrtot)


        if (methInteg.eq.1) then
           ! i.e. methode de MC a faire
           
            V1=0.d0  
            do j=1,maxmes*(maxmes+1)/2
                V1(j)=VC1(j)    
            end do
           
            CALL DMFSD(V1,maxmes,EPS,IER)  !V1 est la partie sup de la chol colonne par colonne
            if (ier.eq.-1) then
                ymarg=9999.d0
                goto 654
            end if
  
            V2=0.d0
            do j=1,maxmes
              do k=1,j
                V2(j,k)=V1(k+j*(j-1)/2)  !partie inferieure de V2 non nulle
              end do
            end do
   
           
           do l=1,nsim
              usim=0.d0
              ysim=0.d0
              do m=1,maxmes
                 SX=1.d0
                 call bgos(SX,0,usim(m),x22,0.d0)
              end do

              ysim=mu+MATMUL(V2,usim)
              do j=1,maxmes
                 niter=0
                 diff=0.d0
                 ier=0
                 
                 ytemp=INV_ISPLINES(ysim(j),splaa,bb,nbzitr,zitr,ier,niter,diff)  
                    if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(j).eq.9999.d0) then
                       ymarg(j)=9999.d0
                    else
                        ymarg(j)=ymarg(j)+ytemp/dble(nsim)
                    end if
              end do
           end do


        else
     
           call gausshermite(gauss,nsim)
      
           do j=1,maxmes
              do l=1,nsim
                 niter=0
                 diff=0.d0
                 ier=0
                 ytemp=mu(j)+sqrt(VC0(j))*gauss(1,l)
                 ytemp2=INV_ISPLINES(ytemp,splaa,bb,nbzitr,zitr,ier,niter,diff) 
                 if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3)) then
                       ymarg(j)=9999.d0
                 else
                    ymarg(j)=ymarg(j)+ytemp2*gauss(2,l)
                 end if
              end do
           end do
        end if



654 continue

  deallocate(zitr,splaa)
  deallocate(ysim,usim,V2)



  return

end subroutine backtransformation







      subroutine dmfsd(a,n,eps,ier)
!
!   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
!   MATRICE = TRANSPOSEE(T)*T
!   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
!            PAR COLONNE DE LA MATRICE A FACTORISER
!   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
!
!   SUBROUTINE APPELE PAR DSINV
!
!   N : DIM. MATRICE
!   EPS : SEUIL DE TOLERANCE
!   IER = 0 PAS D'ERREUR
!   IER = -1 ERREUR
!   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(out)::ier
      double precision,intent(in)::eps
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind

!
!   TEST ON WRONG INPUT PARAMETER N
!        
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
!
!   INITIALIZE DIAGONAL-LOOP
!
      kpiv=0
      do k=1,n
          kpiv=kpiv+k
          ind=kpiv
          lend=k-1
!
!   CALCULATE TOLERANCE
!
          tol=dabs(eps*sngl(A(kpiv)))
!
!   START FACTORIZATION-LOOP OVER K-TH ROW
!
         do i=k,n
            dsum=0.d0
            if (lend.lt.0) goto 2
            if (lend.eq.0) goto 4
            if (lend.gt.0) goto 2
!
!   START INNER LOOP
!
2           do l=1,lend
               lanf=kpiv-l
               lind=ind-l
               dsum=dsum+A(lanf)*A(lind)
            end do

!
!   END OF INNEF LOOP
!
!   TRANSFORM ELEMENT A(IND)
!
4           dsum=A(ind)-dsum
            if (i-k.ne.0) goto 10
            if (i-k.eq.0) goto 5
!   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
!


5           if (sngl(dsum)-tol.le.0) goto 6
            if (sngl(dsum)-tol.gt.0) goto 9
6           if (dsum.le.0) goto 12
            if (dsum.gt.0) goto 7
7           if (ier.le.0) goto 8
            if (ier.gt.0) goto 9
8           ier=k-1
!
!   COMPUTE PIVOT ELEMENT
!
9           dpiv=dsqrt(dsum)
            A(kpiv)=dpiv
            dpiv=1.D0/dpiv
            goto 11
!
!   CALCULATE TERMS IN ROW
!
10          A(ind)=dsum*dpiv
11          ind=ind+i
         end do
      end do

!
!   END OF DIAGONAL-LOOP
!
      return
12    ier=-1
      return

      end subroutine dmfsd



!C ------------------------------------------------------------
!C
!C     EVALUATION OF I-SPLINES and M-splines
!C ------------------------------------------------------------




      SUBROUTINE eval_splines(X00,Ispl,Mspl,splaa,bb,nztr,zi_eval)

      implicit none
      integer ::k,l,i,nztr
      double precision::X00,X0,Ispl,Mspl,som,bb
      double precision,dimension(-1:nztr+2)::zi_eval
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n &
      ,hn,hht,mmeval,mm1eval,mm2eval,imeval,im1eval,im2eval

      double precision,dimension(-1:nztr-1)::splaa


!C ou se trouve la valeur de X


      mmeval=0.d0
      mm1eval=0.d0
      mm2eval=0.d0
      imeval=0.d0
      im1eval=0.d0
      im2eval=0.d0
  
      X0=zi_eval(1)+(zi_eval(nztr)-zi_eval(1))*(1.d0-1.d0/(1.d0+exp(X00)))
      l=0
      do k = 2,nztr
         if ((X0.ge.zi_eval(k-1)).and.(X0.lt.zi_eval(k))) then
            l=k-1
         endif
      end do

      if (X0.eq.zi_eval(nztr)) then
         l=nztr-1
      end if

      ht2 = zi_eval(l+1)-X0
      htm= X0-zi_eval(l-1)
      ht = X0-zi_eval(l)
      ht3 = zi_eval(l+2)-X0
      hht = X0-zi_eval(l-2)
      h = zi_eval(l+1)-zi_eval(l)
      hh= zi_eval(l+1)-zi_eval(l-1)
      hn= zi_eval(l+1)-zi_eval(l-2)
      h2n=zi_eval(l+2)-zi_eval(l-1)
      h2= zi_eval(l+2)-zi_eval(l)
      h3= zi_eval(l+3)-zi_eval(l)

      if (h.eq.0.or.hh.eq.0.or.hn.eq.0.or.h2n.eq.0.or.h2.eq.0.or.h3.eq.0)  then
         Mspl=1.d9
         Ispl=1.d9
         go to 587
      end if

      if (X0.ne.zi_eval(nztr)) then
         mm2eval = (3.d0*ht2*ht2)/(hh*h*hn)
         mm1eval = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
         mmeval  = (3.d0*ht*ht)/(h3*h2*h)
      end if
      if (X0.eq.zi_eval(nztr)) then

         mm2eval = 0.d0
         mm1eval = 0.d0
         mmeval  = 3.d0/h

      end if

      if (mm2eval.lt.0.or.mm1eval.lt.0.or.mmeval.lt.0) then
         Mspl=1.d9
         Ispl=1.d9
         go to 587
      end if

      im2eval=hht*mm2eval/(3.d0)+ h2n*mm1eval/(3.d0)+h3*mmeval/(3.d0)
      im1eval=htm*mm1eval/(3.d0)+h3*mmeval/(3.d0)
      imeval=ht*mmeval/(3.d0)

      som=0.d0
      if (l.gt.1) then
         do i=2,l
            som=som+splaa(i-3)
         end do
      end if

      Ispl=bb+ som +splaa(l-2)*im2eval+splaa(l-1)*im1eval+splaa(l)*imeval
      Mspl= (splaa(l-2)*mm2eval+splaa(l-1)*mm1eval+splaa(l)*mmeval)*      &
            (1.d0-1.d0/((1.d0+exp(X00))**2))*(zi_eval(nztr)-zi_eval(1))

 587   continue


      end subroutine eval_splines





!==================================================================
!
      !  SUBROUTINES INVERSION ISPLINES
!
!==================================================================



!C ------------------------------------------------------------
!C
!C     INVERSE de I-splines
!C
!C     BY USING NEWTON-RAPHSON METHOD
!C ------------------------------------------------------------

      double precision FUNCTION INV_ISPLINES(X00,splaa,bb,nztr,zi_eval,istop,iter,eps)


      implicit none
      integer::iter,istop,nztr
      double precision,dimension(-1:nztr+2)::zi_eval
      double precision::X0,X00,X1,fx0,f1x0,eps,bb,bb1
      double precision,dimension(-1:nztr-1)::splaa

                                 
       !write(*,*)'X00',X00      
      eps=1.d-5                  
      iter=1                   
      X0=1.d10                 
      call eval_splines(X0,fx0,f1x0,splaa,bb,nztr,zi_eval)
                             
      if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
         INV_ISPLINES=fx0        
         istop=3                 
         goto 1234               
      end if                     
                                 
      if (X00.ge.fx0) then       
         INV_ISPLINES=zi_eval(nztr)
         istop=1                 
         goto 1234               
      end if                     
      X0=-1.d10                  
      call eval_splines(X0,fx0,f1x0,splaa,bb,nztr,zi_eval)
      if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
         INV_ISPLINES=fx0        
         istop=3                 
         goto 1234               
      end if                     
      if (X00.le.fx0) then       
         INV_ISPLINES=zi_eval(1)
         istop=1
         goto 1234
      end if
      bb1=bb-X00
      X0=0
      call eval_splines(X0,fx0,f1x0,splaa,bb1,nztr,zi_eval)
      if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
         INV_ISPLINES=fx0
         istop=3
         goto 1234
      end if
      X1=X0-fx0/f1x0
      do while (ABS((X1-X0)/X0).GT.EPS.and.iter.lt.500)
         iter=iter+1
         X0=X1
         call eval_splines(X0,fx0,f1x0,splaa,bb1,nztr,zi_eval)
        if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
            INV_ISPLINES=fx0
            istop=3
            goto 1234
        end if
        X1=X0-fx0/f1x0
      end do
      INV_ISPLINES=zi_eval(1)+(zi_eval(nztr)-zi_eval(1))*exp(X1)/(1.d0+exp(X1))
       
      if (ABS((X1-X0)/X0).le.EPS) then
         istop=1
      else if (iter.ge.500) then
         istop=2
      else
         istop=3
      end if

      eps=ABS((X1-X0)/X0)

1234  continue

      return
      end function INV_ISPLINES




!=============================================================
!
      !SUBROUTINES simulation
!
!=============================================================




!C ******************** BGOS ********************************


      SUBROUTINE BGOS(SX,ID,X1,X2,RO)


!C     ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)

      implicit none
      double precision ::RO,SX
      integer ::ID
      double precision ::F,V1,V2,S,DLS,RO2
      double precision ::X1,X2,UNIRAN
!C     write(*,*)'dans bgos'


 5    CONTINUE

!C     write(*,*)'avant rand :'

!C     X1=RAND()
!C     X2=RAND()

      X1=UNIRAN()
      X2=UNIRAN()

      IF(ID.NE.1) GO TO 10
      F=2.*SQRT(3.)
      X1=(X1-0.5)*F
      X2=(X2-0.5)*F
      GO TO 20
10    CONTINUE
      V1=2.*X1-1
      V2=2.*X2-1
      S=V1*V1+V2*V2
      IF(S.GE.1.) GO TO 5
      DLS=SQRT(-2.*LOG(S)/S)
      X1=V1*DLS
      X2=V2*DLS
20    CONTINUE
      RO2=RO*RO
      IF(ABS(RO).GT.1.E-10) X2=(X1+X2*SQRT(1./RO2-1.))*RO
      X1=X1*SX
      X2=X2*SX

!C      write(*,*) 'X1 ',X1,' X2 ',X2
!C OK, X1 et X2 sont cr??s

!C      write(*,*)'fin bgos'

      RETURN
      END subroutine bgos


!C ------------------- FIN SUBROUTINE BGOS -----------------


!C ------------------------------------------------------

      DOUBLE PRECISION FUNCTION UNIRAN()
!C
!C     Random number generator(RCARRY), adapted from F. James
!C     "A Review of Random Number Generators"
!C      Comp. Phys. Comm. 60(1990), pp. 329-344.
!C
      implicit none
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY, ONE
      PARAMETER ( ONE = 1, TWOM24 = ONE/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS /      &
    0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,      &
    0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
    0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
    0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      UNIRAN = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNIRAN .LT. 0 ) THEN
         UNIRAN = UNIRAN + 1
         CARRY = TWOM24
      ELSE
         CARRY = 0
      ENDIF
      SEEDS(I) = UNIRAN
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )

      return
      END function uniran










! ===================== GAUSSHERMITE NODES =======================


      subroutine gausshermite(gauss,npg)

!
!C
!C     Gauss-Hermite pour f(x)*exp(-x*x/2)/(rac(2*pi))
!C       = somme (f(xg)w(xg))
!C
!C
!C
      implicit none

      double precision, dimension(2,51)::Gauss
      double precision, dimension(51,51)::Wg,Tg
      integer :: npg,i


      DATA ( Wg(I, 5), Tg(I, 5), I = 1, 3) / &
       0.1125741132772071D-01, 0.2856970013872805D+01, &
       0.2220759220056126D+00, 0.1355626179974265D+01, &
       0.5333333333333342D+00, 0.9386691848789097D-16/
      DATA ( Wg(I, 7), Tg(I, 7), I = 1, 4) / &
       0.5482688559722184D-03, 0.3750439717725742D+01, &
       0.3075712396758645D-01, 0.2366759410734542D+01, &
       0.2401231786050126D+00, 0.1154405394739968D+01, &
       0.4571428571428575D+00, 0.2669848554723344D-16/
      DATA ( Wg(I, 9), Tg(I, 9), I = 1, 5) / &
      0.2234584400774664D-04, 0.4512745863399781D+01, &
       0.2789141321231769D-02, 0.3205429002856470D+01, &
       0.4991640676521780D-01, 0.2076847978677829D+01, &
       0.2440975028949394D+00, 0.1023255663789133D+01, &
       0.4063492063492066D+00, 0.0000000000000000D+00/
     DATA ( Wg(I,15), Tg(I,15), I = 1, 8) / &
       0.8589649899633300D-09, 0.6363947888829836D+01, &
       0.5975419597920602D-06, 0.5190093591304780D+01, &
       0.5642146405189029D-04, 0.4196207711269018D+01, &
       0.1567357503549958D-02, 0.3289082424398766D+01, &
       0.1736577449213763D-01, 0.2432436827009758D+01, &
       0.8941779539984458D-01, 0.1606710069028730D+01, &
       0.2324622936097322D+00, 0.7991290683245483D+00, &
       0.3182595182595181D+00, 0.0000000000000000D+00/
      DATA ( Wg(I,20), Tg(I,20), I = 1,10) / &
       0.1257800672437914D-12, 0.7619048541679760D+01, &
       0.2482062362315163D-09, 0.6510590157013660D+01, &
       0.6127490259983006D-07, 0.5578738805893195D+01, &
       0.4402121090230841D-05, 0.4734581334046057D+01, &
       0.1288262799619300D-03, 0.3943967350657311D+01, &
       0.1830103131080496D-02, 0.3189014816553389D+01, &
       0.1399783744710099D-01, 0.2458663611172367D+01, &
       0.6150637206397690D-01, 0.1745247320814126D+01, &
       0.1617393339840001D+00, 0.1042945348802752D+01, &
       0.2607930634495551D+00, 0.3469641570813557D+00/
      DATA ( Wg(I,30), Tg(I,30), I = 1,15) / &
      0.1640807008117853D-20, 0.9706235997359524D+01, &
       0.1585560944966296D-16, 0.8680837722732207D+01, &
       0.1624080129972436D-13, 0.7825051744352813D+01, &
       0.4573425871326147D-11, 0.7055396866960296D+01, &
       0.5178459467189710D-09, 0.6339997686869597D+01, &
       0.2882175154047618D-07, 0.5662381850082873D+01, &
       0.8909088868621158D-06, 0.5012600596486518D+01, &
       0.1657998163067346D-04, 0.4384020365898051D+01, &
       0.1965129439848249D-03, 0.3771894423159236D+01, &
       0.1544707339866097D-02, 0.3172634639420402D+01, &
       0.8295747557723240D-02, 0.2583402100229274D+01, &
       0.3111177018350134D-01, 0.2001858612956431D+01, &
       0.8278683671562172D-01, 0.1426005658374115D+01, &
       0.1580469532090208D+00, 0.8540733517109733D+00, &
       0.2179999718155776D+00, 0.2844387607362094D+00/
      DATA ( Wg(I,40), Tg(I,40), I = 1,20) / &
       0.1461839873869467D-28, 0.1145337784154873D+02, &
       0.4820467940200524D-24, 0.1048156053467427D+02, &
       0.1448609431551587D-20, 0.9673556366934033D+01, &
       0.1122275206827074D-17, 0.8949504543855559D+01, &
       0.3389853443248306D-15, 0.8278940623659475D+01, &
       0.4968088529197761D-13, 0.7646163764541459D+01, &
       0.4037638581695192D-11, 0.7041738406453829D+01, &
       0.1989118526027766D-09, 0.6459423377583766D+01, &
       0.6325897188548972D-08, 0.5894805675372016D+01, &
       0.1360342421574886D-06, 0.5344605445720084D+01, &
       0.2048897436081474D-05, 0.4806287192093873D+01, &
       0.2221177143247582D-04, 0.4277826156362752D+01, &
       0.1770729287992397D-03, 0.3757559776168985D+01, &
       0.1055879016901825D-02, 0.3244088732999869D+01, &
       0.4773544881823334D-02, 0.2736208340465433D+01, &
       0.1653784414256937D-01, 0.2232859218634873D+01, &
       0.4427455520227679D-01, 0.1733090590631720D+01, &
       0.9217657917006089D-01, 0.1236032004799159D+01, &
       0.1499211117635710D+00, 0.7408707252859313D+00, &
       0.1910590096619904D+00, 0.2468328960227240D+00/
      DATA ( Wg(I,50), Tg(I,50), I = 1,25) / &
       0.1034607500576990D-36, 0.1298588445541555D+02, &
       0.9443414659584510D-32, 0.1205301838092448D+02, &
       0.6856280758924735D-28, 0.1127923332148262D+02, &
       0.1206044550761014D-24, 0.1058738174919177D+02, &
       0.7995094477915292D-22, 0.9948035709637500D+01, &
       0.2522482807168144D-19, 0.9346039593575728D+01, &
       0.4368171816201588D-17, 0.8772299579514598D+01, &
       0.4566698246800344D-15, 0.8220815907982127D+01, &
       0.3083828687005300D-13, 0.7687362406712500D+01, &
       0.1414228936126661D-11, 0.7168814837853899D+01, &
       0.4576636712310442D-10, 0.6662775399018720D+01, &
       0.1077060789389039D-08, 0.6167347388659921D+01, &
       0.1888225976835208D-07, 0.5680992291033284D+01, &
       0.2514609880838772D-06, 0.5202434993399912D+01, &
       0.2584937658949391D-05, 0.4730598550228594D+01, &
       0.2078485175734569D-04, 0.4264557843038109D+01, &
       0.1321726328668984D-03, 0.3803505741742012D+01, &
       0.6708280619787080D-03, 0.3346727774732429D+01, &
       0.2738160896935348D-02, 0.2893582727707738D+01, &
       0.9045054154849623D-02, 0.2443487452654017D+01, &
       0.2430481286424306D-01, 0.1995904709795124D+01, &
       0.5334352453170102D-01, 0.1550333214338771D+01, &
       0.9593054035810168D-01, 0.1106299289397183D+01, &
       0.1416854132499443D+00, 0.6633496795082918D+00, &
       0.1721258519924433D+00, 0.2210451816445435D+00/



!      if(npg.ne.5.and.npg.ne.7.and.npg.ne.9.and.npg.ne.15.   &
!          and.npg.ne.20.and.npg.ne.30.and.npg.ne.40.and.npg.ne.50) then
!         write(*,*)'nb pts GH = 5,7,9,15,20,30,40, ou 50'
!         stop
!      end if


!ccccccccccccccccccccccccccccccccccccccccccccccc
      DO I = 1, NPG/2
         GAUSS(1,I) = -Tg(I,NPG)
         GAUSS(2,I) =  Wg(I,NPG)
         GAUSS(1,NPG-I+1) = Tg(I,NPG)
         GAUSS(2,NPG-I+1) = Wg(I,NPG)
      END DO
      IF ( MOD( NPG, 2 ) .EQ. 1 ) THEN
         GAUSS(1, NPG/2 + 1 ) = 0.D0
         GAUSS(2, NPG/2 + 1 ) = Wg( NPG/2 + 1, NPG )
      END IF

!ccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine gausshermite




!      program values
!      
!      implicit none
!      
!      double precision:: val,Ispl,Mspl,bb,invspl
!      double precision, dimension(0:30)::val2
!      integer:: nbzitr,k
!      double precision, dimension(-1:6)::splaa
!      double precision, dimension(-1:9)::zi 
!      double precision::INV_ISPLINES
!      
!      nbzitr=7
!      bb=0.d0
!      splaa(-1:nbzitr-1)=(/1033921.d-5,806659.d-5,1082235.d-5,1769599.d-5,1159447.d-5,1493948.d-5,1316887.d-5,1337304.d-5/)
!      !splaa(-1:nbzitr-1)=(/1033010.d-5,808522.d-5,1081574.d-5,176718.d-5,1161069.d-5,1494364.d-5,1317021.d-5,1337751.d-5/)                                   
!      !splaa(-1:nbzitr-1)=(/1029921.d-5,806659.d-5,1080335.d-5,1771099.d-5,1162447.d-5,1494948.d-5,1318787.d-5,1335804.d-5/)
!      !splaa(-1:nbzitr-1)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0/)       
!      zi(-1:9)=(/0.d0,0.d0,0.d0,10.d0,20.d0,23.d0,26.d0,28.d0,30.d0,30.d0,30.d0/)
!      val2= dble((/0,3,5,7,10,11,13,15,16,17,19,20,21,23,24,26,28,30,31,35,37,40,44,47,51,56,61,67,75,84,100/))
!     
!      do k=0,30
!        val=log(30/(30-dble(k))-1)
!        Ispl=0.d0
!        Mspl=0.d0
!        call eval_splines(val,Ispl,Mspl,splaa,bb,nbzitr,zi)
!        print*,"valeur transf de",k,"=",Ispl
!         !print*,"val2(k)=",val2(k)
!        !invspl=INV_ISPLINES(val2(k),splaa,bb,nbzitr,zi,0,0,0.d0)
!        !print*,"valeur inverse de Ispl =",invspl
!      end do
!
!      
!      stop
!      end program values
      
   

































