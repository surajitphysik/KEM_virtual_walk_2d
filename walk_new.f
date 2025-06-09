c##########################################################################################################################
c############<<<<<< Disorder induced phase transition in kinetic models of opinion dynamics >>>>>> ########################
c################### ###########################<<<<< &  >>>>>>#############################################################
c#################################<<<<<< Virtual walks in spin space >>>>>>>###############################################

         implicit real*8(a-h,o-z)


c**************************************************************************8       
       real*8 , dimension(:),allocatable :: srr
       real*8 , dimension(:),allocatable :: uku,drr,hrr                            !  define iox1 ...iox10, ioy1...ioy10 and remove ioo1...ioo10 also remove ioo
       integer , dimension(:),allocatable :: iox,ioy,ku
       

c       integer, dimension(:),allocatable ::  idox,idoy
c       integer, dimension(:,:),allocatable :: ifo1,ifo2

       integer, dimension(:),allocatable :: iro,ito,mu
c**************************************************************************

c****************************************************************************        
       real*8 , dimension(:,:),allocatable :: avrw1,avrw2


       real*8 , dimension(:,:),allocatable :: sdvm,sdv
          
       real*8 , dimension(:,:),allocatable :: svxy,svnxy
       integer, dimension(:,:),allocatable :: ixy,inxy
c***********************************************************************
c                     RANDOM_SEED INITIALIZATION
c**********************************************************************




            integer :: ivalues(1:8)

             integer,dimension(:), allocatable :: iseed
           


           call date_and_time(values=ivalues)

           call random_seed(size=k)
           allocate(iseed(1:k))
           iseed(:) = ivalues(8)
           call random_seed(put=iseed)






c********************************************************************
c               open input output file
c*******************************************************************
                

         open(102,file='input.dat',status='old') 

          open(787,file='test.dat',status='unknown')
          open(780,file='rw.dat',status='unknown')
          open(120,file='mu_xn.dat',status='unknown')
         open(122,file='mu_x.dat',status='unknown')
         open(121,file='xy_pos_proj.dat',status='unknown')
c         open(123,file='x_neg_proj.dat',status='unknown')
         open(833,file='ratio.dat',status='unknown') 
c          open(622,file='idm.dat',status='unknown')   
          open(922,file='susceptibility.dat',status='unknown')  
         
         
         
         
                                                                  
                                                                       
c*************************************************************************************
c                    initial value
c************************************************************************************* 




         read(102,131)dlta,incd,nl,qtemp,iav,itime
131           format(F5.3,I6,I6,F6.2,I6,I6) 

            ltn=10
          
           ntotal=nl           ! total number of sites(L)
             kt=itime+10

             avtotal=iav
             
              ikt=itime
               ikt1=ikt+2        
                 
                      




c************************************************************************************






                      
c*****************************************Allocate*****************************************************




c        allocate(hr1(1:itime),hr2(1:itime))
        allocate(mu(1:nl),ito(1:nl),srr(1:nl),ioy(1:nl))
        allocate(uku(1:nl),drr(1:nl),hrr(1:nl),iro(1:nl))
        allocate(iox(1:nl),ku(1:nl))           
c        allocate(idox(1:nl),idoy(1:nl))
c        allocate(orp(iav))
c        allocate(ifo1(nl,it1))
c        allocate(ifo2(nl,it2))
c*************************************************************************************************
        allocate(avrw1(kt,kt),avrw2(kt,kt))        
        allocate(sdvm(kt,kt),sdv(kt,kt))        
        allocate(ixy(kt,kt),inxy(kt,kt))
        allocate(svxy(kt,kt),svnxy(kt,kt))


c**********************************************************************************************


   


           do kc=1,kt
           do kp=1,kt
            avrw2(kc,kp)=0.0
            avrw1(kc,kp)=0.0
           svxy(kc,kp)=0.0
           svnxy(kc,kp)=0.0
           sdv(kc,kp)=0.0
           sdvm(kc,kp)=0.0
           ixy(kc,kp)=0.0
           inxy(kc,kp)=0.0
           enddo  
           enddo              




C*********************************************************************************************
            do iit=1,iav                       !do loop for average value calculation
             it=0                              !initialisation of it for every loop
              isum=0






c*********************************************************************************************
c              initialization for random walk                  *******************************
c*********************************************************************************************
            

            do k=1,nl
           
              iox(k)=0
              
               
               ioy(k)=0
               
              
           enddo       

c*********************************************************^^^^^ modify upper portion^^^^^**********************************



            itf=itime+1
            

c*********************************************************************
c            random initial condition iro(m,n)
c********************************************************************* 
          if(incd==1)then

         call random_number(srr)
         

         do m=1,nl
         
         
         

         iro(m)=FLOOR(3.0*srr(m)-1.0)       

         enddo  
          endif
c******************************************************************
c       delta initial condition
c****************************************************************
        if(incd==3)then
         ls1=0
         ls2=0
          ls3=0
        asp=(0.5-(0.5*dlta))*ntotal
         del=ntotal*dlta
          lsp=floor(asp)
         do m=1,nl
c         do n=1,nl

         if (ls1.lt.lsp)then
         ls1=ls1+1

            ku(m)=-1         
          else 
          if (ls2.lt.lsp)then
          ls2=ls2+1
           ku(m)=1
           else
           ku(m)=0
            endif
            endif

            enddo
c            enddo
            
           call random_number(uku)
c           call random_number(yku) 
            do mg=1,nl
c          do ng=1,nl
          iku=FLOOR(nl*uku(mg)+1.0)
c           jku=FLOOR(nl*yku(mg)+1.0)
           iro(mg)=ku(iku)
c            enddo
            enddo


          endif
            
c***********************************************************************
c            polarized initial condition
c**********************************************************************
        if(incd==2)then
        
        do i=1,nl
       

         iro(i)=1

c        enddo
        enddo

          endif        
c********************************************************************************
       do k1=1,nl

       jn=iro(k1)
       if (jn.eq.1)then
       iox(k1)=1
       ioy(k1)=0
       else
       if(jn.eq.-1)then
       iox(k1)=-1
       ioy(k1)=0
       else
       if(jn.eq.0)then
       iox(k1)=0
       ioy(k1)=1
       else
       endif
       endif
       endif
       
       
       enddo

       
       
       
500    continue
c*********************************************************************************
c       mu(ij) INTERACTION PARAMETER  (NEGETIVE WITH PROBABILITY nq)
c*********************************************************************************
          is=0


         nq=ntotal*qtemp



         do m=1,nl

         
         if (is.le.nq)then
         is=is+1

            mu(m)=-1         
          else 
  
          mu(m)=1     
          endif
   

         enddo
         

           
c***********************************************************************
c                  MONTE CARLO STEP START         
c***********************************************************************             

           
             do k=1,nl

              ito(k)=iro(k)

            enddo


            
           call random_number(drr)
           call random_number(hrr)
           call random_number(srr)
          
       

          do m=1,nl
         
            j1=FLOOR(nl*drr(m)+1.0)
        


           j2=FLOOR(nl*hrr(m)+1.0)


           m1=FLOOR(nl*srr(m)+1.0)
           
           jt=it+1
c******************************************           
           

c***********************************************                   

           if(j1==j2)then
           iro(j1)=iro(j1)+1.0*iro(j2)
           else
          iro(j1)=iro(j1)+mu(m1)*iro(j2)
           endif

c***********************************************
           if(iro(j1).gt.1)then
             iro(j1)=1
            else
            
            if(iro(j1).lt.-1)then
             iro(j1)=-1
             else
             endif
             endif 
                    
          
          
         
            enddo


             it=it+1
c/////////////////////////////////////////////////////////////             

c////////////////////////////////////////////////////////////




c////////////////////////////////////////////////////////
              


                 

c***************************************************************************8
c************************RANDOM WALK********************************************                 
c******************************************************************random walk



c******************************************************walk 2***************************************************

                do iu=1,nl
c                do ib=1,nl
                lv1=ito(iu)
                lv2=iro(iu)
                lv=lv1+lv2
               if(lv2.eq.1)then                         
             

c             ioo(iu)=ioo(iu)+1     
               iox(iu)=iox(iu)+1
c               ioy(iu)=ioy(iu)       !define 2 matrix iox(iu) and ioy(iu), nl,remove ioo (done)
              else
               if(lv2.eq.-1)then                         !random walk scheme................................................,,,,,,,,,,,,
c               ioo(iu)=ioo(iu)-1

              iox(iu)=iox(iu)-1
c               ioy(iu)=ioy(iu)   

              else   
              if(lv2.eq.0)then            
c              
c                 ioo(iu)=ioo(iu)+0                
              
               ioy(iu)=ioy(iu)+1
                  
              else

                endif
                endif
                endif
               
               
               

c              enddo
               enddo
                   

c############################################################ walk 1 #######################################################


               ontotal=ntotal
c******************************************************new    ne


      
              if(it.lt.itime)then
              go to 500
               else

                
               endif

             
c****************order parameter***********************             


c**********************************************************************random walk 
c***********************************************************new  new***************************************************
 

c                do ia=1,nl

c                idox(ia)=iox(ia)
c                idoy(ia)=ioy(ia)
c                 enddo                         !define idox,idoy (done)
c**************************************************
c                  itf=it1f
c                 itf1=it1fp
                 !itf=itime+1>itf1, itime>itf
               

!*************************************************************
                do j=1,nl
          
             ka=iox(j)
             kb=ioy(j)
             
             if(ka.ge.0)then	
             ko=ka+1
             koy=kb+1
             ixy(ko,koy)=ixy(ko,koy)+1
             else
             kw=-ka
             koy=kb+1
             inxy(kw,koy)=inxy(kw,koy)+1
             endif

             
             enddo


               do m=1,ikt1
               do ms=1,ikt1

                 gh=ixy(m,ms)      !irwx(m)                          !irxx,irwxm,irwy,irwym>ixy(i,j),inxy(i,j) >i=1,,j=1.itfg
                 gj=ntotal
                 ghm=inxy(m,ms)              
                al=gh/gj                            !al & tl is the probability for different position of random walker
                tl=ghm/gj
                 
c#################################################
                svxy(m,ms)=al                   
                 svnxy(m,ms)=tl
c#################################################   

              enddo
              enddo
              write(787,*)svxy(100,0)
c*************************************************************
c         initialisation of irwx & irwxm
c**************************************************************
           do kgo=1,ikt1
           do igo=1,ikt1
           ixy(kgo,igo)=0
            inxy(kgo,igo)=0
           enddo
           enddo
cllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll


!*************************************************************new new***************************************8
             do m=1,ikt1
             do j=1,ikt1 
          
             
c###################################################
c               s1v(mo,iit)=sv(mo,iit)
c               s1vm(mo,iit)=svm(mo,iit)
c###################################################
              sdv(m,j)=sdv(m,j)+svxy(m,j)
              sdvm(m,j)=sdvm(m,j)+svnxy(m,j)
               enddo
               enddo
               


c              enddo               ! enddo for i


!llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll


            enddo     !end do loop for average 'do iit=1,iav'

ciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii           
c         
c*******************************************************************
c********** configuration average of random walk****************
c********************************************************************

                                                                                               !  do k=1,10
             
c*****************************************************************8            
 
c***************************************************************
c             
c******************************************
            do kr=1,ikt1
            do mr=1,ikt1


                avtotal=iav  

                av1=sdv(kr,mr)
               av2=sdvm(kr,mr)
              avp1=av1/avtotal
               avp2=av2/avtotal
             avrw1(kr,mr)=avp1
             avrw2(kr,mr)=avp2        !- ve rw    

            enddo
            enddo

            

cnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn            +y                       nnnnnnnnnnnnnnnnnnnnnnnnnn
             iktt=ikt+1
             do kap=1,iktt
c             do map=1,itf  
                 
                zs=-(ikt1-kap)
                 mpa=ikt1-kap 
                 
               do jap=1,ikt1                                       
                vzp=(jap-1)                
                   
                 
c               dq=avrw2(mpa)*yrw1(jap)
                dq=avrw2(mpa,jap)
                if(dq.ne.0.0)then 
          
          
               
                
                write(780,*)zs,vzp,dq
                else
                 
                 endif
              enddo 
              
                
                
                write(780,*)
                
              enddo

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
               
                
             do kap=1,ikt1  
                zp=(kap-1) 
             do jap=1,ikt1                                        
                vzp=(jap-1)                
                    
                
                
                aq=avrw1(kap,jap)
                

                 if(aq.ne.0.0)then
              
                 write(780,*)zp,vzp,aq
                 else
                 endif
                 
              enddo 
              
              
              
                
                write(780,*)
                
              enddo
               
!*********************************************               
              iktn=ikt+2
              do ix=1,iktn
              ix1=ix-1
              
c              ixn=ix-ktn
c              ixp=-ixn
              yav=0.0
              yavn=0.0
              dx=0.0
              dxn=0.0
              dy=0.0
              yh=0.0
              yh2=0.0

              do iu=1,iktn
              dx=dx+avrw1(ix,iu)
c              dxn=dxn+avrw2(ixp,iu)
              enddo
              do ig=1,iktn
              dy=dy+avrw1(ig,ix)
              enddo

              if(dx.ne.0.0.or.dy.ne.0.0)then

              do iy=1,iktn
              avy=(avrw1(ix,iy)/dx)*(iy-1)
c              avyn=(avrw2(ixp,iy)/dxn)*iy
              yav=yav+avy
c              yavn=yavn+avyn
              enddo

              do iv=1,iktn
              ivv=iv-1
              yh=yh+((ivv-yav)*(avrw1(ix,iv)/dx))
              yh2=yh2+(((ivv-yav)**2)*(avrw1(ix,iv)/dx))
              enddo
               yh2=yh2**(0.5)
c              write(822,*)ix,yh
              
c***********************************************************              
c                 aen=(jen+ien)             
c                 entropy=log(aen)
                
                write(122,*)ix1,yav
c                write(120,*)ixn,yavn 
                 write(121,*)ix1,dx,dy
c                  write(123,*)ixn,dxn   
                   
                                
              
c***********************************************************8888              
              
              
              
              
              
              
              
              
              
              else
              endif
              
              enddo
             
!**********************************************8
                
              do ig=1,kt
              do kg=1,kt
                 avrw2(ig,kg)=0.0
                  avrw1(ig,kg)=0.0
                  enddo
                  enddo
c             do is=1,itf2
c                 yrw2(is)=0.0
c                  yrw1(is)=0.0
c                  enddo
                 
c                write(922,*)aen,entropy
c433            go to 636
                                                                  !                  enddo                        !enddo for k=1,10

399          continue


c              write(622,*) idm,'give large idm',ka,kb               

c*************************************************************

636            continue               
               end program

c*************************************************************
            






