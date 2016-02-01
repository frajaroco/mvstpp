!
!     Francisco J. Rodriguez-Cortes, July 2015 
!
!     This function provides an kernel 
!     estimator of the temporal mark variogram.  
!

       subroutine gtecore(x,y,txy,n,t,nt,kt,ht,gte)     	   

       implicit real*8(a-h,o-z)

       integer i,j,iv,n,nt,kt
       double precision wij,vij,ht,kernt,gtem,gten,gte,x,y,txy
       double precision hij,tij,xi,yi,ti,pi,two
       dimension x(n),y(n),txy(n),t(nt),gte(nt),gtem(nt),gten(nt)
       
       gtem=0d0
       gten=0d0
      
          two=2d0
		 
       do iv=1,nt
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)      
          do j=1,n
           if (j.ne.i) then
            hij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
            tij=abs(ti-txy(j))
              if (kt.eq.1) then
               kernt=boxkernel((t(iv)-tij)/ht,ht)
                else if (kt.eq.2) then
                 kernt=ekernel((t(iv)-tij)/ht,ht)
                  else if (kt.eq.3) then
                   kernt=qkernel((t(iv)-tij)/ht,ht)
              end if
             if ((kernt.ne.0d0).and.(t(iv).ne.0d0)) then
                    wij=(((hij**two)/two)*kernt)
                    vij=kernt
                    gtem(iv)=gtem(iv)+wij
                    gten(iv)=gten(iv)+vij                    
             end if      
           end if
          end do
          end do 
          if (t(iv).eq.0d0) then
            gte(iv)=0
            else
            gte(iv)=gtem(iv)/gten(iv)   
           end if
          end do
       
        return
        
        end  
