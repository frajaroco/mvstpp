!
!     Francisco J. Rodriguez-Cortes, July 2015
!
!     This function provides an edge-corrected kernel
!     estimator of the spatial mark variogram.
!

       subroutine gspcore(x,y,txy,n,s,ns,ks,hs,gsp)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks
       double precision wij,vij,hs,kerns,gspm,gspn,gsp,x,y,txy
       double precision hij,tij,xi,yi,ti,pi,two
       dimension x(n),y(n),txy(n),s(ns),gsp(ns),gspm(ns),gspn(ns)

       gspm=0d0
       gspn=0d0

          two=2d0
          pi=3.141592654d0

       do iu=1,ns
        do i=1,n
         xi=x(i)
         yi=y(i)
         ti=txy(i)
          do j=1,n
           if (j.ne.i) then
            hij=sqrt(((xi-x(j))**two)+((yi-y(j))**two))
            tij=abs(ti-txy(j))
              if (ks.eq.1) then
               kerns=boxkernel((s(iu)-hij)/hs,hs)
                else if (ks.eq.2) then
                 kerns=ekernel((s(iu)-hij)/hs,hs)
                  else if (ks.eq.3) then
                   kerns=qkernel((s(iu)-hij)/hs,hs)
              end if
             if ((kerns.ne.0d0).and.(s(iu).ne.0d0)) then
                    wij=(((tij**two)/two)*kerns)
                    vij=kerns
                    gspm(iu)=gspm(iu)+wij
                    gspn(iu)=gspn(iu)+vij
             end if
           end if
          end do
          end do
           if (s(iu).eq.0d0) then
            gsp(iu)=0
            else
            gsp(iu)=gspm(iu)/gspn(iu)
           end if
          end do

        return

        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     functions called by :
!     -----------------------------------------
!
!     * boxkernel, ekernel, qkernel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------
!
!     boxkernel
!
!--------------------------------------------------------------------

       function boxkernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x, h

       if (abs(x).le.1d0) then
           boxkernel=1d0/2d0
       else
           boxkernel=0d0
       end if
       boxkernel=boxkernel/h

       return
       end

!--------------------------------------------------------------------
!
!     Epanechnikov kernel
!
!--------------------------------------------------------------------

       function ekernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x

       if (abs(x).le.1d0) then
           ekernel=(3d0/4d0)*(1-x**2)
       else
           ekernel=0d0
       end if
       ekernel=ekernel/h

       return
       end

!--------------------------------------------------------------------
!
!     quartic (biweight) kernel
!
!--------------------------------------------------------------------

       function qkernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x, h

       if (abs(x).le.1d0) then
           qkernel=(15d0/16d0)*(1-x**2)**2
       else
           qkernel=0d0
       end if
       qkernel=qkernel/h

       return
       end

