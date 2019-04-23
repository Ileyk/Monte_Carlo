      subroutine density(x,y,z,rho)

      implicit none

      real*8 x,y,z,rho

      real*8 w,w2,r,r2,h0,h

      real*8 zf

c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f
      w2=x*x+y*y
      w=sqrt(w2)
      r2=w2+z*z
      r=sqrt(r2)

c***** Set up uniform density sphere within the grid
      if(r.gt.1.) then
        rho=0.
      else
        rho=1.
      endif

      return
      end

