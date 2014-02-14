	subroutine loren(gam,bx,by,bz,e,x,y,z,ef1,pxf,pyf,pzf,pf1)

	implicit none

	real*8	gam,bx,by,bz,e,x,y,z
	real*8  gam1,pxf,pyf,pzf,pf1,ef1

	gam1=gam**2/(1.+gam)
	ef1 = gam*(e-bx*x-by*y-bz*z)
	pxf = (1+gam1*bx**2)*x + gam1*bx*(by*y+bz*z) - gam*bx*e
	pyf = (1+gam1*by**2)*y + gam1*by*(bx*x+bz*z) - gam*by*e
	pzf = (1+gam1*bz**2)*z + gam1*bz*(by*y+bx*x) - gam*bz*e
	pf1 = sqrt(pxf**2+pyf**2+pzf**2)

	return
	end
