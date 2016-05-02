#define AM_2d(g,p,d)                    sqrt(g*p/d)
#define XS_2d(d,dtwo,done,ccl,uxr,uxy)  d*((1.0/dtwo)-(1.0/done))*ccl*(1.0-ccl)*(uxr+uyr)
