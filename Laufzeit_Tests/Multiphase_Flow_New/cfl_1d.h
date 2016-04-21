#define AM_1d(g,p,d)                sqrt(g*p/d)
#define XS_1d(d,dtwo,done,ccl,uxr)  d*((1.0/dtwo)-(1.0/done))*ccl*(1.0-ccl)*uxr
