from numpy import sinh,cosh,exp,pi,arange
def dequad_hinf(func,a,width=5e-3,pN=1000,mN=1000):
    ts = width*arange(-mN,pN)
    xs = exp(sinh(ts)*pi/2) + a
    ws = width *pi/2 *cosh(ts)*exp(pi/2*sinh(ts))
    return (ws*func(xs)).sum()