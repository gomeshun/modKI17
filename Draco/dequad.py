from numpy import sinh,cosh,exp,pi,arange,isnan,isinf
def dequad_hinf(func,a,width=5e-3,pN=1000,mN=1000,axis=None):
    '''
    func: func(ndarray_in) = ndarray_out
    axis: define the axis of ndarray_out to use integrate.
    '''
    ts = width*arange(-mN,pN)
    xs = exp(sinh(ts)*pi/2) + a
    ws = width *pi/2 *cosh(ts)*exp(pi/2*sinh(ts))
    fs = func(xs)
    wsfs = ws*fs
    wsfs[isnan(wsfs) & isinf(wsfs)] = 0
    return (ws*fs).sum(axis=axis)