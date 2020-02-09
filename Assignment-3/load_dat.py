from pylab import *
import scipy.special as sp 
data = loadtxt("fitting.dat")

def g(t,A=1.05,B=-0.105):
    y=1.05*sp.jv(2,t)-0.105*t
    return y