#APL - Assignment - 6
#Done by : Ramanan S
#Date : 10-03-2020

#Importing required libraries
import scipy.signal as sp
from pylab import *
import warnings
warnings.filterwarnings('ignore')

#Plotting the magnitude and phase variation of the given function with frequency
H = sp.lti([1,0.5], [1,1,2.5])
w,s,phi = H.bode()
f = figure(figsize=(12,6))

#Magnitude Vs freq.
ax = f.add_subplot(121)
ax.semilogx(w,s)
ax.set_title('Magnitude')

#Phase Vs freq.
ax2 = f.add_subplot(122)
ax2.semilogx(w,phi)
ax2.set_title('Phase')
show()

#Function which solves the laplace equation for a given 'decay rate' and 'frequency'
def solve(d, w):
    M = [1,0,2.25]
    tmp = polymul([1,d],[1,d])
    tmp = polyadd(tmp,w*w)
    a = polymul(M,tmp)
    X = sp.lti([1,d],a)
    #Converts the X in 'laplace' domain to 'time' domain
    t,x=sp.impulse(X,None,linspace(0,100,1001))
    
    return t, x

#Plotting the solution of spring system for different decay rates keeping the frequency constant
f = figure(figsize=(15,5))
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)

#Plot when decay=0.5
t, x = solve(0.5,1.5)
ax.plot(t,x)
ax.set_title('x(t) when decay is 0.5')
ax.set_xlabel('t')
ax.grid()

#Plot when decay=0.05
t, x = solve(0.05,1.5)
ax2.plot(t,x)
ax2.set_title('x(t) when decay is 0.05')
ax2.set_xlabel('t')
ax2.grid()

show()

#Plotting the solution x(t) for various frequencies keeping the decay rate constant
f2 = figure(figsize=(12,15))
ax = f2.add_subplot(111)

#Freq. Range - [1.4,1.6]
for w in arange(1.4,1.6,0.05):
    t,x = solve(0.05,w)
    ax.plot(t, x, label="$w$ = %g rad/s" % (w))
    
ax.set_title('x(t) with varying frequencies')
ax.set_xlabel('t')
ax.set_ylabel('x(t)')
ax.legend()
ax.grid()

#Function which transforms the given Hs into h(t) - (i.e.) time domain
#Function takes the numerator and denominator of Hs function as arguments
def solve_2(num, den):
    num = poly1d(num)
    den = poly1d(den)
    Hs = sp.lti(num, den)
    
    #converting H(s) to h(t)
    t, h = sp.impulse(Hs, None, linspace(0, 100, 1000))
    return t, h

#Solving for the coupled equation which contains both 'x' and 'y'
t1, x = solve_2([1, 0, 2], [1, 0, 3, 0])
t2, y = solve_2([2], [1, 0, 3, 0])

#Plotting 'x(t)' and 'y(t)' in the same graph
plot(t1, x, 'r', label="$x(t)$")
plot(t2, y, 'g', label="$y(t)$")
title("Time evolution of $x(t)$ and $y(t)$ for $0 \leq t \leq 100$. of Coupled spring system ")
xlabel(r"$t \to $")
ylabel(r"$x(t),y(t) \to $")
ylim((-0.5, 2))
legend()
grid()
show()

#Function which solves the transfer function of the given RLC circuit
def RLC_circuit(R=100, L=10**-6, C=10**-6):
    num = poly1d([1])
    den = poly1d([L*C,R*C,1])
    
    Hs = sp.lti(num, den)
    w, mag, phi = Hs.bode()
    return w, mag, phi, Hs

w, mag, phi, H = RLC_circuit()

# plot Magnitude Response
semilogx(w, mag, 'b', label="Magnitude Response")
legend()
title("Magnitude Response Series RLC network")
xlabel(r"$ \log w \to $")
ylabel(r"|H(jw)| $\to $")
grid()
show()

# Plot of phase response
semilogx(w, phi, 'r', label="Phase Response")
legend()
title("Phase response of Series RLC network")
xlabel(r"$ \log w \to $")
ylabel(r"$ \angle H(j\omega)$ $\to $")
grid()
show()

#Defining the time for which values to be noted
t = linspace(0, 90*pow(10, -3), pow(10, 6))
#Defining the input - Vi
vi = cos(t*pow(10, 3))-cos(t*pow(10, 6))

#Convoluting Vi with H(s) and getting the output Vo
t, vo, svec = sp.lsim(H, vi, t)

#As a ideal low-pass filter lets only if freq. <= 10**4 rad/s
vo_ideal = cos(1e3*t)

f3 = figure(figsize=(8,18))
ax1 = f3.add_subplot(311)
ax2 = f3.add_subplot(312)
ax3 = f3.add_subplot(313)

#Plot of the output Vo Vs t - for large values of t
ax1.plot(t, vo, 'r', label="Output voltage $v_0(t)$ for large time")
ax1.legend()
ax1.set_ylim(-2, 2)
ax1.set_title(r"Output voltage $v_0(t)$  of series RLC network at Steady State")
ax1.set_xlabel(r"$ t \to $")
ax1.set_ylabel(r"$ y(t) \to $")
ax1.grid()

#Ideal Plot of Vo Vs the real plot of Vo and time
ax2.plot(t, vo, 'r', label="Output voltage $v_0(t)$ - zoomed in ")
ax2.plot(t, vo_ideal, 'g', label="Ideal Low Pass filter Output with cutoff at $10^4$")
ax2.set_xlim(0.0505, 0.051)
ax2.set_ylim(0.75, 1.1)
ax2.legend()
ax2.set_title(r"Output voltage $v_0(t)$  Vs Ideal Low pass filter Output")
ax2.set_xlabel(r"$ t \to $")
ax2.set_ylabel(r"$ y(t) \to $")
ax2.grid()

#Plot of the output Vo Vs t - for small values of t
ax3.plot(t, vo, 'r', label="Output voltage $v_0(t)$ : $0<t<30\mu sec$")
ax3.legend()
ax3.set_title(r"Output voltage $v_0(t)$ for $0<t<30\mu sec$")
ax3.set_xlim(0, 3e-5)
ax3.set_ylim(-1e-5, 0.3)
ax3.set_xlabel(r"$ t \to $")
ax3.set_ylabel(r"$ v_0(t) \to $")
ax3.grid()
show()