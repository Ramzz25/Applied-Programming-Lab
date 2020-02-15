#Importing required libraries
import numpy as np
from pylab import *
import scipy.integrate as integrate
import warnings
warnings.filterwarnings('ignore')

#defining the functions e^x and cos(cos(x)) which takes scalar/vector inputs
#both the functions exp(x) and cos(x) are defined in the 'pylab' library which take vector inputs
def f1(x):
    return exp(x)
def f2(x):
    return cos(cos(x))

#range of values of x to be plotted
x = linspace(-2*pi,4*pi,500)
p = 2*pi

#plotting the original function e^x
semilogy(x,f1(x),'--r',label='exp(x)')
#plotting the fourier approximation of e^x (as it is periodic with period 2*pi)
semilogy(x,f1(x%p),label='fourier exp(x)')
xlabel('x')
ylabel('f(x)')
legend()
grid()
show()

#plotting the original function cos(cos(x))
semilogy(x,f2(x),'-bo',label='cos(cos(x))')
#plotting the fourier approximation of cos(cos(x))
semilogy(x,f2(x%p),'r',label='fourier cos(cos(x))')
xlabel('x')
ylabel('f(x)')
legend()
grid()
show()

#function to calculate 'a' coefficients of fourier approximation of e^x
def an_f1(n):
    return (1/pi)*integrate.quad(lambda x: (f1(x))*(cos(n*x)),0,2*pi)[0]
#function to calculate 'b' coefficients of fourier approximation of e^x
def bn_f1(n):
    return (1/pi)*integrate.quad(lambda x: (f1(x))*(sin(n*x)),0,2*pi)[0]

#function to calculate 'a' coefficients of fourier approximation of cos(cos(x))
def an_f2(n):
    return (1/pi)*integrate.quad(lambda x: (f2(x))*(cos(n*x)),0,2*pi)[0]
#function to calculate 'b' coefficients of fourier approximation of cos(cos(x))
def bn_f2(n):
    return (1/pi)*integrate.quad(lambda x: (f2(x))*(sin(n*x)),0,2*pi)[0]

#calculating fourier coefficients of e^x and storing them separately in a list
print('Fourier co-efficients of exp(x)')
af1 = []
bf1 = []
f1_vec = [an_f1(0)/2]
#printing the first 50 'a' and 'b' fourier coefficients
print('a'+str(0)+' = '+str(round(an_f1(0)/2,5)))
for i in range(1,51):
    print('a'+str(i)+' = '+str(round(an_f1(i),5)),end=' ')
    print('b'+str(i)+' = '+str(round(bn_f1(i),5)))
    af1.append(abs(an_f1(i)))
    bf1.append(abs(bn_f1(i)))
    f1_vec.append(an_f1(i))
    f1_vec.append(bn_f1(i))

#calculating fourier coefficients of cos(cos(x)) and storing them separately in a list
print('Fourier co-efficients of cos(cos(x))')
af2 = []
bf2 = []
f2_vec = [an_f2(0)/2]
#printing the first 50 'a' and 'b' fourier coefficients
print('a'+str(0)+' = '+str(round(an_f2(0)/2,5)))
for i in range(1,51):
    print('a'+str(i)+' = '+str(round(an_f2(i),7)),end=' ')
    print('b'+str(i)+' = '+str(round(bn_f2(i),7)))
    af2.append(abs(an_f2(i)))
    bf2.append(abs(bn_f2(i)))
    f2_vec.append(an_f2(i))
    f2_vec.append(bn_f2(i))

#plotting the fourier coefficients of e^x in both semilog and log-log plots
f = figure(figsize=(15,5))
f.suptitle('co-efficients of exp(x)')
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)
ax.semilogy(range(1,51),af1,'ro',label='a')
ax.semilogy(range(1,51),bf1,'bo',label='b')
ax2.loglog(range(1,51),af1,'ro',label='a')
ax2.loglog(range(1,51),bf1,'bo',label='b')
ax.title.set_text('semilog graph')
ax2.title.set_text('log-log graph')
ax.grid()
ax2.grid()
ax.legend()
ax2.legend()
show()

#plotting the fourier coefficients of cos(cos(x)) in both semilog and log-log plots
f = figure(figsize=(15,5))
f.suptitle('co-efficients of cos(cos(x))')
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)
ax.semilogy(range(1,51),af2,'ro',label='a')
ax.semilogy(range(1,51),bf2,'bo',label='b')
ax2.loglog(range(1,51),af2,'ro',label='a')
ax2.loglog(range(1,51),bf2,'bo',label='b')
ax.title.set_text('semilog graph')
ax2.title.set_text('log-log graph')
ax.grid()
ax2.grid()
ax.legend()
ax2.legend()
show()

#Filling the A and b matrices for both the functions
#A matrix would be same for both the functions
x=linspace(0,2*pi,401)
x=x[:-1]
b1=f1(x) # f has been written to take a vector
b2=f2(x)
A=zeros((400,51)) # allocate space for A
A[:,0]=1 # col 1 is all ones
for k in range(1,26):
    A[:,2*k-1]=cos(k*x) # cos(kx) column
    A[:,2*k]=sin(k*x) # sin(kx) column
#endfor
c1=lstsq(A,b1)[0]
c2=lstsq(A,b2)[0]# the ’[0]’ is to pull out the
# best fit vector. lstsq returns a list.
#c1,c2 contains the calculated fourier coefficients of the functions e^x and cos(cos(x)) respectively
#coefficients are calculated by minimising the mean square error

#splitting the 'a' and 'b' fourier coefficients present in c1 
a0_f1_ls = c1[0]
an_f1_ls = []
bn_f1_ls = []
for i in range(1,len(c1)):
    if(i%2==1):
        an_f1_ls.append(abs(c1[i]))
    else:
        bn_f1_ls.append(abs(c1[i]))

#splitting the 'a' and 'b' fourier coefficients present in c2
a0_f2_ls = c2[0]
an_f2_ls = []
bn_f2_ls = []
for i in range(1,len(c2)):
    if(i%2==1):
        an_f2_ls.append(abs(c2[i]))
    else:
        bn_f2_ls.append(abs(c2[i]))

#plotting the fourier coefficents of e^x calculated by the integration method and by the least mean square error method
#plots are plotted in both semilog and log-log axes
f = figure(figsize=(15,5))
f.suptitle('co-efficients of exp(x)')
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)
ax.semilogy(range(1,51),af1,'ro',label='a')
ax.semilogy(range(1,51),bf1,'bo',label='b')
ax.semilogy(range(1,26),an_f1_ls,'go',label='a_pred')
ax.semilogy(range(1,26),bn_f1_ls,'yo',label='b_pred')
ax2.loglog(range(1,51),af1,'ro',label='a')
ax2.loglog(range(1,51),bf1,'bo',label='b')
ax2.loglog(range(1,26),an_f1_ls,'go',label='a_pred')
ax2.loglog(range(1,26),bn_f1_ls,'yo',label='b_pred')
ax.title.set_text('semilog graph')
ax2.title.set_text('log-log graph')
ax.grid()
ax2.grid()
ax.legend()
ax2.legend()
show()

#plotting the fourier coefficents of cos(cos(x)) calculated by the integration method and by the least mean square error method
#plots are plotted in both semilog and log-log axes
f = figure(figsize=(15,5))
f.suptitle('co-efficients of cos(cos(x))')
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)
ax.semilogy(range(1,51),af2,'ro',label='a')
ax.semilogy(range(1,51),bf2,'bo',label='b')
ax.semilogy(range(1,26),an_f2_ls,'go',label='a_pred')
ax.semilogy(range(1,26),bn_f2_ls,'yo',label='b_pred')
ax2.loglog(range(1,51),af2,'ro',label='a')
ax2.loglog(range(1,51),bf2,'bo',label='b')
ax2.loglog(range(1,26),an_f2_ls,'go',label='a_pred')
ax2.loglog(range(1,26),bn_f2_ls,'yo',label='b_pred')
ax.title.set_text('semilog graph')
ax2.title.set_text('log-log graph')
ax.grid()
ax2.grid()
ax.legend()
ax2.legend()
show()

#plotting the original function e^x and fourier approximation of the function taking the first 50 coefficients
semilogy(x,f1(x),'--r',label='exp(x)')
semilogy(x,f1(x%p),label='fourier exp(x)')
semilogy(x,np.matmul(A,c1),'go',label='pred_fun')
xlabel('x')
ylabel('f(x)')
legend()
grid()
show()

#plotting the original function cos(cos(x)) and fourier approximation of the function taking the first 50 coefficients
semilogy(x,f2(x),'--r',label='exp(x)')
semilogy(x,f2(x%p),label='fourier exp(x)')
semilogy(x,np.matmul(A,c2),'go',label='pred_fun')
xlabel('x')
ylabel('f(x)')
legend()
grid()
show()

#calculating the deviation of the fourier coefficients from the true value for the function e^x
dev_an_f1 = []
dev_bn_f1 = []
for i in range(25):
    dev_an_f1.append(abs(af1[i]-an_f1_ls[i]))
    dev_bn_f1.append(abs(bf1[i]-bn_f1_ls[i]))

#calculating the deviation of the fourier coefficients from the true value for the function cos(cos(x))
dev_an_f2 = []
dev_bn_f2= []
for i in range(25):
    dev_an_f2.append(abs(af2[i]-an_f2_ls[i]))
    dev_bn_f2.append(abs(bf2[i]-bn_f2_ls[i]))

#plotting the deviation of coefficients vs n for e^x
f = figure(figsize=(15,5))
f.suptitle('Deviation of fourier coefficients of function exp(x)')
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)
ax.plot(range(1,26),dev_an_f1,'-b')
ax2.plot(range(1,26),dev_bn_f1,'-b')
ax.title.set_text('Deviation of an')
ax2.title.set_text('Deviation of bn')
show()
print('For the func. exp(x)')
print('Max deviation of an = '+str(max(dev_an_f1)))
print('Max deviation of bn = '+str(max(dev_bn_f1)))

#plotting the deviation of coefficients vs n for cos(cos(x))
f = figure(figsize=(15,5))
f.suptitle('Deviation of fourier coefficients of function cos(cos(x))')
ax = f.add_subplot(121)
ax2 = f.add_subplot(122)
ax.plot(range(1,26),dev_an_f2,'-b')
ax2.plot(range(1,26),dev_bn_f2,'-b')
ax.title.set_text('Deviation of an')
ax2.title.set_text('Deviation of bn')
show()
print('For the func. cos(cos(x))')
print('Max deviation of an = '+str(max(dev_an_f2)))
print('Max deviation of bn = '+str(max(dev_bn_f2)))