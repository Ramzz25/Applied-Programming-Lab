from pylab import *
import matplotlib.pyplot as plt
import scipy.special as sp
import scipy.linalg as linalg
import time

#Function g(t,A,B) - returns a list with respect to A, B values
def g(t,A=1.05,B=-0.105):
    y=A*sp.jv(2,t)+B*t
    return y

N=101 # no of data points
k=9 # no of sets of data with varying noise
# generate the data points and add noise

t=linspace(0,10,N) # t vector
y=1.05*sp.jv(2,t)-0.105*t # f(t) vector
Y=meshgrid(y,ones(k),indexing='ij')[0] # make k copies
scl=logspace(-3,-1,k) # noise stdev
n=dot(randn(N,k),diag(scl)) # generate k vectors
yy=Y+n # add noise to signal

#writing the values of legend for the plot
arr = []
for i in range(9):
    arr.append('\u03C3'+str(i+1)+' = '+str(round(scl[i],3)))
arr.append('True Value')

#true value of the function without noise as calculated by the function g 
gt = np.array(g(t))[:, None]

#function to solve question 4
def Q4():
    global yy, gt, t
    #gt contains true value of the function without noise
    y = np.concatenate((yy,gt),axis=1)
    #plotting the noisy data and true value of the function
    lines = plot(t, y)
    legend(lines, arr)
    xlabel(r'$t$',size=20)
    ylabel(r'$f(t)+n$',size=20)
    title(r'Figure 0')
    grid(True)
    #savetxt("fitting.dat",c_[t,yy]) # write out matrix to file
    show()
    return

#function to solve question 5
def Q5():
    global yy, gt, t
    data = yy
    #plotting error bar
    errorbar(t[::5],yy[:,8][::5],0.1,fmt='ro',label='err bar')
    #plotting the true value
    plot(t,gt,label='f(t)')
    xlabel(r'$t$',size=15)
    title("Error Bar")
    #adding annotations
    annotate("Noisy Data",(t[0], data[0, 0]), xytext=(40, -40), textcoords="offset points", arrowprops={"arrowstyle": "->"})
    annotate("True Data",(t[0], g(t[0])), xytext=(-20, 35), textcoords="offset points", arrowprops={"arrowstyle": "->"})
    annotate("Noise in data(Red line)",(t[0], data[0, 0]-0.06), xytext=(-20, -60), textcoords="offset points", arrowprops={"arrowstyle": "->"})
    legend()
    grid()
    show()
    return

#function to solve question 6
def Q6():
    global yy, gt, t
    M = np.zeros((101,2))
    p = np.zeros((2,1))
    #initiating M matrix with corresponding values
    for i in range(101):
        M[i,0] = sp.jv(2,t[i]) 
        M[i, 1] = t[i]
    #initiating p matrix with corresponding values
    p[0,0] = 1.0
    p[1,0] = 0.5

    #calculating the ans using matrix multiplication of M and p matrix
    ans = np.matmul(M,p)
    g_ans = np.array(g(t,1.0,0.5))[:,None]

    chk_ans = np.array_equal(ans, g_ans)
    print(chk_ans)
    return M

#function to solve question 7
def Q7():
    global yy, gt, t
    #intialising some values of A and B
    A_chk = linspace(0,2,21)
    B_chk = linspace(-0.2,0,21)
    e = np.zeros((len(A_chk), len(B_chk)))

    #calculating error value for each value of A and B
    for i in range(len(A_chk)):
        for j in range(len(B_chk)):
            e[i, j] = mean(square(g(t,A_chk[i],B_chk[j])-yy[:,0]))

    x, y = np.meshgrid(A_chk, B_chk)
    #plotting contours of error
    co = contour(x, y, e)
    plot([1.05],[-0.105],'ro')
    annotate("True Value\n of Minima", (1.05, -0.105),xytext=(-50, -30), textcoords="offset points")
    xlabel('A')
    ylabel('B')
    clabel(co, co.levels[:5], inline=1, fontsize=8)
    title('Contours of e')
    grid()
    show()

    return

#function to solve question 9
def Q9():
    global yy, gt, t
    M = Q6()
    AB = []
    err_A = []
    err_B = []
    #solving for A and B using different noisy data
    #solving is done in a way such that MSE is minimum
    for i in range(9):
        tmp = linalg.lstsq(M, yy[:,i])
        AB.append(tmp[0])
        err_A.append(square(tmp[0][0]-1.05))
        err_B.append(square(tmp[0][1]+0.105))
    #plotting error in A and in normal scale
    plot(scl,err_A,'--ro',label='Error in A')
    plot(scl,err_B,'--go',label='Error in B')
    xlabel('\u03C3')
    ylabel('Error')
    grid()
    legend()
    show()

    #plotting error in A and in log-log scale
    x = logspace(-3,-1,9)
    loglog(x,err_A,'ro')
    loglog(x,err_B,'go')
    errorbar(x,err_A,0.0005,fmt='ro',label='errbar_A')
    errorbar(x,err_B,0.0005,fmt='go',label='errbar_B')
    legend()
    grid()
    show()
    return

if(__name__=="__main__"):
    Q4()
    Q5()
    Q6()
    Q7()
    Q9()