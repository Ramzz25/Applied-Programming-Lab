#Assignment-2 EE2703 APL
#Done by : Ramanan S
#Date : 31-01-2020

#Importing required libraries
import sys 
import numpy as np
from math import sin, cos

#checking if the given netlist has ac components or dc components
is_ac = False
#Function to read and return the data in the given filename
def create_list(filename):
    global is_ac
    try:
        #opening the file
        f = open(filename)
    except:
        print("Incorrect File Name")
        return [False]

    #reading the file
    data = f.read()
    f.close()
        
    l = data.splitlines()
    arr = []
    data = []
    #checking if its a proper netlist file
    if(('.circuit' in l) and ('.end' in l) ):
        for i in l:
            arr.append(i.split(' '))
        #reversing the contents
        arr_n = arr[::-1]
        #variables to help loop the through the read data
        chk = True
        start = False
        i = 0
        while(chk):
            if(arr_n[i][0] == '.ac'):
                is_ac = True    #checking if the circuit is ac
                start = True
            if(arr_n[i][0] == '.end'):
                start = True
            if(arr_n[i][0] == '.circuit'):
                chk = False
                break
            if(start):
                #collecting the required data inside the file
                data.append([])
                temp = list(reversed(arr_n[i]))
                end = len(temp)
                for l in range(len(temp)):
                    if(arr_n[i][l]=='' or arr_n[i][l]==' '):
                        end = l
                        break
                for j in range(len(temp)-end,len(temp)):
                    data[-1].append(temp[j])
                data[-1] = list(reversed(data[-1]))
            i += 1
        #Info from the file is extracted and stored in the list 'data'
        return [True, data]
    else:
        print('Doesn\'t contain circuit info.' )
        return [False]

#function to parse the values of different elements
def par_value(l):
    pos = l[-1].find('e')               #finding if 'e' is present in the value to be parsed
    if(pos!=-1):
        v1 = float(l[-1][:pos])            #if 'e' is present the value is 'x*10^y'
        v2 = float(l[-1][pos+1:])
        return v1*(10**v2)
    else:
        return float(l[-1])

#class defined to hold the characteristics of different components like resistor, capacior, inductor
class component:
    def __init__(self,l):
        self.name = l[0]
        self.nodes = []
        for i in range(1,3):
            self.nodes.append(d[l[i]])
        self.value = l[3]

#class defined to hold characteristics of independent sources (both ac and dc)
class source:
    def __init__(self,l):
        self.name = l[0]
        self.nodes = []
        for i in range(1,3):
            self.nodes.append(d[l[i]])
        self.type = l[3]
        if(self.type=='ac'):
            self.value = (par_value([l[4]]))/2
            self.phase = (float(l[5])*(np.pi))/(180)
        else:
            self.value = par_value([l[4]])
            self.phase = 'na'

fname = sys.argv[1]
l = create_list(fname)
d = {'GND':0, 0:'GND'}          #d is a dictionary to hold the no. of nodes
if(l[0]):                       #check if the file contains data in a proper manner
    data = list(l[1])
    if(is_ac):                  #if the circuit is ac, this finds frequency and omega of the ac source
        freq = par_value([data[0][2]])
        w = 2*(np.pi)*freq
        data = data[2:]
    else:
        data = data[1:]
    
    count = 1               #numbering the different type of nodes for calculation purposes
    for i in data:
        if(i[1] not in d):
            d[i[1]] = count
            d[count] = i[1]
            count += 1
        if(i[2] not in d):
            d[i[2]] = count
            d[count] = i[2]
            count += 1
    
    comp_list = []                  #list to hold details of all components in the circuit like R,L,C
    source_list = []                #list to hold details of all independent in the circuit sources
    
    #separting sources and components for further processing
    for j in range(len(data)):
        if(data[j][0][0]=='V' or data[j][0][0]=='I'):
            if('ac' not in data[j] and 'dc' not in data[j] ):
                t = data[j][-1]
                del data[j][-1]
                data[j].append('dc')
                data[j].append(t)
            source_list.append(source(data[j]))
        else:  
            data[j][-1] = par_value(data[j])
            comp_list.append(component(data[j]))

    mat_len = len(d)//2         
    nodes = mat_len                 #mat_len contains the dimensions of the 'M' and 'b' matrices to be formed
    vcs = {}                        #storing the indices of current thro volt sources in the 'M' matrix
    for i in source_list:
        if(i.name[0]=='V'):
            mat_len += 1
            vcs[i.name] = mat_len-1
            d[mat_len-1] = i.name
    
    M = np.zeros((mat_len,mat_len),dtype=complex)
    b = np.zeros(mat_len,dtype=complex)

    M[d['GND']][d['GND']] += complex(1,0)       #Eqn. Vgnd = 0
    for j in source_list:
        if(j.name[0]=='V'):                     #Eqn. for independent voltage sources both ac and dc
            if(j.type=='ac'):
                b[vcs[j.name]] += complex(j.value*cos(j.phase),j.value*sin(j.phase))
            else:
                b[vcs[j.name]] += complex(j.value,0)
            M[vcs[j.name]][j.nodes[0]] -= complex(1,0)
            M[vcs[j.name]][j.nodes[1]] += complex(1,0)
            M[j.nodes[0]][vcs[j.name]] += complex(1,0)
            M[j.nodes[1]][vcs[j.name]] -= complex(1,0)
        elif(j.name[0]=='I'):                   #Eqn. for independent current sources both ac and dc
            if(j.type=='ac'):
                b[j.nodes[0]] += complex(j.value*cos(j.phase),j.value*sin(j.phase))
                b[j.nodes[1]] -= complex(j.value*cos(j.phase),j.value*sin(j.phase))
            else:
                b[j.nodes[0]] += complex(j.value,0)
                b[j.nodes[1]] -= complex(j.value,0)

    for j in comp_list:
        if(j.name[0]=='R'):
            #Equation for the nodes having resistors
            M[j.nodes[0]][j.nodes[0]] += complex((1/j.value),0)
            M[j.nodes[1]][j.nodes[1]] += complex((1/j.value),0)
            M[j.nodes[0]][j.nodes[1]] -= complex((1/j.value),0)
            M[j.nodes[1]][j.nodes[0]] -= complex((1/j.value),0)
        elif(j.name[0]=='L'):
            #Equation for the nodes having inductors
            M[j.nodes[0]][j.nodes[0]] += complex(0,-1/(j.value*w))
            M[j.nodes[1]][j.nodes[1]] += complex(0,-1/(j.value*w))
            M[j.nodes[0]][j.nodes[1]] -= complex(0,-1/(j.value*w))
            M[j.nodes[1]][j.nodes[0]] -= complex(0,-1/(j.value*w))
        elif(j.name[0]=='C'):
            #Equation for the nodes having capacitors
            M[j.nodes[0]][j.nodes[0]] += complex(0,w*j.value)
            M[j.nodes[1]][j.nodes[1]] += complex(0,w*j.value)
            M[j.nodes[0]][j.nodes[1]] -= complex(0,w*j.value)
            M[j.nodes[1]][j.nodes[0]] -= complex(0,w*j.value)

    #finding the x matrix
    x = np.linalg.solve(M, b)
    #printing the solution - (x matrix)
    for i in range(len(x)):
        if(i<nodes):
            print('Voltage at node '+ str(d[i]) + ' : ',end=' ')
            print(round(x[i], 6))
        else:
            print('Curent through voltage source ' + str(d[i]) + ' : ',end=' ')
            print(round(x[i], 6))