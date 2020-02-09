#Assignment-1 EE2073
#Date : 26-01-2020
#Done by : Ramanan S

#Have to give argument in the terminal 
import sys
filename = sys.argv[1]                                      #gets the name of the file to be checked
try:
    f = open(filename)
except:
    print("Incorrect File Name")

l = f.read()                                                #reads the file, if present
f.close()
    
l = l.splitlines()                                          #Split the lines using '\n' as delimiter
arr = []
data = []
if(('.circuit' in l) and ('.end' in l) ):                   #check if '.circuit' and '.end' is present in the file or not
    for i in l:
        arr.append(i.split(' '))
    arr_n = arr[::-1]                                       #reveres the list of lines
    chk = True
    start = False
    i = 0
    while(chk):
        if(arr_n[i][0] == '.end'):
            start = True
            i += 1
            continue
        if(arr_n[i][0] == '.circuit'):
            chk = False
            break

        if(start):
            data.append([])
            temp = list(reversed(arr_n[i]))                 #printing lines in a reversed manner
            end = len(temp)
            for l in range(len(temp)):
                if(arr_n[i][l]=='' or arr_n[i][l]==' ' or arr_n[i][l]=='#'):
                    end = l
                    break
            for j in range(len(temp)-end,len(temp)):
                print(temp[j],end =' ')
                data[-1].append(temp[j])
            print()
        i += 1
    #Info from the file is extracted and stored in the list 'data'
else:
    print('Doesn\'t contain circuit info.' )