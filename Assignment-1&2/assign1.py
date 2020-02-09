#Have to give argument in the terminal 

def create_list(filename):
    try:
        f = open(filename)
    except:
        print("Incorrect File Name")
        return [False]

    data = f.read()
    f.close()
        
    l = data.splitlines()
    arr = []
    data = []
    if(('.circuit' in l) and ('.end' in l) ):
        for i in l:
            arr.append(i.split(' '))

        arr_n = arr[::-1]
        chk = True
        start = False
        i = 0
        while(chk):
            if(arr_n[i][0] == '.ac' or arr_n[i][0] == '.end'):
                start = True
            if(arr_n[i][0] == '.circuit'):
                chk = False
                break
            if(start):
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