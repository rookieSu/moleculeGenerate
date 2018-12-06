#!/usr/bin/python

atomDict = {}
bondAdj = []
atomNumber = 0
bondNumber = 0
coordinates = []
replace = {}

def read_file(filename):
    global atomNumber, bondNumber, atomDict, bondAdj, coordinates
    with open(filename,'r') as f:
        lines = f.readlines()
        atomNumber, bondNumber = eval(lines[3].split()[0]), eval(lines[3].split()[1])
        for i in range(0,atomNumber):
           x, y, z = eval(lines[i+4].split()[0]), eval(lines[i+4].split()[1]), eval(lines[i+4].split()[2])
           atomDict[i+1] = lines[i+4].split()[3]
           coordinates.append([x,y,z])
           temp = []
           for j in range(0, atomNumber):
               temp.append(0)
           bondAdj.append(temp)
        for i in range(atomNumber+4,len(lines)-1):
            row, column = eval(lines[i].split()[0]), eval(lines[i].split()[1])
            if eval(lines[i].split()[2]) == 2:
                bondAdj[row-1][column-1] = 2 
                bondAdj[column-1][row-1] = 2
            else: 
                bondAdj[row-1][column-1] = 1
                bondAdj[column-1][row-1] = 1  
    f.close()

def is_replace():
    global replace, bondAdj
    bondCounter = [sum(bond) for bond in bondAdj]
    print(bondCounter)
    for (count, index) in zip(bondCounter,range(0, atomNumber)):
        if count == 5:
            pass
        elif count == 4:
            if atomDict[index+1] == 'N':
                replace[index+1] = 1
            else:replace[index+1] = 0
        elif count == 2:
            if atomDict[index+1] == 'N':
                replace[index+1] = 1
            elif atomDict[index+1] == 'C':
                replace[index+1] = 2
            else:replace[index+1] = 0
        elif count == 3:
            if atomDict[index+1] == 'C':
                replace[index+1] = 1
            else:replace[index+1] = 0
        elif count == 1:
            if atomDict[index+1] == 'O':
                replace[index+1] = 1
            elif atomDict[index+1] == 'N':
                replace[index+1] = 2
            elif atomDict[index+1] == 'C':
                replace[index+1] = 3
            else:replace[index+1] = 0
        else:
            replace[index+1] = 0
                
if __name__ == "__main__":
    read_file('p35_173.mol')
    is_replace()
    print(replace)

