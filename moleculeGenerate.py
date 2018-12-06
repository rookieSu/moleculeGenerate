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
           bondAdj.append([0 for x in range(0, atomNumber)])
        for i in range(atomNumber+4,len(lines)-1):
            row, column = eval(lines[i].split()[0])-1, eval(lines[i].split()[1])-1
            if eval(lines[i].split()[2]) == 2:
                bondAdj[row][column] = 2 
                bondAdj[column][row] = 2
            else: 
                bondAdj[row][column] = 1
                bondAdj[column][row] = 1  
    f.close()

def is_replace():
    global replace, bondAdj
    bondCounter = [sum(bond) for bond in bondAdj]
    for (count, index) in zip(bondCounter,range(1, atomNumber+1)):
        if count == 5:
            pass
        elif count == 4:
            if atomDict[index] == 'N':
                replace[index] = 1
            else:replace[index] = 0
        elif count == 2:
            if atomDict[index] == 'N':
                replace[index] = 1
            elif atomDict[index] == 'C':
                replace[index] = 2
            else:replace[index] = 0
        elif count == 3:
            if atomDict[index] == 'C':
                replace[index] = 1
            else:replace[index] = 0
        elif count == 1:
            if atomDict[index] == 'O':
                replace[index] = 1
            elif atomDict[index] == 'N':
                replace[index] = 2
            elif atomDict[index] == 'C':
                replace[index] = 3
            else:replace[index] = 0
        else:
            replace[index] = 0
                
if __name__ == "__main__":
    read_file('ben.mol')
    is_replace()
    print(replace)

