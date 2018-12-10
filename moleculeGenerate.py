#!/usr/bin/python
import os
import copy
atomDict = {}
bondAdj = []
atomNumber = 0
bondNumber = 0
coordinates = []
rePosition = {}
def read_molfile(molname):
    global atomNumber, bondNumber, atomDict, bondAdj, coordinates
    with open(molname,'r') as f:
        lines = f.readlines()
        atomNumber, bondNumber = eval(lines[3].split()[0]), eval(lines[3].split()[1])
        for i in range(0,atomNumber):
           x, y, z = eval(lines[i+4].split()[0]), eval(lines[i+4].split()[1]), eval(lines[i+4].split()[2])
           atomDict[i+1] = lines[i+4].split()[3]
           coordinates.append([x,y,z])
           bondAdj.append([0 for x in range(0, atomNumber)])
        for i in range(atomNumber+4,len(lines)-1):
            row, column = eval(lines[i].split()[0])-1, eval(lines[i].split()[1])-1
            bondValue = eval(lines[i].split()[2])
            bondAdj[row][column], bondAdj[column][row] = bondValue, bondValue
    f.close()

def is_replace():
    global rePosition, bondAdj
    bondCounter = [sum(bond) for bond in bondAdj]
    for (count, index) in zip(bondCounter,range(1, atomNumber+1)):
        if count == 5:
            pass
        elif count == 4:
            if atomDict[index] == 'N':
                rePosition[index] = 1
            else:rePosition[index] = 0
        elif count == 2:
            if atomDict[index] == 'N':
                rePosition[index] = 1
            elif atomDict[index] == 'C':
                rePosition[index] = 2
            else:rePosition[index] = 0
        elif count == 3:
            if atomDict[index] == 'C':
                rePosition[index] = 1
            else:rePosition[index] = 0
        elif count == 1:
            if atomDict[index] == 'O':
                rePosition[index] = 1
            elif atomDict[index] == 'N':
                rePosition[index] = 2
            elif atomDict[index] == 'C':
                rePosition[index] = 3
            else:rePosition[index] = 0
        else:
            rePosition[index] = 0
    print('max relpace number is: {}'.format(sum(rePosition.values())))

def read_sub(subname):
    subSmiles = []
    with open(subname,'r') as f:
        lines = f.readlines()
        for index, item in enumerate(lines):
            subSmiles.append('I'+item.strip())
    f.close()
    replace(subSmiles)

def replace(subSmiles):
    reFileName = []
    count = 1
    if os.path.exists('mol'):
        os.system("rm mol/*")
    else:
        os.system("mkdir mol")
    for index, item in enumerate(subSmiles):
        os.system("obabel -:\"{}\" -O {}.mol --gen3D".format(item,index))
        reFileName.append("{}.mol".format(index))
    for name in reFileName:
        with open(name, 'r') as f:
            lines = f.readlines()
            subAtom, subBond = eval(lines[3].split()[0]), eval(lines[3].split()[1])
            subAtomInfo, subBondInfo = lines[5:4+subAtom], lines[4+subAtom:4+subAtom+subBond]
            with open('p35n.mol', 'r') as molFile:
                molLines = molFile.readlines()
                for key, value in zip(rePosition.keys(), rePosition.values()):
                    lines = copy.deepcopy(molLines)
                    if value > 0:
                        lines[4+atomNumber:4+atomNumber] = subAtomInfo
                        results = []
                        for index, item in enumerate(subBondInfo):
                            if index == 0:
                                item = item.split()
                                if eval(item[0]) > eval(item[1]):
                                    item[0] = str(eval(item[0])-1)
                                    item[1] = str(key)
                                else:
                                    item[0] = str(key)
                                    item[1] = str(atomNumber+eval(item[1])-1)
                            else:
                                item = item.split()
                                item[0] = str(atomNumber+eval(item[0])-1)
                                item[1] = str(atomNumber+eval(item[1])-1)
                            results.append(" "+" ".join(item)+"\n")
                        changeAtomBond = lines[3].split()
                        changeAtomBond[0], changeAtomBond[1] = str(atomNumber+subAtom-1), str(bondNumber+subBond)
                        lines[3]="  ".join(changeAtomBond)+"\n"
                        lines[-1:-1] = results
                        with open('conf{}.mol'.format(count),'w') as conf:
                            conf.write("".join(lines))
                        conf.close()
                        os.system("obabel conf{}.mol -O conf{}.mol --gen3D".format(count,count))
                        os.system("mv conf{}.mol mol".format(count))
                        count+=1
            molFile.close()
        f.close()

def inchi_filter():
    os.chdir(os.getcwd()+'/mol')
    root = os.getcwd()
    filenames = os.walk(root).__next__()[2]
    currentInchi = ""
    inchiList = []
    for index, name in enumerate(filenames):
        if name.split(".")[1] == "mol":
            currentInchi = os.popen("obabel -i mol {} -o inchi".format(name)).read()
            if index != 0:
                if currentInchi in inchiList:
                    os.system("rm {}".format(name))
                    print("rm {} !!!".format(name))
                else:
                    inchiList.append(currentInchi)
            else:
                inchiList.append(currentInchi)

def smiles_filter():
    root = os.getcwd()
    filenames = os.walk(root).__next__()[2]
    currentSmiles = ""
    smilesList = []
    for index, name in enumerate(filenames):
        if name.split(".")[1] == "mol":
            currentSmiles = os.popen("obabel -i mol {} -o smi".format(name)).read()
            if index != 0:
                if currentSmiles in smilesList:
                    os.system("rm {}".format(name))
                    print("rm {}!!!".format(name))
                else:
                    smilesList.append(currentSmiles)
            else:
                smilesList.append(currentSmiles)

if __name__ == "__main__":
    read_molfile('p35n.mol')
    is_replace()
    read_sub('sub.txt')
    inchi_filter()
    smiles_filter()
