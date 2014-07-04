__author__='pier-luc'
__date__ = "2014-06-05"


import Error


def checkFiles(listOfFile):
    for i in listOfFile:
        if not(isValid(i)):
            Error.error("File"+str(i)+"can't be opened")

def is_valid(file_name):
    valid = 1
    try:
        a = open(file_name, 'r')
        a.close()
    except:
        valid = 0
    return valid

def countLines(file_name):
    lines = 0
    try:
        for i in open(file):
            lines += 1
    except:
        lines = 0
    if lines == 0:
        Error.error("File "+str(file)+" is empty")
    return lines
