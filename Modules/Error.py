__author__='pier-luc'


import sys

"""
Source:
https://docs.python.org/2/tutorial/errors.html
"""

class Error(Exception):
    """
    Base exception class
    """
    pass

class LogicError(Error):
    """
    Class for Logic error. Its your job to be clear in the message.
    """
    def __init__(self, message):
        self.msg = message

    def __str__(self):
        return repr(self.msg)







#TODO: removed these functions
def warning(statement):
    sys.stderr.write("Warning:"+str(statement)+"\n")

def error(statement):
    """

    :rtype : object
    """
    sys.stderr.write("Error:"+str(statement)+"\n")
    sys.exit(1)