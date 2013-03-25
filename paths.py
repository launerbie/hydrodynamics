import os

''' Class that contains some useful directories. Must be updated 
manually when directory names change.'''

class Directories(object):
    def __init__(self):
        self.root = os.path.dirname(__file__ )
        self.movies = self.root+"/movies"
        self.hydroresults  = self.root+"hydroresults"
        self.images = self.root+"images"
        self.bodies = self.root+"bodies"
        self.plots = self.root+"plots"

paths = Directories()


