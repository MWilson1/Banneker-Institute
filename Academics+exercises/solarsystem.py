#!/usr/bin/env python

'''Some fun with object-oriented programming: make a file called ``solarsystem.py`` that defines ``Star``, ``Planet``, and ``System`` objects that enable the provided notebook file to work as demonstrated.  If you finish this, feel free to add more functionality to these objects! '''

class Star(object):
    def __init__(self, Name='Sun', Mass=1.00, Radius=1.00):
        self.Name = Name
        self.Mass = Mass
        self.Radius = Radius
        print('Sol: '+str(self.Mass)+' M_'+self.Name+', '+str(self.Radius)+' R_'+self.Name)
        #return
