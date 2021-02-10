from numpy import *
################################################################################
class  BCTab:
    pass
#===============================================================================
    def __init__(self,tab=None):
        self.CH = {}
        self.NH = {}
        self.CHnu = {}
        self.NHnu = {}
        self.CHHz = {}
        self.NHHz = {}
        self.NH = {}
        self.N = {}
        self.H = {}
        self.C = {}
        self.num = {}
        if tab:
            self.readTab(tab)
#===============================================================================    
    def readTab(self,tab):
        lines = map(str.strip,open(tab).readlines()[1:])
        Tab1 = {x.split()[0]:float(x.split()[1]) for x in lines }
        Tab2 = {x.split()[0]:float(x.split()[2]) for x in lines }
        try:
            self.num = {x.split()[0]:int(x.split()[3]) for x in lines }
        except IndexError:
            pass
        bds = sorted(Tab1.keys())
        if all(array(Tab2.values())>9):
            self.N = {x:Tab1[x] for x in bds}
            self.H = {x:Tab2[x] for x in bds}
            self.NH = {x:[Tab1[x],Tab2[x]] for x in bds}
            self.NHHz = {x:[Tab1[x]*60,Tab2[x]*600] for x in bds}
            for x in bds:
                if x[0] in ['G','A']:
                    self.NHnu[x] = {'N1':Tab1[x],'H1':Tab2[x]} 
                if x[0] in ['U','C']:
                    self.NHnu[x] = {'N3':Tab1[x],'H3':Tab2[x]}
        elif all(array(Tab2.values())<8.6) and all(6.4<array(Tab2.values())) and\
           all(array(Tab1.values())<148) and all(134<array(Tab1.values())):
            self.C = {x:Tab1[x] for x in bds}
            self.H = {x:Tab2[x] for x in bds}
            self.CH = {x:[Tab1[x],Tab2[x]] for x in bds}
            self.CHHz = {x:[Tab1[x]*150,Tab2[x]*600] for x in bds}
            for x in bds:
                if x[0] in ['G','A']:
                    self.CHnu[x] = {'C8':Tab1[x],'H8':Tab2[x]} 
                if x[0] in ['U','C']:
                    self.CHnu[x] = {'C6':Tab1[x],'H6':Tab2[x]}
            
        elif all(array(Tab2.values())<6.4) and all(4.8<array(Tab2.values())) and\
           all(array(Tab1.values())<98) and all(94<array(Tab1.values())):
            self.C = {x:Tab1[x] for x in bds}
            self.H = {x:Tab2[x] for x in bds}
            self.CH = {x:[Tab1[x],Tab2[x]] for x in bds}
            self.CHHz = {x:[Tab1[x]*150,Tab2[x]*600] for x in bds}
            for x in bds:
                self.CHnu[x] = {"C1'":Tab1[x],"H1'":Tab2[x]}
        elif all(array(Tab2.values())<8.4) and all(6.0<array(Tab2.values())) and\
           all(148<array(Tab1.values())) and all(array(Tab1.values())<158):
            self.C = {x:Tab1[x] for x in bds}
            self.H = {x:Tab2[x] for x in bds}
            self.CH = {x:[Tab1[x],Tab2[x]] for x in bds}
            self.CHHz = {x:[Tab1[x]*150,Tab2[x]*600] for x in bds}
            for x in bds:
                self.CHnu[x] = {"C2":Tab1[x],"H2":Tab2[x]}


      
        

            

