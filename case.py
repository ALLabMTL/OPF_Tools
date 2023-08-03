import json
import numpy as np
import cvxpy as cp
import pandas as pd

class Case:
    '''Defines a Case object containing all relevant data from a MatPower case.
    Principal data entries are:
    mva -> base MVA for the case
    bus -> bus information in the network
    gen -> Generator information and location
    branch -> Connection information on branch connections
    genCost -> Cost coefficients for generators
    '''
    def __init__(self, dct=None):
        '''Takes a dictionary obtained from .txt file'''
        if dct is None:
            self.mva = 100
            self.bus = None
            self.gen = None
            self.branch = None
            self.gencost = None
            self.N = 0

        if isinstance(dct, str):
            dct = self.getDct(dct)

        if isinstance(dct, dict):
            self.mva = dct['baseMVA']
            self.bus = dct['bus']
            self.gen = dct['gen']
            self.branch = dct['branch']
            self.gencost = dct['gencost']
            self.N = len(self.bus)
            if not isinstance(self.gen[0], list):
                self.gen = [self.gen]
            if not isinstance(self.gencost[0], list):
                self.gencost = [self.gencost]

    def getDct(self, filename):
        "Load json dictionary from file"
        with open(filename, 'r') as file:
            dct = json.loads(file.read())
            file.close()
        return dct     


class UsefulCase(Case):
    '''Case object specifically geared towards solving the OPF problem
    Fields corespond to the format used in the tutorials of the ELE8350 class

    adj -> Adjecency matrix of the network contains admittance of lines
    smax -> Matrix containing the maximum apparent power in lines
    genData -> Maximum and minimum active and reactive generation
    loadData -> Load numbers
    cost -> Cost coefficients for the generators

    '''
    def __init__(self, dct=None):
        '''Intatiate the object
        Input: dct -> Either loaded json dictionary from .txt file or the filename to load
        
        If dct = None, intatiates empty instance'''
        super().__init__(dct)
        if dct is not None:
            self.adj = self.getAdj()
            self.smax = self.getsmax()
            self.genData = self.getGenData()
            self.loadData = self.getLoadData()
            self.cost = self.getCost()
            self.vlim = self.getVlim()


    def getAdj(self):
        '''Generates adjacency matrix for the case
        
        Element (i,j) is admittance of line between buses i -> j
        Zero when no line exists between two buses
        '''
        # adjecency matrix
        ans = np.zeros((self.N, self.N), dtype=complex)
        for line in self.branch:
            i, j = line[0] - 1, line[1] - 1
            r, x = line[2], line[3]
            gam = 1/(r + 1.0j*x)
            if i<self.N and j < self.N:
                # symmetric matrix
                ans[i, j] = gam
                ans[j, i] = gam
        return ans
            
    def getsmax(self):
        '''Generates maximum apparent power matrix for the case
        
        Element (i,j) is maximum apperent power of line between buses i -> j
        Zero when no line exists between two buses'''
        ans = np.zeros((self.N, self.N))
        for line in self.branch:
            i, j = line[0] - 1, line[1] - 1
            val = line[5]/self.mva
            if i< self.N and j< self.N:
                ans[i, j] = val
                ans[j, i] = val
        return ans
    
    def getCost(self):
        '''Generate generation cost matrix
        Ordered as c_2, c_1, c_0'''
        ans = np.zeros((self.N, 3))
        for i in range(len(self.gen)):
            bus = self.gen[i][0] - 1
            if bus < self.N:
                ans[bus, 0] = self.gencost[i][4]
                ans[bus, 1] = self.gencost[i][5]
                ans[bus, 2] = self.gencost[i][6]
        return ans
    
    def getGenData(self):
        '''Generate bounds on generation
        All values will be zero if no generator is present'''
        ans = np.zeros((self.N, 4))
        if not isinstance(self.gen[0], list):
            self.gen = [self.gen]
        for line in self.gen:
            bus = line[0] - 1
            ans[bus, 2] = line[3]/self.mva
            ans[bus, 3] = line[4]/self.mva
            ans[bus, 0] = line[8]/self.mva
            ans[bus, 1] = line[9]/self.mva
        return ans
    
    def getLoadData(self):
        '''Load information
        Ordered as active power demand, reactive power demand'''
        ans = np.zeros((self.N, 2))
        for line in self.bus:
            bus = line[0] - 1
            ans[bus, 0] = line[2]/self.mva
            ans[bus, 1] = line[3]/self.mva
        return ans
    
    def getVlim(self):
        '''Generate voltage limits
        '''
        ans = np.zeros((self.N, 2))
        for line in self.bus:
            bus = line[0]-1
            ans[bus, 0] = line[11]
            ans[bus, 1] = line[12]
        return ans

    def getLines(self):
        ans = []
        for i in range(len(self.branch)):
            entry = (self.branch[i][0]-1, self.branch[i][1]-1)
            ans.append(entry)
        return ans
    
    def displayData(self):
        Load_labels=["P_d [MVa]","Q_d [MVar]"]
        Load_df=pd.DataFrame(self.loadData, columns=Load_labels)
        display(Load_df)

        Gen_labels = ['P_max [Mva]', 'P_min [Mva]', 'Q_max [MVar]', 'Q_min [MVar]']
        Gen_df=pd.DataFrame(self.genData, columns=Gen_labels)
        display(Gen_df)

        Lines_labels = ["from bus", "to bus", "R (p.u.)", "X (p.u.)","S_max (MVA)"]
        costs_df = pd.DataFrame(np.asarray(self.branch)[:,(0,1,2,3,5)], columns=Lines_labels)
        display(costs_df)

        Costs_labels=["c2 [$/MW^2]", "c1 [$/MW]", "c0 [$]"]
        costs_df = pd.DataFrame(self.cost, columns=Costs_labels)
        display(costs_df)
    

def loadCase(filename):
    with open(filename) as file:
        dct = json.loads(file.read())
        file.close()
    return UsefulCase(dct)