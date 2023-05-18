class RunResult:
    '''Object contains relevant information on optimal values obtained from optimization'''
    def __init__(self):
        '''Instantiate the object'''
        self.W = None
        self.p = None
        self.q = None
        self.loss = 0

    def setAll(self, p, q, W, loss):
        '''Set all parameters including loss'''
        self.W = W
        self.p = p
        self.q = q
        self.loss = loss

    def setFromVars(self, p, q, W):
        '''Set parameters without loss'''
        self.W = W.value
        self.p = p.value
        self.q = q.value

    def setLoss(self, loss):
        '''Set optimal losses'''
        self.loss = loss