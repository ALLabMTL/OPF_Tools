import networkx as nx

class RunResult:
    '''Object contains relevant information on optimal values obtained from optimization'''
    def __init__(self):
        '''Instantiate the object'''
        self.W = None
        self.p = None
        self.q = None
        self.loss = 0
        self.network = None
        self.chordal_extension = None
        self.ordering = None
        self.solve_time = 0
        self.compilation_time = 0
        self.relaxation = None
        self.number_of_cliques = 0
        self.fill_in = 0
        self.linking_constraints = 0
        self.mean_size_of_cliques = 0

    def setAll(self, p, q, W, loss, network, solve_time, compilation_time, relaxation):
        '''Set all parameters including loss'''
        self.W = W
        self.p = p
        self.q = q
        self.loss = loss
        self.network = network
        self.solve_time = solve_time
        self.compilation_time = compilation_time
        self.relaxation = relaxation

    def setFromVars(self, p, q, W):
        '''Set parameters without loss'''
        self.W = W.value
        self.p = p.value
        self.q = q.value

    def setLoss(self, loss):
        '''Set optimal losses'''
        self.loss = loss

    def setAll_chordal(self, p, q, W, loss, network, chordal_extension, ordering, solve_time, compilation_time, relaxation, number_of_cliques, fill_in, linking_constraints, mean_size_of_cliques):
        '''Set all parameters including loss'''
        self.W = W
        self.p = p
        self.q = q
        self.loss = loss
        self.network = network
        self.chordal_extension = chordal_extension
        self.ordering = ordering
        self.solve_time = solve_time
        self.compilation_time = compilation_time
        self.relaxation = relaxation
        self.number_of_cliques = number_of_cliques
        self.fill_in = fill_in
        self.linking_constraints = linking_constraints
        self.mean_size_of_cliques = mean_size_of_cliques

    def is_radial(self):
        return nx.is_tree(self.network)