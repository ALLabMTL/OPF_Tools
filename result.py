import networkx as nx
from OPF_Tools.chordalification import Relaxation

class RunResult:
    '''Object contains relevant information on optimal values obtained from optimization'''
    def __init__(self):
        '''Instantiate the object'''
        self.status = None
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
        self.nb_of_sdp_coeffs = 0
        self.nb_of_independent_sdp_coeffs = 0

    def setAll(self, status, loss, network, solve_time, compilation_time, relaxation):
        '''Set all parameters including loss'''
        self.status = status
        self.loss = loss
        self.network = network
        self.solve_time = solve_time
        self.compilation_time = compilation_time
        self.relaxation = relaxation
    
    def setStatus(self, status):
        '''Set status of optimization'''
        self.status = status

    def set_p_q(self, p, q):
        '''Set parameters without loss'''
        self.p = p
        self.q = q
    
    def set_W(self, W):
        '''Set parameters without loss'''
        self.W = W

    def setLoss(self, loss):
        '''Set optimal losses'''
        self.loss = loss

    def setAll_chordal(self, status, loss, network, chordal_extension, ordering, solve_time, compilation_time, relaxation, number_of_cliques, fill_in, linking_constraints, mean_size_of_cliques, nb_of_sdp_coeffs, nb_of_independent_sdp_coeffs):
        '''Set all parameters including loss'''
        self.status = status
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
        self.nb_of_sdp_coeffs = nb_of_sdp_coeffs
        self.nb_of_independent_sdp_coeffs = nb_of_independent_sdp_coeffs

    def is_radial(self):
        return nx.is_tree(self.network)