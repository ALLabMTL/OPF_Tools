import numpy as np
import cvxpy as cp
import networkx as nx
from OPF_Tools.result import RunResult
from OPF_Tools.chordalification import Relaxation, getChordalExtension


def runOPF(case, relaxation_type='SDR', verb=False, solver=None):
    '''Find the solution to the OPF specified.
    The solver uses a sparse representation of the problem

    Relaxation defines the type of relaxation used. These are 
    'SDR' (default)
    'Chordal_MFI'
    'Chordal_AMD'
    'Chordal_MD'
    'Chordal_MCS_M'
    'SOCR'
    'TCR'
    'STCR'

    verb invokes the verbose option in cvxpy

    Returns:
        RunResult instance with the optimal values
        -> None if optimizer does not converge 
    '''

    try:
        relaxation = Relaxation[relaxation_type]
    except KeyError:
        print(f"'{relaxation_type}' is not a valid Relaxation")
        return None
    linking_constraints = 0

    # load data from case
    n = case.N
    baseMVA = case.mva
    Load_data = case.loadData
    Gen_data = case.genData
    Costs_data = case.cost
    Y = case.adj
    Sij = case.smax
    v_lim = case.vlim
    lines = case.getLines()

    # Network graph
    Network = nx.Graph()
    Network.add_nodes_from(range(n))
    Network.add_edges_from(lines)

    #Voltage matrix
    if relaxation in Relaxation.Non_Chordal_Relaxations:
        W = cp.Variable((n,n), hermitian=True)
    
    if relaxation in Relaxation.Chordal_Relaxations:
        chordal_extension, ordering, cliques = getChordalExtension(Network,relaxation)
        N_cliques = len(cliques)
        fill_in = chordal_extension.number_of_edges() - Network.number_of_edges()
        mean_size_of_cliques = np.mean([len(clique) for clique in cliques])
        
        W_K = {clique: cp.Variable((len(clique), len(clique)), hermitian=True) for clique in cliques}

    # power transfer variables
    pij = cp.Variable((len(lines)))
    pji = cp.Variable((len(lines)))
    qij = cp.Variable((len(lines)))
    qji = cp.Variable((len(lines)))

    # generation variable
    pi_g = cp.Variable((n))
    qi_g = cp.Variable((n))

    # Define constraints
    constraints=[] 

    # Calculate the sum of all inbound power flows to each bus
    for i in range(n) :  
        psum = 0
        qsum = 0
        for line in range(len(lines)):
            start, end = lines[line]
            if start == i:
                psum += pij[line]
                qsum += qij[line]
            elif end == i:
                psum += pji[line]
                qsum += qji[line]

        # Sum pij = pi
        constraints += [psum == pi_g[i]-Load_data[i,0]]
        # Sum qij = qi 
        constraints += [qsum == qi_g[i]-Load_data[i,1]]

        #Constraints on active and reactive generation (min-max)
        constraints+=[pi_g[i] >= Gen_data[i,1], pi_g[i] <= Gen_data[i,0]]
        constraints+=[qi_g[i] >= Gen_data[i,3], qi_g[i] <= Gen_data[i,2]]

        # Voltage limits
        if relaxation in Relaxation.Non_Chordal_Relaxations:
            constraints+=[cp.real(W[i,i])>= (v_lim[i,1])**2, cp.real(W[i,i]) <= (v_lim[i,0])**2]
        
        if relaxation in Relaxation.Chordal_Relaxations:
            for clique in cliques:
                for i in clique:
                    pos = clique.index(i)
                    constraints += [v_lim[i,1]**2 <= cp.real(W_K[clique][pos,pos]), cp.real(W_K[clique][pos,pos]) <= v_lim[i,0]**2]
                    break

    # Power flow equations (sparse representation)
    for line in range(len(lines)):
        i, j = lines[line]
            
        #Powerflow
        if relaxation in Relaxation.Non_Chordal_Relaxations:
            constraints+=[pij[line] + 1j*qij[line]==(W[i,i]-W[i,j])*np.conjugate(Y[i,j])] 
            constraints+=[pji[line] + 1j*qji[line]==(W[j,j]-W[j,i])*np.conjugate(Y[j,i])] 
        
        if relaxation in Relaxation.Chordal_Relaxations:
            for clique in cliques:
                if i in clique and j in clique:
                    pos_i = clique.index(i)
                    pos_j = clique.index(j)
                    constraints+=[pij[line] + 1j*qij[line]==(W_K[clique][pos_i,pos_i]-W_K[clique][pos_i,pos_j])*np.conjugate(Y[i,j])]
                    constraints+=[pji[line] + 1j*qji[line]==(W_K[clique][pos_j,pos_j]-W_K[clique][pos_j,pos_i])*np.conjugate(Y[j,i])]
                    break
            
        if not Sij[i,j] == 0:
        #Apparent power capacity S_bar
            constraints+=[cp.square(pij[line])+cp.square(qij[line])<=cp.square(Sij[i,j])]
            constraints+=[cp.square(pji[line])+cp.square(qji[line])<=cp.square(Sij[j,i])]
               
    # SDR problem
    if relaxation == Relaxation.SDR:
        constraints += [W >> 0]

    # SOCR problem
    if relaxation == Relaxation.SOCR:
        for line in range(len(lines)):
            i, j = lines[line]
            constraints+=[cp.norm(cp.hstack([2*W[i,j],(W[i,i]-W[j,j])])) <= cp.real(W[i,i]+W[j,j])]
            constraints+=[cp.norm(cp.hstack([2*W[j,i],(W[j,j]-W[i,i])])) <= cp.real(W[j,j]+W[i,i])]

    if relaxation in Relaxation.Chordal_Relaxations:

        # Computing the clique tree
        from itertools import combinations
        clique_graph = nx.Graph()
        clique_graph.add_nodes_from(cliques, type="clique")
        for edge in combinations(cliques, 2):
            set_edge_0 = set(edge[0])
            set_edge_1 = set(edge[1])
            if not set_edge_0.isdisjoint(set_edge_1):
                sepset = tuple(sorted(set_edge_0.intersection(set_edge_1)))
                clique_graph.add_edge(edge[0], edge[1], weight=len(sepset), sepset=sepset)
        clique_tree = nx.maximum_spanning_tree(clique_graph)

        # Adding the linking constraints from the clique tree
        for edge in clique_tree.edges(data=True):
            clique_1 = edge[0]
            clique_2 = edge[1]
            sepset = edge[2]["sepset"]
            
            pos_sepset_clique_1 = [i for i in range(len(clique_1)) if clique_1[i] in sepset]
            pos_sepset_clique_2 = [i for i in range(len(clique_2)) if clique_2[i] in sepset]
            
            constraints += [W_K[clique_1][np.ix_(pos_sepset_clique_1,pos_sepset_clique_1)] == W_K[clique_2][np.ix_(pos_sepset_clique_2,pos_sepset_clique_2)]]
            linking_constraints += len(sepset)**2
        for clique in cliques:
            constraints += [W_K[clique] >> 0]

    # TCR problem
    if relaxation == Relaxation.TCR:
        v = cp.Variable((n), complex=True)
        constraints += [cp.real(W[0,0]) <= (v_lim[0,0]+v_lim[0,1])*cp.real(v[0])-v_lim[0,0]*v_lim[0,1]]
        constraints += [cp.imag(v[0]) == 0]
        for i in range(n):
            for j in range(n):
                line1 = cp.hstack([1, cp.conj(v[i]), cp.conj(v[j])])
                line2 = cp.hstack([v[i], W[i,i], W[i,j]])
                line3 = cp.hstack([v[j], cp.conj(W[i,j]), W[j,j]])
                mat = cp.vstack([line1, line2, line3])
                constraints += [mat >> 0]
    
    # STCR problem
    if relaxation == Relaxation.STCR:
        for i in range(n):
            for j in range(n):
                line1 = cp.hstack([W[0,0], W[0,i], W[0,j]])
                line2 = cp.hstack([cp.conj(W[0,i]), W[i,i], W[i,j]])
                line3 = cp.hstack([cp.conj(W[0,j]), cp.conj(W[i,j]), W[j,j]])
                mat = cp.vstack([line1, line2, line3])
                constraints += [mat >> 0]

    # Define costs
    Costs = 0
    for i in range(n): 
        c0 = Costs_data[i][2]
        c1 = Costs_data[i][1]
        c2 = Costs_data[i][0]
        if c1 > 0.0: # Bus has a generator installed
            Costs += c0+c1*pi_g[i]*baseMVA+c2*cp.square(pi_g[i]*baseMVA)
    
    prob = cp.Problem(cp.Minimize(Costs),constraints)
    try:
        if solver is not None:
            prob.solve(verbose=verb, solver=solver)
        else:
            prob.solve(verbose=verb)
        loss = prob.value
    except cp.SolverError:
        loss = np.NaN
    status = prob.status
    ans = RunResult()
    if status == 'optimal':
        ans.set_p_q(pi_g*baseMVA, qi_g*baseMVA)
        if relaxation in Relaxation.Non_Chordal_Relaxations:
            ans.set_W(W)
    if relaxation in Relaxation.Non_Chordal_Relaxations:
        ans.setAll(status, loss, Network, prob._solve_time, prob._compilation_time, relaxation)
    if relaxation in Relaxation.Chordal_Relaxations:
        ans.setAll_chordal(status, loss, Network, chordal_extension, ordering, prob._solve_time, prob._compilation_time, relaxation, N_cliques, fill_in, linking_constraints, mean_size_of_cliques)
    return ans