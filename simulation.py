import numpy as np
import cvxpy as cp
import networkx as nx
from OPF_Tools.result import RunResult
from OPF_Tools.chordalification import Relaxation, getChordalExtension


def runOPF(case, relaxation_type='SDR', verb=False, solver=None, perturb_loads=(0,0)):
    '''Find the solution to the OPF specified.
    The solver uses a sparse representation of the problem

    Relaxation defines the type of relaxation used. These are:
    'SDR' (default)    -> Semidefinite relaxation
    'Chordal_MFI'      -> Chordal relaxation using the Minimum Fill-In ordering
    'Chordal_AMD'      -> Chordal relaxation using the Approximate Minimum Degree ordering
    'Chordal_MD'       -> Chordal relaxation using the Minimum Degree ordering
    'Chordal_MCS_M'    -> Chordal relaxation using the Maximum Cardinality Search ordering
    'SOCR'             -> Second Order Cone relaxation
    'TCR'              -> Tight and Cheap relaxation
    'STCR'             -> Strong Tight and Cheap relaxation

    verb invokes the verbose option in cvxpy
    
    perturb_loads is a tuple (std_dev, seed) that perturbates the loads with a normal distribution centered at 0 with a standard deviation proportional to the original value

    Returns:
        RunResult instance with the optimal values
        -> RunResult.status is None if the solver does not converge
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
    
    if perturb_loads[0]:
        n = len(Load_data)
        np.random.seed(perturb_loads[1])
        sum_active = np.sum(Load_data[:,0])
        sum_reactive = np.sum(Load_data[:,1])
        Load_data += np.random.normal(0, np.abs(perturb_loads[0]*Load_data), (n,2))
        Load_data = np.clip(Load_data, 0, None)
        sum_active_new = np.sum(Load_data[:,0])
        sum_reactive_new = np.sum(Load_data[:,1])
        diff_active = sum_active_new - sum_active
        diff_reactive = sum_reactive_new - sum_reactive
        print(f'A total of {diff_active:.2f} MW and {diff_reactive:.2f} MVar were added to the system')

    # Network graph
    Network = nx.Graph()
    Network.add_nodes_from(range(n))
    Network.add_edges_from(lines)
    
    print('Setting up the variables')

    #Voltage matrix
    if relaxation in Relaxation.Non_Chordal_Relaxations:
        W = cp.Variable((n,n), hermitian=True)
    
    if relaxation in Relaxation.Chordal_Relaxations:
        chordal_extension, ordering, cliques = getChordalExtension(Network,relaxation)
        N_cliques = len(cliques)
        fill_in = chordal_extension.number_of_edges() - Network.number_of_edges()
        mean_size_of_cliques = np.mean([len(clique) for clique in cliques])
        
        W_K = {clique: cp.Variable((len(clique), len(clique)), hermitian=True) for clique in cliques}
        
        # Creating dictionaries that map the nodes and edges to their respective cliques
        
        clique_sets = [set(clique) for clique in cliques]
        
        node_to_cliques = [[] for _ in range(n)]
        for idx, clique in enumerate(clique_sets):
            for node in clique:
                node_to_cliques[node].append(idx)

        edge_to_clique = {}
        for edge in lines:
            edge_set = set(edge)
            found = False
            for node in edge:
                for clique_idx in node_to_cliques[node]:
                    if edge_set.issubset(clique_sets[clique_idx]):
                        edge_to_clique[edge] = clique_idx
                        found = True
                        break
                if found:
                    break

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
    print('Setting up voltage limits and power balance constraints')
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
        
        # Voltage limits
        if relaxation in Relaxation.Non_Chordal_Relaxations:
            constraints+=[cp.real(W[i,i])>= (v_lim[i,1])**2, cp.real(W[i,i]) <= (v_lim[i,0])**2]
        if relaxation in Relaxation.Chordal_Relaxations:
            clique = cliques[node_to_cliques[i][0]]
            pos = clique.index(i)
            constraints += [v_lim[i,1]**2 <= cp.real(W_K[clique][pos,pos]), cp.real(W_K[clique][pos,pos]) <= v_lim[i,0]**2]

        # Constraints on active and reactive generation (min-max)
        if Gen_data[i, 0] is not None:
            constraints += [pi_g[i] <= Gen_data[i, 0]]
        if Gen_data[i, 1] is not None:
            constraints += [pi_g[i] >= Gen_data[i, 1]]
        if Gen_data[i, 2] is not None:
            constraints += [qi_g[i] <= Gen_data[i, 2]]
        if Gen_data[i, 3] is not None:
            constraints += [qi_g[i] >= Gen_data[i, 3]]

    # Power flow equations (sparse representation)
    print('Setting up the power flow equations')
    for line in range(len(lines)):
        i, j = lines[line]
        
        #Powerflow
        if relaxation in Relaxation.Non_Chordal_Relaxations:
            constraints+=[pij[line] + 1j*qij[line]==(W[i,i]-W[i,j])*np.conjugate(Y[i,j])] 
            constraints+=[pji[line] + 1j*qji[line]==(W[j,j]-W[j,i])*np.conjugate(Y[j,i])] 
        
        if relaxation in Relaxation.Chordal_Relaxations:
            clique = cliques[edge_to_clique[(i,j)]]
            pos_i = clique.index(i)
            pos_j = clique.index(j)
            constraints+=[pij[line] + 1j*qij[line]==(W_K[clique][pos_i,pos_i]-W_K[clique][pos_i,pos_j])*np.conjugate(Y[i,j])]
            constraints+=[pji[line] + 1j*qji[line]==(W_K[clique][pos_j,pos_j]-W_K[clique][pos_j,pos_i])*np.conjugate(Y[j,i])]  
            
        if not Sij[i,j] == 0:
        #Apparent power capacity S_bar
            constraints+=[cp.square(pij[line])+cp.square(qij[line])<=cp.square(Sij[i,j])]
            constraints+=[cp.square(pji[line])+cp.square(qji[line])<=cp.square(Sij[j,i])]
    
    print('Setting up the relaxation')
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
    print('Solving the problem')
    try:
        if solver is not None:
            if solver == 'MOSEK':
                # import os
                # prob.solve(warm_start = True, solver='MOSEK',mosek_params = {"MSK_IPAR_NUM_THREADS":os.cpu_count()},verbose=verb,ignore_dpp = True)
                prob.solve(verbose=verb, solver=solver)
            else:
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