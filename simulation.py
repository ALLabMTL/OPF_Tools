import numpy as np
import cvxpy as cp

from OPF_Tools.result import *

def runOPF(case, type=0, verb=False, solver=None, max_iters=10000):
    '''Find the solution to the OPF specified.
    The solver uses a sparse representation of the problem

    Type defines the type of relaxation used. These are 
    0 (default) -> SDR
    1 -> SOCR
    2 -> TCR
    3 -> STCR

    verb invokes the verbose option in cvxpy

    Returns:
        RunResult instance with the optimal values
        -> None if optimizer does not converge 
    '''
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

    #Voltage matrix
    V = cp.Variable((n,n), hermitian=True) 

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
        constraints+=[cp.real(V[i,i])>= (v_lim[i,1])**2, cp.real(V[i,i]) <= (v_lim[i,0])**2]
        constraints += [cp.imag(V[i,i]) == 0]

    # Power flow equations (sparse representation)
    for line in range(len(lines)):
        i, j = lines[line]
            
        #Powerflow
        constraints+=[pij[line] + 1j*qij[line]==(V[i,i]-V[i,j])*np.conjugate(Y[i,j])] 
        constraints+=[pji[line] + 1j*qji[line]==(V[j,j]-V[j,i])*np.conjugate(Y[j,i])] 
        
        if not Sij[i,j] == 0:
        #Apparent power capacity S_bar
            constraints+=[cp.square(pij[line])+cp.square(qij[line])<=cp.square(Sij[i,j])]
            constraints+=[cp.square(pji[line])+cp.square(qji[line])<=cp.square(Sij[j,i])]
            
           
    # SDR problem
    if type == 0:
        constraints += [V >> 0.]

    # SOCR problem
    if type == 1:
        for line in range(len(lines)):
            i, j = lines[line]
            constraints+=[cp.norm(cp.hstack([2*V[i,j],(V[i,i]-V[j,j])])) <= cp.real(V[i,i]+V[j,j])]
            constraints+=[cp.norm(cp.hstack([2*V[j,i],(V[j,j]-V[i,i])])) <= cp.real(V[j,j]+V[i,i])]

    # TCR problem
    if type == 2:
        v = cp.Variable((n), complex=True)
        constraints += [cp.real(V[0,0]) <= (v_lim[0,0]+v_lim[0,1])*cp.real(v[0])-v_lim[0,0]*v_lim[0,1]]
        constraints += [cp.imag(v[0]) == 0]
        for i in range(n):
            for j in range(n):
                line1 = cp.hstack([1, cp.conj(v[i]), cp.conj(v[j])])
                line2 = cp.hstack([v[i], V[i,i], V[i,j]])
                line3 = cp.hstack([v[j], cp.conj(V[i,j]), V[j,j]])
                mat = cp.vstack([line1, line2, line3])
                constraints += [mat >> 0.]
    
    # STCR problem
    if type == 3:
        for i in range(n):
            for j in range(n):
                line1 = cp.hstack([V[0,0], V[0,i], V[0,j]])
                line2 = cp.hstack([cp.conj(V[0,i]), V[i,i], V[i,j]])
                line3 = cp.hstack([cp.conj(V[0,j]), cp.conj(V[i,j]), V[j,j]])
                mat = cp.vstack([line1, line2, line3])
                constraints += [mat >> 0.]

    # Define costs
    Costs = 0
    for i in range(n): 
        c0 = Costs_data[i][2]
        c1 = Costs_data[i][1]
        c2 = Costs_data[i][0]
        if c1 > 0.0: # Bus has a generator installed
            Costs += c0+c1*pi_g[i]*baseMVA+c2*cp.square(pi_g[i]*baseMVA)
    
    prob = cp.Problem(cp.Minimize(Costs),constraints)
    if solver is not None:
        prob.solve(verbose=verb, solver=solver, max_iters=max_iters)
    else:
        prob.solve(verbose=verb)
    ans = RunResult()
    try:
        ans.setAll(pi_g.value*baseMVA, qi_g*baseMVA, V.value, prob.value)
    except TypeError:
        return None
    return ans