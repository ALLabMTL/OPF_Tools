# OPF_Tools

This library enables the loading and solving of matpower cases in Python. 

# Imports

```python
from OPF_Tools import *
```


# Usage
There are 3 basic steps to utilization. First we load a case either from a file or from a dictionary. This case can then be solved using the runOPF function. The result is stored in a RunResult object that permits access to the important optimal values. 

## Option 1
```python
case = loadCase('cases/case9.json')

result = runOPF(case)

print('Optimal result is : {}'.format(result.loss))
```

## Option 2
```python
case = Case('cases/case9.json')

# Run using TCR relaxation and verbose option
result = runOPF(case,'TCR',True)

# relaxations = ['SDR','Chordal_MFI','Chordal_AMD','Chordal_MD','Chordal_MCS_M','SOCR','TCR','STCR']

print('Optimal generation is: ')
print(result.p)
```

## Chordal Relaxation



# Case
A case lists all required information for a standard MatPower test case. It can be loaded directly from a json file or entered manually. All json files in the "cases" folder are compatible with the Case format. The fields in a Case object share the names of the standard MatPower mpc object.

## UsefulCase
UsefulCase inherits from the Case object. This object defines a few values that are useful for solving the OPF with a cost minimization as objective. This includes all the fields presented in the standard MatPower case and:
- An adjacency matrix representing admittance between two nodes (zero if no line connects the two)
- A maximum apparent power matrix 
- A cost matrix that contains the c0, c1 and c2 coefficients for generators
- Upper and lower bounds on active and reactive power generation
- Active and reactive power demands at all nodes
- Limits on voltage at each node



# Result
The RunResult object stores the run result of an opf instance.


## Format
Four elements are saved in a result. The (complex) square voltage matrix "W", the active power generation "p", the reactive power generation "q" and the optimal value loss.

It is possible to set these parametrs using either a numpy array or a cvxpy variable.

## Chordal Relaxation

If a chordal relaxation has been used, the result stores : 
- the active power generation in ```result.p```
- the reactive power generation in ```result.q```
- the objective value in ```result.loss```
- A networkx graph representation of the original network in ```result.network```
- The chordal extension in ```result.chordal_extension```
- The ordering of the nodes used to compute the chordal extension in ```result.ordering```
- The solve time in ```result.solve_time```
- The compilaiton time in ```result.compilation_time```
- Information about the clique decomposition :
```python
from OPF_Tools import *
import networkx as nx

case = loadCase('OPF_Tools/cases/case30.json')
result = runOPF(case,'Chordal_AMD', verb = False, solver = 'MOSEK')

print('Results for case30.json with Chordal_AMD relaxation:')
print('\n')

print(f'Optimal result is : {result.loss:.3f}')
print(f'Solve time is : {result.solve_time:.2f} seconds')
print(f'Compilation time is : {result.compilation_time:.2f} seconds')
print('\n')

print(f'Number of cliques : {result.number_of_cliques}')
print(f'Fill in : {result.fill_in}')
print(f'Linking constraints : {result.linking_constraints}')
print(f'Mean size of cliques : {result.mean_size_of_cliques}')

nx.draw(result.network, with_labels=True)
nx.draw(result.chordal_extension, with_labels=True)
```

```console
Results for case30.json with Chordal_AMD relaxation:


Optimal result is : 580.666
Solve time is : 0.09 seconds
Compilation time is : 1.83 seconds


Number of cliques : 26
Fill in : 14
Linking constraints : 93
Mean size of cliques : 2.96154
```
