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
result = runOPF(case, 'TCR', True)

print('Optimal generation is: ')
print(result.p)
```

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
