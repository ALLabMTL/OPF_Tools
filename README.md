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
dct = loadCase('cases/case9.json')
case = Case(dct)

result = runOPF(case)

print('Optimal result is : {}'.format(result.loss))
```

## Option 2
```python
case = Case('cases/case9.json')

# Run using TCR relaxation and verbose option
result = runOPF(case, 2, True)

print('Optimal generation is: ')
print(result.p)
```
