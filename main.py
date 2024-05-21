from simulation import runOPF
from case import loadCase

instance_name = 'case30'
relaxation = 'Chordal_AMD'
solver = 'MOSEK'
Verbose = False

case = loadCase('cases/'+instance_name+'.json')

ans = runOPF(case, relaxation, Verbose, solver)

print('Case: {}'.format(instance_name))
print('Optimal result is : {:.3f}'.format(ans.loss))
print('Solve time is : {:.2f} seconds'.format(ans.solve_time))
print('Compilation time is : {:.2f} seconds'.format(ans.compilation_time))
print()

# relaxations = ['SDR','Chordal_MFI','Chordal_AMD','Chordal_MD','Chordal_MCS_M','SOCR','TCR','STCR']
# cases = ['case3','case4_dist','case4gs','case5','case6ww','case9','case9Q','case9target','case10ba','case10tp2','case12da','case14','case15da','case15nbr','case16am','case16ci','case17me','case18','case18nbr','case22','case24_ieee_rts','case28da','case30','case30pwl','case30Q','case33bw','case33mg','case34sa','case38si','case39','case51ga','case51he','case57','case60nordic','case69','case70da','case74ds','case85','case89pegase','case94pi','case118','case118zh','case136ma','case141','case145','case300','case1354pegase','case1888rte','case1951rte','case2383wp','case2736sp','case2737sop','case2746wop','case2746wp','case2848rte','case2868rte','case2869pegase','case3012wp','case3120sp','case3375wp','case6468rte','case6470rte','case6495rte','case6515rte','case8387pegase','case9241pegase','case13659pegase','contab_ACTIVSg10k','contab_ACTIVSg200','contab_ACTIVSg500','contab_ACTIVSg2000','scenarios_ACTIVSg200','scenarios_ACTIVSg2000','testcase','case_ACTIVSg10k','case_ACTIVSg25k','case_ACTIVSg70k','case_ACTIVSg200','case_ACTIVSg500','case_ACTIVSg2000','case_ieee30','case_RTS_GMLC','case_SyntheticUSA']

