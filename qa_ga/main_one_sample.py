from utils.RT import ReverseTranslation
from utils.extract_peptides import get_peptides

# peptides = get_peptides()
peptide = "MFVFLVLLPLVSSQCVNLTT"
project = ReverseTranslation(peptide, expression="e_coli_316407", constant_f=0.1, constant_gc=1, constant_r=0.1,
                             constant_epsilon=5, rou_t=0.5, big_num=50, verbose=True)
# arr = [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
#        0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
#        0, 1, 0, 0, 0, 1, 0]
# val = project.evaluate(arr)
# print(val)

res = project.compare()
# print(res)
