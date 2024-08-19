import sys
from utils.RT import ReverseTranslation
from utils.extract_peptides import get_peptides
import matplotlib.pyplot as plt
import time
import os

start_time = time.time()
# sel_all = True
sel_all = False
num_peptides = 16
peptides = get_peptides(sel_all=sel_all, num_peptides=num_peptides)
energy_values = []

output_dir = "output"
if sel_all:
    working_dir = os.path.join(output_dir,
                               "AllPeptides_" + time.strftime("%Y-%m-%d-%H-%M-%S"))

else:
    working_dir = os.path.join(output_dir,
                               "NumOfPeptides_" + str(num_peptides) + "_" + time.strftime("%Y-%m-%d-%H-%M-%S"))
os.makedirs(working_dir)
os.chdir(working_dir)

with open('output.out', 'w') as f:
    original_stdout = sys.stdout
    sys.stdout = f

    try:
        for peptide in peptides:
            project = ReverseTranslation(peptide, expression="e_coli_316407", constant_f=0.1, constant_gc=1,
                                         constant_r=0.1,
                                         constant_epsilon=5, rou_t=0.5, big_num=50, verbose=False)
            ga_energy, qa_energy = project.compare()
            energy_values.append((ga_energy, qa_energy))
        end_time = time.time()

        print("Energy Values:")
        for energy in energy_values:
            print(energy)
        print(f"Total Time taken: {end_time - start_time} seconds")
    finally:
        sys.stdout = original_stdout

ga_energies, qa_energies = zip(*energy_values)

# Create scatter plot
plt.scatter(ga_energies, qa_energies, label='Energy values')
plt.plot([min(ga_energies), max(ga_energies)], [min(ga_energies), max(ga_energies)], 'r--', label='y=x')

# Add labels and title
plt.xlabel('GA Energy')
plt.ylabel('QA Energy')
plt.title('GA Energy vs QA Energy')
plt.legend()

# Save plot to file
plt.savefig('res.png')
