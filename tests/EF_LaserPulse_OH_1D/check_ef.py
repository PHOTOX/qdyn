from os.path import exists
from numpy import genfromtxt, max, abs

# Exact factorization files to be compared
# todo: add gd-tdpes to tests
files = ['el_coefficients_ef', 'gi-tdpes', 'nuclear_density', 'nuclear_phase']
# files = ['el_coefficients_ef', 'gd-tdpes', 'gi-tdpes', 'nuclear_density', 'nuclear_phase']


result = True  # if mismatch between result and reference appear, it's overwritten with False
output = open('results.dat', mode='w')

for file in files:
    # read calculated data
    data_file = file + '.dat'
    if not exists(data_file):
        print(f'ERROR: {data_file:s} not found!')
        exit(1)
    data = genfromtxt(data_file)

    # read reference
    ref_file = file + '_ref.dat'
    if not exists(ref_file):
        print(f'ERROR: {ref_file:s} not found!')
        exit(1)
    ref = genfromtxt(ref_file)

    # calculate maximum difference in the files (including x axis which shouldn't change)
    max_difference = max(abs(data - ref))

    # compare max difference to threshold
    if max_difference > 1e-10:
        result = False
        output.write(f"Mismatch found in '{data_file}'! (max difference: {max_difference:.2e})\n")

# close output file
output.close()

# print test result
print(result)