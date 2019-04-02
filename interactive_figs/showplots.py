import pickle
import matplotlib.pyplot as plt
from glob import glob

files = glob('matplotlib_*')
files = sorted(files)
for filename in files:
    p = pickle.load(open(filename, 'rb'))
    plt.show()

