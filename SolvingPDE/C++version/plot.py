from matplotlib import pyplot as plt

x = []
y_approx = []
y_exact = []
with open('data_for_plots.txt', 'r') as inp:
    x = [float(s) for s in inp.readline().strip().split()]
    y_approx = [float(s) for s in inp.readline().strip().split()]
    y_exact = [float(s) for s in inp.readline().strip().split()]

plt.plot(x, y_approx, linestyle='dashed', color='black', label='Approximate solution')
plt.plot(x, y_exact, color='red', label='Exact solution')
plt.savefig(fname='fig.png')
plt.legend()
