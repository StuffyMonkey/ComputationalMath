from matplotlib import pyplot as plt
h = 0.05
N = round(1/h)

def get_matrix(filename):
    with open(filename, 'r') as inp:
        lst = inp.readline().split()
        data = []
        while(lst):
            data.append([float(s) for s in lst])
            lst = inp.readline().split()
        return data

domain = [i*h for i in range(0, N+1)]

fig = plt.figure()
ax = plt.axes(projection='3d')

Z1 = get_matrix('approximate_solution.txt')
Z2 = get_matrix('exact_solution.txt')

# black - original solution, plasma - numerical
ax.contour3D(domain, domain, Z1, 150, cmap='bone')
# colorful - numerical solution
ax.contour3D(domain, domain, Z2, 250, cmap='plasma')
# ax.plot_surface(T, X, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.view_init(30, 110)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.savefig(fname='3dplot.png')