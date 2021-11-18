import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

# borrowed fancy arrow class
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

coords = [[-2,-2,-2], [2,-2,-2], [2,2,-2], [-2,2,-2],
          [-2,-2,2], [2,-2,2], [2,2,2] , [-2,2,2], 
          [-2,-2,6], [2,-2,6], [2,2,6] , [-2,2,6]]
#U = np.loadtxt("U.txt")
# result from FEM by using 10*F for the global force vect to see more clearly the results
U = np.loadtxt("U.txt") 
print(U)

X, Y = np.meshgrid(np.arange(-2, 3), np.arange(-2, 3))
Z = 0*X-2

fig = plt.figure(figsize=(14,15))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, alpha=0.3)  # the horizontal plane
ax.plot_surface(X, Y, Z+4, alpha=0.3)  # the horizontal plane
ax.plot_surface(X, Y, Z+8, alpha=0.3)  # the horizontal plane

ax.plot_surface(Z, Y, X, alpha=0.3)
ax.plot_surface(Z+4, Y, X, alpha=0.3)

ax.plot_surface(Z, Y, X+4, alpha=0.3)
ax.plot_surface(Z+4, Y, X+4, alpha=0.3)


for node in np.arange(12):
    print("node = ",node)
    
    for i in np.arange(3):
        print(coords[node][i], ",", U[3*node+i])
        
    ax.scatter(coords[node][0],coords[node][1],coords[node][2],color='blue' ) # plot the point (2,3,4) on the figure
    ax.scatter(coords[node][0]+ U[3*node+0],coords[node][1]+U[3*node+1],coords[node][2]+U[3*node+2],color='red' ) # plot the point (2,3,4) on the figure
    arw = Arrow3D([coords[node][0],coords[node][0]+ U[3*node+0]],[coords[node][1],coords[node][1]+U[3*node+1]],[coords[node][2],coords[node][2]+U[3*node+2]], arrowstyle="->", color="purple", lw = 1, mutation_scale=10)
    ax.add_artist(arw)
    
    

ax.view_init(45, 30)
plt.show()
