
# note: each time you change mesh need to update the nodelist, solution list, and the N over here in this code 
# 21 feb: Import all text files K, F and nodeList and update N

import sys

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch


        
        
# =================================
N = int(sys.argv[1])
nb_nodes = (N+1)**3
fileSaveName = "result" + str(N) +".png"
# =================================


# ==============================
#         MATRIX K 
# ==============================
with open("K.txt") as K_in:
    K_lines = []
    for Ki in K_in: 
        Ki_float = np.array(Ki.rstrip(' \n').split(' '))
        K_lines.append(Ki_float)
    KK = np.asarray(K_lines, dtype = np.float64, order = 'C')



# print("KK ======== \n")
# print(KK)
# print("Check KK size \n")
# print(KK.shape)


print("nb_nodes: \n")
print(nb_nodes)

# ==============================
#         VECTOR F
# ==============================


with open("F.txt") as F_in:
    F_lines = []

    for line in F_in:
        F_lines.append(float(line.rstrip('\n')))
        
    FF = np.array(F_lines)

            
# print("FF ======== \n")     
# print(FF)
# print("Check FF size \n")
# print(FF.shape)



# ==============================
#           NODE LIST
# ==============================


with open("nodeList.txt") as nodes_in:
    nodes_lines = []
    i = 0;
    coords = []
    for line in nodes_in:
        
        #print(line)
        if (line != "--node-position--\n") and (line != "\n"):
            if i == 3 :
                #print("Coords:",coords)
                nodes_lines.append(coords)
                coords = []
                i = 0
            
            #print("i = ",i)

            coords.append(int(line.rstrip('\n')))
            

            i += 1
            #nodes_lines.append(int(line.rstrip('\n')))

            
# print(nodes_lines)

# =================================
#   END OF DAT EXTRACTION 
# =================================

# =================================
#   SOLVE BY PARTITION
# =================================

UU =  [1] * 3*nb_nodes
count = 0
for ii in np.arange(3*nb_nodes):
#     print("ii: ", ii ,"\n")
    ii_mod = ii-count;
#     print("ii_mod: ",ii_mod ,"\n")
    if (KK[ii_mod][ii_mod] == 1):
        
        
        
#         remove row and column
#       size of KK and FF chlange at eachl ite so need use ii_mod
        KK = np.delete(KK, (ii_mod), axis=0)
        KK = np.delete(KK, (ii_mod), axis=1)

        FF = np.delete(FF, (ii_mod), axis=0)        
        
        UU[ii] = 0
#         print("Check size of K and F:\n", KK.shape, "\n", FF.shape)
        count += 1
#         print("Boundary found, count = ", count, "\n")
#         print("======================================")
        
        
UU_res = np.linalg.solve(KK, FF)
        
        
count = 0;
for ii in np.arange(3*nb_nodes):
    if UU[ii] == 0:
        count += 1
    else:
        UU[ii] = UU_res[ii-count];

#UU_res = KK\FF;

# print("Solution, U:\n", UU)
# print(len(UU))


# =================================
#   PLOTS
# =================================

X, Y = np.meshgrid(np.arange(-2, 12), np.arange(-2, 12))
Z = 0*X-2

fig = plt.figure(figsize=(14,15))
ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z, alpha=0.3)  # the horizontal plane
# ax.plot_surface(X, Y, Z+4, alpha=0.3)  # the horizontal plane
# ax.plot_surface(X, Y, Z+8, alpha=0.3)  # the horizontal plane

# ax.plot_surface(Z, Y, X, alpha=0.3)
# ax.plot_surface(Z+4, Y, X, alpha=0.3)

# ax.plot_surface(Z, Y, X+4, alpha=0.3)
# ax.plot_surface(Z+4, Y, X+4, alpha=0.3)
X = []; Y = []; Z = [];

for node in np.arange(nb_nodes-1):
    #print("node = ",node)
    #X.append(nodes_lines[node][0]+ soln_lines[3*node+0]); Y.append(nodes_lines[node][1]+soln_lines[3*node+1]); Z.append(nodes_lines[node][2]+soln_lines[3*node+2]); 

    #for i in np.arange(3):
        #print(nodes_lines[node][i], ",", soln_lines[3*node+i])
        
    ax.scatter(nodes_lines[node][0],nodes_lines[node][1],nodes_lines[node][2],color='black' ,linewidths=1) # plot the point (2,3,4) on the figure
    ax.scatter(nodes_lines[node][0]+ UU[3*node+0],nodes_lines[node][1]+UU[3*node+1],nodes_lines[node][2]+UU[3*node+2],color='red' ,linewidths=1) # plot the point (2,3,4) on the figure
    #arw = Arrow3D([nodes_lines[node][0],nodes_lines[node][0]+ soln_lines[3*node+0]],[nodes_lines[node][1],nodes_lines[node][1]+soln_lines[3*node+1]],[nodes_lines[node][2],nodes_lines[node][2]+soln_lines[3*node+2]], arrowstyle="->", color="purple", lw = 1, mutation_scale=10)
    
    # GLOBAL
    #ax.scatter(soln_lines[3*node+0],soln_lines[3*node+1],soln_lines[3*node+2],color='red' ) # plot the point (2,3,4) on the figure
    #arw = Arrow3D([nodes_lines[node][0], soln_lines[3*node+0]],[nodes_lines[node][1],soln_lines[3*node+1]],[nodes_lines[node][2],soln_lines[3*node+2]], arrowstyle="->", color="purple", lw = 1, mutation_scale=10)
    #ax.add_artist(arw)
    ax.plot([nodes_lines[node][0],nodes_lines[node][0]+ UU[3*node+0]],[nodes_lines[node][1],nodes_lines[node][1]+UU[3*node+1]],[nodes_lines[node][2],nodes_lines[node][2]+UU[3*node+2]], 'b-');
    



#ax.contour(X,Y,Z, color='black')
ax.view_init(0,90)

# top view ----------
# ax.view_init(90,45)

ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
ax.set_zlim([-1,N+1])
# plt.savefig('result.png')
plt.savefig(fileSaveName)
plt.show()
