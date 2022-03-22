
# note: each time you change mesh need to update the nodelist, solution list, and the N over here in this code 
# 21 feb: Import all text files K, F and nodeList and update N

import sys

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch


        
        
# =================================
# N = 3
N = int(sys.argv[1])
nb_nodes = (N+1)**3
# fileSaveName = "result" + str(N) +".png"
# =================================


# ==============================
#         MATRIX K 
# ==============================
KK_iters_list = []
with open("K.txt") as K_in:
    K_lines = []
    for Ki in K_in: 
        if (Ki != "===\n"):
            Ki_float = np.array(Ki.rstrip(' \n').split(' '))
            K_lines.append(Ki_float)
        else:
            # create the matrix KK with collected lines, put it in the list of matrices KK_iters_list, then clear the temporary box K_lines
            KK = np.asarray(K_lines, dtype = np.float64, order = 'C')
            KK_iters_list.append(KK)
            K_lines = []
            
    # KK = np.asarray(K_lines, dtype = np.float64, order = 'C')


for KK in KK_iters_list:
    print("KK ======== \n")
    print(KK)
# print("Check KK size \n")
# print(KK.shape)


# print("nb_nodes: \n")
# print(nb_nodes)

# ==============================
#         VECTOR F
# ==============================

FF_iters_list = []
with open("F.txt") as F_in:
    F_lines = []

    for line in F_in:
        if (line != "===\n"):
            F_lines.append(float(line.rstrip('\n')))
        else:
            FF = np.array(F_lines)
            FF_iters_list.append(FF)
            F_lines = []
        

for FF in FF_iters_list:        
    print("FF ======== \n")     
    print(FF)
# print("Check FF size \n")
# print(FF.shape)
print(len(FF_iters_list), len(KK_iters_list))
iters = len(FF_iters_list)

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

            
            #print("i = ",i)

            coords.append(int(line.rstrip('\n')))
            i += 1

            if i == 3 :
                #print("Coords:",coords)
                nodes_lines.append(coords)
                coords = []
                i = 0
            

           
            #nodes_lines.append(int(line.rstrip('\n')))

            
# print(nodes_lines)
elem_iters_list = []

# with open("elemList.txt") as elem_in:
#     elem_lines = []
#     for elem in elem_in: 
#         print(elem)
#         if (elem != "===\n"):
#             elem_lines.append((int(elem.rstrip('\n'))))
#             # elem_lines.append(elem)
#         else:
#             # create the matrix KK with collected lines, put it in the list of matrices KK_iters_list, then clear the temporary box K_lines
#             elem_iters_list.append(elem_lines)
#             elem_lines = []
            
#     # KK = np.asarray(K_lines, dtype = np.float64, order = 'C')

# print(elem_iters_list)
# for elem_lines in elem_iters_list:
#     print("Checking elements")
#     print("Elements ======== \n")
#     print(elem_lines)
# =================================
#   END OF DATA EXTRACTION 
# =================================

# =================================
#   SOLVE BY PARTITION
# =================================
for i in np.arange(len(FF_iters_list)):
# for i in np.arange(1):
    print("Running through FF iteration list\n")
    KK = KK_iters_list[i]
    FF = FF_iters_list[i]
    # elements = elem_iters_list[i]
    fileSaveName = "standAloneVerticalHole" + str(N) + "-" + str(i) +".png"
    
    UU =  [1] * 3*nb_nodes
    count = 0
    for ii in np.arange(3*nb_nodes):
        # print("ii: ", ii ,"\n")
        # print("ii//3: ", ii//3 ,"\n")
        # print("Length of nodeslines:", len(nodes_lines), "\n")
        ii_mod = ii-count;
    #     print("ii_mod: ",ii_mod ,"\n")
        # if (KK[ii_mod][ii_mod] == 1):
        if (nodes_lines[ii//3][2] == 0) or (KK[ii_mod][ii_mod] == 0):

            
            
            
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
            
    if (np.linalg.det(KK) == 0):
        print("Matrix is non-singular")
        continue        
    else:
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
    # for element in elements:
    #     print("Element: ", element)
    #     level = ((element-1)-(element-1)%(N*N))/(N*N)
    #     j = (element-1)%N
    #     print('level',level, 'j',j, 'i', ((element-1)-(element-1)%N-level*N*N)/N, '\n')
        
    X, Y = np.meshgrid(np.arange(0, N+1), np.arange(0, N+1))
    Z = 0*X

    fig = plt.figure(figsize=(14,15))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, alpha=0.3)  # the horizontal plane
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
            
        # we plot only if its the firt layer or it is a node that moved. 
        if (UU[3*node+0] != 0):
            
            ax.scatter(nodes_lines[node][0],nodes_lines[node][1],nodes_lines[node][2],color='black' ,linewidths=1) # plot the point (2,3,4) on the figure
            ax.scatter(nodes_lines[node][0]+ UU[3*node+0],nodes_lines[node][1]+UU[3*node+1],nodes_lines[node][2]+UU[3*node+2],color='red' ,linewidths=1) # plot the point (2,3,4) on the figure
            #arw = Arrow3D([nodes_lines[node][0],nodes_lines[node][0]+ soln_lines[3*node+0]],[nodes_lines[node][1],nodes_lines[node][1]+soln_lines[3*node+1]],[nodes_lines[node][2],nodes_lines[node][2]+soln_lines[3*node+2]], arrowstyle="->", color="purple", lw = 1, mutation_scale=10)
        
        # GLOBAL
        #ax.scatter(soln_lines[3*node+0],soln_lines[3*node+1],soln_lines[3*node+2],color='red' ) # plot the point (2,3,4) on the figure
        #arw = Arrow3D([nodes_lines[node][0], soln_lines[3*node+0]],[nodes_lines[node][1],soln_lines[3*node+1]],[nodes_lines[node][2],soln_lines[3*node+2]], arrowstyle="->", color="purple", lw = 1, mutation_scale=10)
        #ax.add_artist(arw)
        
        ax.plot([nodes_lines[node][0],nodes_lines[node][0]+ UU[3*node+0]],[nodes_lines[node][1],nodes_lines[node][1]+UU[3*node+1]],[nodes_lines[node][2],nodes_lines[node][2]+UU[3*node+2]], 'b-');
        



    #ax.contour(X,Y,Z, color='black')
    ax.view_init(20,80)
    # ax.view_init(80,0)

    # top view ----------
    # ax.view_init(90,45)

    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$')
    ax.set_zlim([-1,N+1])
    # plt.savefig('result.png')
    print("gonna save the stupid file\n")
    plt.savefig(fileSaveName)
    # plt.show()
