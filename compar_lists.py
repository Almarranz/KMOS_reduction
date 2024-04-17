"""
compare_lists.py
Almarranz 2021
Compare coordinates of two lists and extract the rows whose coordinates columns distance are smaller than variable 'distancia'.
If more than two points in the second list are closer than 'distancia' to a point in the first list,the closest one is choosen 
Input parametres: list1, list2, col_index x1, col_index y1, col_index x2, col_index y2, max distance for a positive match
It computes euclidean distances
Both lists must be numpy arrays
"""
from scipy.spatial import distance
import numpy as np
import time
def compare_lists(list1,list2,ls1_x,ls1_y,ls2_x,ls2_y,distancia):
    ls1_com = np.empty((0,list1.shape[1]))
    ls2_com = np.empty((0,list2.shape[1]))
    coord1 =np.array([list1[:,ls1_x],list1[:,ls1_y]]).T
    coord2 =np.array([list2[:,ls2_x],list2[:,ls2_y]]).T
    # coord1 =np.array([list1[:,0],list1[:,1]]).T
    # coord2 =np.array([list2[:,0],list2[:,1]]).T
    tic1 = time.perf_counter()
    dist_A = distance.cdist(coord1[:,0:2],coord2[:,0:2], 'euclidean')
    toc1 = time.perf_counter()
    print('distance matrix tooks %.2f seconds'%(toc1-tic1))
    # dist   = distance.cdist(coord1[:,0:2],coord2[:,0:2], 'euclidean')
    dist = dist_A # for some reason when doing this python also modify dist_A
    # sys.exit()
    # dis the is n_elem(l1)*n_elem[l2] matrix. Each row "i" is the distance between
    # the i-th element in l1 with all the elements in l2
    l1 =[]
    l2 =[]
    # New aproach:
    #     -Flat and sort the distance matrix
    #     -Find the position of the first sorted value in the flatted matrix (row, col)
    #      -If this value is a minimun in its row and column, store both point as 
    #         a match, mark as inf all values in these rows and columns and move to
    #         the nex value
    #      -If not, move to the next sorted distance
    #     -Repeat till the value of the sorted distance is bigger than min distance
    marked = np.zeros((dist.shape[0],dist.shape[1]))
    dist_fs = np.sort(np.reshape(dist_A,dist_A.shape[0]*dist_A.shape[1]))
    tic = time.perf_counter()
    for i in range(len(np.where(dist_fs<distancia)[0])):
        inds = np.argwhere(dist == dist_fs[i])
        for j in range(len(inds)):
            if np.isfinite(dist[inds[j][0],inds[j][1]]) and np.isfinite(dist[inds[j][0],inds[j][1]]):
                if (dist[inds[j][0],inds[j][1]] <= min(dist[inds[j][0]])) and (dist[inds[j][0],inds[j][1]] <= min(dist[:,inds[j][1]])):
                    dist[inds[j][0]] = float('inf') 
                    dist[:,inds[j][1]] = float('inf') 
                    l1.append(inds[j][0])
                    l2.append(inds[j][1])
                    # print(inds[j][0],inds[j][1])
    return list1[l1],list2[l2], l1, l2 # Common lists and their indices
