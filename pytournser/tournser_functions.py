import numpy as np
from pytourn import run_tournser

def tournser_edges(vertex_values, edge_list, filtration='vertexMax', approx=False, approx_val=100000, count_only=False):

    return run_tournser(vertex_values, edge_list, filtration, approx, approx_val, count_only, False, 'null')



def tournser(vertex_values, adjacency_matrix, filtration='vertexMax', approx=False, approx_val=100000, count_only=False):

    return tournser_edges(vertex_values, [[i[0],i[1],adjacency_matrix[i[0],i[1]]] for i in np.transpose(np.nonzero(adjacency_matrix))], filtration, approx, approx_val, count_only)



def tournser_file(in_file, filtration='vertexMax', approx=False, approx_val=100000, count_only=False):

    return run_tournser(np.array([]), np.array([]), filtration, approx, approx_val, count_only, True, in_file)


#TODO implement print and print-dist