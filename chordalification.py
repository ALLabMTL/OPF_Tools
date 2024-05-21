import networkx as nx
import matplotlib.pyplot as plt
from cvxopt import spmatrix, amd, matrix
import chompack as cp
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
from enum import Flag, auto

class Relaxation(Flag):
    SDR = auto()
    SOCR = auto()
    TCR = auto()
    STCR = auto()
    Chordal_MFI = auto()
    Chordal_AMD = auto()
    Chordal_MD = auto()
    Chordal_MCS_M = auto()
    
    Non_Chordal_Relaxations = SDR | SOCR | TCR | STCR
    Chordal_Relaxations = Chordal_MFI | Chordal_AMD | Chordal_MD | Chordal_MCS_M

def AMD_Cholesky(Network):
    """Chordal extension using Approximate Minimum Degree ordering and Cholesky factorization.
    
    Args:
        Network (nx.Graph): Original graph.
    
    Returns:
        tuple: Extended chordal graph, elimination ordering and cliques.
    """
    p = amd.order                                                 # Approximate Minimum Degree ordering
    Network.add_edges_from([(i,i) for i in Network.nodes])        # Add self loops for symbolic factorization
    Sparsity_pattern = nx.to_scipy_sparse_array(Network)
    coo = Sparsity_pattern.tocoo()
    
    # Symbolic Cholesky factorization
    SP_matrix = cp.symbolic(spmatrix(coo.data.tolist(), coo.row.tolist(), coo.col.tolist(), size=Sparsity_pattern.shape),p)
    SP_filled = SP_matrix.sparsity_pattern(reordered=False)         # Sparsity pattern with filled in values
    SP_csr = csr_matrix((list(SP_filled.V), (list(SP_filled.I), list(SP_filled.J))), shape=SP_filled.size)
    chordal_extension = nx.from_scipy_sparse_array(SP_csr)
    chordal_extension.remove_edges_from(nx.selfloop_edges(chordal_extension))
    cliques = SP_matrix.cliques(reordered = False)
    cliques = [tuple(sorted(i)) for i in cliques]
    return chordal_extension, SP_matrix.p, cliques


def elimination_game(Network, heuristic = 'MD'):
    """Elimination game algorithm to find a chordal extension of a graph.

    Args:
        Network (nx.Graph): Original graph.
        heuristic (str, optional): Heuristic used to select the node to eliminate at each step.
        'MD' for Minimum Degree, 'MFI' for Minimum Fill-In. Defaults to 'MD'.

    Returns:
        _type_: Extended chordal graph and elimination ordering.
    """
    n = Network.number_of_nodes()
    Gk = Network.copy()
    Gplus = Network.copy()
    Ordering = dict()
    if heuristic == 'MD':
        for k in range(n):
            degrees = dict(nx.degree(Gk))
            vertex_min_degree = min(degrees, key=degrees.get)
            Ordering[vertex_min_degree] = k
            missing_edges = missing_edges_for_clique(Gk, vertex_min_degree)
            Gk.remove_node(vertex_min_degree)
            if missing_edges is None:
                continue
            Gk.add_edges_from(missing_edges)
            Gplus.add_edges_from(missing_edges)
        return Gplus, Ordering
    if heuristic == 'MFI':
        for k in range(n):
            number_of_missing_edges_per_node = dict()
            for node in Gk:
                missing_edges = missing_edges_for_clique(Gk, node)
                if missing_edges is None:
                    number_of_missing_edges_per_node[node] = 0
                else:
                    number_of_missing_edges_per_node[node] = len(missing_edges)
            vertex_min_missing_edges = min(number_of_missing_edges_per_node, key=number_of_missing_edges_per_node.get)
            Ordering[vertex_min_missing_edges] = k
            missing_edges = missing_edges_for_clique(Gk, vertex_min_missing_edges)
            Gk.remove_node(vertex_min_missing_edges)
            if missing_edges is None:
                continue
            Gk.add_edges_from(missing_edges)
            Gplus.add_edges_from(missing_edges)
        return Gplus, Ordering

def missing_edges_for_clique(graph, node):
    """
    Find the set of edges needed to make the neighborhood of a node a clique.

    Parameters:
        graph (nx.Graph): Input graph.
        node: Node whose neighborhood needs to form a clique.

    Returns:
        set: Set of edges needed to make the neighborhood of the node a clique.
    """
    neighborhood = set(graph.neighbors(node))
    
    if is_subclique(graph, neighborhood):
        return None
    
    missing_edges = set()

    # Check for missing edges to form a clique
    for u in neighborhood:
        for v in neighborhood:
            if u < v and not graph.has_edge(u, v):
                missing_edges.add((u, v))

    return missing_edges

def is_subclique(G,nodelist):
    """Tests if the subgraph induced by nodelist is a clique.

    Args:
        G (nx.Graph): Input graph.
        nodelist (): set of nodes.

    Returns:
        boolean: True if the subgraph induced by nodelist is a clique.
    """
    H = G.subgraph(nodelist)
    n = len(nodelist)
    return H.size() == n*(n-1)/2

def getChordalExtension(Network, relaxation):
    if relaxation in Relaxation.Non_Chordal_Relaxations:
        return None
    if relaxation == Relaxation.Chordal_AMD:
        return AMD_Cholesky(Network)
    if relaxation == Relaxation.Chordal_MD:
        chordal_extension, ordering = elimination_game(Network, 'MD')
        cliques = nx.chordal_graph_cliques(chordal_extension)
    if relaxation == Relaxation.Chordal_MFI:
        chordal_extension, ordering = elimination_game(Network, 'MFI')
        cliques = nx.chordal_graph_cliques(chordal_extension)
    if relaxation == Relaxation.Chordal_MCS_M:
        chordal_extension, ordering = nx.complete_to_chordal_graph(Network)
        cliques = nx.chordal_graph_cliques(chordal_extension)
    cliques = [tuple(sorted(i)) for i in cliques]
    return chordal_extension, ordering, cliques


