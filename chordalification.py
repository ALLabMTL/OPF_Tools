import networkx as nx
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
    
    W_Relaxations = SDR | TCR | STCR
    Chordal_Relaxations = Chordal_MFI | Chordal_AMD | Chordal_MD | Chordal_MCS_M
    Non_Chordal_Relaxations = SDR | SOCR | TCR | STCR

def AMD_Cholesky(Network):
    """Chordal extension using Approximate Minimum Degree ordering and Cholesky factorization.
    
    Args:
        Network (nx.Graph): Original graph.
    
    Returns:
        tuple: Extended chordal graph, elimination ordering and cliques.
    """
    from cvxopt import spmatrix, amd, matrix
    import chompack as cp
    from scipy.sparse import csr_matrix
    p = amd.order                                                           # Approximate Minimum Degree ordering
    Network_copy = Network.copy()
    Network_copy.add_edges_from([(i,i) for i in Network_copy.nodes])        # Add self loops for symbolic factorization
    Sparsity_pattern = nx.to_scipy_sparse_array(Network_copy)
    coo = Sparsity_pattern.tocoo()
    
    # Symbolic Cholesky factorization
    SP_matrix = cp.symbolic(spmatrix(coo.data.tolist(), coo.row.tolist(), coo.col.tolist(), size=Sparsity_pattern.shape),p)
    SP_filled = SP_matrix.sparsity_pattern(reordered=False)                 # Sparsity pattern with filled in values
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
        number_of_missing_edges_per_node = dict()
        missing_edges_per_node = dict()
        for node in Gk:
                missing_edges = missing_edges_for_clique(Gk, node)
                missing_edges_per_node[node] = missing_edges
                if missing_edges is None:
                    number_of_missing_edges_per_node[node] = 0
                else:
                    number_of_missing_edges_per_node[node] = len(missing_edges)
        for k in range(n):
            vertex_min_missing_edges = min(number_of_missing_edges_per_node, key=number_of_missing_edges_per_node.get)      
            Ordering[vertex_min_missing_edges] = k
            missing_edges = missing_edges_per_node[vertex_min_missing_edges]
            if missing_edges is None:
                nodes_affected = set(Gk.neighbors(vertex_min_missing_edges))
            else:
                nodes_affected = nx.single_source_shortest_path_length(Gk, vertex_min_missing_edges, cutoff=2)
                nodes_affected.pop(vertex_min_missing_edges)
                nodes_affected = set(nodes_affected.keys())
                Gk.add_edges_from(missing_edges)
                Gplus.add_edges_from(missing_edges)
            Gk.remove_node(vertex_min_missing_edges)
            number_of_missing_edges_per_node.pop(vertex_min_missing_edges)
            missing_edges_per_node.pop(vertex_min_missing_edges)
            for node in nodes_affected:
                missing_edges = missing_edges_for_clique(Gk, node)
                missing_edges_per_node[node] = missing_edges
                if missing_edges is None:
                    number_of_missing_edges_per_node[node] = 0
                else:
                    number_of_missing_edges_per_node[node] = len(missing_edges)
        return Gplus, Ordering
    
    
def get_clique_graph(cliques):
    """Construct a clique graph from a list of cliques.

    Args:
        cliques (list): List of cliques.

    Returns:
        nx.Graph: Clique graph.
    """
    from itertools import combinations
    clique_graph = nx.Graph()
    clique_graph.add_nodes_from(cliques, type="clique")
    for edge in combinations(cliques, 2):
        set_edge_0 = set(edge[0])
        set_edge_1 = set(edge[1])
        if not set_edge_0.isdisjoint(set_edge_1):
            sepset = tuple(sorted(set_edge_0.intersection(set_edge_1)))
            clique_graph.add_edge(edge[0], edge[1], weight=len(sepset), sepset=sepset)
    return clique_graph

def get_clique_tree(cliques, merge_strategy = 'No Merge', chordal_extension = None):
    clique_graph = get_clique_graph(cliques)
    if merge_strategy == 'No Merge':
        return cliques, nx.maximum_spanning_tree(clique_graph)
    try:
        merge_strat = MERGING_STRATEGY[merge_strategy]
    except KeyError as err:
        msg = f"{merge_strategy} is not a valid merging strategy."
        raise ValueError(msg) from err
    
    merged_clique_tree = merge_strat(clique_graph, chordal_extension)
    return list(merged_clique_tree.nodes), merged_clique_tree
    
def Parent_Child_Merge(clique_graph, chordal_extension):
    clique_tree = nx.maximum_spanning_tree(clique_graph)
    """ 
    TODO
    Sun and Andersen - Decomposition in conic optimization with partially separable structure (2014)
    """
    clique_tree_merged = clique_tree.copy()
    return clique_tree_merged

def Clique_Graph_Merge(clique_graph, chordal_extension):
    G = clique_graph.copy()
    edge_list = list(G.edges)
    edge_to_weight = dict()
    for edge in edge_list:
        weight = len(edge[0])**3 + len(edge[1])**3 - len(set(edge[0]).union(set(edge[1])))**3
        edge_to_weight[frozenset(edge)] = weight
    while any(weight > 0 for weight in edge_to_weight.values()):
        for edge, weight in edge_to_weight.items():
            if weight > 0:
                for node in edge:
                    for neighbor in G.neighbors(node):
                        F = frozenset((node, neighbor))
                        if not F == edge:
                            edge_to_weight.pop(F)
                edge_to_weight.pop(edge)
                edge_tuple = tuple(edge)
                clique_union = tuple(set(edge_tuple[0]).union(set(edge_tuple[1])))
                nx.contracted_edge(G, edge_tuple, self_loops=False, copy=False)
                nx.relabel_nodes(G,{edge_tuple[0]:clique_union}, copy=False)
                for new_edge in G.edges(clique_union):
                    edge_to_weight[frozenset(new_edge)] = len(new_edge[0])**3 + len(new_edge[1])**3 - len(set(new_edge[0]).union(set(new_edge[1])))**3
                    set_edge_0 = set(new_edge[0])
                    set_edge_1 = set(new_edge[1])
                    sepset = tuple(sorted(set_edge_0.intersection(set_edge_1)))
                    nx.set_edge_attributes(G, {new_edge: {'weight': len(sepset), 'sepset': sepset}})
                # Update chordal extension
                make_subgraph_complete(chordal_extension, clique_union)
                break
    return nx.maximum_spanning_tree(G)
    

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

def make_subgraph_complete(G, subgraph_nodes):
    for i, u in enumerate(subgraph_nodes):
        for v in subgraph_nodes[i+1:]:
            if not G.has_edge(u, v):
                G.add_edge(u, v)
    return G

def getChordalExtension(Network, relaxation, merge_strategy):
    if relaxation in Relaxation.Non_Chordal_Relaxations:
        return None
    if relaxation == Relaxation.Chordal_AMD:
        chordal_extension, ordering, cliques = AMD_Cholesky(Network)
    if relaxation == Relaxation.Chordal_MD:
        chordal_extension, ordering = elimination_game(Network, 'MD')
        cliques = nx.chordal_graph_cliques(chordal_extension)
    if relaxation == Relaxation.Chordal_MFI:
        chordal_extension, ordering = elimination_game(Network, 'MFI')
        cliques = nx.chordal_graph_cliques(chordal_extension)
    if relaxation == Relaxation.Chordal_MCS_M:
        chordal_extension, ordering = nx.complete_to_chordal_graph(Network)
        cliques = nx.chordal_graph_cliques(chordal_extension)
    
    ordered_cliques = [tuple(sorted(i)) for i in cliques]
    cliques, clique_tree = get_clique_tree(ordered_cliques, merge_strategy, chordal_extension)
    
    cliques = [tuple(sorted(i)) for i in cliques]
    nx.relabel_nodes(clique_tree, {i:tuple(sorted(i)) for i in clique_tree.nodes}, copy=False)
    
    return chordal_extension, ordering, cliques, clique_tree

MERGING_STRATEGY = {
    "ParentChild": Parent_Child_Merge,
    "CliqueGraph": Clique_Graph_Merge,
}


