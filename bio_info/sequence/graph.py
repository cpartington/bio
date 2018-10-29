import random

from bio_info.util.graph import Graph


def overlap_graph(pattern_list):
    # TODO change to Graph form
    adjacencies = dict()

    for i in range(len(pattern_list)):
        # Setup
        pattern = pattern_list[i]
        adjacencies[pattern] = list()
        suffix = pattern[1:]
        other_patterns = pattern_list[:i] + pattern_list[i+1:]
        # Check other patterns for matching prefixes
        for p in other_patterns:
            prefix = p[:-1]
            if suffix == prefix:
                # Add to adjacency list
                adjacencies.get(pattern).append(p)

    return adjacencies


def de_bruijn(pattern_list):
    """
    Builds a De Bruijn graph from a given set of patterns.

    :param pattern_list: a list of DNA sequences

    :return: a Graph object with one edge for each DNA sequence
             in :param pattern_list
    """
    g = Graph()
    label_count = dict()

    # Add the k-mers
    for i in range(len(pattern_list)):
        prefix = g.add_node(pattern_list[i][:-1])
        suffix = g.add_node(pattern_list[i][1:])
        g.add_edge(prefix, suffix, pattern_list[i])

        # Add to label count dictionary
        if prefix.label in label_count:
            label_count.get(prefix.label).append(prefix)
        else:
            label_count[prefix.label] = [prefix]
        if suffix.label in label_count:
            label_count.get(suffix.label).append(suffix)
        else:
            label_count[suffix.label] = [suffix]

    # Use label counts
    for label in label_count.keys():
        label_nodes = label_count.get(label)
        if len(label_nodes) > 1:
            g.merge_nodes(label_nodes, label)

    return g


def eulerian_path(input_graph):
    """
    Find a Eulerian path in a graph, if it exists. If it starts and ends
    at the same node, it's a Eulerian cycle.

    :param graph: a de Bruijn graph of type Graph

    :return: the edges in :param graph forming the path
    """
    if type(input_graph) != Graph:
        raise TypeError("param graph must be of type Graph")

    graph = input_graph.copy()

    in_edges = dict()
    out_edges = dict()
    # Initialize dicts
    for node in graph.nodes:
        in_edges[node] = 0
    # Fill dicts
    for node in graph.nodes:
        out_edges[node] = len(node.edges)
        for edge in node.edges:
            in_edges[edge.dest_node] += 1
    # Check for unbalanced nodes
    unbalanced_nodes = [n for n in graph.nodes if in_edges[n] != out_edges[n]]
    if len(unbalanced_nodes) == 2:
        if in_edges[unbalanced_nodes[0]] == 0:
            start = unbalanced_nodes[0]
        else:
            start = unbalanced_nodes[1]
        return _eulerian_path_(graph, start, list())
    elif len(unbalanced_nodes) == 0:
        start = graph.nodes[random.randint(0, len(graph.nodes))]
        return _eulerian_path_(graph, start, list())
    else:
        return


def _eulerian_path_(graph, node, cycle):
    cycle.append(node)
    if len(node.edges) == 0:
        return cycle

    while len(node.edges) > 0:
        # Get destination node of first edge
        tmp_edge = node.edges[0]
        graph.remove_edge(tmp_edge)
        # Recursive call
        sub_cycle = _eulerian_path_(graph, tmp_edge.dest_node, list())
        cycle = cycle[:1] + sub_cycle + cycle[1:]

    return cycle
