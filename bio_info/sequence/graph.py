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


def eulerian_cycle(graph):
    """
    Find a Eulerian cycle in a graph, if it exists.

    :param graph: a de Bruijn graph of type Graph

    :return: the edges in :param graph forming the cycle
    """
    if type(graph) != Graph:
        raise TypeError("param graph must be of type Graph")

    unexplored_edge_nodes = list()  # list of nodes with unexplored edges
    visited_edges = list()  # list of visited edges
    start_node = graph.nodes[random.randint(0, len(graph.nodes) - 1)]  # random start node

    # Start building the cycle
    while len(visited_edges) < len(graph.edges):
        edges = [e for e in start_node.edges if e not in visited_edges]

        if len(edges) == 0:
            # Pick new starting node
            if len(unexplored_edge_nodes) == 0:
                # There is no Eulerian cycle
                return
            else:
                # Start with previously explored node with additional unexplored edges
                start_node = unexplored_edge_nodes.pop()
                node_visited_edges = [e for e in start_node.edges if e in visited_edges]
                index = visited_edges.index(node_visited_edges[0])
                visited_edges = visited_edges[index:] + visited_edges[:index]

        elif len(edges) == 1:
            if start_node in unexplored_edge_nodes:
                # Remove from list of nodes with unexplored edges
                unexplored_edge_nodes.remove(start_node)
            visited_edges.append(edges[0])
            start_node = edges[0].dest_node

        elif len(edges) > 1:
            if start_node not in unexplored_edge_nodes:
                # Add to list of nodes with unexplored edges
                unexplored_edge_nodes.append(start_node)
            edge = edges[random.randint(0, len(edges) - 1)]
            visited_edges.append(edge)
            start_node = edge.dest_node

    return visited_edges
