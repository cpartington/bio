import random

try:
    from graph_tool.all import Graph as GTGraph
    from graph_tool.all import graph_draw
    graphing = True
except ImportError:
    graphing = False
    GTGraph = None
    graph_draw = None


class Graph:
    def __init__(self):
        self.nodes = list()
        self.edges = list()

    def add_node(self, label):
        node = Node(label)
        self.nodes.append(node)
        return node

    def add_edge(self, from_node, dest_node, label=None):
        """
        Adds an edge to the graph and its origin node.

        :param from_node: the origin of the edge
        :param dest_node: the destination of the edge
        :param label: the label for the edge
        """
        edge = Edge(from_node, dest_node, label)
        self.edges.append(edge)
        from_node.edges.append(edge)

    def remove_node(self, node):
        """
        Removes a node and its edges from the graph.

        :param node: a Node object

        :return: True if no edges had to removed when removing
                 the node
        """
        # TODO add error checking
        no_change = True
        for edge in self.edges:
            if edge.dest_node == node or edge.from_node == node:
                self.edges.remove(edge)
                no_change = False
        self.nodes.remove(node)
        return no_change

    def remove_nodes(self, node_list):
        """
        Removes a list of nodes and their edges from the graph.

        :param node_list: a list of Node objects

        :return: True if no edges had to be removed when removing
                 any of the nodes
        """
        no_change = True
        for node in node_list:
            safely_removed = self.remove_node(node)
            if safely_removed is False:
                no_change = False
        return no_change

    def merge_nodes(self, node_list, label=None):
        """
        Merges a list of nodes into a single node and updates the
        graph nodes and edges.

        :param node_list: a list of Node objects
        :param label: the label for the new node
        """
        master_node = self.add_node(label)
        for edge in self.edges:
            if edge.from_node in node_list:
                edge.from_node = master_node
                master_node.edges.append(edge)
            if edge.dest_node in node_list:
                edge.dest_node = master_node
        safe_removal = self.remove_nodes(node_list)
        assert(safe_removal is True)

    def draw(self, output_file, node_labels=True, edge_labels=False):
        """
        Uses the graph-tool library to generate a graphical version
        of the graph and saves it to a file. Installation:
        https://graph-tool.skewed.de/

        :param output_file: the name to save the image as
        :param node_labels: if True, include node labels in graph image
        :param edge_labels: if True, include edge labels in graph image
        """
        if graphing is False:
            print("graph-tool library necessary for this function")
            print("https://graph-tool.skewed.de/")
            return
        graph = GTGraph()
        v_prop = None
        e_prop = None
        if node_labels:
            v_prop = graph.new_vertex_property("string")
        if edge_labels:
            e_prop = graph.new_edge_property("string")
        vertices = dict()

        for node in self.nodes:
            v = graph.add_vertex()
            if node_labels:
                v_prop[v] = node.label
            vertices[node.label] = v

        for edge in self.edges:
            e = graph.add_edge(vertices[edge.from_node.label],
                               vertices[edge.dest_node.label])
            if edge_labels:
                e_prop[e] = edge.label

        graph_draw(graph, vertex_text=v_prop, edge_text=e_prop, vertex_font_size=18,
                   output_size=(500, 500), output=output_file)


class Node:
    def __init__(self, label=None):
        self.label = label
        self.edges = list()


class Edge:
    def __init__(self, from_node, dest_node, label=None):
        self.label = label
        self.from_node = from_node
        self.dest_node = dest_node


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
    unexplored_edge_nodes = list()  # list of nodes with unexplored edges
    visited_edges = list()  # list of visited edges
    start_node = graph.nodes[random.randint(0, len(graph.nodes) - 1)]  # random start node
    while len(visited_edges) < len(graph.edges):
        edges = [e for e in start_node.edges if e not in visited_edges]
        if len(edges) == 0:
            # Pick new starting node
            if len(unexplored_edge_nodes) == 0:
                # There is no Eulerian cycle
                return
            else:
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
