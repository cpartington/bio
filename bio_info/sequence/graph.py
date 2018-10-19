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
