import traceback

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

    """ Graph Operations """

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
        from_node.out_degree += 1
        dest_node.in_degree += 1
        return edge

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
                edge.dest_node.in_degree -= 1
                edge.from_node.out_degree -= 1
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
                master_node.out_degree += 1
            if edge.dest_node in node_list:
                edge.dest_node = master_node
                master_node.in_degree += 1
        safe_removal = self.remove_nodes(node_list)
        return safe_removal

    def merge_edges(self, edge_list, merge_nodes=False):
        """
        Merges a list of edges by merging the associated nodes
        if desired and then updating the edge label.

        :param edge_list: a list of edges to merge
        :param merge_nodes

        :return the new merged node
        """
        # TODO finish function
        new_label = []
        if merge_nodes:
            nodes = list()
        for edge in edge_list:
            new_label.append(edge.label[0])
            if merge_nodes:
                nodes.append(edge.from_node)
        new_label.append(edge_list[-1].label[1:])
        if merge_nodes:
            nodes.append(edge_list[-1].dest_node)
            # Add new master node
            # Remove any edges pointing within the merged nodes
        # Add one master edge pointing to itself
        new_edge = self.add_edge(edge_list[0].from_node, edge_list[-1].dest_node, "".join(new_label))
        for edge in edge_list:
            self.remove_edge(edge)
        return new_edge

    def merge_equal_edges(self):
        for node in self.nodes:
            i = 0
            while i < len(node.edges):
                equal_edges = [e for e in node.edges if e.dest_node == node.edges[i].dest_node]
                for edge in equal_edges[1:]:
                    self.remove_edge(edge)
                i += 1

    def remove_edge(self, edge):
        edge.from_node.edges.remove(edge)
        edge.dest_node.in_degree -= 1
        edge.from_node.out_degree -= 1
        self.edges.remove(edge)

    def copy(self):
        g = Graph()
        nodes = dict()
        for node in self.nodes:
            new_node = g.add_node(node.label)
            nodes[id(node)] = new_node
        for edge in self.edges:
            from_node = nodes[id(edge.from_node)]
            dest_node = nodes[id(edge.dest_node)]
            g.add_edge(from_node, dest_node, edge.label)
        return g

    def combine_simple_nodes(self):
        extra_edges = list()
        paths = [path for path in self.maximal_nonbranching_paths() if len(path) > 1]
        for path in paths:
            from_node = path[0].from_node
            for edge in path:
                # Update edge to be merged k-mer
                edge.dest_node.edges[0].label = "".join([edge.label[0], edge.dest_node.edges[0].label])
                # Get new node label
                new_label = "".join([from_node.label[0], edge.dest_node.label])
                # Merge source node and destination node
                new_node = self.merge_nodes([from_node, edge.dest_node], new_label)
                extra_edges.append(edge)
                from_node = new_node
        for edge in extra_edges:
            self.remove_edge(edge)

    def remove_tips(self):
        pass

    def remove_bubbles(self):
        pass

    """ Graph Utilization """

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
                if len(node.label) > 10:
                    label = "{}...{}".format(node.label[:3], node.label[-3:])
                    v_prop[v] = label
                else:
                    v_prop[v] = node.label
            vertices[node.label] = v

        for edge in self.edges:
            e = graph.add_edge(vertices[edge.from_node.label],
                               vertices[edge.dest_node.label])
            if edge_labels:
                if len(edge.label) > 15:
                    label = "{}...{}".format(edge.label[:4], edge.label[-4:])
                    e_prop[e] = label
                else:
                    e_prop[e] = edge.label

        graph_draw(graph, vertex_text=v_prop, edge_text=e_prop, vertex_font_size=18,
                   output_size=(5000, 5000), output=output_file, fit_view=True)

    def maximal_nonbranching_paths(self):
        """
        Finds maximal non-branching paths in the graph. A non-branching
        path is a path where the middle nodes have an in-degree and out-
        degree of one.

        :return: the list of all maximal non-branching paths through the
                 graph, in the form of a list of lists of Nodes
        """
        paths = []
        for node in self.nodes:
            if node.in_degree != 1 or node.out_degree != 1:
                if node.out_degree > 0:
                    for edge in node.edges:
                        non_b_path = [edge]
                        new_node = edge.dest_node
                        while new_node.in_degree == new_node.out_degree == 1:
                            non_b_path += [new_node.edges[0]]
                            new_node = new_node.edges[0].dest_node
                        paths += [non_b_path]
        return paths

    def contigs(self):
        """
        Generate contigs formed by the graph.

        :return: the list of created contigs
        """
        paths = self.maximal_nonbranching_paths()
        contigs = list()
        for path in paths:
            if len(path) > 1:
                contig = list()
                for edge in path:
                    contig.append(edge.label[0])
                contig.append(path[-1].label[1:])
                contigs.append("".join(contig))
            else:
                contigs.append(path[0].label)
        return contigs


class Node:
    def __init__(self, label=None):
        self.label = label
        self.edges = list()
        self.in_degree = 0
        self.out_degree = 0


class Edge:
    def __init__(self, from_node, dest_node, label=None):
        self.label = label
        self.from_node = from_node
        self.dest_node = dest_node