#!/usr/bin/python
# -*- coding: utf-8 -*-

########################################################################################################################
                                    # Problem 1
########################################################################################################################

#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import Counter, defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
from random import randrange
import snap

__author__ = 'dghosh'

COLLAB_GRAPH_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/ca-GrQc.txt"
NUM_NODES = 5242
NUM_EDGES = 14484


def problem_1_1():
    """
    Generate a random graph from both the Erd˝os-R´enyi (i.e., G(n, m)) and Small-World models and
    read in the collaboration network. Delete all of the self-edges in the collaboration network (there
    should be 14,484 total edges remaining).

    Plot the degree distribution of all three networks in the same plot on a log-log scale. In other words,
    generate a plot with the horizontal axis representing node degrees and the vertical axis representing
    the proportion of nodes with a given degree (by “log-log scale” we mean that both the horizontal
    and vertical axis must be in logarithmic scale). In one to two sentences, describe one key difference
    between the degree distribution of the collaboration network and the degree distributions of the
    random graph models.

    Hint: Use SNAP’s GenRndGnm function to generate from the Erd˝os-R´enyi model. However, you
    need to write a routine to generate a random instance of the Small-World model as described above.
    """

    # 1. Generate a random graph from Erd˝os-R´enyi (i.e., G(n, m))
    erdos_graph = snap.GenRndGnm(snap.PUNGraph, NUM_NODES, NUM_EDGES)  # 5242 nodes, 14484 undirected graph
    print "1. # of nodes in erdos graph is %d" % erdos_graph.GetNodes()

    # 2. Small world model
    small_wold = generate_small_world(NUM_NODES)

    # 3. read collab graph
    collab_graph = snap.LoadEdgeList(snap.PUNGraph, COLLAB_GRAPH_PATH, 0, 1)

    # 4. delete self edges from collab graph
    snap.DelSelfEdges(collab_graph)
    print "3. # of edges in collab graph (w/o self edges) are %d" % collab_graph.GetEdges()

    # 5. Plot the degree distribution of all three networks
    plot_graph([erdos_graph, small_wold, collab_graph])

    return [erdos_graph, small_wold, collab_graph]


def problem_1_2a(networks):

    pk_erdos = expected_degree_distribution(networks[0])
    pk_small_world = expected_degree_distribution(networks[1])
    pk_collab = expected_degree_distribution(networks[2])

    max_expected_degree = max(max(pk_erdos), max(pk_small_world), max(pk_collab))
    x = range(1, max_expected_degree+1)

    expected_deg_erdos = [pk_erdos[d] for d in xrange(1, max_expected_degree + 1)]
    expected_deg_small_world = [pk_small_world[d] for d in xrange(1, max_expected_degree + 1)]
    expected_deg_collab = [pk_collab[d] for d in xrange(1, max_expected_degree +1)]

    print "1.2a) Expected degree calculations.."
    print "* Erdos network : %f" % sum([m * n for m, n in zip(x, expected_deg_erdos)])
    print "* Small world network: %f" % sum([m * n for m, n in zip(x, expected_deg_small_world)])
    print "* Collab network : %f" % sum([m * n for m, n in zip(x, expected_deg_collab)])

    qk_erdos = excess_degree_distribution(networks[0])
    qk_small_world = excess_degree_distribution(networks[1])
    qk_collab = excess_degree_distribution(networks[2])

    max_excess_deg = max(max(qk_erdos), max(qk_small_world), max(qk_collab))
    x = range(1, max_excess_deg + 1)

    y_erdos = [qk_erdos[d] for d in xrange(1, max_excess_deg + 1)]
    y_small_world = [qk_small_world[d] for d in xrange(1, max_excess_deg + 1)]
    y_collab = [qk_collab[d] for d in xrange(1, max_excess_deg + 1)]

    print "1.2a) Expected excess degree calculations.."
    print "* Erdos network : %f" % sum([m * n for m, n in zip(x, y_erdos)])
    print "* Small world network: %f" % sum([m * n for m, n in zip(x, y_small_world)])
    print "* Collab network : %f" % sum([m * n for m, n in zip(x, y_collab)])

    plot_excess_deg_dist([[x, y_erdos], [x, y_small_world], [x, y_collab]])


def problem_1_2b():
    pass


def problem_1_3_2(networks):
    """
    :param networks: list: Array of networks - erdos, small world, collab in that order
    """
    coeffs = []
    for idx, network in enumerate(networks):
        coeffs.append([])
        for node in network.Nodes():
            ki = node.GetDeg()
            if ki < 2:
                coeffs[idx].append(0)  # ci = 0 if ki < 2
            else:
                neighbors = [edge for edge in node.GetOutEdges()]
                ei = 0  # ei is the number of edges between the neighbors of vi or node in this case
                for neighbor_1, neighbor_2 in combinations(neighbors, 2):
                    if network.IsEdge(neighbor_1, neighbor_2) and neighbor_1 != node.GetId() and \
                                    neighbor_2 != node.GetId():
                        ei += 1
                coeffs[idx].append(2 * ei / (1.0*ki*(ki-1)))

    print "1.3a) Clustering coefficient for Erdos network is %f" % (sum(coeffs[0]) / (1.0*NUM_NODES))
    print "1.3b) Clustering coefficient for small-world network is %f" % (sum(coeffs[1]) / (1.0*NUM_NODES))
    print "1.3c) Clustering coefficient for collab network is %f" % (sum(coeffs[2]) / (1.0*NUM_NODES))
    print "1.3d) Collab network has highest clustering coefficient = %f" % (sum(coeffs[2]) / (1.0*NUM_NODES))
    return coeffs


def expected_degree_distribution(graph):
    degrees = Counter([n.GetDeg() for n in graph.Nodes()])
    sum_degrees = float(sum(degrees.values()))
    pk = defaultdict(int)
    pk.update({x: (y / sum_degrees) for x,y in degrees.items()})
    return pk


def excess_degree_distribution(graph):
    qk = defaultdict(int)
    deg_dists = []
    for edge in graph.Edges():
        node_1, node_2 = edge.GetId()
        deg_dists.append(graph.GetNI(node_1).GetDeg())
        deg_dists.append(graph.GetNI(node_2).GetDeg())

    deg_dists = Counter(deg_dists)
    sum_q_prime = float(sum(deg_dists.values()))
    qk.update({x: (y/sum_q_prime) for x, y in deg_dists.items()})
    return qk


def generate_small_world(num_nodes):
    """
        Generate an instance from this model as follows: begin with
        n = 5242 nodes arranged as a ring, i.e., imagine the nodes form a circle and each node is
        connected to its two direct neighbors (e.g., node 399 is connected to nodes 398 and 400),
        giving us 5242 edges. Next, connect each node to the neighbors of its neighbors (e.g., node
        399 is also connected to nodes 397 and 401). This gives us another 5242 edges. Finally,
        randomly select 4000 pairs of nodes not yet connected and add an edge between them. In
        total, this will make m = 5242 · 2 + 4000 = 14484 edges. (Write code to construct instances
        of this model, i.e., do not call a SNAP function as above).
    """
    graph = snap.TUNGraph.New()
    for node in xrange(num_nodes):
        graph.AddNode(node)

    for node in xrange(1, num_nodes-1):  # [1 - 5241)
        graph.AddEdge(node, node-1)
        graph.AddEdge(node, node+1)
    graph.AddEdge(0, num_nodes-1)

    # connect each node to the neighbors of its neighbors
    for node in xrange(2, num_nodes-2):  # [2 - 5240)
        graph.AddEdge(node, node-2)  #
        graph.AddEdge(node, node+2)
    graph.AddEdge(1, num_nodes-1)
    graph.AddEdge(0, num_nodes-2)

    # randomly select 4000 pairs of nodes not yet connected and add an edge between them
    pair_count = 0
    while pair_count < 4000:
        node_1 = randrange(0, num_nodes)
        node_2 = randrange(0, num_nodes)
        if node_1 != node_2 and not graph.IsEdge(node_1, node_2):
            graph.AddEdge(node_1, node_2)
            pair_count += 1

    print "2. # of edges in small world graph is %d" % graph.GetEdges()
    return graph


def plot_graph(networks):
    plt.title("Log-log degree distribution of erdos, small-world and citation network")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("log(node_degree)")
    plt.ylabel("log(nodes_of_given_degree)")
    plt.ylim([0.0001, 1])
    plt.grid(True)
    labels = ["erdos", "small-world", "collab"]
    colors = ['red',  'blue',  'green']

    degree_dict = {}
    y_data = []
    for idx, graph in enumerate(networks):
        if idx == 2:
            # collab n/w
            deg_cnt_vector = snap.TIntPrV()
            snap.GetDegCnt(graph, deg_cnt_vector)
            x_data = map(lambda x: x.GetVal1(), deg_cnt_vector)
            y_data = map(lambda x: x.GetVal2() / (1.0*NUM_NODES), deg_cnt_vector)
        else:
            for node in xrange(NUM_NODES):
                degree = graph.GetNI(node).GetDeg()
                if degree in degree_dict:
                    degree_dict[degree] += 1
                else:
                    degree_dict[degree] = 1
            x_data = degree_dict.keys()
            y_data = map(lambda x: x / (1.0*NUM_NODES), degree_dict.values())

        plt.plot(x_data, y_data, linestyle="dashed", marker="o", color=colors[idx], label=labels[idx])

    plt.legend()
    plt.show()
    plt.clf()


def plot_excess_deg_dist(network_qk_list):
    plt.title("Log-log excess degree distribution for erdos, small-world and citation network")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("log(excess degree distribution)")
    plt.ylabel("log(proportion of nodes)")
    plt.ylim([0.0001, 1])
    plt.grid(True)
    labels = ["erdos", "small-world", "collab"]
    colors = ['red',  'blue',  'green']

    for idx, network_qk in enumerate(network_qk_list):
        plt.plot(network_qk_list[idx][0], network_qk_list[idx][1], linestyle="dashed", marker="o",
                 color=colors[idx], label=labels[idx])

    plt.legend()
    plt.show()
    plt.clf()

if __name__ == '__main__':
    print "1.1 Degree Distribution [10 points]\n"
    graphs = problem_1_1()
    print "1.2 Excess Degree Distribution [15 points]\n"
    print "1.2 (a) [10 points]\n"
    problem_1_2a(graphs)
    print "1.2 (b) [5 points]\n"
    problem_1_2b()
    print "1.3 Clustering Coefficient [10 points]\n"
    problem_1_3_2(graphs)

########################################################################################################################
                                    # Problem 2
########################################################################################################################
#!/usr/bin/python
# -*- coding: utf-8 -*-

import operator
import snap

__author__ = 'dghosh'

IMDB_ACTOR_EDGE_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/imdb_actor_edges.tsv"
IMDB_ACTOR_KEY_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/imdb_actors_key.tsv"


def init():
    graph = snap.LoadEdgeList(snap.PUNGraph, IMDB_ACTOR_EDGE_PATH, 0, 1)
    print "IMDB graph nodes count = %d" % graph.GetNodes()
    print "IMDB graph edge count = %d" % graph.GetEdges()
    wcc = snap.GetMxWcc(graph)
    actor_name_dict, imdb_db = generate_imdb_db()
    return wcc, actor_name_dict, imdb_db


def problem_2_1(graph, actor_name_dict, imdb_db):
    actor_degree_dict = calc_degree_centrality(graph)
    sorted_actor = sorted(actor_degree_dict.items(), key=operator.itemgetter(1))

    for idx, item in enumerate(reversed(sorted_actor[-20:])):
        actor_id, degree = item
        actor_name = actor_name_dict[actor_id]
        movie_count = imdb_db[actor_id][1]
        genre = imdb_db[actor_id][2]
        print ' #%d - %s:   %s    %s    %s' % (idx+1, actor_name, degree, movie_count, genre)


def problem_2_2(graph, actor_name_dict, imdb_db):
    betweenness_dict = calc_betweenness_centrality(graph)
    sorted_actor = sorted(betweenness_dict.items(), key=operator.itemgetter(1))
    for idx, item in enumerate(reversed(sorted_actor[-20:])):
        actor_id, centrality = item
        actor_name = actor_name_dict[actor_id]
        movie_count = imdb_db[actor_id][1]
        print ' #%d - %s:   %s    %s' % (idx+1, actor_name, centrality, movie_count)


def problem_2_3(graph, actor_name_dict, imdb_db):
    closeness_dict = calc_closeness_centrality(graph)
    sorted_actor = sorted(closeness_dict.items(), key=operator.itemgetter(1))
    for idx, item in enumerate(reversed(sorted_actor[-20:])):
        actor_id, centrality = item
        actor_name = actor_name_dict[actor_id]
        movie_count = imdb_db[actor_id][1]
        genre = imdb_db[actor_id][2]
        print ' #%d - %s:   %s    %s    %s' % (idx+1, actor_name, centrality, movie_count, genre)


def calc_degree_centrality(imdb_graph):
    actor_degree_dict = {}
    for actor in imdb_graph.Nodes():
        actor_id = actor.GetId()
        actor_degree_dict[actor_id] = snap.GetDegreeCentr(imdb_graph, actor_id)
    return actor_degree_dict


def calc_betweenness_centrality(imdb_graph):
    nodes = snap.TIntFltH()
    edges = snap.TIntPrFltH()
    snap.GetBetweennessCentr(imdb_graph, nodes, edges, 1.0)
    betweenness_dict = {node_id: nodes[node_id] for node_id in nodes}
    return betweenness_dict


def calc_closeness_centrality(imdb_graph):
    closeness_dict = {}
    for actor in imdb_graph.Nodes():
        actor_id = actor.GetId()
        closeness_dict[actor_id] = snap.GetClosenessCentr(imdb_graph, actor.GetId())
    return closeness_dict


def generate_imdb_db():
    actor_name_dict = {}
    imdb_db = {}
    count = 0
    with open(IMDB_ACTOR_KEY_PATH) as f:
        for line in f:
            tokens = line.split("\t")
            if not tokens[0].isdigit():
                continue
            actor_id = int(tokens[0])
            actor_name_dict[actor_id] = str.strip(tokens[1])  # store id vs name e.g. {15629: "Rudder, Michael (I)"}
            imdb_db[actor_id] = tokens[1:]
            count += 1
            if count < 3:
                print "Actor id: %d, token[1]: %s\n" % (actor_id, tokens[1])
                print "line: %s\n" % tokens[1:]
    return actor_name_dict, imdb_db


if __name__ == '__main__':
    WCC, id_name_dict, imdb = init()
    print "2.1 Degree centrality [10 points]\n"
    problem_2_1(WCC, id_name_dict, imdb)
    print "2.2 Betweenness centrality [10 points]\n"
    problem_2_2(WCC, id_name_dict, imdb)
    print "2.3 Closeness centrality [10 points]\n"
    problem_2_3(WCC, id_name_dict, imdb)

########################################################################################################################
                                    # Problem 3
########################################################################################################################
#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
import numpy as np
import snap

__author__ = 'dghosh'


def simulation():
    hT = 10  # height of tree, T is 10
    b = 2  # every non-leaf node has exactly 2 children
    k = 5  # every node has 5 outgoing edges / 5 neighbors
    alphas = np.arange(0.1, 11, 0.1)
    N = int((math.pow(b, hT+1) - 1) / (b - 1))  # a complete binary tree with non-leaf nodes being virtual
    E = N-1  # no of edges in complete binary tree with N nodes = (N-1), otherwise there will be cycle

    graph_db = {}
    for alpha in alphas:
        btree, leaf_nodes, virtual_nodes = create_complete_binary_tree(b, hT)
        edge_prob_table = generate_edge_probability_table(btree, leaf_nodes, virtual_nodes)
        graph = populate_network(btree, edge_prob_table)
        for sample in xrange(1000):
            # pick (s,t) randomly s != t
            # s, t = two randomly picked nodes in graph
            search_success = do_decentralized_search(graph)  # pass (s,t)
            graph_db[alpha] = (graph, search_success)

    path_dict = {}
    for alpha, graph_data in graph_db.iteritems():
        # pick 1000 pair of nodes from graph_data[0] i.e. the graph
        # compute average path length for successful searches
        mean_path_length = 0.0
        path_dict[alpha] = mean_path_length

    # 1. plot alpha vs mean_path_length
    # 2. plot alpha vs search success probability


def create_random_network(btree):
    """
        Use leaf nodes of btree to generate a random network with virtual nodes (i.e. the set of all non-leaf nodes)
    """
    graph = snap.PNGraph.New()
    return graph


def create_complete_binary_tree(b, height):
    btree = snap.GenTree(snap.PNGraph, b, height, True, False)
    leaf_nodes = set([node.GetId() for node in btree.Nodes() if node.GetOutDeg() == 0])
    virtual_nodes = {}
    for node in btree.Nodes():
        if node.GetOutDeg() > 0:
            virtual_nodes[node.GetId()] = node.GetOutEdges()
    print "Count of leaf nodes = %d" % len(leaf_nodes)  # snap.CntDegNodes(btree, 1)
    print "Count of virtual nodes = %d" % len(virtual_nodes.keys())  # snap.CntDegNodes(btree, 3)
    print "Node count of  tree  = %d" % btree.GetNodes()
    print "Edge count  = %d" % btree.GetEdges()
    return btree, leaf_nodes, virtual_nodes


def generate_edge_probability_table(btree, leaves, virtual_nodes, alpha, b):
    """
    * Create table based on b^(-alpha * h(v,w))
    :param graph:
    :param leaves:
    :param virtual_nodes:
    :param alpha:
    :param b:
    :return:
    """
    # create leaf_nodes X leaf_nodes sized 2D matrix filled with 0s
    # no self edge, thus diagonal should remain 0 i.e. (0,0) or (1,1) etc.
    # each row represents edge probabilities of node indexed by row-id to rest of the leaves
    table = np.zeros((len(leaves), len(leaves)), dtype=np.float)

    root_id = snap.GetTreeRootNId(btree)

    for r_idx, v in enumerate(leaves):
        for c_idx, w in enumerate(leaves):
            table[r_idx, c_idx] = calc_edge_prob(alpha, b, v, w, btree)
    return table


def calc_edge_prob(alpha, b, v, w, btree):
    """

    :param alpha:
    :param b:
    :param v:
    :param w:
    :param btree:
    :return:
    """


def do_decentralized_search(graph, s, t):
    search_success = []
    # 0.1 calculate h(s,t)
    # 1.  find all neighbors of s
    # 2.  for each neighbor u calculate h(u,t)
    # 3.  find min(h(u,t)) and store u
    # 4.  if two or more u's have same value as min(h(u,t))
    # 5.  pick the first u in the neighbors list (i.e. break ties arbitrarily)
    # 6.  if u == t:
    #       search completed successfully
    #       search_success.append(u)
    # 7.  if h(s,t) > h(u,t):
    # 8.    s = u
    # 8.1   do_decentralized_search(graph, s, t)  # now s = u
    # 9.  else if h(s,t) <= h(u,t):
    #       search failed

    return search_success

if __name__ == '__main__':
    simulation()

########################################################################################################################