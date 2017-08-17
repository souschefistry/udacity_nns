############################ Problem 1 ############################
#!/usr/bin/python
# -*- coding: utf-8 -*-

from itertools import combinations
import matplotlib.pyplot as plt
from random import randint, shuffle
import snap
import numpy as np

__author__ = 'dghosh'

US_POWER_GRID_GRAPH_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/USpowergrid_n4941.txt"


def problem_1_1(verbose=False):
    """
    0.000486
    Average clustering coefficient for 100 samples = 0.000486
    """
    # 1. read power grid graph
    power_grid_graph = snap.LoadEdgeList(snap.PUNGraph, US_POWER_GRID_GRAPH_PATH, 0, 1)
    coeffs = {}
    for idx in xrange(100):
        random_graph = generate_graph_via_stub(power_grid_graph, verbose)
        while snap.CntSelfEdges(random_graph) > 0:
            if verbose:
                print "Generated network has self-edge, rerunning stub generation.."
            # if the network is not simple (has self-loops or multi-edges)
            random_graph = generate_graph_via_stub(power_grid_graph)
        # at this point random_graph should be simple (no self-loops or multi-edges)
        # coeffs[random_graph] = calc_coeff([random_graph], idx)
        coeffs[random_graph] = snap.GetClustCf(random_graph, -1)
    print sum(coeffs.values())
    print "Average clustering coefficient for 100 samples = %f" % (sum(coeffs.values()) / 100.0)


def generate_graph_via_stub(power_grid_graph, verbose=False):
    """
    A random network is generated from the stub-matching algorithm as follows:
    1. First, calculate the degree sequence of the real world network by creating a vector −→k =
    {k1, k2, ...kn} where ki is the degree of node i and n is the number of nodes in the graph.

    2. Then create an array v and fill it by writing the index i exactly ki times for each element in
    the degree sequence. Each of these i’s represents a stub of the i-th node in the graph.

    3. Next, randomly permute the elements of v (using, e.g., random.shuffle())

    4. Finally, create a random network by connecting adjacent pairs in v with an edge. In other
    words, connect up the nodes corresponding to the first and second elements, the third and
    fourth, and so on. Formally, v is of length 2m and the m undirected edges of the graph are
    (v1, v2),(v3, v4), . . . ,(v2m−1, v2m).
    """
    graph = snap.TUNGraph.New()  # this will store our stub generated graph
    for node in power_grid_graph.Nodes():
        graph.AddNode(node.GetId())

    # calculate size of vector, vector is 0 indexed while node-ids are 1 indexed in graph
    DegToCntV = snap.TIntPrV()
    snap.GetDegCnt(power_grid_graph, DegToCntV)
    vector_size = 0
    for item in DegToCntV:
        vector_size += item.GetVal2() * item.GetVal1()  # no of nodes with degree d * value of d
    vector = snap.TIntV(vector_size)
    vector_idx = 0

    if verbose:
        print "Vector size: %d" % vector_size
    # compute node-id | degree table
    result_degree = snap.TIntV()
    snap.GetDegSeqV(power_grid_graph, result_degree)
    if verbose:
        for i in range(0, result_degree.Len()):
            print "Node %s has degree %s" % (i, result_degree[i])

    # fill vector
    for node in power_grid_graph.Nodes():
        if verbose:
            print "ID: %d -- deg: %d" % (node.GetId(), result_degree[node.GetId()-1])
        for each_edge in xrange(result_degree[node.GetId()-1]):
            # fill it by writing the index i exactly ki times for each element
            vector[vector_idx] = node.GetId()
            vector_idx += 1
    if verbose:
        for i in range(0, vector.Len()):
            print i, vector[i]

    if verbose:
        print "permute by random.shuffle(x), shuffle the sequence x in place"

    # permute by random.shuffle(x), shuffle the sequence x in place
    shuffle(vector)

    if verbose:
        for i in range(0, vector.Len()-1):
            print i, vector[i]

    # Finally, create a random network by connecting adjacent pairs in v with an edge.
    for i in range(0, vector.Len()-1):
        graph.AddEdge(vector[i], vector[i+1])

    return graph


def normalize_graph(graph):
    """
    One import note is that this algorithm can sometimes create improper (i.e., nonsimple)
    networks with self-loops or multiple edges between two nodes; if this happens, then you
    must reject this sampled network and try again. In other words, you must make sure the sampled
    network is simple (no self-loops or multi-edges) after you construct it and draw a new sample if the
    constructed network is non-simple.
    :param graph:
    :return:
    """
    return snap.DelSelfEdges(graph)


def calc_coeff(networks, sample_id, verbose=False):
    """
    :param networks:
    """
    coeffs = []
    num_nodes = networks[0].GetNodes()
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

    if verbose:
        print "Clustering coefficient network (id=%d) is %f\n" % (sample_id, sum(coeffs[0]) / (1.0*num_nodes))
    return sum(coeffs[0]) / (1.0*num_nodes)


def problem_1_2(verbose=False):
    """
    Simulating the configuration model through rewiring
    """
    # 1. read power grid graph
    power_grid_graph = snap.LoadEdgeList(snap.PUNGraph, US_POWER_GRID_GRAPH_PATH, 0, 1)
    max_iterations = 10000
    step_size = 100
    rewire_count = 0
    # ccoeff_start = calc_coeff([power_grid_graph], 1)
    ccoeff_start = snap.GetClustCf(power_grid_graph, -1)
    if verbose:
        print "Starting clustering coefficient : %f" % ccoeff_start

    coeffs = []

    while True:
        success = rewire_edges(power_grid_graph)
        if success:
            rewire_count += 1
            if rewire_count % step_size == 0:
                # coeff = calc_coeff([power_grid_graph], 1)
                coeff = snap.GetClustCf(power_grid_graph, -1)
                coeffs.append(coeff)

            if rewire_count > max_iterations:
                break

    if verbose:
        print "Average clustering coeff: %f" % (sum(coeffs) / (1.0*max_iterations))

    plot_coeffs(coeffs)


def rewire_edges(graph, verbose=False):
    rewire_success = False
    node_a, node_b = select_random_edge(graph)
    node_c, node_d = select_random_edge(graph)

    e1_1 = (node_a, node_d)
    e2_1 = (node_b, node_c)

    if is_simple_network(e1_1, graph) and is_simple_network(e2_1, graph):
        if verbose:
            print "rewiring success!"
        graph.DelEdge(node_a, node_b)
        graph.DelEdge(node_c, node_d)
        graph.AddEdge(node_a, node_d)
        graph.AddEdge(node_b, node_c)
        rewire_success = True
    return rewire_success


def is_simple_network(edge, graph):
    """
    :param edge: (node_a, node_b)
    :param graph:
    :return: boolean True if no self-loop / multi-edges
    """
    if edge[0] == edge[1]:
        return False
    if graph.IsEdge(edge[0], edge[1]):
        return False
    return True


def select_random_edge(graph, verbose=True):
    node_a = randint(1, graph.GetNodes())
    degree_a = graph.GetNI(node_a).GetDeg()
    rand_edge_id = randint(0, degree_a-1)
    node_b = graph.GetNI(node_a).GetOutNId(rand_edge_id)
    return node_a, node_b


def plot_coeffs(coeffs):
    # plot
    iters = np.arange(1, len(coeffs)+1,1)
    plt.plot(iters, coeffs)
    plt.xlabel("iterations")
    plt.ylabel("average clustering coefficient")
    plt.title("Average clustering coefficient change with rewiring model")
    plt.savefig("rewire_graph.png")
    plt.show()


if __name__ == '__main__':
    print "The Configuration Model [25 points – Will, Zhedi, Jessica, Ben]\n"
    print "1.1 Simulating the configuration model through stub-matching [10 points]\n"
    problem_1_1()
    print "1.2 Simulating the configuration model through rewiring [15 points]\n"
    problem_1_2(verbose=True)

############################ Problem 2 ############################

#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
import snap
from random import randint

EPINIONS_DATA_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/epinions-signed.txt"


def problem_2_1(verbose=True):
    # 1. read graph
    epinions_graph = snap.LoadEdgeList(snap.PUNGraph, EPINIONS_DATA_PATH, 0, 1)
    if verbose:
        print epinions_graph.GetNodes()
        print epinions_graph.GetEdges()

    # get rid of self edges
    snap.DelSelfEdges(epinions_graph)

    # 2. built adjacency list edge -> sign
    # (node_1, node_2) -> sign
    # (node_m, node_p) -> sign
    adjacency_list, pos_edges, neg_edges = create_adjacency_list()
    if verbose:
        print len(adjacency_list.keys())

    print "# of positive edges = %d" % pos_edges
    print "# of negative edges = %d" % neg_edges

    triads = {"3-": 0, "2-": 0, "1-": 0, "3+": 0}  # t0, t1, t2, t3 respectively

    # 3. traverse graph and build triads
    for node in epinions_graph.Nodes():
        neighbors = snap.TIntV()
        for edge_id in node.GetOutEdges():
            neighbors.Add(edge_id)

        for n1_idx, edge_id in enumerate(neighbors):  # n1 -> first neighbor of current node
            n1 = neighbors[n1_idx]

            for n2_idx in xrange(n1_idx+1, len(neighbors)):  # n2 -> second neighbor of current node
                n2 = neighbors[n2_idx]

                if n1 == n2:  # ignore self loops
                    continue
                if not epinions_graph.IsEdge(n1, n2):  # go to next neighbor if n1 is not connected to n2
                    continue

                # at this point we should have a triad
                e1_sign = adjacency_list[(node.GetId(), n1)]
                e2_sign = adjacency_list[(node.GetId(), n2)]
                e3_sign = adjacency_list[(n1, n2)]

                neg_edge_count = 0
                if e1_sign == "-1":
                    neg_edge_count += 1
                if e2_sign == "-1":
                    neg_edge_count += 1
                if e3_sign == "-1":
                    neg_edge_count += 1

                if neg_edge_count == 3:
                    triads["3-"] += 1
                elif neg_edge_count == 2:
                    triads["2-"] += 1
                elif neg_edge_count == 1:
                    triads["1-"] += 1
                elif neg_edge_count == 0:
                    triads["3+"] += 1

    # reset triad count as we triple count
    triads = {ttype: count/3 for ttype, count in triads.iteritems()}
    print triads

    print "a. [5 points] Calculate the count and fraction of triads of each type in the Epinions network."
    print "t0: %d (%f), t1: %d (%f), t2: %d (%f), t3: %d (%f)" \
          % (triads["3-"], triads["3-"] / (1.0*sum(triads.values())),
             triads["2-"], triads["2-"] / (1.0*sum(triads.values())),
             triads["1-"], triads["1-"] / (1.0*sum(triads.values())),
             triads["3+"], triads["3+"] / (1.0*sum(triads.values())))
    print "b. [5 points] Calculate the fraction of positive and negative edges in the graph. " \
          "Calculate the probability of each type of triad."
    print "Fraction of +edges = %f, fraction of -edges = %f" % ((pos_edges / (1.0*(pos_edges+neg_edges))),
                                                                (neg_edges / (1.0*(pos_edges+neg_edges))))

    p = pos_edges / (1.0*(pos_edges+neg_edges))  # fraction of positive edges in graph
    if verbose:
        print p
    print "Probability of 3- edges = %f" % (math.pow((1-p), 3))
    print "Probability of 3+ edges = %f" % (math.pow(p, 3))
    print "Probability of 2+ 1- edges = %f" % (3*math.pow(p, 2)*(1-p))
    print "Probability of 1+ 2- edges = %f" % (3*p*math.pow(1-p, 2))
    print "c. [5 points] Compare the probabilities from part (b) with fractions calculated in part (a)."
    print "Q. Which type of triads do you see more in data as compared to the random baseline values?"
    print "Ans. Type 3 triads with all +ve edges are overwhelmingly more than the rest almost by a factor of 10."
    print "Q. Which type of triads do you see less?"
    print "Ans. Triads with all three -ve edges are the least in this network."
    print "Provide an explanation for this observation."
    print "Ans. In closely knit social communities or the ones with user identity is known e.g. Facebook  down-voting" \
          "flagging negatively is generally avoided because of social cost associated; this can cause users to view" \
          "down voters negatively."


def create_adjacency_list():
    adjacency_list = {}
    pos_edges, neg_edges = 0, 0
    with open(EPINIONS_DATA_PATH) as f:
        for line in f:
            tokens = line.split()
            if not tokens[0].isdigit():
                continue
            v1, v2, sign = tokens
            adjacency_list[(int(v1), int(v2))] = sign
            adjacency_list[(int(v2), int(v1))] = sign
            if sign == "1":
                pos_edges += 1
            elif sign == "-1":
                neg_edges += 1
    return adjacency_list, pos_edges, neg_edges


def problem_2_2():
    pass


def problem_2_3():
    pass


def problem_2_4(verbose=False):
    node_count = 10
    simulations = 1000000
    graph_count = 100
    balanced_count = 0

    for gid in xrange(graph_count):
        adjacency_list = dict()     # (node_k, node_p) -> sign
        graph = snap.TUNGraph.New()

        for nid in xrange(node_count):
            graph.AddNode(nid)

        for n1 in xrange(node_count):
            for n2 in xrange(n1, node_count):
                if n1 == n2:    # ignore self-loops
                    continue

                # flip coin and add edge with sign +1 / -1
                if randint(0, 1) == 0:
                    adjacency_list[(n1, n2)] = 1
                    adjacency_list[(n2, n1)] = 1
                else:
                    adjacency_list[(n1, n2)] = -1
                    adjacency_list[(n2, n1)] = -1

                graph.AddEdge(n1, n2)
        if verbose:
            print adjacency_list

        sid = 0     # simulation id
        balanced = False
        while sid < simulations:
            balanced = is_balanced_graph(graph, adjacency_list)
            if balanced:
                balanced_count += 1
                break

            # pick triad and check if balanced
            triad = select_triad(node_count)    # (node, neighbor_1, neighbor_2)

            sid += 1
            if is_balanced_triad(triad, adjacency_list):
                continue

            # randomly flip edge if unbalanced triad found
            flip_edge(triad, adjacency_list)

    print "2.4) What fraction of the networks end up balanced?"
    print "Out of %d graphs balanced count is = %d, hence percentage of balanced graphs = %f" % \
          (graph_count, balanced_count, graph_count*100/(1.0*balanced_count))


def select_triad(node_count):
    random_node_id = randint(0, node_count-1)  # pick a random node
    neighbor_1_id, neighbor_2_id = random_node_id, random_node_id  # initialize
    while neighbor_1_id == random_node_id:
        random_node_id = randint(0, node_count-1)

    # keep selecting till n2 is different than n1 or random node
    while neighbor_2_id == neighbor_1_id or neighbor_2_id == random_node_id:
        neighbor_2_id = randint(0, node_count-1)

    # triad = (nid, n1, n2), nid -> start node id, n1 -> 1st neighbor's id, n2 -> 2nd neighbor's id
    return random_node_id, neighbor_1_id, neighbor_2_id


def is_balanced_triad(triad, adjacency_list):
    # triad = (nid, n1, n2), nid -> start node id, n1 -> 1st neighbor's id, n2 -> 2nd neighbor's id
    # Based on the number of positive edges we label triads with odd number of +1s as balanced,
    # and triads with even positive edges as unbalanced.
    # sum(signs) = +3 (+1+1+1) or -1 (+1-2) == balanced
    # sum(signs) = +1 (+1+1-1) or -3 (-1-1-1) == unbalanced
    signs = list()
    signs.append(adjacency_list[(triad[0], triad[1])])
    signs.append(adjacency_list[(triad[0], triad[2])])
    signs.append(adjacency_list[(triad[1], triad[2])])

    if sum(signs) == -3 or sum(signs) == 1:
        return False
    return True


def is_balanced_graph(graph, adjacency_list):
    for node in graph.Nodes():
        nid = node.GetId()
        for node1 in graph.Nodes():
            n1 = node1.GetId()
            if n1 == nid:
                continue
            for node2 in graph.Nodes():
                n2 = node2.GetId()
                if n2 == nid or n2 == n1:
                    continue
                if not is_balanced_triad((nid, n1, n2), adjacency_list):
                    return False
    return True


def flip_edge(triad, adjacency_list):
    # triad = (nid, n1, n2), nid -> start node id, n1 -> 1st neighbor's id, n2 -> 2nd neighbor's id
    selector = randint(0, 2)
    if selector == 0:
        adjacency_list[(triad[0], triad[1])] *= -1
        adjacency_list[(triad[1], triad[0])] *= -1
    elif selector == 1:
        adjacency_list[(triad[0], triad[2])] *= -1
        adjacency_list[(triad[2], triad[0])] *= -1
    elif selector == 2:
        adjacency_list[(triad[1], triad[2])] *= -1
        adjacency_list[(triad[2], triad[1])] *= -1

if __name__ == '__main__':
    print "2.1 Signed Triads in Epinions [15 points]\n"
    problem_2_1()
    print "2.2 Balance in a Random Signed Network [8 points]\n"
    problem_2_2()
    print "2.3 A Dynamic Process Model of Balance [5 points]\n"
    problem_2_3()
    print "2.4 Simulation [8 points]\n"
    problem_2_4()

############################ Problem 3 ############################

#!/usr/bin/python
# -*- coding: utf-8 -*-

import operator
import snap
import matplotlib.pyplot as plt

GRAPH_1_DATA_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/graph1.txt"
GRAPH_2_DATA_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/graph2.txt"
CURRENT_ALTERNATING_VOTE = "A"


def init():
    # 1. read graphs
    graph_1 = snap.LoadEdgeList(snap.PUNGraph, GRAPH_1_DATA_PATH, 0, 1)
    graph_2 = snap.LoadEdgeList(snap.PUNGraph, GRAPH_2_DATA_PATH, 0, 1)
    return graph_1, graph_2


def problem_3_1(graph_1, graph_2, verbose=True):
    global CURRENT_ALTERNATING_VOTE
    simulations = 10
    if verbose:
        print graph_1.GetNodes()
        print graph_2.GetNodes()

    voter_support_1 = ["-"] * graph_1.GetNodes()  # table indexed by position for "A" / "B" / "-" (undecided) in graph_1
    voter_support_2 = ["-"] * graph_2.GetNodes()  # table indexed by position for "A" / "B" / "-" (undecided) in graph_2

    # assign candidates to voters
    for node in graph_1.Nodes():
        if node.GetId() % 10 in [0, 1, 2, 3]:
            voter_support_1[node.GetId()] = "A"
        elif node.GetId() % 10 in [4, 5, 6, 7]:
            voter_support_1[node.GetId()] = "B"

    if verbose:
        print "A's supporter in graph_1= %d (percentage = %f)" % \
              (voter_support_1.count("A"), voter_support_1.count("A") / (1.0*len(voter_support_1)))
        print "B's supporter in graph_1= %d (percentage = %f)" % \
              (voter_support_1.count("B"), voter_support_1.count("B") / (1.0*len(voter_support_1)))
        print "Undecided voters in graph_1= %d (percentage = %f)" % \
              (voter_support_1.count("-"), voter_support_1.count("-") / (1.0*len(voter_support_1)))

    run_simulations(simulations, graph_1, voter_support_1, 1, False)

    # for sim in range(simulations):
    #     voter_support_1 = run_simulations_v2(graph_1, voter_support_1, 1, False)

    for node in graph_2.Nodes():
        if node.GetId() % 10 in [0, 1, 2, 3]:
            voter_support_2[node.GetId()] = "A"
        elif node.GetId() % 10 in [4, 5, 6, 7]:
            voter_support_2[node.GetId()] = "B"

    if verbose:
        print "A's supporter in graph_2= %d (percentage = %f)" % \
              (voter_support_2.count("A"), voter_support_2.count("A") / (1.0*len(voter_support_2)))
        print "B's supporter in graph_2= %d (percentage = %f)" % \
              (voter_support_2.count("B"), voter_support_2.count("B") / (1.0*len(voter_support_2)))
        print "Undecided voters in graph_2= %d (percentage = %f)" \
              % (voter_support_2.count("-"), voter_support_2.count("-") / (1.0*len(voter_support_2)))

    # for sim in range(simulations):
    #     voter_support_1 = run_simulations_v2(graph_1, voter_support_1, 1, False)

    run_simulations(simulations, graph_2, voter_support_2, 2, False)


def run_simulations_v2(graph, voter_support, gid, verbose):
    global CURRENT_ALTERNATING_VOTE
    for node in graph.Nodes():
        if voter_support[node.GetId()] == "":    # this is an undecided voter
            a_supporters, b_supporters = 0, 0
            for edge_id in node.GetOutEdges():
                # find out your friend's affiliation and store
                if voter_support[edge_id] == "A":
                    a_supporters += 1
                if voter_support[edge_id] == "B":
                    b_supporters += 1

            # print a_supporters, b_supporters
            if a_supporters > b_supporters:
                if verbose:
                    print "current alternate = %s" % CURRENT_ALTERNATING_VOTE
                    print "I am going with A"
                voter_support[node.GetId()] = "A"

            if b_supporters > a_supporters:
                if verbose:
                    print "current alternate = %s" % CURRENT_ALTERNATING_VOTE
                    print "I am going with B"
                voter_support[node.GetId()] = "B"

            if a_supporters == b_supporters:
                if CURRENT_ALTERNATING_VOTE == "A":
                    CURRENT_ALTERNATING_VOTE = "B"
                else:
                    CURRENT_ALTERNATING_VOTE = "A"
                voter_support[node.GetId()] = CURRENT_ALTERNATING_VOTE

    if verbose:
        print "New count of A's supporters in graph %d = %d" % (gid, voter_support.count("A") - 4000)
        print "New count of B's supporters in graph %d = %d" % (gid, voter_support.count("B") - 4000)

    if voter_support.count("A") > voter_support.count("B"):
        print "A wins graph %d by %d votes" % (gid, voter_support.count("A")-voter_support.count("B"))
    elif voter_support.count("A") < voter_support.count("B"):
        print "B wins graph %d by %d votes" % (gid, voter_support.count("B")-voter_support.count("A"))
    return voter_support


def run_simulations(simulations, graph, voter_support, gid, verbose):
    global CURRENT_ALTERNATING_VOTE
    for sid in range(simulations):
        for node in graph.Nodes():
            if voter_support[node.GetId()] == "-":    # this is an undecided voter
                a_supporters, b_supporters = 0, 0
                for edge_id in node.GetOutEdges():
                    # find out your friend's affiliation and store
                    if voter_support[edge_id] == "A":
                        a_supporters += 1
                    if voter_support[edge_id] == "B":
                        b_supporters += 1

                if a_supporters > b_supporters:
                    if verbose:
                        print "current alternate = %s" % CURRENT_ALTERNATING_VOTE
                        print "I am going with A"
                    voter_support[node.GetId()] = "A"

                if b_supporters > a_supporters:
                    if verbose:
                        print "current alternate = %s" % CURRENT_ALTERNATING_VOTE
                        print "I am going with B"
                    voter_support[node.GetId()] = "B"

                if a_supporters == b_supporters:
                    if CURRENT_ALTERNATING_VOTE == "A":
                        CURRENT_ALTERNATING_VOTE = "B"
                    else:
                        CURRENT_ALTERNATING_VOTE = "A"
                    voter_support[node.GetId()] = CURRENT_ALTERNATING_VOTE

    if verbose:
        print "New count of A's supporters in graph %d = %d" % (gid, voter_support.count("A") - 4000)
        print "New count of B's supporters in graph %d = %d" % (gid, voter_support.count("B") - 4000)

    print "A=%d, B=%d, -=%d" % (voter_support.count("A"), voter_support.count("B"), voter_support.count("-"))


def problem_3_2(graph_1, graph_2, verbose):
    voter_support_1, voter_support_2 = {}, {}
    vote_count_g1, vote_count_g2 = [], []
    simulations = 100
    min_spent_1, min_spent_2 = [], []

    for budget in range(0, 10000, 1000):
        if verbose:
            print "processing graph 1 with budget = %d" % budget
        voter_support_1 = ["-"] * graph_1.GetNodes()
        voter_support_2 = ["-"] * graph_2.GetNodes()

        voter_support_1 = assign_candidates(graph_1, voter_support_1)
        turned_voters, voter_support_1 = run_ad_campaign(voter_support_1, budget)
        voter_support_1 = run_simulations_type_2(simulations, graph_1, voter_support_1, turned_voters, 1, False)
        vote_count_g1.append(voter_support_1.count("A") - voter_support_1.count("B"))
        if voter_support_1.count("A") > voter_support_1.count("B"):
            min_spent_1.append(budget)

        if verbose:
            print "processing graph 2 with budget = %d" % budget
        voter_support_2 = assign_candidates(graph_1, voter_support_2)
        turned_voters, voter_support_2 = run_ad_campaign(voter_support_2, budget)
        voter_support_2 = run_simulations_type_2(simulations, graph_2, voter_support_2, turned_voters, 2, False)
        vote_count_g2.append(voter_support_2.count("A") - voter_support_2.count("B"))
        if voter_support_2.count("A") > voter_support_2.count("B"):
            min_spent_2.append(budget)

    print "Minimum spending needed = %d, %d" % (min(min_spent_1), min(min_spent_2))

    budgets = range(0, 10000, 1000)
    plt.title("Campaign budget vs vote margin")
    labels = ["graph-1", "graph-2"]
    colors = ['red', 'blue']
    plt.xlabel('Amount spent for ad campaign')
    plt.ylabel('Wining margin for candidate A')
    plt.grid(True)
    plt.plot(budgets, vote_count_g1, linestyle="dashed", marker="o", color=colors[0], label=labels[0])
    plt.plot(budgets, vote_count_g2, linestyle="dashed", marker="o", color=colors[1], label=labels[1])
    plt.legend(loc='upper left', shadow=True)
    plt.savefig("3.2.png")
    plt.show()


def assign_candidates(graph, voter_support, verbose=False):
    # assign candidates to voters
    for node in graph.Nodes():
        if node.GetId() % 10 in [0, 1, 2, 3]:
            voter_support[node.GetId()] = "A"
        elif node.GetId() % 10 in [4, 5, 6, 7]:
            voter_support[node.GetId()] = "B"

    if verbose:
        print "A's supporter in graph= %d (percentage = %f)" % \
              (voter_support.count("A"), voter_support.count("A") / (1.0*len(voter_support)))
        print "B's supporter in graph= %d (percentage = %f)" % \
              (voter_support.count("B"), voter_support.count("B") / (1.0*len(voter_support)))
        print "Undecided voters in graph= %d (percentage = %f)" % \
              (voter_support.count("-"), voter_support.count("-") / (1.0*len(voter_support)))
    return voter_support


def run_ad_campaign(voter_support, budget, verbose=False):
    step = 0
    turned_voters = []
    old = budget
    while budget > 0:
        for idx in xrange(3000 + step*10, 3010 + step*10):
            turned_voters.append(idx)
            voter_support[idx] = "A"
        step += 1
        budget -= 1000
    if verbose:
        print "budget %d, A count %d, turned voters list %s" % (old, voter_support.count("A"), turned_voters)
    return turned_voters, voter_support


def run_simulations_type_2(simulations, graph, voter_support, turned_voters, gid, verbose):
    global CURRENT_ALTERNATING_VOTE
    loyalties = {}
    for sid in range(simulations):
        for node in graph.Nodes():
            if voter_support[node.GetId()] == "-" and node.GetId() not in turned_voters:    # this is an undecided voter
                a_supporters, b_supporters = 0, 0
                for edge_id in node.GetOutEdges():
                    # find out your friend's affiliation and store
                    if voter_support[edge_id] == "A":
                        a_supporters += 1
                    if voter_support[edge_id] == "B":
                        b_supporters += 1

                    voter_support[node.GetId()] = count_votes_and_break_ties(a_supporters, b_supporters)

    if verbose:
        print "New count of A's supporters in graph %d = %d" % (gid, voter_support.count("A") - 4000)
        print "New count of B's supporters in graph %d = %d" % (gid, voter_support.count("B") - 4000)
        print "A=%d, B=%d, -=%d" % (voter_support.count("A"), voter_support.count("B"), voter_support.count("-"))
        print loyalties

    return voter_support


def problem_3_3(graph_1, graph_2, verbose=False):
    voter_support_1, voter_support_2 = {}, {}
    vote_count_g1, vote_count_g2 = [], []
    simulations = 10

    whales_g1 = tag_whales(graph_1)
    whales_g2 = tag_whales(graph_2)
    min_spent_1, min_spent_2 = [], []

    for budget in range(0,10000,1000):
        voter_support_1 = ["-"] * graph_1.GetNodes()
        voter_support_2 = ["-"] * graph_2.GetNodes()

        voter_support_1 = assign_candidates(graph_1, voter_support_1)
        # how many top rollers did we persuade?
        invited_whales, voter_support_1 = persuade_whales(whales_g1, voter_support_1, budget)
        voter_support_1 = run_simulations_type_3(simulations, graph_1, voter_support_1, invited_whales, verbose=False)
        vote_count_g1.append(voter_support_1.count("A") - voter_support_1.count("B"))
        if voter_support_1.count("A") > voter_support_1.count("B"):
            min_spent_1.append(budget)

        voter_support_2 = assign_candidates(graph_2, voter_support_2)
        # how many top rollers did we persuade?
        invited_whales, voter_support_2 = persuade_whales(whales_g2, voter_support_2, budget)
        voter_support_2 = run_simulations_type_3(simulations, graph_2, voter_support_2, invited_whales, verbose=False)
        vote_count_g2.append(voter_support_2.count("A") - voter_support_2.count("B"))
        if voter_support_2.count("A") > voter_support_2.count("B"):
            min_spent_2.append(budget)

    print "Minimum spending to dine high rollers: %d, %d" % (min(min_spent_1), min(min_spent_2))
    budgets = range(0, 10000, 1000)
    plt.title("Campaign budget vs vote margin with high rollers")
    labels = ["graph-1", "graph-2"]
    colors = ['red', 'blue']
    plt.xlabel('Amount spent for high rollers')
    plt.ylabel('Wining margin for candidate A')
    plt.grid(True)
    plt.plot(budgets, vote_count_g1, linestyle="dashed", marker="o", color=colors[0], label=labels[0])
    plt.plot(budgets, vote_count_g2, linestyle="dashed", marker="o", color=colors[1], label=labels[1])
    plt.legend(loc='upper left', shadow=True)
    plt.savefig("3.3.png")
    plt.show()


def tag_whales(graph):
    rollers = {}    # node-id vs degree count sorted
    for node in graph.Nodes():
        rollers[node.GetId()] = node.GetOutDeg()
    sorted_x = sorted(rollers.items(), key=operator.itemgetter(1), reverse=True)
    return sorted_x


def persuade_whales(whales, voter_support, budget):
    # everyone that comes to your dinner is instantly persuaded to vote for candidate A
    idx = 0
    invited_whales = []
    while budget > 0:   # while you have money
        whale_id = whales[idx][0]  # select roller with top-k friends
        voter_support[whale_id] = "A"   # persuade roller to vote for candidate A
        invited_whales.append(whale_id)
        idx += 1
        budget -= 1000
    return invited_whales, voter_support


def run_simulations_type_3(simulations, graph, voter_support, invited_whales, verbose=False):
    global CURRENT_ALTERNATING_VOTE
    for sid in range(simulations):
        for node in graph.Nodes():
            if node.GetId() not in invited_whales and voter_support[node.GetId()] == "-":    # undecided voter who is not a whale
                a_supporters, b_supporters = 0, 0
                for edge_id in node.GetOutEdges():
                    # find out your friend's affiliation and store
                    if voter_support[edge_id] == "A":
                        a_supporters += 1
                    if voter_support[edge_id] == "B":
                        b_supporters += 1

                voter_support[node.GetId()] = count_votes_and_break_ties(a_supporters, b_supporters)
    if verbose:
        print "New count of A's supporters in graph = %d" % (voter_support.count("A") - 4000)
        print "New count of B's supporters in graph = %d" % (voter_support.count("B") - 4000)
        print "A=%d, B=%d, -=%d" % (voter_support.count("A"), voter_support.count("B"), voter_support.count("-"))

    return voter_support


def count_votes_and_break_ties(a_supporters, b_supporters, verbose=False):
    global CURRENT_ALTERNATING_VOTE
    loyalty = ""
    if a_supporters > b_supporters:
        if verbose:
            print "current alternate = %s" % CURRENT_ALTERNATING_VOTE
            print "I am going with A"
        loyalty = "A"

    if b_supporters > a_supporters:
        if verbose:
            print "current alternate = %s" % CURRENT_ALTERNATING_VOTE
            print "I am going with B"
        loyalty = "B"

    if b_supporters == a_supporters:
        if CURRENT_ALTERNATING_VOTE == "A":
            CURRENT_ALTERNATING_VOTE = "B"
        else:
            CURRENT_ALTERNATING_VOTE = "A"
        loyalty = CURRENT_ALTERNATING_VOTE
    return loyalty


def problem_3_4(graph_1, graph_2):
    DegToCntV_g1 = snap.TIntPrV()
    snap.GetDegCnt(graph_1, DegToCntV_g1)
    degrees_g1 = []
    nodes_g1 = []

    for item in DegToCntV_g1:
        nodes_g1.append(item.GetVal2() / float(10000))
        degrees_g1.append(item.GetVal1())

    DegToCntV_g2 = snap.TIntPrV()
    snap.GetDegCnt(graph_2, DegToCntV_g2)
    degrees_g2 = []
    nodes_g2 = []

    for item in DegToCntV_g2:
        nodes_g2.append(item.GetVal2() / float(10000))
        degrees_g2.append(item.GetVal1())

    plt.title("Degree distribution plot")
    labels = ["graph-1", "graph-2"]
    colors = ['red', 'blue']
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('degree')
    plt.ylabel('number of nodes')
    plt.grid(True)
    plt.plot(degrees_g1, nodes_g1, color=colors[0], label=labels[0])
    plt.plot(degrees_g2, nodes_g2, color=colors[1], label=labels[1])
    plt.legend(loc='upper left', shadow=True)
    plt.savefig("3.4.png")
    plt.show()


if __name__ == '__main__':
    g1, g2 = init()
    print "3.1 Basic setup and forecasting [8 points]\n"
    problem_3_1(g1, g2)
    print "3.2 TV Advertising [8 points]\n"
    problem_3_2(g1, g2, verbose=False)
    print "3.3 Wining and Dining the High Rollers [8 points]\n"
    problem_3_3(g1, g2)
    print "3.4 Analysis [6 points]\n"
    problem_3_4(g1, g2)
