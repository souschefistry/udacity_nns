#!/usr/bin/python
# -*- coding: utf-8 -*-

import snap
import matplotlib.pyplot as plt

__author__ = 'dghosh'

HISTOGRAM_PATH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/thresholds.txt"


def problem_1_c():
    histogram_data, y_axes = [], []
    with open(HISTOGRAM_PATH, "rb") as f:
        for line in f:
            tokens = line.split("\t")[0]
            histogram_data.append(int(tokens))
    y_axes = [sum(histogram_data[0:k+1]) for k in range(len(histogram_data)-1)]
    _plot_histogram(y_axes)

    final_rioters = 0

    # find out final # of rioters
    for idx in range(0, len(y_axes)-1):
        if idx == 0 or y_axes[idx-1] >= idx:
            final_rioters += histogram_data[idx]
        else:
            break

    print "Final # of rioters = %d" % final_rioters


def _plot_histogram(y_axes):
    plt.hist(y_axes, len(y_axes), cumulative=True, color="red")
    plt.xlabel("value")
    plt.ylabel("cumulative sum")
    plt.grid(True)
    plt.title("Cumulative histogram of rioters")
    plt.savefig("cumulative-histogram.png")
    plt.show()

if __name__ == "__main__":
    print "Mob Psychology and Thresholds [10 points – Will]"
    print "(c) [4 points] The cumulative histogram"
    problem_1_c()

############### Problem 2 ###############
# coding=utf-8

import math
from collections import Counter
import numpy as np
from random import uniform
from matplotlib import pyplot as plt

__author__ = 'dghosh'

SAMPLES = 100000
ALPHA = 2
XMIN = 1


def problem_2_2():
    """
        Plot the empirical distribution on a log-log scale
        by rst rounding each sample to the nearest integer and then plotting the empirical distribution
        over these rounded values. Also plot the true probability density function for the power law (this
        will help you verify that you generated the data correctly).
    :return:
    """
    empirical_dataset = create_dataset()
    hist_x, hist_y = create_empirical_hist(empirical_dataset)
    plt.title("Empirical distribution and true PDF graph on log-log scale")
    plt.loglog(hist_x, hist_y, "r--", basex=10, basey=10, label='empirical')
    plt.grid(True)
    lx = sorted(hist_x)
    ly = [1/(x**2) for x in lx]
    plt.loglog(lx, ly, "g.", basex=10, basey=10, label='true pdf')
    plt.xlabel("data points")
    plt.ylabel("probability distribution")
    plt.savefig("hw3-2_2-empirical.png")
    plt.show()


def create_dataset():
    dataset = []
    for i in xrange(SAMPLES):
        data_point = round(1/(1 - uniform(0, 1)))
        dataset.append(data_point)
    return dataset


# generating empirical histogram
def create_empirical_hist(dataset):
    counter = Counter(dataset)
    x_axes, y_axes = zip(*(counter.items()))
    y_axes = [idx/(1.0 * SAMPLES) for idx in y_axes]
    return np.array(x_axes), np.array(y_axes)


# 2.3 Least squares estimation with the PDF [5 points]
def problem_2_3():
    """
    2.3 Least squares estimation with the PDF [5 points]
    :return:
    """
    empirical_dataset = create_dataset()
    x_axes_1, y_axes_1 = create_empirical_hist(empirical_dataset)
    arr = np.array([np.log(x_axes_1), np.ones(len(x_axes_1))]).T
    lse = np.linalg.lstsq(arr, np.log(y_axes_1))
    print "least-squares estimate of alpha=%f" % -lse[0][0]

    # plot as in 2.2
    plt.title("(2.3) Least squares estimation with the PDF and alpha, alpha'")
    plt.loglog(x_axes_1, y_axes_1, "r--", basex=10, basey=10, label='empirical')
    plt.grid(True)
    lx = sorted(x_axes_1)
    ly = [1/(x**2) for x in lx]
    plt.loglog(lx, ly, "g.", basex=10, basey=10, label='true pdf')
    # calc pdf
    pdf = [0] * len(empirical_dataset)
    for idx in xrange(1, len(empirical_dataset)):
        val = (ALPHA - 1) / (1.0 * (float(XMIN) * math.pow(float(idx) / XMIN, -ALPHA)))
        pdf[idx - 1] = val
    # plot ls,pdf and ls',pdf
    plt.loglog(range(1, len(empirical_dataset)), pdf, "y--", basex=10, basey=10, label='alpha-ls-pdf')

    # In 1{2 sentences, briely explain how you can improve the accuracy of the estimate by ignoring
    # some of your data. Use this procedure to compute a new estimate 0
    # ls,pdf.

    data = filter(lambda val: val > float(1) / 1000, empirical_dataset)
    x_axes_2, y_axes_2 = create_empirical_hist(data)
    arr_prime = np.array([np.log(x_axes_2), np.ones(len(x_axes_2))]).T
    lse_prime = np.linalg.lstsq(arr_prime, np.log(y_axes_2))
    print "least-squares estimate of alpha-prime=%f" % -lse_prime[0][0]
    plt.loglog(range(1, len(empirical_dataset)), lse_prime[0][:], "b.", basex=10, basey=10, label='alpha-prime-ls-pdf')
    plt.xlabel("data points")
    plt.ylabel("probability distribution")
    plt.savefig("hw3-2.3.png")
    plt.show()


# 2.4 Least squares estimate from the CCDF [5 points]
def problem_2_4():
    """
    2.4 Least squares estimate from the CCDF [5 points]
    :return:
    """
    empirical_dataset = create_dataset()
    x_axes_1, y_axes_1 = create_empirical_hist(empirical_dataset)
    x_axes_2, x_axes_2 = create_ccdf_hist(empirical_dataset)
    arr = np.array([np.log(x_axes_2), np.ones(len(x_axes_2))]).T
    lse_ccdf = np.linalg.lstsq(arr, np.log(x_axes_2))
    alpha_ccdf = (1-lse_ccdf[0][0])
    print "least-squares estimate of alpha-ccdf=%f" % alpha_ccdf
    plt.title("(2.4) Least squares estimate from the CCDF with alpha-ccdf")
    plt.loglog(x_axes_1, y_axes_1, "r--", basex=10, basey=10, label='empirical')
    plt.grid(True)
    lx = sorted(x_axes_1)
    ly = [1/(x**2) for x in lx]
    plt.loglog(lx, ly, "g.", basex=10, basey=10, label='true pdf')
    # calc pdf for ccdf
    pdf_ccdf = [0] * len(empirical_dataset)
    for idx in xrange(1, len(empirical_dataset)):
        val = (alpha_ccdf - 1) / (1.0 * (float(XMIN) * math.pow(float(idx) / XMIN, -alpha_ccdf)))
        pdf_ccdf[idx - 1] = val
    plt.loglog(range(1, len(empirical_dataset)), pdf_ccdf, "b--", basex=10, basey=10, label="alpha-ccdf")
    plt.xlabel("data points")
    plt.ylabel("probability distribution")
    plt.savefig("hw3-2.4.png")
    plt.show()


# Calculate the empirical CCDF from the data.
def create_ccdf_hist(dataset):
    dataset_counter = Counter(dataset)
    ccdf_sum = 0.0
    ccdf = {}
    for x in sorted(dataset_counter.keys(), reverse=True):
        ccdf_sum += dataset_counter[x]
        ccdf[x] = ccdf_sum/(1.0*SAMPLES)
    x_axes, y_axes = zip(*(ccdf.items()))
    return np.array(x_axes), np.array(y_axes)


# 2.5 Maximum likelihood estimation [5 points]
def problem_2_5():
    empirical_dataset = create_dataset()
    x_axes_1, y_axes_1 = create_empirical_hist(empirical_dataset)
    mle_sum = sum([np.log(data) for data in empirical_dataset])
    alpha_mle = (1+SAMPLES/(1.0*mle_sum))
    print "Derive the maximum likelihood estimate of alpha =%f" % alpha_mle
    plt.title("(2.5) Maximum likelihood estimation with alpha-mle")
    plt.loglog(x_axes_1, y_axes_1, "r--", basex=10, basey=10, label='empirical')
    plt.grid(True)
    lx = sorted(x_axes_1)
    ly = [1/(x**2) for x in lx]
    plt.loglog(lx, ly, "g.", basex=10, basey=10, label='true pdf')
    # calc pdf
    pdf_mle = [0] * len(empirical_dataset)
    for idx in xrange(1, len(empirical_dataset)):
        val = (alpha_mle - 1) / (1.0 * (float(XMIN) * math.pow(float(idx) / XMIN, -alpha_mle)))
        pdf_mle[idx - 1] = val
    plt.loglog(range(1, len(empirical_dataset)), pdf_mle, "b--",  basex=10, basey=10, label='alpha-mle')
    plt.xlabel("data points")
    plt.ylabel("probability distribution")
    plt.savefig("hw3-2.5.png")
    plt.show()


# 2.6 Comparison [5 points]
def problem_2_6():
    alpha_ls_pdfs, alpha_prime_ls_pdfs, alpha_ls_ccdfs, alpha_mles = [], [], [], []
    for run in xrange(100):
        empirical_dataset = create_dataset()
        # 1. compute alpha-ls-pdf
        x_axes_1, y_axes_1 = create_empirical_hist(empirical_dataset)
        arr = np.array([np.log(x_axes_1), np.ones(len(x_axes_1))]).T
        lse = np.linalg.lstsq(arr, np.log(y_axes_1))
        alpha_ls_pdfs.append(-lse[0][0])
        # print "least-squares estimate of alpha-ls-pdf=%f" % -lse[0][0]
        # 2. compute alpha-prime-ls-pdf
        # Let's ignore values < 100
        data = filter(lambda val: val < 100, empirical_dataset)
        x_axes_2, y_axes_2 = create_empirical_hist(data)
        arr_prime = np.array([np.log(x_axes_2), np.ones(len(x_axes_2))]).T
        lse_prime = np.linalg.lstsq(arr_prime, np.log(y_axes_2))
        alpha_prime_ls_pdfs.append(-lse_prime[0][0])
        # print "least-squares estimate of alpha-prime=%f" % -lse_prime[0][0]
        # 3. compute alpha-ls-ccdf as in problem 2.4
        x_axes_3, y_axes_3 = create_ccdf_hist(empirical_dataset)
        arr = np.array([np.log(x_axes_3), np.ones(len(x_axes_3))]).T
        lse_ccdf = np.linalg.lstsq(arr, np.log(y_axes_3))
        alpha_ls_ccdfs.append((1-lse_ccdf[0][0]))
        # print "least-squares estimate of alpha-ccdf=%f" % (1-lse_ccdf[0][0])
        # 4. compute alpha-mle as in 2.5
        mle_sum = sum([np.log(data) for data in empirical_dataset])
        alpha_mle = (1+SAMPLES/(1.0 * mle_sum))
        alpha_mles.append(alpha_mle)

    print "alpha-ls-pdf mean = %f, standard-deviation = %f" % (np.mean(alpha_ls_pdfs), np.std(alpha_ls_pdfs))
    print "alpha-ls-prime-pdf mean = %f, standard-deviation = %f" % (np.mean(alpha_prime_ls_pdfs), np.std(alpha_prime_ls_pdfs))
    print "alpha-ls-ccdf mean = %f, standard-deviation = %f" % (np.mean(alpha_ls_ccdfs), np.std(alpha_ls_ccdfs))
    print "alpha-mle mean = %f, standard-deviation = %f" % (np.mean(alpha_mles), np.std(alpha_mles))


if __name__ == "__main__":
    print "Empirical Power Laws [30 points – Austin]"
    print "2.2 Sampling [5 points]\n"
    problem_2_2()
    print "2.3 Least squares estimation with the PDF [5 points]\n"
    problem_2_3()
    print "2.4 Least squares estimate from the CCDF [5 points]"
    problem_2_4()
    print "2.5 Maximum likelihood estimation [5 points]"
    problem_2_5()
    print "2.6 Comparison [5 points]"
    problem_2_6()

############### Problem 3 ###############
#!/usr/bin/python
# -*- coding: utf-8 -*-

import operator
import snap
import numpy as np
from scipy.stats import chi2_contingency, mannwhitneyu
from random import randint, random

__author__ = 'dghosh'

ACTORS_GRAPH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/imdb_actor_edges.tsv"
SIR_ERDOS_GRAPH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/SIR_erdos_renyi.txt"
SIR_PREF_GRAPH = "/Users/dghosh/Documents/education/autumn2016/cs224w/resources/SIR_preferential_attachment.txt"
SIMULATIONS = 100
BETA = 0.05
DELTA = 0.5
GRAPHS_STR_LIST = ["Erdos", "PA", "Actors"]


def init(verbose=False):
    actors_graph = snap.LoadEdgeList(snap.PUNGraph, ACTORS_GRAPH, 0, 1)
    sir_erdos_graph = snap.LoadEdgeList(snap.PUNGraph, SIR_ERDOS_GRAPH, 0, 1)
    sir_pref_graph = snap.LoadEdgeList(snap.PUNGraph, SIR_PREF_GRAPH, 0, 1)
    if verbose:
        print actors_graph.GetNodes(), sir_erdos_graph.GetNodes(), sir_pref_graph.GetNodes()
        print actors_graph.GetEdges(), sir_erdos_graph.GetEdges(), sir_pref_graph.GetEdges()
    return sir_erdos_graph, sir_pref_graph, actors_graph


# 3.1 [10 points]

def problem_3_1_a(graphs, verbose=False):
    # epidemics = []
    # for gid, graph in enumerate(graphs):
    #     S, I, R, epidemic_count = run_epidemic_simulation_v4(graphs[gid], gid)
    #     if epidemic_count >= 50:
    #         print "Epidemic in %s with perc=%f" % (GRAPHS_STR_LIST[gid], epidemic_count)
    #         epidemics.append(epidemic_count)

    # [*] Proportion of simulations that infected at least 50 perc of the erdos network = 74
    # [*] Proportion of simulations that infected at least 50 perc of the PA network = 66
    # [*] Proportion of simulations that infected at least 50 perc of the actors network = 46
    epidemics = [74, 66, 46]

    chi2, p_val, dof, expected = chi2_contingency([[epidemics[0], 100-epidemics[0]], [epidemics[1], 100-epidemics[1]]])
    print "Erdos vs PA: chi2 = %f, P-value = %f" % (chi2, p_val)
    chi2, p_val, dof, expected = chi2_contingency([[epidemics[1], 100-epidemics[1]], [epidemics[2], 100-epidemics[2]]])
    print "PA vs Actors: chi2 = %f, P-value = %f" % (chi2, p_val)
    chi2, p_val, dof, expected = chi2_contingency([[epidemics[2], 100-epidemics[2]], [epidemics[0], 100-epidemics[0]]])
    print "Actors vs Erdos: chi2 = %f, P-value = %f" % (chi2, p_val)

    # Erdos vs PA: chi2 = 1.166667, P-value = 0.280087
    # PA vs Actors: chi2 = 7.325487, P-value = 0.006798
    # Actors vs Erdos: chi2 = 15.187500, P-value = 0.000097


def problem_3_1_b(graphs, verbose=False):
    r_counts = []
    for gid, graph in enumerate(graphs):
        S, I, R, r_count_arr = run_epidemic_simulation_3b(graphs[gid], gid)
        r_counts.append(r_count_arr)
        print "[*] %s mean infection = %f" % (GRAPHS_STR_LIST[gid], np.mean(r_count_arr))

    # [*] Proportion of simulations that infected at least 50 perc of the erdos network = 74
    # [*] Proportion of simulations that infected at least 50 perc of the PA network = 66
    # [*] Proportion of simulations that infected at least 50 perc of the actors network = 46
    # epidemics = [74, 66, 46]

    # find len(R) / node_count for each simulation when epidemic happened
    percs_erdos = [(r_count / (1.0 * graphs[0].GetNodes())) for r_count in r_counts[0]]
    percs_pa = [(r_count / (1.0 * graphs[1].GetNodes())) for r_count in r_counts[1]]
    percs_actors = [(r_count / (1.0 * graphs[2].GetNodes())) for r_count in r_counts[2]]
    print percs_erdos, percs_pa, percs_actors

    statistic, p_val = mannwhitneyu(percs_pa, percs_erdos)
    print "PA vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_pa, percs_actors)
    print "PA vs Actors: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_actors, percs_erdos)
    print "Actors vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)


def problem_3_1_c(graphs, verbose=False):
    r_counts = []
    for gid, graph in enumerate(graphs):
        S, I, R, r_count_arr = run_epidemic_simulation_3c_v2(graph, gid)
        r_counts.append(r_count_arr)
        print "[*]%s mean infection = %f" % (GRAPHS_STR_LIST[gid], np.mean(r_count_arr))

    # [*] Proportion of simulations that infected at least 50 perc of the erdos network = 74
    # [*] Proportion of simulations that infected at least 50 perc of the PA network = 66
    # [*] Proportion of simulations that infected at least 50 perc of the actors network = 46
    # epidemics = [74, 66, 46]

    # find len(R) / node_count for each simulation when epidemic happened
    percs_erdos = [(r_count / (1.0 * graphs[0].GetNodes())) for r_count in r_counts[0]]
    percs_pa = [(r_count / (1.0 * graphs[1].GetNodes())) for r_count in r_counts[1]]
    percs_actors = [(r_count / (1.0 * graphs[2].GetNodes())) for r_count in r_counts[2]]

    statistic, p_val = mannwhitneyu(percs_pa, percs_erdos)
    print "PA vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_pa, percs_actors)
    print "PA vs Actors: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_actors, percs_erdos)
    print "Actors vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)


# util methods


def run_epidemic_simulation_v4(graph, gid):
    epidemic_count = 0
    mark_for_infection, mark_for_recovery = False, False
    S, I, R = set(), set(), set()

    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        infected_node_id = randint(1, graph.GetNodes())
        I.add(infected_node_id)
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        S.discard(infected_node_id)

        while len(I) > 0:
            for u in graph.Nodes():
                if u.GetId() in S:
                    for neighbor in u.GetOutEdges():
                        if neighbor in I:
                            if decide_to_infect():
                                mark_for_infection = True
                                break
                elif u.GetId() in I:
                    if decide_to_recover():
                        mark_for_recovery = True

                # now update all at once
                if mark_for_infection:
                    S.discard(u.GetId())
                    I.add(u.GetId())
                    mark_for_infection = False
                if mark_for_recovery:
                    I.discard(u.GetId())
                    R.add(u.GetId())
                    mark_for_recovery = False
        # import pdb; pdb.set_trace()
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        if len(R) > (graph.GetNodes() / 2.0):
            print "[*] Epidemic started..."
            epidemic_count += 1
    print "[*] Proportion of simulations that infected at least 50 perc of %s network = %d" % (GRAPHS_STR_LIST[gid], epidemic_count)
    return S, I, R, epidemic_count


def run_epidemic_simulation_3b(graph, gid):
    print "3.2 (b) Compute the mean proportion of individuals infected in each network, " \
          "conditioned on the infection becoming an epidemic."
    """
    :param graph: snap.PUNGraph
    :param gid: int
    :return: list of len(R) for simulations with epidemic
    """
    epidemic_count = 0
    mark_for_infection, mark_for_recovery = False, False
    S, I, R = set(), set(), set()
    r_counts = []

    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        infected_node_id = randint(1, graph.GetNodes())
        I.add(infected_node_id)
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        S.discard(infected_node_id)

        while len(I) > 0:
            for u in graph.Nodes():
                if u.GetId() in S:
                    for neighbor in u.GetOutEdges():
                        if neighbor in I:
                            if decide_to_infect():
                                mark_for_infection = True
                                break
                elif u.GetId() in I:
                    if decide_to_recover():
                        mark_for_recovery = True

                # now update all at once
                if mark_for_infection:
                    S.discard(u.GetId())
                    I.add(u.GetId())
                    mark_for_infection = False
                if mark_for_recovery:
                    I.discard(u.GetId())
                    R.add(u.GetId())
                    mark_for_recovery = False
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        if len(R) >= (graph.GetNodes() / 2.0):
            print "[*] Epidemic started..."
            epidemic_count += 1
            r_counts.append(len(R))
    print "[*] Mean proportion of individuals infected in %s network, conditioned on the infection becoming an " \
          "epidemic = %d" % (GRAPHS_STR_LIST[gid], (sum(r_counts) / (1.0 * epidemic_count)))
    return S, I, R, r_counts


def run_epidemic_simulation_3c(graph, gid):
    """
    :param graph: snap.PUNGraph
    :param gid: int
    :return: list of len(R) for simulations with epidemic
    """
    mark_for_infection, mark_for_recovery = {}, {}  # node -> whether to infect, node -> whether to recover
    S, I, R = set(), set(), set()
    r_counts = []

    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        mark_for_infection, mark_for_recovery = {}, {}
        infected_node_id = randint(1, graph.GetNodes())
        I.add(infected_node_id)
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        S.discard(infected_node_id)

        while len(I) > 0:
            for u in graph.Nodes():
                if u.GetId() in S:
                    for neighbor in u.GetOutEdges():
                        if neighbor in I:
                            if decide_to_infect():
                                mark_for_infection[u.GetId()] = True
                                break
                elif u.GetId() in I:
                    if decide_to_recover():
                        mark_for_recovery[u.GetId()] = True

            # now update all at once
            for node_id, infection_decision in mark_for_infection.iteritems():
                S.discard(node_id)
                I.add(node_id)
                mark_for_infection[node_id] = False
            for node_id, infection_decision in mark_for_recovery.iteritems():
                I.discard(node_id)
                R.add(node_id)
                mark_for_recovery[node_id] = False
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        r_counts.append(len(R))
    return S, I, R, r_counts


def run_epidemic_simulation_3c_v2(graph, gid):
    """
    :param graph: snap.PUNGraph
    :param gid: int
    :return: list of len(R) for simulations with epidemic
    """
    S, I, R = set(), set(), set()
    r_counts = []

    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        infected_node_id = randint(1, graph.GetNodes())
        I.add(infected_node_id)
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        S.discard(infected_node_id)
        print "[*]recreated infection set.."

        while len(I) > 0:
            for u in graph.Nodes():
                # print "[*]processing node u=%d.." % u.GetId()
                if u.GetId() in S:
                    print "[*] node u=%d is in S" % u.GetId()
                    for neighbor in u.GetOutEdges():
                        # print "[*]processing neighbor=%d.." % neighbor
                        if neighbor in I:
                            print "[*] neighbor=%d is in I" % neighbor
                            if decide_to_infect():
                                print "[*] proceed to infect neighbor.."
                                S.discard(u.GetId())
                                I.add(u.GetId())
                                break
                elif u.GetId() in I:
                    print "[*] u=%d is in I" % u.GetId()
                    if decide_to_recover():
                        print "[*] proceed to recover node u=%d.." % u.GetId()
                        I.discard(u.GetId())
                        R.add(u.GetId())
                else:
                    continue
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        r_counts.append(len(R))
    return S, I, R, r_counts


# instead of selecting a random starting node, infect the node with the highest degree.
def problem_3_2(graphs, verbose=False):
    r_counts = []
    for gid, graph in enumerate(graphs):
        S, I, R, r_count_arr = run_epidemic_simulation_3_2(graph, gid)
        r_counts.append(r_count_arr)
        print "[*]%s mean infection = %f" % (GRAPHS_STR_LIST[gid], np.mean(r_count_arr))

    # find len(R) / node_count for each simulation when epidemic happened
    percs_erdos = [(r_count / (1.0 * graphs[0].GetNodes())) for r_count in r_counts[0]]
    percs_pa = [(r_count / (1.0 * graphs[1].GetNodes())) for r_count in r_counts[1]]
    percs_actors = [(r_count / (1.0 * graphs[2].GetNodes())) for r_count in r_counts[2]]

    statistic, p_val = mannwhitneyu(percs_pa, percs_erdos)
    print "PA vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_pa, percs_actors)
    print "PA vs Actors: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_actors, percs_erdos)
    print "Actors vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)


def run_epidemic_simulation_3_2(graph, gid):
    S, I, R = set(), set(), set()
    r_counts = []

    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        infected_node_id = snap.GetMxDegNId(graph)
        I.add(infected_node_id)
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        S.discard(infected_node_id)
        print "[*]recreated infection set.."

        while len(I) > 0:
            for u in graph.Nodes():
                # print "[*]processing node u=%d.." % u.GetId()
                if u.GetId() in S:
                    print "[*] node u=%d is in S" % u.GetId()
                    for neighbor in u.GetOutEdges():
                        # print "[*]processing neighbor=%d.." % neighbor
                        if neighbor in I:
                            print "[*] neighbor=%d is in I" % neighbor
                            if decide_to_infect():
                                print "[*] proceed to infect neighbor.."
                                S.discard(u.GetId())
                                I.add(u.GetId())
                                break
                elif u.GetId() in I:
                    print "[*] u=%d is in I" % u.GetId()
                    if decide_to_recover():
                        print "[*] proceed to recover node u=%d.." % u.GetId()
                        I.discard(u.GetId())
                        R.add(u.GetId())
                else:
                    continue
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        r_counts.append(len(R))
    return S, I, R, r_counts


# instead initialize the infected set to be 10 random
def problem_3_4_a(graphs, verbose=False):
    r_counts = []
    for gid, graph in enumerate(graphs):
        S, I, R, r_count_arr = run_epidemic_simulation_3_4a(graph, gid)
        r_counts.append(r_count_arr)
        print "[*] %s mean infection (10 random nodes) = %f" % (GRAPHS_STR_LIST[gid], np.mean(r_count_arr))

    print "3.4) Part a: Mean values, U statistic and P-values for infections with 10 random nodes"

    # find len(R) / node_count for each simulation when epidemic happened
    percs_erdos = [(r_count / (1.0 * graphs[0].GetNodes())) for r_count in r_counts[0]]
    percs_pa = [(r_count / (1.0 * graphs[1].GetNodes())) for r_count in r_counts[1]]
    percs_actors = [(r_count / (1.0 * graphs[2].GetNodes())) for r_count in r_counts[2]]

    statistic, p_val = mannwhitneyu(percs_pa, percs_erdos)
    print "PA vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_pa, percs_actors)
    print "PA vs Actors: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_actors, percs_erdos)
    print "Actors vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)


def run_epidemic_simulation_3_4a(graph, gid):
    S, I, R = set(), set(), set()
    r_counts = []

    # part A: start with a set of 10 random infected nodes
    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        for idx in range(10):
            I.add(randint(1, graph.GetNodes()))
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        for node_id in I:
            S.discard(node_id)
        print "[*]recreated infection set.."

        while len(I) > 0:
            for u in graph.Nodes():
                # print "[*]processing node u=%d.." % u.GetId()
                if u.GetId() in S:
                    print "[*] node u=%d is in S" % u.GetId()
                    for neighbor in u.GetOutEdges():
                        # print "[*]processing neighbor=%d.." % neighbor
                        if neighbor in I:
                            print "[*] neighbor=%d is in I" % neighbor
                            if decide_to_infect():
                                print "[*] proceed to infect neighbor.."
                                S.discard(u.GetId())
                                I.add(u.GetId())
                                break
                elif u.GetId() in I:
                    print "[*] u=%d is in I" % u.GetId()
                    if decide_to_recover():
                        print "[*] proceed to recover node u=%d.." % u.GetId()
                        I.discard(u.GetId())
                        R.add(u.GetId())
                else:
                    continue
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        r_counts.append(len(R))
    return S, I, R, r_counts


# instead initialize the top 10 highest degree nodes
def problem_3_4_b(graphs, verbose=False):
    r_counts = []
    for gid, graph in enumerate(graphs):
        S, I, R, r_count_arr = run_epidemic_simulation_3_4b(graph, gid)
        r_counts.append(r_count_arr)
        print "[*]%s mean infection = %f" % (GRAPHS_STR_LIST[gid], np.mean(r_count_arr))

    print "3.4) Part b: Mean values, U statistic and P-values for infections with top 10 nodes by degree"

    # find len(R) / node_count for each simulation when epidemic happened
    percs_erdos = [(r_count / (1.0 * graphs[0].GetNodes())) for r_count in r_counts[0]]
    percs_pa = [(r_count / (1.0 * graphs[1].GetNodes())) for r_count in r_counts[1]]
    percs_actors = [(r_count / (1.0 * graphs[2].GetNodes())) for r_count in r_counts[2]]

    statistic, p_val = mannwhitneyu(percs_pa, percs_erdos)
    print "PA vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_pa, percs_actors)
    print "PA vs Actors: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)
    statistic, p_val = mannwhitneyu(percs_actors, percs_erdos)
    print "Actors vs Erdos: Mann-Whitney U statistic = %f, P-value = %f" % (statistic, p_val)


def run_epidemic_simulation_3_4b(graph, gid):
    S, I, R = set(), set(), set()
    r_counts = []

    node_deg = {node.GetId(): node.GetOutDeg() for node in graph.Nodes()}
    top_10 = sorted(node_deg.items(), key=operator.itemgetter(1), reverse=True)[:10]    # take top 10 degree nodes

    # part A: instead initialize the top 10 highest degree nodes
    for simulation in xrange(SIMULATIONS):
        print "[*]simulation run=%d" % simulation

        # reset infected set
        S, I, R = set(), set(), set()
        for node_id, deg in enumerate(top_10):
            I.add(node_id)
        _ = [S.add(x.GetId()) for x in graph.Nodes()]
        for node_id in I:
            S.discard(node_id)
        print "[*]recreated infection set.."

        while len(I) > 0:
            for u in graph.Nodes():
                # print "[*]processing node u=%d.." % u.GetId()
                if u.GetId() in S:
                    print "[*] node u=%d is in S" % u.GetId()
                    for neighbor in u.GetOutEdges():
                        # print "[*]processing neighbor=%d.." % neighbor
                        if neighbor in I:
                            print "[*] neighbor=%d is in I" % neighbor
                            if decide_to_infect():
                                print "[*] proceed to infect neighbor.."
                                S.discard(u.GetId())
                                I.add(u.GetId())
                                break
                elif u.GetId() in I:
                    print "[*] u=%d is in I" % u.GetId()
                    if decide_to_recover():
                        print "[*] proceed to recover node u=%d.." % u.GetId()
                        I.discard(u.GetId())
                        R.add(u.GetId())
                else:
                    continue
        print "[*] Completed simulation run=%d, graph=%s, perc=%f" % (simulation, GRAPHS_STR_LIST[gid], 100 * (len(R) / (1.0 * graph.GetNodes())))
        r_counts.append(len(R))
    return S, I, R, r_counts


def infect_node(u, S, I):
    if decide_to_infect():
        S.discard(u.GetId())
        I.add(u.GetId())
    return S, I


def recover_node(u, I, R):
    if decide_to_recover():
        I.discard(u.GetId())
        R.add(u.GetId())
    return I, R


def decide_to_infect():
    if random() <= BETA:
        return True
    return False


def decide_to_recover():
    if random() <= DELTA:
        return True
    return False

if __name__ == "__main__":
    networks = init()
    print "3.1 [10 points]"
    problem_3_1_a(networks, True)
    print "3.2 [10 points]"
    problem_3_1_b(networks, True)
    problem_3_1_c(networks, True)
    print "3.2 [10 points] Initialize with max degree node"
    problem_3_2(networks, True)
    print "3.4 (A) Initialize with to 10 degree nodes [10 points]"
    problem_3_4_a(networks, True)
    print "3.4 (B) Initialize with 10 random nodes [10 points]"
    problem_3_4_b(networks, True)
