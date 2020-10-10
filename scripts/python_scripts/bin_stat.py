import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from collections import defaultdict
import csv

# count how many QPs are held for each bin

dir_analysis_data = "SMP_Repository_analysis_tables/merged_csvs.csv"
stat_dir = "SMP_statistical_count/"

# bin of data based on the value of metrics in SMP_Repository
norm_bin = [(0, 5), (5, 10), (10, 100), (100, 500), (500, 1500), (1500, 2500), \
    (2500, 4500), (4500, 10000), (10000, 100000), (100000, np.inf)]
constrs_bin = [(0, 5), (5, 10), (10, 50), (50, 100), (100, 200), (200, 500), (500, 750), \
    (750, 1000), (1000, 1500), (1500, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, np.inf)]
numerical_range_bin = [(0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 5.0), (5.0, 10.0), (10.0, 100.0), (100.0, np.inf)]

header = ["bin", "count"]

def find_correct_bin(bin_lst, val):
    """
    find the bin to store the given value
    """
    start = 0
    end = len(bin_lst)

    while start != end:
        mid = (start + end) // 2
        if bin_lst[mid][0] <= val and val < bin_lst[mid][1]:
            return bin_lst[mid]
        elif bin_lst[mid][0] > val:
            end = mid
        elif bin_lst[mid][1] < val:
            start = mid
        elif bin_lst[mid][1] == val:
            return bin_lst[mid + 1]

def count_by_bin(dir_analysis_data):
    """
    count how many QPs in each bin
    """
    norm_stat, constrs_stat, numerical_range_stat = defaultdict(lambda: 0), defaultdict(lambda: 0), defaultdict(lambda: 0)
    df = pd.read_csv(dir_analysis_data)

    norm_lst = df["QP norm (order = infty)"].values

    ineq = df["Number of inequality constraints"].values
    i = 0
    while i < len(ineq):
        if ineq[i] == "\\":
            ineq[i] = 0
        i += 1
    ineq = np.int32(ineq)
    eq = np.array(df["Number of equality constraints"].values)
    i = 0
    while i < len(eq):
        if eq[i] == "\\":
            eq[i] = 0
        i += 1
    eq = np.int32(eq)
    constr_lst = ineq + eq

    qp_min = df["QP min"].values
    qp_max = df["QP max"].values
    numerical_range_lst = np.abs(qp_max / qp_min)

    for norm_val in norm_lst:
        norm_stat[find_correct_bin(norm_bin, norm_val)] += 1
    
    for constr_count in constr_lst:
        constrs_stat[find_correct_bin(constrs_bin, constr_count)] += 1
    
    for numerical_range in numerical_range_lst:
        numerical_range_stat[find_correct_bin(numerical_range_bin, numerical_range)] += 1
    
    plt.figure(figsize=(10, 10))
    norm_lst = np.array(list(norm_stat.values()))
    norm_keys = sorted(list(norm_stat.keys()), key=lambda x: x[0])
    sorted_indices = sorted(range(len(list(norm_stat.keys()))), key=lambda k: list(norm_stat.keys())[k][0])

    idx = list(range(norm_lst.shape[0]))
    plt.bar(idx, np.array(list(norm_stat.values()))[sorted_indices])
    plt.xticks(idx, norm_keys, rotation=30)
    plt.title("norm stat")
    plt.xlabel("norm bin")
    plt.ylabel("count")
    plt.savefig(stat_dir + "norm_stat.png")

    plt.figure(figsize=(10, 10))
    constr_lst = np.array(list(constrs_stat.values()))
    constr_keys = sorted(list(constrs_stat.keys()), key=lambda x: x[0])
    sorted_indices = sorted(range(len(list(constrs_stat.keys()))), key=lambda k: list(constrs_stat.keys())[k][0])

    idx = list(range(constr_lst.shape[0]))
    plt.bar(idx, np.array(list(constrs_stat.values()))[sorted_indices])
    plt.xticks(idx, constr_keys, rotation=30)
    plt.title("constr number stat")
    plt.xlabel("constr number bin")
    plt.ylabel("count")
    plt.savefig(stat_dir + "constrs_stat.png")
    
    plt.figure(figsize=(10, 10))
    numerical_range_lst = np.array(list(numerical_range_stat.values()))
    numerical_range_keys = sorted(list(numerical_range_stat.keys()), key=lambda x: x[0])
    sorted_indices = sorted(range(len(list(numerical_range_stat.keys()))), key=lambda k: list(numerical_range_stat.keys())[k][0])

    idx = list(range(numerical_range_lst.shape[0]))
    plt.bar(idx, np.array(list(numerical_range_stat.values()))[sorted_indices])
    plt.xticks(idx, numerical_range_keys, rotation=30)
    plt.title("numerical range stat")
    plt.xlabel("numerical range bin")
    plt.ylabel("count")
    plt.savefig(stat_dir + "numerical_range_stat.png")

    # fout1 = open()

if __name__ == "__main__":
    if not os.path.exists(stat_dir):
        os.makedirs(stat_dir)
    count_by_bin(dir_analysis_data)