import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from collections import defaultdict 

dir_params = "test_smp_tune_csvs/"
dir_analysis_data = "test_smp_analysis_tables/merged_csvs.csv"

plot_path = "test_smp_metrics_plots/"

# the status of successful cases
OPTIMAL = 1
SOLUTION_PRESENT = [OPTIMAL]
# the max timming for success
MAX_VAL = 1600
# run_time for failed cases
MAX_TIMING = 1e8

# all tuned parameters
mode_lst = ("max_iter", "stop_tol", "diag_perturb")

# all nasoq tools
tool_lst = ("nasoq-fixed", "nasoq-tuned", "nasoq-custom")

# all metrics tested
metric_lst = ("norm", "nnz", "constr", "numerical_range")

# all values of parameters to be tested against
max_iter_lst = sorted([0, 1, 2, 3, 4, 5, 10, 20])
stop_tol_lst = sorted([-13, -15, -16, -17])
diag_perturb_lst = sorted([-6, -7, -8, -9, -10, -11, -12])

# eps tested
eps_lst = sorted([-3, -6])
eps_realval_lst = sorted([1e-3, 1e-6])

# gmean of reference tool (nasoq-fixed) here
# gmean_ref_without = defaultdict(lambda: np.zeros(6))
# gmean_ref_with = defaultdict(lambda: np.zeros(6))

# bin of data based on the value of metrics in SMP_Repository
norm_bin = [(0, 5), (5, 10), (10, 100), (100, 500), (500, 1500), (1500, 2500), \
    (2500, 4500), (4500, 10000), (10000, 100000), (100000, np.inf)]
nnz_bin = [(0, 100), (100, 500), (500, 1000), (1000, 5000), (10000, 20000), \
    (20000, 30000), (30000, 45000), (45000, 60000), (60000, 150000), (150000, np.inf)]
constrs_bin = [(0, 5), (5, 10), (10, 50), (50, 100), (100, 200), (200, 500), (500, 750), \
    (750, 1000), (1000, 1500), (1500, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, np.inf)]
numerical_range_bin = [(0, 1.5), (1.5, 2), (2, 2,5), (2.5, 3), (3, 5), (5, 10), (10, 100), (100, np.inf)]

# geometric mean of reference tool (nasoq-fixed here)
ref_gmens_with = {}
ref_gmens_without = {}

def find_gmean(param, eps):
    """
    find where to insert and extract ref-gmean (nasoq-fixed usually)
    """
    return len(mode_lst) * eps_lst.index(eps) + mode_lst.index(param)

def compute_numerical_range(path="SMP_Repository_analysis_tables/merged_csvs.csv"):
    """
    compute numerical range of each QP and write to a txt
    """
    df_analysis = pd.read_csv(path)
    qp_min = df_analysis["QP min"].values
    qp_max = df_analysis["QP max"].values
    numerical_range = np.abs(qp_max / qp_min)
    
    with open("SMP_numerical_range.txt", "w") as fp:
        i = 0
        for n in numerical_range:
            fp.write(str(i) + ": " + str(n) + '\n')
            i += 1

def geom_mean_ignore_failure(run_times, shift=10.):
        """Compute the shifted geometric mean using formula from
        http://plato.asu.edu/ftp/shgeom.html
        NB. Use logarithms to avoid numeric overflows
        """
        # return np.exp(np.sum(np.log(np.maximum(1, run_times + shift))/len(run_times))) - shift

        # compute the sum of running time for successful cases only
        # and also the number of successful cases only
        sum_time = 0.0
        length_time = 0
        for time in run_times:
            if time <= MAX_VAL:
                sum_time += np.log(np.maximum(1, time + shift))
                length_time += 1
        
        # if the solver succeed at no cases, return 0
        if sum_time == 0.0:
            return sum_time
        
        # return the geometric mean
        return np.exp(sum_time / length_time) - shift
    
def geom_mean_with_failure(run_times, shift=10.):
    """Compute the shifted geometric mean using formula from
    http://plato.asu.edu/ftp/shgeom.html
    NB. Use logarithms to avoid numeric overflows
    """
    return np.exp(np.sum(np.log(np.maximum(1, run_times + shift))/len(run_times))) - shift

def compute_frate_by_bin(stat_dict):
    """
    find the number of failed cases by the number of problems solved by
    the solver under such a condition do not have the optimal status
    """
    for key, value in stat_dict.items():
        n_problems = len(value)
        failed_statuses = np.logical_and([np.array(value) != s
                                    for s in SOLUTION_PRESENT],True)
        n_failed_problems = np.sum(failed_statuses)
        # the failure rate for the tool at the tuned parameters and eps
        failure_rate = n_failed_problems / n_problems
        stat_dict[key] = failure_rate
    return stat_dict

def find_correct_bin(bin_lst, val):
    """
    find the bin to store the given value
    """
    for bin_t in bin_lst:
        if val >= bin_t[0] and val < bin_t[1]:
            return bin_t

def compute_gmean_by_bin(stat_dict, with_failure=True):
    """
    compute the geometric mean for each bin
    """
    result_stat = {}
    for key, value in stat_dict.items():
        if with_failure:
            result_stat[key] = geom_mean_with_failure(np.array(value), 1.0)
        else:
            result_stat[key] = geom_mean_ignore_failure(np.array(value), 1.0)
    return result_stat

def compute_speedup(gmean, metric, bin_t, param, param_val, eps, tool="nasoq-fixed", with_failure=True):
    """
    compute the speedup of the tool at the given condition (find gmean of reference tool at
    given conditions)
    """
    # format:
    # ref_gmens_with[eps][param][param_val][metric][bin_t] = nasoq_gmen_with_failure
    # ref_gmens_without[eps][param][param_val][metric][bin_t] = nasoq_gmen_without_failure
    if tool == "nasoq-fixed":
        # if nasoq-fixed failed on all given QPs under the circumstance, speedup = 0
        if not with_failure: # without failure
            if gmean == 0.0: speedup = 0
            else: speedup = 1
            ref_gmens_without[eps][param][param_val][metric][bin_t] = gmean
        else:
            speedup = 1
            ref_gmens_with[eps][param][param_val][metric][bin_t] = gmean
    else:
        if not with_failure: # without failure
            if gmean == 0.0:
                speedup = np.inf
            elif ref_gmens_without[eps][param][param_val][metric][bin_t] == 0.0:
                speedup = 0.0
            else:
                speedup = ref_gmens_without[eps][param][param_val][metric][bin_t] / gmean
        else:
            speedup = ref_gmens_with[eps][param][param_val][metric][bin_t] / gmean
    return speedup

def plot_speedup(param="diag_perturb", tool="nasoq-fixed", eps=-3):
    """
    plot the speedup (with and without failure cases) for the
    given parameter, tool and required accuracy
    """
    # read out the data needed and transfer to problem name to lowercase letter for merging
    df_analysis = pd.read_csv(dir_analysis_data)
    df_analysis["Problem name"] = df_analysis["Problem name"].str.lower()

    # choose the list of tested values by the parameter tuned
    if param == "max_iter":
        lst = max_iter_lst
    elif param == "diag_perturb":
        lst = diag_perturb_lst
    elif param == "stop_tol":
        lst = stop_tol_lst

    # list all csv files and make graphs
    perfcsv_for_param = os.listdir(dir_params + param)
    fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(15, 15))

    for csv_name in perfcsv_for_param:
        # find the specified csv
        if param in csv_name and tool in csv_name and "eps{}".format(eps) in csv_name and csv_name.endswith(".csv"):
            # merge the csv for plotting like SQL
            # easy to extract data
            df_param = pd.read_csv(dir_params + param + "/" + csv_name)
            df_param["Problem Name"] = df_param["Problem Name"].str.lower()
            df_merged = df_analysis.merge(df_param, left_on="Problem name", right_on="Problem Name", how="left")
            
            # read out the data needed for plotting (metrics required)
            status_lst = df_merged["Status"].values
            time_lst = df_merged["Time (s)"].values
            n_problem = len(df_merged)
            for i in range(n_problem):
                if status_lst[i] not in SOLUTION_PRESENT:
                    time_lst[i] = MAX_TIMING
    
            norm_lst = df_merged["QP norm (order = infty)"].values
            nnz_lst = df_merged["QP nonzeros"].values
            
            ineq = df_merged["Number of inequality constraints"].values
            i = 0
            while i < len(ineq):
                if ineq[i] == "\\":
                    ineq[i] = 0
                i += 1
            ineq = np.int32(ineq)
            eq = np.array(df_merged["Number of equality constraints"].values)
            i = 0
            while i < len(eq):
                if eq[i] == "\\":
                    eq[i] = 0
                i += 1
            eq = np.int32(eq)
            constr_lst = ineq + eq

            qp_min = df_analysis["QP min"].values
            qp_max = df_analysis["QP max"].values
            num_range_lst = np.abs(qp_max / qp_min)

            # build dict for statistics of each bin
            norm_stat, nnz_stat, constr_stat, num_range_stat = defaultdict(list), \
                defaultdict(list), defaultdict(list), defaultdict(list)
            
            # store the statistics needed (status for each qp at each metric) into proper bins
            i = 0
            while i < len(df_merged):
                norm_stat[find_correct_bin(norm_bin, norm_lst[i])].append(time_lst[i])
                nnz_stat[find_correct_bin(nnz_bin, nnz_lst[i])].append(time_lst[i])
                constr_stat[find_correct_bin(constrs_bin, constr_lst[i])].append(time_lst[i])
                num_range_stat[find_correct_bin(numerical_range_bin, num_range_lst[i])].append(time_lst[i])

                i += 1
            
            # find the value of the tuned parameter for the case to write into csv
            # one value of tuned parameter for the given tool at given eps
            # (the tuned parameter value for the given csv data from archived output of NASOQ-BIN)
            param_val = lst[0]
            for val in lst:
                if param + str(val) in csv_name and \
                    param + str(val) + '0' not in csv_name:
                    param_val = val
                    break
            
            # compute the gmean in each bin (with and without failure)
            norm_stat_with = compute_gmean_by_bin(norm_stat, True)
            nnz_stat_with = compute_gmean_by_bin(nnz_stat, True)
            constr_stat_with = compute_gmean_by_bin(constr_stat, True)
            num_range_stat_with = compute_gmean_by_bin(num_range_stat, True)

            norm_stat_without = compute_gmean_by_bin(norm_stat, False)
            nnz_stat_without = compute_gmean_by_bin(nnz_stat, False)
            constr_stat_without = compute_gmean_by_bin(constr_stat, False)
            num_range_stat_without = compute_gmean_by_bin(num_range_stat, False)

            # compute speedup in each bin
            # compute_speedup(gmean, metric, bin_t, param, param_val, eps, tool="nasoq-fixed", with_failure=True)
            # metric_lst = ["norm", "nnz", "constr", "numerical_range"]
            norm_array_with = [compute_speedup(norm_stat_with[bin_t], "norm", bin_t, param, param_val, eps, tool, True) for \
                bin_t in norm_bin if bin_t in norm_stat_with]
            nnz_array_with = [compute_speedup(nnz_stat_with[bin_t], "nnz", bin_t, param, param_val, eps, tool, True) for \
                bin_t in nnz_bin if bin_t in nnz_stat_with]
            constr_array_with = [compute_speedup(constr_stat_with[bin_t], "constr", bin_t, param, param_val, eps, tool, True) for \
                bin_t in constrs_bin if bin_t in constr_stat_with]
            num_range_array_with = [compute_speedup(num_range_stat_with[bin_t], "numerical_range", bin_t, param, param_val, eps, tool, True) for \
                bin_t in numerical_range_bin if bin_t in num_range_stat_with]

            norm_array_without = [compute_speedup(norm_stat_without[bin_t], "norm", bin_t, param, param_val, eps, tool, False) for \
                bin_t in norm_bin if bin_t in norm_stat_without]
            nnz_array_without = [compute_speedup(nnz_stat_without[bin_t], "nnz", bin_t, param, param_val, eps, tool, False) for \
                bin_t in nnz_bin if bin_t in nnz_stat_without]
            constr_array_without = [compute_speedup(constr_stat_without[bin_t], "constr", bin_t, param, param_val, eps, tool, False) for \
                bin_t in constrs_bin if bin_t in constr_stat_without]
            num_range_array_without = [compute_speedup(num_range_stat_without[bin_t], "numerical_range", bin_t, param, param_val, eps, tool, False) for \
                bin_t in numerical_range_bin if bin_t in num_range_stat_without]

            # plot the performance in each subplot and save
            axs[0, 0].plot(np.arange(len(norm_array_without)), norm_array_without, label=param + "={}".format(param_val))
            axs[0, 1].plot(np.arange(len(norm_array_with)), norm_array_with, label=param + "={}".format(param_val))
            axs[1, 0].plot(np.arange(len(nnz_array_without)), nnz_array_without, label=param + "={}".format(param_val))
            axs[1, 1].plot(np.arange(len(nnz_array_with)), nnz_array_with, label=param + "={}".format(param_val))

            axs[2, 0].plot(np.arange(len(constr_array_without)), constr_array_without, label=param + "={}".format(param_val))
            axs[2, 1].plot(np.arange(len(constr_array_with)), constr_array_with, label=param + "={}".format(param_val))
            axs[3, 0].plot(np.arange(len(num_range_array_without)), num_range_array_without, label=param + "={}".format(param_val))
            axs[3, 1].plot(np.arange(len(num_range_array_with)), num_range_array_with, label=param + "={}".format(param_val))

            axs[0, 0].legend(); axs[0, 1].legend(); axs[1, 0].legend(); axs[1, 1].legend()
            axs[2, 0].legend(); axs[2, 1].legend(); axs[3, 0].legend(); axs[3, 1].legend()

            axs[0, 0].set_title("{} speedup vs norm \n of QP (eps = {}) over {} (without failure)".format(tool, eps, param));
            axs[0, 1].set_title("{} speedup vs norm \n of QP (eps = {}) over {} (with failure)".format(tool, eps, param));
            axs[1, 0].set_title("{} speedup vs nnz \n of QP (eps = {}) over {} (without failure)".format(tool, eps, param));
            axs[1, 1].set_title("{} speedup vs nnz \n of QP (eps = {}) over {} (with failure)".format(tool, eps, param));

            axs[2, 0].set_title("{} speedup vs number of constraints \n of QP (eps = {}) over {} (without failure)".format(tool, eps, param));
            axs[2, 1].set_title("{} speedup vs number of constraints \n of QP (eps = {}) over {} (with failure)".format(tool, eps, param));
            axs[3, 0].set_title("{} speedup vs numerical range \n of QP (eps = {}) over {} (without failure)".format(tool, eps, param));
            axs[3, 1].set_title("{} speedup vs numerical range \n of QP (eps = {}) over {} (with failure)".format(tool, eps, param));

            axs[0, 0].set(xlabel="norm bin", ylabel="speedup"); axs[0, 1].set(xlabel="norm bin", ylabel="speedup")
            axs[1, 0].set(xlabel="nnz bin", ylabel="speedup"); axs[1, 1].set(xlabel="nnz bin", ylabel="speedup")
            axs[2, 0].set(xlabel="num of constraints bin", ylabel="speedup"); axs[2, 1].set(xlabel="num of constraints bin", ylabel="speedup");
            axs[3, 0].set(xlabel="numerical range bin", ylabel="speedup"); axs[3, 1].set(xlabel="numerical range bin", ylabel="speedup");

            plt.sca(axs[0, 0]); plt.xticks(np.arange(len(norm_array_without)), sorted(list(norm_stat_without.keys()), key=lambda t: t[0]))
            plt.sca(axs[0, 1]); plt.xticks(np.arange(len(norm_array_with)), sorted(list(norm_stat_with.keys()), key=lambda t: t[0]))
            plt.sca(axs[1, 0]); plt.xticks(np.arange(len(nnz_array_without)), sorted(list(nnz_stat_without.keys()), key=lambda t: t[0]))
            plt.sca(axs[1, 1]); plt.xticks(np.arange(len(nnz_array_with)), sorted(list(nnz_stat_with.keys()), key=lambda t: t[0]))

            plt.sca(axs[2, 0]); plt.xticks(np.arange(len(constr_array_without)), sorted(list(constr_stat_without.keys()), key=lambda t: t[0]))
            plt.sca(axs[2, 1]); plt.xticks(np.arange(len(constr_array_with)), sorted(list(constr_stat_with.keys()), key=lambda t: t[0]))
            plt.sca(axs[3, 0]); plt.xticks(np.arange(len(num_range_array_without)), sorted(list(num_range_stat_without.keys()), key=lambda t: t[0]))
            plt.sca(axs[3, 1]); plt.xticks(np.arange(len(num_range_array_with)), sorted(list(num_range_stat_with.keys()), key=lambda t: t[0]))

            
            print("csv name: {}".format(csv_name), norm_array_with)
            print("csv name: {}".format(csv_name), norm_array_without)

    plt.tight_layout()
    plt.savefig(plot_path + "speedup/" + "speedup_{}_{}_eps{}.png".format(tool, param, eps))

def plot_failure_rate(param="diag_perturb", tool="nasoq-custom", eps=-6):
    """
    plot the failure rate for the given parameter, tool and required accuracy
    """
    # read out the data needed and transfer to problem name to lowercase letter for merging
    df_analysis = pd.read_csv(dir_analysis_data)
    df_analysis["Problem name"] = df_analysis["Problem name"].str.lower()

    # choose the list of tested values by the parameter tuned
    if param == "max_iter":
        lst = max_iter_lst
    elif param == "diag_perturb":
        lst = diag_perturb_lst
    elif param == "stop_tol":
        lst = stop_tol_lst

    perfcsv_for_param = os.listdir(dir_params + param)
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 15))

    for csv_name in perfcsv_for_param:
        # find the specified csv
        if param in csv_name and tool in csv_name and "eps{}".format(eps) in csv_name:
            # merge the csv for plotting like SQL
            # easy to extract data
            df_param = pd.read_csv(dir_params + param + "/" + csv_name)
            df_param["Problem Name"] = df_param["Problem Name"].str.lower()
            df_merged = df_analysis.merge(df_param, left_on="Problem name", right_on="Problem Name", how="left")
            
            # read out the data needed for plotting (metrics required)
            status_lst = df_merged["Status"].values
            norm_lst = df_merged["QP norm (order = infty)"].values
            nnz_lst = df_merged["QP nonzeros"].values
            
            ineq = df_merged["Number of inequality constraints"].values
            i = 0
            while i < len(ineq):
                if ineq[i] == "\\":
                    ineq[i] = 0
                i += 1
            ineq = np.int32(ineq)
            eq = np.array(df_merged["Number of equality constraints"].values)
            i = 0
            while i < len(eq):
                if eq[i] == "\\":
                    eq[i] = 0
                i += 1
            eq = np.int32(eq)
            constr_lst = ineq + eq

            qp_min = df_analysis["QP min"].values
            qp_max = df_analysis["QP max"].values
            num_range_lst = np.abs(qp_max / qp_min)

            # store the statistics needed (status for each qp at each metric) into proper bins
            norm_stat, nnz_stat, constr_stat, num_range_stat = defaultdict(list), \
                defaultdict(list), defaultdict(list), defaultdict(list)
            
            # store the information needed (status for each qp at each metric) into proper bins
            i = 0
            while i < len(df_merged):
                norm_stat[find_correct_bin(norm_bin, norm_lst[i])].append(status_lst[i])
                nnz_stat[find_correct_bin(nnz_bin, nnz_lst[i])].append(status_lst[i])
                constr_stat[find_correct_bin(constrs_bin, constr_lst[i])].append(status_lst[i])
                num_range_stat[find_correct_bin(numerical_range_bin, num_range_lst[i])].append(status_lst[i])

                i += 1
            
            # compute the failure rate in each bin
            norm_frate_stat = compute_frate_by_bin(norm_stat)
            nnz_frate_stat = compute_frate_by_bin(nnz_stat)
            constr_frate_stat = compute_frate_by_bin(constr_stat)
            num_range_frate_stat = compute_frate_by_bin(num_range_stat)

            # read out the failure rate
            norm_frate_array = [norm_frate_stat[bin_t] for bin_t in norm_bin if bin_t in norm_frate_stat]
            nnz_frate_array = [nnz_frate_stat[bin_t] for bin_t in nnz_bin if bin_t in nnz_frate_stat]
            constr_frate_array = [constr_frate_stat[bin_t] for bin_t in constrs_bin if bin_t in constr_frate_stat]
            num_range_frate_array = [num_range_frate_stat[bin_t] for bin_t in numerical_range_bin if bin_t in num_range_frate_stat]

            # find the value of the tuned parameter for the case to write into csv
            # one value of tuned parameter for the given tool at given eps
            # (the tuned parameter value for the given csv data from archived output of NASOQ-BIN)
            param_val = lst[0]
            for val in lst:
                if param + str(val) in csv_name and \
                    param + str(val) + '0' not in csv_name:
                    param_val = val
                    break

            # plot the performance in each subplot and save
            axs[0][0].plot(np.arange(len(norm_frate_array)), norm_frate_array, label=param + "={}".format(param_val))
            axs[0][1].plot(np.arange(len(nnz_frate_array)), nnz_frate_array, label=param + "={}".format(param_val))
            axs[1][0].plot(np.arange(len(constr_frate_array)), constr_frate_array, label=param + "={}".format(param_val))
            axs[1][1].plot(np.arange(len(num_range_frate_array)), num_range_frate_array, label=param + "={}".format(param_val))
            axs[0][0].legend(); axs[0][1].legend(); axs[1][0].legend(); axs[1][1].legend()
            axs[0][0].set_title("{} failure rate vs norm \n of QP (eps = {}) over {}".format(tool, eps, param));
            axs[0][1].set_title("{} failure rate vs nnz \n of QP (eps = {}) over {}".format(tool, eps, param));
            axs[1][0].set_title("{} failure rate vs number of constraints \n of QP (eps = {}) over {}".format(tool, eps, param));
            axs[1][1].set_title("{} failure rate vs numerical range \n of QP (eps = {}) over {}".format(tool, eps, param));
            axs[0][0].set(xlabel="norm bin", ylabel="failure rate"); axs[0][1].set(xlabel="nnz bin", ylabel="failure rate")
            axs[1][0].set(xlabel="num of constraints bin", ylabel="failure rate"); axs[1][1].set(xlabel="numerical range bin", ylabel="failure rate")
            plt.sca(axs[0, 0]); plt.xticks(np.arange(len(norm_frate_array)), sorted(list(norm_frate_stat.keys()), key=lambda t: t[0]))
            plt.sca(axs[0, 1]); plt.xticks(np.arange(len(nnz_frate_array)), sorted(list(nnz_frate_stat.keys()), key=lambda t: t[0]))
            plt.sca(axs[1, 0]); plt.xticks(np.arange(len(constr_frate_array)), sorted(list(constr_frate_stat.keys()), key=lambda t: t[0]))
            plt.sca(axs[1, 1]); plt.xticks(np.arange(len(num_range_frate_array)), sorted(list(num_range_frate_stat.keys()), key=lambda t: t[0]))

    plt.tight_layout()
    plt.savefig(plot_path + "failure_rate/" + "failure_rate_{}_{}_eps{}.png".format(tool, param, eps))

def main():
    # build directory for plots
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
        os.makedirs(plot_path + "failure_rate/")
        os.makedirs(plot_path + "speedup/")
    
    # store geometric mean of reference tool at a specific condition
    for eps in eps_lst:
        ref_gmens_with[eps] = {}
        ref_gmens_without[eps] = {}
        for param in mode_lst:
            ref_gmens_with[eps][param] = {}
            ref_gmens_without[eps][param] = {}
            if param == "max_iter":
                lst = max_iter_lst
            elif param == "diag_perturb":
                lst = diag_perturb_lst
            elif param == "stop_tol":
                lst = stop_tol_lst
            for item in lst:
                ref_gmens_with[eps][param][item] = {}
                ref_gmens_without[eps][param][item] = {}
                for metric in metric_lst:
                    ref_gmens_with[eps][param][item][metric] = {}
                    ref_gmens_without[eps][param][item][metric] = {}
                    if metric == "norm":
                        for bin_t in norm_bin:
                            ref_gmens_with[eps][param][item][metric][bin_t] = 0.0
                            ref_gmens_without[eps][param][item][metric][bin_t] = 0.0
                    elif metric == "nnz":
                        for bin_t in nnz_bin:
                            ref_gmens_with[eps][param][item][metric][bin_t] = 0.0
                            ref_gmens_without[eps][param][item][metric][bin_t] = 0.0
                    elif metric == "constr":
                        for bin_t in constrs_bin:
                            ref_gmens_with[eps][param][item][metric][bin_t] = 0.0
                            ref_gmens_without[eps][param][item][metric][bin_t] = 0.0
                    elif metric == "numerical_range":
                        ref_gmens_with[eps][param][item][metric][bin_t] = 0.0
                        ref_gmens_without[eps][param][item][metric][bin_t] = 0.0

    # compute the failure rate
    for param in mode_lst:
        for tool in tool_lst:
            for eps in eps_lst:
                plot_failure_rate(param, tool, eps)
    
    # compute the speedup
    # CAUTION: in tool_lst, "nasoq-fixed" must locate at the head for reference purpose
    for tool in tool_lst:
        for param in mode_lst:
            for eps in eps_lst:
                plot_speedup(param, tool, eps)

if __name__ == "__main__":
    # plot_failure_rate()
    # compute_numerical_range()
    main()
