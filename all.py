import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import shutil
import glob

"""
this Python script compute all metrics based on
the performance of solvers and then generate plots on them
"""

# time for computing geometric mean if the solver does not converge
MAX_TIMING = 1e8
# the status if a solver converges at the solver
SOLUTION_PRESENT = [1, "solved"]

def compute_failure_rate(status):
    """
    compute the failure rate for each solver
    """
    tools_failure_rate = {}

    for tool in status.keys():
        # failure rate = 1- (# successful cases / # all cases)
        tools_failure_rate[tool] = 1 - (float(status[tool].count(1) + status[tool].count("solved")) / float(len(status[tool])))
    return tools_failure_rate

def compute_speedup(status, time_stat):
    """
    compute the speedup for each solver
    """
    def geom_mean(t, shift=1.0):
        """Compute the shifted geometric mean using formula from
        http://plato.asu.edu/ftp/shgeom.html
        NB. Use logarithms to avoid numeric overflows
        """
        return np.exp(np.sum(np.log(np.maximum(1, t + shift))/len(t))) - shift

    time_stat_copy = {}
    tools_speedup = {}
    ref = None

    for tool in status.keys():
        time_stat_copy[tool] = []
        # choose nasoq fixed as the ref while computing the speedup
        if "nasoq-fixed" in tool:
            ref = tool

        i = 0
        while i < len(status[tool]):
            time_stat_copy[tool].append(time_stat[tool][i])
            if status[tool][i] != 1 and status[tool][i] != "solved":
                time_stat_copy[tool][i] = MAX_TIMING
            i += 1
        # compute the geometric mean of running time for each solver
        gmean = geom_mean(np.array(time_stat_copy[tool]), 1.0)
        tools_speedup[tool] = gmean
    
    # if nasoq fixed is not used, choose the first variate searched as ref
    if not ref:
        ref = list(status.key())[0]
    
    ref_gmean = tools_speedup[ref]
    # compute the speedup for each solver
    # speedup = geometric mean of time(ref) / geometric mean of time(another tool)
    for tool in status.keys():
        tools_speedup[tool] = ref_gmean / tools_speedup[tool]
    return tools_speedup

def compute_performance_profile(status, time_stat):
    """
    compute the performance profile for each tool
    a balance between speed and convergence
    """
    time_stat_copy = {}
    for solver in status.keys():
        n_problems = len(status[solver])

        time_stat_copy[solver] = []
        for idx in range(n_problems):
            time_stat_copy[solver].append(time_stat[solver][idx])
            if status[solver][idx] not in SOLUTION_PRESENT:
                time_stat_copy[solver][idx] = MAX_TIMING
        
    r = {}  # Dictionary of relative times for each solver/problem
    for s in status.keys():
        r[s] = np.zeros(n_problems)

    # Iterate over all problems to find best timing between solvers
    for p in range(n_problems):
        # Get minimum time
        min_time = np.min([time_stat_copy[s][p] for s in status.keys()])

        # Normalize t for minimum time
        if min_time == 0:
            print("division by zero")
        for s in status.keys():
            r[s][p] = time_stat_copy[s][p] / min_time
    
    # Compute curve for all solvers
    n_tau = 1000
    tau_vec = np.logspace(0, 4, n_tau)
    rho = {'tau': tau_vec}  # Dictionary of all the curves

    for s in status.keys():
        rho[s] = np.zeros(n_tau)
        for tau_idx in range(n_tau):
            count_problems = 0  # Count number of problems with t[p, s] <= tau
            for p in range(n_problems):
                if r[s][p] <= tau_vec[tau_idx]:
                    count_problems += 1
            rho[s][tau_idx] = count_problems / n_problems

    return rho

def plot_histogram(data, title, tol):
    """
    plot the histogram plots for speedup and failure rate
    """
    plt.figure(figsize=(10, 10))
    x = np.arange(len(list(data.keys())))

    ordered_keys = sorted(data.keys(), key=lambda x: x.lower())
    values = [data[ordered_keys[i]] for i in range(len(ordered_keys))]
    value_series = pd.Series(values)
    for i in range(len(ordered_keys)):
        ordered_keys[i] = ordered_keys[i][:-4]
    
    # Plot the figure.
    plt.figure(figsize=(6, 4))
    ax = value_series.plot(kind='bar')
    ax.set_title('{} of all tools for tolerence {}'.format(title, tol))
    ax.set_xlabel('tools')
    ax.set_ylabel(title)
    ax.set_xticklabels(ordered_keys)
    ax.set_ylim(bottom=0)

    for tick in ax.get_xticklabels():
        tick.set_rotation(30)
    
    def add_value_labels(ax, title, spacing=0.0):
        """
        from: https://stackoverflow.com/questions/28931224/adding-value-labels-on-a-matplotlib-bar-chart
        a general approach to add labels on bar chart
        Add labels to the end of each bar in a bar chart.

        Arguments:
            ax (matplotlib.axes.Axes): The matplotlib object containing the axes
                of the plot to annotate.
            spacing (int): The distance between the labels and the bars.
        """

        # For each bar: Place a label
        for rect in ax.patches:
            # Get X and Y placement of label from rect.
            y_value = rect.get_height()
            x_value = rect.get_x() + rect.get_width() / 2

            # Number of points between bar and label. Change to your liking.
            space = spacing
            # Vertical alignment for positive values
            va = 'bottom'

            # If value of bar is negative: Place label below bar
            if y_value < 0:
                # Invert space to place label below
                space *= -1
                # Vertically align label at top
                va = 'top'

            # Use Y value as label and format number with one decimal place
            if title == "failure_rate":
                label = "{:.2f}%".format(y_value * 100)
            elif title == "speedup":
                label = "{:.2f}x".format(y_value)

            # Create annotation
            ax.annotate(
                label,                      # Use `label` as label
                (x_value, y_value),         # Place label at end of the bar
                xytext=(0, space),          # Vertically shift label by `space`
                textcoords="offset points", # Interpret `xytext` as offset in points
                ha='center',                # Horizontally center label
                va=va)                      # Vertically align label differently for
                                            # positive and negative values.

    # add the value on each bar
    add_value_labels(ax, title)   

    plt.tight_layout()
    plt.savefig("all_plots/" + title + "_tol{}".format(tol), dpi=300)

def plot_performance_profiles(rho, tol=''):
    """
    Plot performance profiles in matplotlib for specified problems and solvers
    """
    plt.figure(figsize=(10, 10))
    ordered_keys = sorted(rho.keys(), key=lambda x: x.lower())
    for solver in rho.keys():
        if solver == "tau": continue
        plt.plot(rho["tau"], rho[solver], label=solver[:-4])
    plt.xlim(1., 10000.)
    plt.ylim(0., 1.)
    plt.xlabel(r'Performance ratio $\tau$')
    plt.ylabel('Ratio of problems solved (tol = %s)' %tol)
    plt.xscale('log')
    plt.legend()
    plt.title("performance_profile_eps{}".format(tol))
    results_file = "all_plots/" + "performance_profile_tol{}.png".format(tol)
    plt.savefig(results_file, dpi=300)


def merge_data(path):
    """
    since the script makes performance data output csv by class
    (each csv in the raw output is for a group of QPs in one folder),
    this function merges all output csv into one by eps and tool
    and store the merged table in the subfolder "merge"
    """
    if os.path.exists(path + "/merge"):
        shutil.rmtree(path + "/merge")
    if not os.path.exists(path + "/merge"):
        os.mkdir(path + "/merge")

    p1 = "-e-3"
    p2 = "-e-6"

    solvers = ["nasoq-fixed", "nasoq-tuned", "nasoq-custom", "gurobi", "mosek", "osqp", "osqp-polished"]

    # merge csv tables by eps and solver
    for solver in solvers:
        for p in [p1, p2]:
            li = []
            for f in glob.glob(path + "/*.csv"):
                if solver in f and p in f:
                    if solver == "osqp" and "osqp-polished" in f:
                        continue
                    df = pd.read_csv(f, index_col=False, header=0, na_values="N/A", na_filter=False)
                    li.append(df)
            if li:
                frame = pd.concat(li, ignore_index=True)
                frame.to_csv(path + "/merge/" + solver + p + ".csv")


def read_dir(path):
    """
    read the merged data, compute the metrics (performance profile, speedup
    and failure rate) and make plots on them
    """
    if os.path.exists("all_plots/"):
        shutil.rmtree("all_plots/")
    os.makedirs("all_plots/")

    time_3 = {}
    status_3 = {}
    time_6 = {}
    status_6 = {}

    for f in os.listdir(path):
        if ".csv" in f:
            try:
                df = pd.read_csv(path + f)
            except pd.errors.EmptyDataError:
                print(path + f + " is empty, a seg fault may occur")
                continue
            f_snippet = f[:-4]
             
            if "-3" in f:
                time_3[f_snippet] = df["Time (s)"].values.tolist()
                status_3[f_snippet] = df["Status"].values.tolist()
            elif "-6" in f:
                time_6[f_snippet] = df["Time (s)"].values.tolist()
                status_6[f_snippet] = df["Status"].values.tolist()

            if "-3" in f:
                time_3[f_snippet] = df["Time (s)"].values.tolist()
                status_3[f_snippet] = df["Status"].values.tolist()
            elif "-6" in f:
                time_6[f_snippet] = df["Time (s)"].values.tolist()
                status_6[f_snippet] = df["Status"].values.tolist()
    
    tools_failure_rate_3 = compute_failure_rate(status_3)
    tools_failure_rate_6 = compute_failure_rate(status_6)
    tools_speedup_3 = compute_speedup(status_3, time_3)
    tools_speedup_6 = compute_speedup(status_6, time_6)
    tools_performance_profile_3 = compute_performance_profile(status_3, time_3)
    tools_performance_profile_6 = compute_performance_profile(status_6, time_6)

    plot_performance_profiles(tools_performance_profile_3, -3)
    plot_performance_profiles(tools_performance_profile_6, -6)
    plot_histogram(tools_failure_rate_3, "failure_rate", -3)
    plot_histogram(tools_failure_rate_6, "failure_rate", -6)
    plot_histogram(tools_speedup_3, "speedup", -3)
    plot_histogram(tools_speedup_6, "speedup", -6)

if __name__ == "__main__":
    path = sys.argv[1]
    if not path.endswith("/"):
        path = path + "/"
    merge_data(path)
    read_dir(path + "/merge/")
