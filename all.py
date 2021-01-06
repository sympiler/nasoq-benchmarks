import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import shutil

MAX_TIMING = 1e8
SOLUTION_PRESENT = [1, "solved"]

def compute_failure_rate(status):
    tools_failure_rate = {}

    for tool in status.keys():
        tools_failure_rate[tool] = 1 - (float(status[tool].count(1) + status[tool].count("solved")) / float(len(status[tool])))
    return tools_failure_rate

def compute_speedup(status, time_stat):
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
        if "nasoq-fixed" in tool:
            ref = tool

        i = 0
        while i < len(status[tool]):
            time_stat_copy[tool].append(time_stat[tool][i])
            if status[tool][i] != 1 or status[tool][i] != "solved":
                time_stat_copy[tool][i] = MAX_TIMING
            i += 1
        gmean = geom_mean(np.array(time_stat_copy[tool]), 1.0)
        tools_speedup[tool] = gmean
    
    if not ref:
        ref = list(status.key())[0]
    
    ref_gmean = tools_speedup[ref]
    for tool in status.keys():
        tools_speedup[tool] = ref_gmean / tools_speedup[tool]
    return tools_speedup

def compute_performance_profile(status, time_stat):
    time_stat_copy = {}
    for solver in status.keys():
        n_problems = len(status[solver])
        # t[solver] = df["Time (s)"].values
        # status[solver] = df["Status"].values

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

    # Store final pandas dataframe
    # df_performance_profiles = pd.DataFrame(rho)
    # df_performance_profiles.to_csv(plots_path + "diff_variates_performance_profile_eps{}.csv".format(eps), index=False)
    return rho

def plot_histogram(data, title, tol):
    plt.figure(figsize=(10, 10))
    x = np.arange(len(list(data.keys())))
    # width = 0.25

    ordered_keys = sorted(data.keys(), key=lambda x: x.lower())
    values = [data[ordered_keys[i]] for i in range(len(ordered_keys))]
    value_series = pd.Series(values)
    
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
    
    def add_value_labels(ax, spacing=0.0):
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
            label = "{:.3f}".format(y_value)

            # Create annotation
            ax.annotate(
                label,                      # Use `label` as label
                (x_value, y_value),         # Place label at end of the bar
                xytext=(0, space),          # Vertically shift label by `space`
                textcoords="offset points", # Interpret `xytext` as offset in points
                ha='center',                # Horizontally center label
                va=va)                      # Vertically align label differently for
                                            # positive and negative values.


    # plt.bar(x, height=values, width=width)
    # Add some text for labels, title and custom x-axis tick labels, etc.
    # plt.xticks(x, ordered_keys, rotation=30)
    # plt.ylabel(title)
    # plt.title('{} of all tools for tolerence {}'.format(title, tol))

    # ax.set_ylabel(title)
    # ax.set_title('{} of all tools for tolerence {}'.format(title, tol))
    # plt.ylim(bottom=0)

    add_value_labels(ax)   

    # plt.grid()
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
        plt.plot(rho["tau"], rho[solver], label=solver)
    plt.xlim(1., 10000.)
    plt.ylim(0., 1.)
    plt.xlabel(r'Performance ratio $\tau$')
    plt.ylabel('Ratio of problems solved (tol = %s)' %tol)
    plt.xscale('log')
    plt.legend()
    plt.title("performance_profile_eps{}".format(tol))
    results_file = "all_plots/" + "performance_profile_tol{}.png".format(tol)
    plt.savefig(results_file, dpi=300)

def read_dir():
    path = sys.argv[1]
    if not path.endswith("/"):
        path = path + "/"

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

    # print(tools_failure_rate_3)
    # print(tools_failure_rate_6)
    # print(tools_speedup_3)
    # print(tools_speedup_6)
    # print(tools_performance_profile_3)
    # print(tools_performance_profile_6)
    plot_performance_profiles(tools_performance_profile_3, -3)
    plot_performance_profiles(tools_performance_profile_6, -6)
    plot_histogram(tools_failure_rate_3, "failure_rate", -3)
    plot_histogram(tools_failure_rate_6, "failure_rate", -6)
    plot_histogram(tools_speedup_3, "speedup", -3)
    plot_histogram(tools_speedup_6, "speedup", -6)

if __name__ == "__main__":
    read_dir()
