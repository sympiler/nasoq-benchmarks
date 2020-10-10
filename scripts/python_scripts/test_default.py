import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

MAX_TIMING = 1e8
SOLUTION_PRESENT = [1]

path1 = "test_smp_test_settings_csvs/"
path2 = "test_smp_tune_csvs/"
plots_path = "test_smp_opt_test_plots/"

# recommended settings:
# diag_perturb = -8, max_iter = 3

def geom_mean(t, shift=10.):
    """Compute the shifted geometric mean using formula from
    http://plato.asu.edu/ftp/shgeom.html
    NB. Use logarithms to avoid numeric overflows
    """
    return np.exp(np.sum(np.log(np.maximum(1, t + shift))/len(t))) - shift

def compute_speedup1():
    """
    compare the speedup between Nasoq-custom with recommended settings
    and two other NASOQ variates with the recommended settings
    """
    data = os.listdir(path1)
    custom_data = []
    for datum in data:
        if "nasoq-custom" in datum:
            custom_data.append(datum)

    ref_gmean = {}
    
    for custom_datum in custom_data:
        if "eps-3" in custom_datum:
            df = pd.read_csv(path1 + custom_datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean1 = geom_mean(time, 1.0)
            ref_gmean["eps-3"] = gmean1
        
        elif "eps-6" in custom_datum:
            df = pd.read_csv(path1 + custom_datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean2 = geom_mean(time, 1.0)
            ref_gmean["eps-6"] = gmean2

    
    speedup = {"nasoq-custom_eps-3": 1.0, "nasoq-custom_eps-6": 1.0}
    for datum in data:
        if "nasoq-custom" not in datum:
            df = pd.read_csv(path1 + datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean = geom_mean(time, 1.0)

            if "eps-3" in datum:
                if "nasoq-fixed" in datum:
                    speedup["nasoq-fixed_eps-3"] = ref_gmean["eps-3"] / gmean
                elif "nasoq-tuned" in datum:
                    speedup["nasoq-tuned_eps-3"] = ref_gmean["eps-3"] / gmean
            
            elif "eps-6" in datum:
                if "nasoq-fixed" in datum:
                    speedup["nasoq-fixed_eps-6"] = ref_gmean["eps-6"] / gmean
                elif "nasoq-tuned" in datum:
                    speedup["nasoq-tuned_eps-6"] = ref_gmean["eps-6"] / gmean
        
    eps3_speedup = [speedup["nasoq-fixed_eps-3"], speedup["nasoq-tuned_eps-3"], speedup["nasoq-custom_eps-3"]]
    eps6_speedup = [speedup["nasoq-fixed_eps-6"], speedup["nasoq-tuned_eps-6"], speedup["nasoq-custom_eps-6"]]

    plt.figure()
    labels = ['nasoq-fixed', 'nasoq-tuned', 'nasoq-custom']
    x = np.arange(3)
    width = 0.35

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, eps3_speedup, width, label='eps-3')
    rects2 = ax.bar(x + width/2, eps6_speedup, width, label='eps-6')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('speedup')
    ax.set_title('speedup of three tools with recommended settings')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    plt.savefig(plots_path + "different_solvers.png")

def compute_speedup2():
    """
    compare the speedup between Nasoq-custom with recommended settings
    and other NASOQ custom with settings previously set
    """
    data = os.listdir(path1)
    custom_data = []
    for datum in data:
        if "nasoq-custom" in datum:
            custom_data.append(datum)
    data1 = os.listdir(path2 + "diag_perturb/")
    data2 = os.listdir(path2 + "max_iter/")

    ref_gmean = {}
    for custom_datum in custom_data:
        if "eps-3" in custom_datum:
            df = pd.read_csv(path1 + custom_datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean1 = geom_mean(time, 1.0)
            ref_gmean["eps-3"] = gmean1
        
        elif "eps-6" in custom_datum:
            df = pd.read_csv(path1 + custom_datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean2 = geom_mean(time, 1.0)
            ref_gmean["eps-6"] = gmean2
    
    diag_perturb_lst = ["-12", "-11", "-10", "-9", "-8", "-7", "-6"]
    max_iter_lst = ["0", "1", "2", "3", "4", "5", "10", "20"]

    speedup1 = {"tested-3": 1.0, "tested-6": 1.0}
    labels1 = ["optimized nasoq-custom", "diag_perturb -12", "diag_perturb -11", "diag_perturb -10", "diag_perturb -9", \
        "diag_perturb -8", "diag_perturb -7", "diag_perturb -6"]

    for datum in data1:
        if "nasoq-custom" in datum:
            df = pd.read_csv(path2 + "diag_perturb/" + datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean = geom_mean(time, 1.0)

            if "eps-3" in datum:
                for diag_perturb in diag_perturb_lst:
                    if diag_perturb in datum and diag_perturb + "0" not in datum:
                        speedup1["diag_perturb" + diag_perturb + "eps-3"] = ref_gmean["eps-3"] / gmean
            
            elif "eps-6" in datum:
                for diag_perturb in diag_perturb_lst:
                    if diag_perturb in datum and diag_perturb + "0" not in datum:
                        speedup1["diag_perturb" + diag_perturb + "eps-6"] = ref_gmean["eps-6"] / gmean
    
    speedup2 = {"tested-3": 1.0, "tested-6": 1.0}
    labels2 = ["optimized nasoq-custom", "max_iter 0", "max_iter 1", "max_iter 2", "max_iter 3", \
        "max_iter 4", "max_iter 5", "max_iter 10", "max_iter 20"]
    
    for datum in data2:
        if "nasoq-custom" in datum:
            df = pd.read_csv(path2 + "max_iter/" + datum)
            status = df["Status"].values
            time = df["Time (s)"].values

            i = 0
            while i < len(status):
                if status[i] != 1:
                    time[i] = MAX_TIMING
                i += 1
            gmean = geom_mean(time, 1.0)

            if "eps-3" in datum:
                for max_iter in max_iter_lst:
                    if max_iter in datum and max_iter + "0" not in datum:
                        speedup2["max_iter" + max_iter + "eps-3"] = ref_gmean["eps-3"] / gmean
            
            elif "eps-6" in datum:
                for max_iter in max_iter_lst:
                    if max_iter in datum and max_iter + "0" not in datum:
                        speedup2["max_iter" + max_iter + "eps-6"] = ref_gmean["eps-6"] / gmean
    
    eps3_speedup1 = [speedup1["tested-3"], speedup1["diag_perturb-12eps-3"], speedup1["diag_perturb-11eps-3"], speedup1["diag_perturb-10eps-3"], \
        speedup1["diag_perturb-9eps-3"], speedup1["diag_perturb-8eps-3"], speedup1["diag_perturb-7eps-3"], speedup1["diag_perturb-6eps-3"]]

    eps6_speedup1 = [speedup1["tested-6"], speedup1["diag_perturb-12eps-6"], speedup1["diag_perturb-11eps-6"], speedup1["diag_perturb-10eps-6"], \
        speedup1["diag_perturb-9eps-6"], speedup1["diag_perturb-8eps-6"], speedup1["diag_perturb-7eps-6"], speedup1["diag_perturb-6eps-6"]]
    
    eps3_speedup2 = [speedup2["tested-3"], speedup2["max_iter0eps-3"], speedup2["max_iter1eps-3"], speedup2["max_iter2eps-3"], speedup2["max_iter3eps-3"], \
        speedup2["max_iter4eps-3"], speedup2["max_iter5eps-3"], speedup2["max_iter10eps-3"], speedup2["max_iter20eps-3"]]
    
    eps6_speedup2 = [speedup2["tested-6"], speedup2["max_iter0eps-6"], speedup2["max_iter1eps-6"], speedup2["max_iter2eps-6"], speedup2["max_iter3eps-6"], \
        speedup2["max_iter4eps-6"], speedup2["max_iter5eps-6"], speedup2["max_iter10eps-6"], speedup2["max_iter20eps-6"]]
    
    plt.figure()
    x = np.arange(len(labels1))
    width = 0.35

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, eps3_speedup1, width, label='eps-3')
    rects2 = ax.bar(x + width/2, eps6_speedup1, width, label='eps-6')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('speedup')
    ax.set_title('speedup of nasoq-custom (optimized vs non-optimized) for diag_perturb')
    ax.set_xticks(x)
    ax.set_xticklabels(labels1)
    ax.legend()
    for tick in ax.get_xticklabels():
        tick.set_rotation(20)

    plt.savefig(plots_path + "different_diag_perturb.png")

    plt.figure()
    x = np.arange(len(labels2))
    width = 0.35

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, eps3_speedup2, width, label='eps-3')
    rects2 = ax.bar(x + width/2, eps6_speedup2, width, label='eps-6')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('speedup')
    ax.set_title('speedup of nasoq-custom (optimized vs non-optimized) for max_iter')
    ax.set_xticks(x)
    ax.set_xticklabels(labels2)
    ax.legend()
    for tick in ax.get_xticklabels():
        tick.set_rotation(20)

    plt.savefig(plots_path + "different_max_iter.png")

def computer_performance_profile1():
    """
    plot the performance profile for NASOQ-CUSTOM optimized version
    and other NASOQ variates with the same settings
    """
    data = os.listdir(path1)
    solvers = ["nasoq-fixed", "nasoq-tuned", "nasoq-custom"]
    
    def compute_helper(eps):
        t = {}
        status = {}
        file_name = ""
        for solver in solvers:
            for datum in data:
                if solver in datum and eps in datum:
                    file_name = datum
                    break
            df = pd.read_csv(path1 + file_name)
            n_problems = len(df)
            t[solver] = df["Time (s)"].values
            status[solver] = df["Status"].values

            for idx in range(n_problems):
                if status[solver][idx] not in SOLUTION_PRESENT:
                    t[solver][idx] = MAX_TIMING
            
        r = {}  # Dictionary of relative times for each solver/problem
        for s in solvers:
            r[s] = np.zeros(n_problems)

        # Iterate over all problems to find best timing between solvers
        for p in range(n_problems):
            # Get minimum time
            min_time = np.min([t[s][p] for s in solvers])

            # Normalize t for minimum time
            if min_time == 0:
                print("division by zero")
            for s in solvers:
                r[s][p] = t[s][p] / min_time
        
        # Compute curve for all solvers
        n_tau = 1000
        tau_vec = np.logspace(0, 4, n_tau)
        rho = {'tau': tau_vec}  # Dictionary of all the curves

        for s in solvers:
            rho[s] = np.zeros(n_tau)
            for tau_idx in range(n_tau):
                count_problems = 0  # Count number of problems with t[p, s] <= tau
                for p in range(n_problems):
                    if r[s][p] <= tau_vec[tau_idx]:
                        count_problems += 1
                rho[s][tau_idx] = count_problems / n_problems

        # Store final pandas dataframe
        df_performance_profiles = pd.DataFrame(rho)
        df_performance_profiles.to_csv(plots_path + "diff_variates_performance_profile_eps{}.csv".format(eps), index=False)
    
    def plot_performance_profiles(tol=''):
        """
        Plot performance profiles in matplotlib for specified problems and solvers
        """
        df = pd.read_csv(plots_path + "diff_variates_performance_profile_eps{}.csv".format(tol))
        plt.figure()
        for solver in solvers:
            plt.plot(df["tau"], df[solver], label=solver)
        plt.xlim(1., 10000.)
        plt.ylim(0., 1.)
        plt.xlabel(r'Performance ratio $\tau$')
        plt.ylabel('Ratio of problems solved (tol = %s)' %tol)
        plt.xscale('log')
        plt.legend()
        plt.title("diff_variates_performance_profile_eps{}".format(tol))
        results_file = plots_path + "diff_variates_performance_profile_eps{}.png".format(tol)
        plt.savefig(results_file, dpi=100)
        plt.savefig(results_file, dpi=300)
    
    for eps in ["-3", "-6"]:
        compute_helper(eps)
        plot_performance_profiles(eps)

def computer_performance_profile2():
    """
    plot the performance profile for NASOQ-CUSTOM optimized version
    and non-optimized version
    """
    data1 = os.listdir(path1)
    data2 = os.listdir(path2)

    solvers1 = ["optimized nasoq-custom", "diag_perturb-12", "diag_perturb-11", "diag_perturb-10", "diag_perturb-9", \
        "diag_perturb-8", "diag_perturb-7", "diag_perturb-6"]
    solvers2 = ["optimized nasoq-custom", "max_iter0", "max_iter1", "max_iter2", "max_iter3", \
        "max_iter4", "max_iter5", "max_iter10", "max_iter20"]
    
    def compute_helper(eps, solvers):
        t = {}
        status = {}
        file_name = ""

        param = ""
        if solvers[1] == solvers1[1]:
            param = "diag_perturb/"
        else:
            param = "max_iter/"
        data2 = os.listdir(path2 + param)

        df = pd.read_csv(path1 + "nasoq-custom-eps{}.csv".format(eps))
        n_problems = len(df)
        t["optimized nasoq-custom"] = df["Time (s)"].values
        status["optimized nasoq-custom"] = df["Status"].values
        for idx in range(n_problems):
            if status["optimized nasoq-custom"][idx] not in SOLUTION_PRESENT:
                t[solver][idx] = MAX_TIMING
        
        for solver in solvers[1:]:
            for datum in data2:
                if solver in datum and "eps" + eps in datum:
                    file_name = datum
                    break

            df = pd.read_csv(path2 + param + file_name)
            n_problems = len(df)
            t[solver] = df["Time (s)"].values
            status[solver] = df["Status"].values

            for idx in range(n_problems):
                if status[solver][idx] not in SOLUTION_PRESENT:
                    t[solver][idx] = MAX_TIMING
            
        r = {}  # Dictionary of relative times for each solver/problem
        for s in solvers:
            r[s] = np.zeros(n_problems)

        # Iterate over all problems to find best timing between solvers
        for p in range(n_problems):
            # Get minimum time
            min_time = np.min([t[s][p] for s in solvers])

            # Normalize t for minimum time
            if min_time == 0:
                print("division by zero")
            for s in solvers:
                r[s][p] = t[s][p] / min_time
        
        # Compute curve for all solvers
        n_tau = 1000
        tau_vec = np.logspace(0, 4, n_tau)
        rho = {'tau': tau_vec}  # Dictionary of all the curves

        for s in solvers:
            rho[s] = np.zeros(n_tau)
            for tau_idx in range(n_tau):
                count_problems = 0  # Count number of problems with t[p, s] <= tau
                for p in range(n_problems):
                    if r[s][p] <= tau_vec[tau_idx]:
                        count_problems += 1
                rho[s][tau_idx] = count_problems / n_problems

        # Store final pandas dataframe
        df_performance_profiles = pd.DataFrame(rho)
        df_performance_profiles.to_csv(plots_path + "{}_performance_profile_eps{}.csv".format(param[:-1], eps), index=False)
    
    def plot_performance_profiles(solvers, tol=''):
        """
        Plot performance profiles in matplotlib for specified problems and solvers
        """
        param = ""
        if solvers[1] == solvers1[1]:
            param = "diag_perturb"
        else:
            param = "max_iter"
    
        df = pd.read_csv(plots_path + "{}_performance_profile_eps{}.csv".format(param, tol))
        plt.figure()
        for solver in solvers:
            plt.plot(df["tau"], df[solver], label=solver)
        plt.xlim(1., 10000.)
        plt.ylim(0., 1.)
        plt.xlabel(r'Performance ratio $\tau$')
        plt.ylabel('Ratio of problems solved (tol = %s)' %tol)
        plt.xscale('log')
        plt.legend()
        plt.title("{}_performance_profile_eps{}".format(param, tol))
        results_file = plots_path + "{}_performance_profile_eps{}.png".format(param, tol)
        plt.savefig(results_file, dpi=100)
        plt.savefig(results_file, dpi=300)
    
    for eps in ["-3", "-6"]:
        for solvers in [solvers1, solvers2]:
            compute_helper(eps, solvers)
            plot_performance_profiles(solvers, eps)

if __name__ == "__main__":
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)
    compute_speedup1()
    compute_speedup2()
    computer_performance_profile1()
    computer_performance_profile2()
