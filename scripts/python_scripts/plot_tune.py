import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# dir = "test_smp_tune_plots_data/"
# des_dir = "test_smp_tune_plots/"

dir = "SMP_Repository_tune_plots_data/"
des_dir = "SMP_Repository_tune_plots/"

param_lst = ["max_iter", "stop_tol", "diag_perturb"]

# all values of parameters to be tested against
max_iter_lst = [0, 1, 2, 3, 4, 5, 10, 20]
stop_tol_lst = [-13, -15, -16, -17]
diag_perturb_lst = [-6, -7, -8, -9, -10, -11, -12]

eps_lst = [-3, -6]

def plot_failure_rate(eps=-3, param="max_iter"):
    """
    plot failure rates vs tuned parameters
    """
    second_dir = "failure_rates/"

    # load the data from the given directory for failure rate
    data_lst = os.listdir(dir + second_dir)
    data_loaded = {}
    for data in data_lst:
        if "eps{}".format(eps) in data and param in data:
            if "nasoq-custom" in data:
                data_loaded["nasoq-custom"] = pd.read_csv(dir + second_dir + data)
            elif "nasoq-fixed" in data:
                data_loaded["nasoq-fixed"] = pd.read_csv(dir + second_dir + data)
            elif "nasoq-tuned" in data:
                data_loaded["nasoq-tuned"] = pd.read_csv(dir + second_dir + data)
    
    # load the value list for the given tuned parameter
    if param == "max_iter":
        lst = max_iter_lst
    elif param == "stop_tol":
        lst = stop_tol_lst
    elif param == "diag_perturb":
        lst = diag_perturb_lst
    
    # load the failure rate data by ascending order of parameter values from the computed dataset
    custom_failure_rate = data_loaded["nasoq-custom"]["failure rate * 100"].values
    custom_paramvals = data_loaded["nasoq-custom"][param].values
    custom_failure_rate = custom_failure_rate[np.argsort(custom_paramvals)]

    fixed_failure_rate = data_loaded["nasoq-fixed"]["failure rate * 100"].values
    fixed_paramvals = data_loaded["nasoq-fixed"][param].values
    fixed_failure_rate = fixed_failure_rate[np.argsort(fixed_paramvals)]

    tuned_failure_rate = data_loaded["nasoq-tuned"]["failure rate * 100"].values
    tuned_paramvals = data_loaded["nasoq-tuned"][param].values
    tuned_failure_rate = tuned_failure_rate[np.argsort(tuned_paramvals)]
    
    # plot the graph
    plt.figure()
    plt.title("failure rate vs {} (eps = {})".format(param, eps))
    plt.plot(sorted(lst), custom_failure_rate, label="nasoq-custom")
    plt.plot(sorted(lst), fixed_failure_rate, label="nasoq-fixed")
    plt.plot(sorted(lst), tuned_failure_rate, label="nasoq-tuned")
    plt.xlabel(param)
    plt.ylabel("failure_rate * 100")
    plt.legend()
    plt.savefig(des_dir + second_dir + "failure_rate_" + "{}_eps{}.png".format(param, eps))

def plot_speedup(eps=-3, param="max_iter"):
    """
    plot speedup vs tuned parameters
    """
    second_dir = "speedup/"

    # load the data from the given directory for speedup
    data_lst = os.listdir(dir + second_dir)
    data_loaded = {}
    for data in data_lst:
        if "eps{}".format(eps) in data and param in data:
            if "nasoq-custom" in data:
                data_loaded["nasoq-custom"] = pd.read_csv(dir + second_dir + data)
            elif "nasoq-fixed" in data:
                data_loaded["nasoq-fixed"] = pd.read_csv(dir + second_dir + data)
            elif "nasoq-tuned" in data:
                data_loaded["nasoq-tuned"] = pd.read_csv(dir + second_dir + data)
    
    # load the value list for the given tuned parameter
    if param == "max_iter":
        lst = max_iter_lst
    elif param == "stop_tol":
        lst = stop_tol_lst
    elif param == "diag_perturb":
        lst = diag_perturb_lst
    
    # load the speedup data by ascending order of parameter values from the computed dataset
    custom_speedup_ignore_failure = data_loaded["nasoq-custom"]["speedup (ignore failure)"].values
    custom_speedup_with_failure = data_loaded["nasoq-custom"]["speedup (with failure)"].values
    custom_paramvals = data_loaded["nasoq-custom"][param].values
    custom_speedup_ignore_failure = custom_speedup_ignore_failure[np.argsort(custom_paramvals)]
    custom_speedup_with_failure = custom_speedup_with_failure[np.argsort(custom_paramvals)]

    fixed_speedup_ignore_failure = data_loaded["nasoq-fixed"]["speedup (ignore failure)"].values
    fixed_speedup_with_failure = data_loaded["nasoq-fixed"]["speedup (with failure)"].values
    fixed_paramvals = data_loaded["nasoq-fixed"][param].values
    fixed_speedup_ignore_failure = fixed_speedup_ignore_failure[np.argsort(fixed_paramvals)]
    fixed_speedup_with_failure = fixed_speedup_with_failure[np.argsort(fixed_paramvals)]

    tuned_speedup_ignore_failure = data_loaded["nasoq-tuned"]["speedup (ignore failure)"].values
    tuned_speedup_with_failure = data_loaded["nasoq-tuned"]["speedup (with failure)"].values
    tuned_paramvals = data_loaded["nasoq-tuned"][param].values
    tuned_speedup_ignore_failure = tuned_speedup_ignore_failure[np.argsort(tuned_paramvals)]
    tuned_speedup_with_failure = tuned_speedup_with_failure[np.argsort(tuned_paramvals)]
    
    # plot the graph
    plt.figure()
    fig, axs = plt.subplots(2)
    # plt.title("failure rate vs {} (eps = {})".format(param, eps))
    axs[0].set_title('speedup (ignore failure) vs {} (eps = {})'.format(param, eps))
    axs[0].plot(sorted(lst), custom_speedup_ignore_failure, label="nasoq-custom")
    axs[0].plot(sorted(lst), fixed_speedup_ignore_failure, label="nasoq-fixed")
    axs[0].plot(sorted(lst), tuned_speedup_ignore_failure, label="nasoq-tuned")
    axs[0].set(xlabel=param, ylabel="speedup (ignore failure)")

    axs[1].set_title('speedup (with failure) vs {} (eps = {})'.format(param, eps))
    axs[1].plot(sorted(lst), custom_speedup_with_failure, label="nasoq-custom")
    axs[1].plot(sorted(lst), fixed_speedup_with_failure, label="nasoq-fixed")
    axs[1].plot(sorted(lst), tuned_speedup_with_failure, label="nasoq-tuned")
    axs[1].set(xlabel=param, ylabel="speedup (with failure)")

    axs[0].legend(); axs[1].legend()
    plt.tight_layout()
    plt.savefig(des_dir + second_dir + "speedup_" + "{}_eps{}.png".format(param, eps))

def main():
    if not os.path.exists(des_dir):
        os.makedirs(des_dir)
        os.makedirs(des_dir + "failure_rates/")
        os.makedirs(des_dir + "speedup/")

    for eps in eps_lst:
        for param in param_lst:
            plot_failure_rate(eps, param)
            plot_speedup(eps, param)

    # for param in param_lst:
    #     plot_failure_rate(-3, param)
    #     plot_speedup(-3, param)

if __name__ == "__main__":
    main()