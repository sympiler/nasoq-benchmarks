import pandas as pd
import numpy as np
import csv
import os

# the name of the folder and subfolder for output
# dir = "test_smp_tune_csvs/"
# des_dir = "test_smp_tune_plots_data/"

dir = "SMP_Repository_tune_csvs/"
des_dir = "SMP_Repository_tune_plots_data/"

# the max timing for success
MAX_VAL = 1600

# run_time for failed cases
MAX_TIMING = 1e8

# the status of successful cases
OPTIMAL = 1
SOLUTION_PRESENT = [OPTIMAL]

# all tuned parameters
mode_lst = ["max_iter", "stop_tol", "diag_perturb"]

# all values of parameters to be tested against
max_iter_lst = [0, 1, 2, 3, 4, 5, 10, 20]
stop_tol_lst = [-13, -15, -16, -17]
diag_perturb_lst = [-6, -7, -8, -9, -10, -11, -12]

# eps tested
eps_lst = [-3, -6]
eps_realval_lst = [1e-3, 1e-6]

# the header of the output csv for failure rate and speedup ("mode" is just a placeholder)
header_failure_rate = ["tool", "eps", "failure rate * 100", "mode"]
header_speedup = ["tool", "eps", "gmean (ignore failure)", "speedup (ignore failure)", \
    "gmean (with failure)", "speedup (with failure)", "mode"]

def compute_failure_rates(mode="max_iter"):
    """
    the function to compute the failure rate of the tool for
    all tested values of the given parameter under all eps
    and store outputs in csvs
    """

    # the subfolder to store the output of this function
    second_des_dir = "failure_rates/"

    def write_csv(tool, eps, writer, csv_name, mode):
        """
        write the output of failure rate to csv for the tuned
        parameter input
        """
        # read out the data (performance of tool at the value
        # of the given output csv of NASOQ-BIN)
        df = pd.read_csv(dir + mode + "/" + csv_name)
        # the number of problems
        n_problems = len(df)
        
        # find the number of failed cases by the number of problems solved by
        # the solver under such a condition do not have the optimal status
        failed_statuses = np.logical_and([df['Status'].values != s
                                    for s in SOLUTION_PRESENT],True)
        n_failed_problems = np.sum(failed_statuses)
        # the failure rate for the tool at the tuned parameters and eps
        failure_rate = n_failed_problems / n_problems

        # choose the list of tested values by the parameter tuned
        if mode == "max_iter":
            lst = max_iter_lst
        elif mode == "diag_perturb":
            lst = diag_perturb_lst
        elif mode == "stop_tol":
            lst = stop_tol_lst

        # write the failure rate for the given tuned parameter under given eps into corresponding csv
        i = 0
        while i < len(lst):
            if mode + str(lst[i]) in csv_name and \
                mode + str(lst[i]) + "0" not in csv_name:
                row = {header_failure_rate[0]: tool, header_failure_rate[1]: eps, \
                    header_failure_rate[2]: failure_rate * 100, mode: lst[i]}
                writer.writerow(row)
                break
            i += 1

    # generate the csv files to write (for the given tested parameter)
    # the csv has the failure rate over all tested values of the tested parameter
    # and also all the values of the tested parameter
    fout1 = open(des_dir + second_des_dir + "nasoq-fixed_eps-3_{}_failrate.csv".format(mode), "w")
    fout2 = open(des_dir + second_des_dir + "nasoq-fixed_eps-6_{}_failrate.csv".format(mode), "w")
    fout3 = open(des_dir + second_des_dir + "nasoq-tuned_eps-3_{}_failrate.csv".format(mode), "w")
    fout4 = open(des_dir + second_des_dir + "nasoq-tuned_eps-6_{}_failrate.csv".format(mode), "w")
    fout5 = open(des_dir + second_des_dir + "nasoq-custom_eps-3_{}_failrate.csv".format(mode), "w")
    fout6 = open(des_dir + second_des_dir + "nasoq-custom_eps-6_{}_failrate.csv".format(mode), "w")

    header_failure_rate_mode = header_failure_rate.copy()
    header_failure_rate_mode[-1] = mode

    writer1 = csv.DictWriter(fout1, fieldnames=header_failure_rate_mode)
    writer1.writeheader()
    writer2 = csv.DictWriter(fout2, fieldnames=header_failure_rate_mode)
    writer2.writeheader()
    writer3 = csv.DictWriter(fout3, fieldnames=header_failure_rate_mode)
    writer3.writeheader()
    writer4 = csv.DictWriter(fout4, fieldnames=header_failure_rate_mode)
    writer4.writeheader()
    writer5 = csv.DictWriter(fout5, fieldnames=header_failure_rate_mode)
    writer5.writeheader()
    writer6 = csv.DictWriter(fout6, fieldnames=header_failure_rate_mode)
    writer6.writeheader()

    # compute the failure rate for all tools under all cases
    for csv_name in os.listdir(dir + mode):
        if "eps-3" in csv_name:
            if "nasoq-fixed" in csv_name:
                write_csv("nasoq-fixed", 1e-3, writer1, csv_name, mode)
            elif "nasoq-tuned" in csv_name:
                write_csv("nasoq-tuned", 1e-3, writer3, csv_name, mode)
            elif "nasoq-custom" in csv_name:
                write_csv("nasoq-custom", 1e-3, writer5, csv_name, mode)
        elif "eps-6" in csv_name:
            if "nasoq-fixed" in csv_name:
                write_csv("nasoq-fixed", 1e-6, writer2, csv_name, mode)
            elif "nasoq-tuned" in csv_name:
                write_csv("nasoq-tuned", 1e-6, writer4, csv_name, mode)
            elif "nasoq-custom" in csv_name:
                write_csv("nasoq-custom", 1e-6, writer6, csv_name, mode)
                
    fout1.close()
    fout2.close()
    fout3.close()
    fout4.close()
    fout5.close()
    fout6.close()

def compute_speedup(mode):
    """
    the function to compute the speedup of the tool for
    all tested values of the given parameter under all eps
    and store outputs in csvs
    """

    # the subfolder to store the output of this function
    second_des_dir = "speedup/"

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
    
    def write_csv(tool, eps, writer, csv_name, mode):
        """
        write the output of speedup to csv for the tuned
        parameter input
        """

        # pick the list of tested values by the parameter tuned
        if mode == "max_iter":
            lst = max_iter_lst
        elif mode == "diag_perturb":
            lst = diag_perturb_lst
        elif mode == "stop_tol":
            lst = stop_tol_lst
        
        # find the value of the tuned parameter for the case to write into csv
        # one value of tuned parameter for the given tool at given eps
        # (the tuned parameter value for the given csv data from archived output of NASOQ-BIN)
        param = lst[0]
        for val in lst:
            if mode + str(val) in csv_name and \
                mode + str(val) + '0' not in csv_name:
                param = val
                break

        # read out the data (performance of tool at the value
        # of the given output csv of NASOQ-BIN)
        df = pd.read_csv(dir + mode + "/" + csv_name)
        n_problems = len(df)
        # compute the geometric mean by running time and status (fail or succeed)
        run_times = df['Time (s)'].values
        statuses = df['Status'].values
        # if a case fails, label its runnint time to MAX_TIMING
        for i in range(n_problems):
            if statuses[i] not in SOLUTION_PRESENT:
                run_times[i] = MAX_TIMING
        # compute the geometric mean
        gmean_ignore_failure = geom_mean_ignore_failure(run_times, 1)
        gmean_with_failure = geom_mean_with_failure(run_times, 1)

        # speedup for the reference method (nasoq-fixed)
        if tool == "nasoq-fixed":
            # if nasoq-fixed failed on all given QPs under the circumstance, speedup = 0
            if gmean_ignore_failure == 0.0: speedup_ignore_failure = 0
            else: speedup_ignore_failure = 1

            speedup_with_failure = 1
        # compute speedup for other tools
        # read out the gmean of nasoq-fixed at first since nasoq-fixed is the reference
        # and its gmean will be used to compute speedup
        else:
            if eps == 1e-3:
                nf_data = pd.read_csv(des_dir + second_des_dir + "nasoq-fixed_eps-3_{}_speedup.csv".format(mode))
            elif eps == 1e-6:
                nf_data = pd.read_csv(des_dir + second_des_dir + "nasoq-fixed_eps-6_{}_speedup.csv".format(mode))

            # compute the speedup for other methods
            nf_param = nf_data[mode].values
            nf_gmean_index = np.where(nf_param == param)[0][0]
            nf_gmean_ignore_failure = nf_data["gmean (ignore failure)"].values[nf_gmean_index]
            nf_gmean_with_failure = nf_data["gmean (with failure)"].values[nf_gmean_index]

            speedup_with_failure = nf_gmean_with_failure / gmean_with_failure

            # if the solver fails on all cases, let its speed be infinity
            # since speedup_of_this solver = gmean_ref / gmean_this_solver
            if gmean_ignore_failure == 0.0:
                speedup_ignore_failure = np.inf
            # if this solver succeeds on some cases but nasoq-fixed fails all,
            # let the speedup be 0 (not be able to use nasoq-fixed as reference
            # for this tuned parameter value)
            elif nf_gmean_ignore_failure == 0.0:
                speedup_ignore_failure = 0
            # normal case for speedup computation (a ratio between gmean_ref and gmean_this_tool)
            else:
                speedup_ignore_failure = nf_gmean_ignore_failure / gmean_ignore_failure

         # write the geometric mean and speedup for the given tuned parameter
         # under given eps into corresponding csv
        row = {header_speedup[0]: tool, header_speedup[1]: eps, header_speedup[2]: gmean_ignore_failure, \
            header_speedup[3]: speedup_ignore_failure, header_speedup[4]: gmean_with_failure, \
                header_speedup[5]: speedup_with_failure, mode: param}
        writer.writerow(row)

    # generate the csv files to write (for the given tested parameter)
    # the csv has the speedup over all tested values of the tested parameter
    # and also all the values of the tested parameter
    fout1 = open(des_dir + second_des_dir + "nasoq-fixed_eps-3_{}_speedup.csv".format(mode), "w")
    fout2 = open(des_dir + second_des_dir + "nasoq-fixed_eps-6_{}_speedup.csv".format(mode), "w")
    fout3 = open(des_dir + second_des_dir + "nasoq-tuned_eps-3_{}_speedup.csv".format(mode), "w")
    fout4 = open(des_dir + second_des_dir + "nasoq-tuned_eps-6_{}_speedup.csv".format(mode), "w")
    fout5 = open(des_dir + second_des_dir + "nasoq-custom_eps-3_{}_speedup.csv".format(mode), "w")
    fout6 = open(des_dir + second_des_dir + "nasoq-custom_eps-6_{}_speedup.csv".format(mode), "w")

    header_speedup_mode = header_speedup.copy()
    header_speedup_mode[-1] = mode

    writer1 = csv.DictWriter(fout1, fieldnames=header_speedup_mode)
    writer1.writeheader()
    writer2 = csv.DictWriter(fout2, fieldnames=header_speedup_mode)
    writer2.writeheader()
    writer3 = csv.DictWriter(fout3, fieldnames=header_speedup_mode)
    writer3.writeheader()
    writer4 = csv.DictWriter(fout4, fieldnames=header_speedup_mode)
    writer4.writeheader()
    writer5 = csv.DictWriter(fout5, fieldnames=header_speedup_mode)
    writer5.writeheader()
    writer6 = csv.DictWriter(fout6, fieldnames=header_speedup_mode)
    writer6.writeheader()

    # use nasoq-fixed as reference, so compute the gmean and speedup
    # about nasoq-fixed under all cases at first
    csv_names = os.listdir(dir + mode)
    nofixed_csv_names = csv_names.copy()
    for csv_name in csv_names:
        if "eps-3" in csv_name and "nasoq-fixed" in csv_name:
            write_csv("nasoq-fixed", 1e-3, writer1, csv_name, mode)
            nofixed_csv_names.remove(csv_name)
        elif "eps-6" in csv_name and "nasoq-fixed" in csv_name:
            write_csv("nasoq-fixed", 1e-6, writer2, csv_name, mode)
            nofixed_csv_names.remove(csv_name)
    
    fout1.close()
    fout2.close()
    
    # compute the speedup for other tools under all cases, with reference to nasoq-fixed
    for csv_name in nofixed_csv_names:
        if "eps-3" in csv_name:
            if "nasoq-tuned" in csv_name:
                write_csv("nasoq-tuned", 1e-3, writer3, csv_name, mode)
            elif "nasoq-custom" in csv_name:
                write_csv("nasoq-custom", 1e-3, writer5, csv_name, mode)
        elif "eps-6" in csv_name:
            if "nasoq-tuned" in csv_name:
                write_csv("nasoq-tuned", 1e-6, writer4, csv_name, mode)
            elif "nasoq-custom" in csv_name:
                write_csv("nasoq-custom", 1e-6, writer6, csv_name, mode)
    
    fout3.close()
    fout4.close()
    fout5.close()
    fout6.close()

def main():
    """
    generate the analysis data for benchmarks of each tool
    under all given required accuracies (eps) for all tuned
    parameters and their tested values listed above
    """

    # create folder and subfolder for the analysis result
    if not os.path.exists(des_dir):
        os.makedirs(des_dir)
        os.makedirs(des_dir + "failure_rates/")
        os.makedirs(des_dir + "speedup/")
    
    # make outputs of failure rate
    for mode in mode_lst:
        compute_failure_rates(mode)
    
    # make outputs of speedup
    for mode in mode_lst:
        compute_speedup(mode)

if __name__ == "__main__":
    main()