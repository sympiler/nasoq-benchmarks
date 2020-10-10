import os
import pandas as pd
import numpy as np
import csv

path = "SMP_Repository_tune_csvs/"
dir_analysis_data = "SMP_Repository_analysis_tables/merged_csvs.csv"
final_dir = "classfication_SMP_Repository/"

# all values of parameters to be tested against
max_iter_lst = sorted([0, 1, 2, 3, 4, 5, 10, 20])
diag_perturb_lst = sorted([-6, -7, -8, -9, -10, -11, -12])

header1 = ["problem name", "eps", "norm", "number of constraints", "numerical range"]
header2 = ["problem name", "eps", "norm", "number of constraints", "numerical range", "status for each value"]
# header3 = ["problem name", "eps", "norm", "number of constraints", "numerical range"]

def classify_dataset(eps, path):
    """
    classify the dataset for each eps for max_iter and diag_perturb
    for nasoq-custom based on three standards:
        1. always fail
        2. sometimes fail
        3. never fail
    """
    df_analysis = pd.read_csv(dir_analysis_data)
    # row = df_analysis.loc[df_analysis['Problem name'] == "test05_0"]
    # print(type(row["equality min"].values[0]))

    max_iter_df_array = []
    for val in max_iter_lst:
        max_iter_df_array.append(pd.read_csv(path + "max_iter" + "/" + "nasoq-custom-eps{}_max_iter{}.csv".format(eps, val)))
    
    diag_perturb_df_array = []
    for val in diag_perturb_lst:
        diag_perturb_df_array.append(pd.read_csv(path + "diag_perturb" + "/" + "nasoq-custom-eps{}_diag_perturb{}.csv".format(eps, val)))

    fout1 = open(final_dir + "max_iter_eps{}_always_fail.csv".format(eps), "w")
    fout2 = open(final_dir + "max_iter_eps{}_sometimes_fail.csv".format(eps), "w")
    fout3 = open(final_dir + "max_iter_eps{}_never_fail.csv".format(eps), "w")
    fout4 = open(final_dir + "diag_perturb_eps{}_always_fail.csv".format(eps), "w")
    fout5 = open(final_dir + "diag_perturb_eps{}_sometimes_fail.csv".format(eps), "w")
    fout6 = open(final_dir + "diag_perturb_eps{}_never_fail.csv".format(eps), "w")

    writer1 = csv.DictWriter(fout1, fieldnames=header1); writer1.writeheader()
    writer2 = csv.DictWriter(fout2, fieldnames=header2); writer2.writeheader()
    writer3 = csv.DictWriter(fout3, fieldnames=header1); writer3.writeheader()
    writer4 = csv.DictWriter(fout4, fieldnames=header1); writer4.writeheader()
    writer5 = csv.DictWriter(fout5, fieldnames=header2); writer5.writeheader()
    writer6 = csv.DictWriter(fout6, fieldnames=header1); writer6.writeheader()

    for i in range(len(max_iter_df_array[0])):
        status = []
        info = df_analysis.loc[df_analysis['Problem name'].str.lower() == max_iter_df_array[0]["Problem Name"].values[i]]
        if len(info) > 1:
            info = info.iloc[[0]]
        p_name = info["Problem name"].values[0]
        p_norm = np.float(info["QP norm (order = infty)"])
        p_num_constrs = 0
        if info["Number of equality constraints"].values[0] != "\\":
            p_num_constrs += np.int(info["Number of equality constraints"].values[0])
        if info["Number of inequality constraints"].values[0] != "\\":
            p_num_constrs += np.int(info["Number of inequality constraints"].values[0])
        p_min = np.float(info["QP min"].values[0])
        p_max = np.float(info["QP max"].values[0])
        p_numerical_range = np.abs(p_max / p_min)

        for j in range(len(max_iter_df_array)):
            status.append(np.int(max_iter_df_array[j]["Status"].values[i]))

        convergence = []
        for s in status:
            if s == 1: convergence.append(1)
            else: convergence.append(0)

        if convergence.count(1) == len(convergence):
            row = {header1[0]: p_name, header1[1]: eps, header1[2]: p_norm, header1[3]: p_num_constrs, header1[4]: p_numerical_range}
            writer3.writerow(row)
        elif convergence.count(1) == 0:
            row = {header1[0]: p_name, header1[1]: eps, header1[2]: p_norm, header1[3]: p_num_constrs, header1[4]: p_numerical_range}
            writer1.writerow(row)
        else:
            row = {header2[0]: p_name, header2[1]: eps, header2[2]: p_norm, header2[3]: p_num_constrs, header2[4]: p_numerical_range, header2[5]: status}
            writer2.writerow(row)
        
        status.clear()
        for j in range(len(diag_perturb_df_array)):
            status.append(np.int(diag_perturb_df_array[j]["Status"].values[i]))
        
        convergence.clear()
        for s in status:
            if s == 1: convergence.append(1)
            else: convergence.append(0)
        
        if convergence.count(1) == len(convergence):
            row = {header1[0]: p_name, header1[1]: eps, header1[2]: p_norm, header1[3]: p_num_constrs, header1[4]: p_numerical_range}
            writer6.writerow(row)
        elif convergence.count(1) == 0:
            row = {header1[0]: p_name, header1[1]: eps, header1[2]: p_norm, header1[3]: p_num_constrs, header1[4]: p_numerical_range}
            writer4.writerow(row)
        else:
            row = {header2[0]: p_name, header2[1]: eps, header2[2]: p_norm, header2[3]: p_num_constrs, header2[4]: p_numerical_range, header2[5]: status}
            writer5.writerow(row)

    fout1.close()
    fout2.close()
    fout3.close()
    fout4.close()
    fout5.close()
    fout6.close()

def sort_csv(final_dir):
    """
    sort each row by norm, number of constraints, numerical range
    """
    files = os.listdir(final_dir)
    real_files = []
    for f in files:
        if f.endswith(".csv"):
            real_files.append(f)
    
    for f in real_files:
        df = pd.read_csv(final_dir + f)
        df = df.sort_values(by=["norm", "number of constraints", "numerical range"], ascending=True)
        df.index = np.arange(len(df))
        df.to_csv(final_dir + f)


if __name__ == "__main__":
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)

    classify_dataset(-3, path)
    classify_dataset(-6, path)
    sort_csv(final_dir)