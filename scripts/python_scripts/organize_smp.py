import pandas as pd
import numpy as np
import os
import csv

path = "SMP_Repository_analysis_tables/"
suffix = "_QPs_analysis.csv"

header = ["Problem Name", "Number of problems", "Number of variables", "Number of constraints", "nnz", "Source", "Application"]

def main(path):
    probs = os.listdir(path)
    probs.remove("merged_csvs.csv")
    fp = open("smp_info.csv", "w")
    writer = csv.DictWriter(fp, fieldnames=header)
    writer.writeheader()
    # print(probs)

    for prob in probs:
        if prob.endswith(suffix):
            prob_name = prob[: prob.find(suffix)]

            df = pd.read_csv(path + prob)
            num_probs = len(df)
            # print(df.head())
            # print(df.columns.values.tolist())
            # print(df["Number of variables"])

            num_vars = np.unique(df["Number of variables"].values)
            num_vars = (np.min(num_vars), np.max(num_vars))

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
            num_constrs = np.unique(eq + ineq)
            num_constrs = (np.min(num_constrs), np.max(num_constrs))

            nnz = np.unique(df["QP nonzeros"].values)
            nnz = (np.min(nnz), np.max(nnz))

            source = list(set(df["Problem source"].values))[0]

            application = list(set(df["Application"].values))[0]

            row = {header[0]: prob_name, header[1]: num_probs, header[2]: num_vars, header[3]: num_constrs, \
                header[4]: nnz, header[5]: source, header[6]: application}
            writer.writerow(row)
    fp.close()

if __name__ == "__main__":
    main(path)