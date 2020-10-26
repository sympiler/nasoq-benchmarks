import os
import numpy as np
import pandas as pd

analysis_path = "SMP_Repository_analysis_tables/merged_csvs.csv"
perf_path = "SMP_Repository_log/"

max_iter_labels = {"max_iter0": 0, "max_iter1": 1, "max_iter2": 2, "max_iter3": 3, "max_iter4": 4, "max_iter5": 5, "max_iter10": 10, "max_iter20": 20}
diag_perturb_labels = {"diag_perturb-10": -10, "diag_perturb-9": -9, "diag_perturb-8": -8, "diag_perturb-7": -7, "diag_perturb-6": -6}
eps_labels = {"eps-3": -3, "eps-6": -6}

def compute_metrics(analysis_path):
    """
    compute and process the metrics of each QP
    """
    def compute_numerical_range(cell1, cell2):
        """
        compute numerical range: abs(max / min)
        """
        result = 0.0
        if cell2 != 0.0:
            result = np.abs(cell1 / cell2)
        return np.float64(result)

    # read out the analysis data and process the data
    df = pd.read_csv(analysis_path)
    df["Problem name"] = df["Problem name"].str.lower()

    df.loc[df["Number of inequality constraints"] == "\\", "Number of inequality constraints"] = 0
    df.loc[df["inequality norm (order = infty)"] == "\\", "inequality norm (order = infty)"] = 0
    df.loc[df["inequality nonzeros"] == "\\", "inequality nonzeros"] = 0
    df.loc[df["inequality max"] == "\\", "inequality max"] = 0
    df.loc[df["inequality min"] == "\\", "inequality min"] = 0

    df.loc[df["Number of equality constraints"] == "\\", "Number of equality constraints"] = 0
    df.loc[df["equality norm (order = infty)"] == "\\", "equality norm (order = infty)"] = 0
    df.loc[df["equality nonzeros"] == "\\", "equality nonzeros"] = 0
    df.loc[df["equality max"] == "\\", "equality max"] = 0
    df.loc[df["equality min"] == "\\", "equality min"] = 0

    # compute the number of constraints: equality constr # + inequality constr #
    df["number of constraints"] = np.int32(df["Number of equality constraints"].values) + np.int32(df["Number of inequality constraints"].values)
    df["QP max"] = np.float64(df["QP max"])
    df["inequality max"] = np.float64(df["inequality max"])
    df["equality max"] = np.float64(df["equality max"])
    df["QP min"] = np.float64(df["QP min"])
    df["inequality min"] = np.float64(df["inequality min"])
    df["equality min"] = np.float64(df["equality min"])

    # compute the numerical range
    df["objective numerical range (abs(max / eps-3))"] = np.abs(df["QP max"]) / 1e-3
    df["inequality numerical range (abs(max / eps-3))"] = np.abs(df["inequality max"]) / 1e-3
    df["equality numerical range (abs(max / eps-3))"] = np.abs(df["equality max"]) / 1e-3

    df["objective numerical range (abs(max / eps-6))"] = np.abs(df["QP max"]) / 1e-6
    df["inequality numerical range (abs(max / eps-6))"] = np.abs(df["inequality max"]) / 1e-6
    df["equality numerical range (abs(max / eps-6))"] = np.abs(df["equality max"]) / 1e-6

    df["objective numerical range (abs(max / min))"] = np.abs(df["QP max"] / df["QP min"])
    df["inequality numerical range (abs(max / min))"] = df.apply(lambda x: compute_numerical_range(x['inequality max'], x['inequality min']), axis=1)
    df["equality numerical range (abs(max / min))"] = df.apply(lambda x: compute_numerical_range(x['equality max'], x['equality min']), axis=1)

    return df[["Problem name", "QP max", "QP min", "Number of inequality constraints", "inequality norm (order = infty)", "inequality nonzeros", \
        "inequality max", "inequality min", "Number of equality constraints", "equality norm (order = infty)", \
            "equality nonzeros", "equality max", "equality min", "number of constraints", "objective numerical range (abs(max / eps-3))", \
                "inequality numerical range (abs(max / eps-3))", "equality numerical range (abs(max / eps-3))", \
                    "objective numerical range (abs(max / eps-6))", "inequality numerical range (abs(max / eps-6))", \
                        "equality numerical range (abs(max / eps-6))", "objective numerical range (abs(max / min))", \
                            "inequality numerical range (abs(max / min))", "equality numerical range (abs(max / min))"]]

def get_status(perf_path, eps="eps-3"):
    """
    get the a dataframe of status for each diag_perturb and max_iter
    and select rows with any failure
    """
    # print(eps)
    diag_perturb_data = os.listdir(perf_path + "diag_perturb/")
    max_iter_data = os.listdir(perf_path + "max_iter/")

    diag_perturb_df_dict = {}
    max_iter_df_dict = {}

    # filter the status from diag_perturb dataframes
    for label in diag_perturb_labels.keys():
        for datum in diag_perturb_data:
            if datum.endswith(".csv") and label in datum and label + '0' not in datum and eps in datum:
                # print(datum)
                diag_perturb_df_dict[label] = pd.read_csv(perf_path + "diag_perturb/" + datum)
                diag_perturb_df_dict[label] = diag_perturb_df_dict[label][["Problem Name", "Status"]]
                diag_perturb_df_dict[label].rename(columns={'Problem Name': 'Problem name', 'Status': '{} status'.format(label)}, inplace=True)
                # print(diag_perturb_df_dict[label]['{} status'.format(label)][0])
                diag_perturb_df_dict[label]['{} status'.format(label)] = np.int32(diag_perturb_df_dict[label]['{} status'.format(label)])
    
    # filter the status from max_iter dataframes
    for label in max_iter_labels.keys():
        for datum in max_iter_data:
            if datum.endswith(".csv") and label in datum and label + '0' not in datum and eps in datum:
                # print(datum)
                max_iter_df_dict[label] = pd.read_csv(perf_path + "max_iter/" + datum)
                max_iter_df_dict[label] = max_iter_df_dict[label][["Problem Name", "Status"]]
                max_iter_df_dict[label].rename(columns={'Problem Name': 'Problem name', 'Status': '{} status'.format(label)}, inplace=True)
                max_iter_df_dict[label]['{} status'.format(label)] = np.int32(max_iter_df_dict[label]['{} status'.format(label)])
    
    df = diag_perturb_df_dict["diag_perturb-10"]
    del diag_perturb_df_dict["diag_perturb-10"]

    diag_perturb_data_lst = list(diag_perturb_df_dict.values())
    max_iter_data_lst = list(max_iter_df_dict.values())

    # merge all datasets of each status
    for df_m in diag_perturb_data_lst:
        df = df.merge(df_m, on="Problem name", how="inner")
    for df_m in max_iter_data_lst:
        df = df.merge(df_m, on="Problem name", how="inner")
    
    # drop the row which has no failure
    df.drop(df[(df['max_iter0 status'] == 1) & (df['max_iter1 status'] == 1) & (df['max_iter2 status'] == 1) & (df['max_iter3 status'] == 1) & \
        (df['max_iter4 status'] == 1) & (df['max_iter5 status'] == 1) & (df['max_iter10 status'] == 1) & (df['max_iter20 status'] == 1) & \
            (df['diag_perturb-10 status'] == 1) & (df['diag_perturb-9 status'] == 1) & (df['diag_perturb-8 status'] == 1) & \
                (df['diag_perturb-7 status'] == 1) & (df['diag_perturb-6 status'] == 1)].index, inplace = True)
    
    return df

def df_inner_join(df_analysis, df1_status, df2_status):
    """
    merge dataframes with status and metrics
    """
    df1_status = df1_status.merge(df_analysis, on="Problem name", how="inner")
    df2_status = df2_status.merge(df_analysis, on="Problem name", how="inner")
    df1_status.to_csv("status_eps-3.csv")
    df2_status.to_csv("status_eps-6.csv")
    return df1_status, df2_status

if __name__ == "__main__":
    df_analysis = compute_metrics(analysis_path)
    df_analysis.to_csv("out.csv")
    # df.drop(df[(df['objective numerical range (abs(max / min))'] > 2) | (df["equality norm (order = infty)"] == 0)].index, inplace=True)
    # print(len(df))
    # df.to_csv("out.csv")
    # get_status(perf_path)
    # for label in max_iter_labels.keys():
    #     print(label)
    
    # for label in diag_perturb_labels.keys():
    #     print(label)

    df1_status = get_status(perf_path, "eps-3")
    df2_status = get_status(perf_path, "eps-6")
    df1_status.to_csv("eps-3.csv")
    df1_status.to_csv("eps-6.csv")
    df_inner_join(df_analysis, df1_status, df2_status)
