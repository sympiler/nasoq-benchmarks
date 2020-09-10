from subprocess import call
import sys
import os

def main():
    if len(sys.argv) != 3:
        call(["echo", "<dataset>", "<build folder>"])
        sys.exit(1)
    
    if os.path.exists("../../setting_plots"):
        call(["rm", "-rf", "../../setting_plots"])
    os.makedirs("../../setting_plots")
    
    dataset, build_folder = sys.argv[1], sys.argv[2]

    suggested_eps = map(str, [-3, -6])
    suggested_max_iter = map(str, [0, 1, 2, 3, 4, 5, 10, 20])
    suggested_stop_tol = map(str, [-13, -15, -16, -17])
    suggested_diag_perturb = [-6, -7, -8, -9, -10, -11, -12]

    for eps in suggested_eps:
        for max_iter in suggested_max_iter:
            for stop_tol in suggested_stop_tol:
                for diag_perturb in suggested_diag_perturb:
                    os.makedirs("../../setting_plots/eps{}_max_iter{}_stop_tol{}_diag_perturb{}".format(eps, max_iter, stop_tol, diag_perturb))

                    call(["echo", "Running NASOQ-Fixed ..."])
                    call(["bash", "../NASOQ_bench.sh", "../../{}/nasoq/NASOQ-BIN".format(build_folder), dataset, eps, \
                        "-p {} -r {} -t {}".format(diag_perturb, max_iter, stop_tol), ">", "../../logs/nasoq-fixed-e{}.csv".format(eps)])

                    call(["echo", "Running NASOQ-Tuned ..."])
                    call(["bash","../NASOQ_bench.sh", "../../{}/nasoq/NASOQ-BIN".format(build_folder), dataset, eps, \
                        "-p {} -r {} -t {} -v tuned".format(diag_perturb, max_iter, stop_tol), ">", "../../logs/nasoq-tuned-e{}.csv".format(eps)])

                    call(["echo", "Running customized NASOQ ..."])
                    call(["bash", "../NASOQ_bench.sh", "../../{}/nasoq/NASOQ-BIN".format(build_folder), dataset, eps, \
                        "-p {} -r {} -t {} -v predet".format(diag_perturb, max_iter, stop_tol), ">", "../../logs/nasoq-custom-e{}.csv".format(eps)])

                    call(["python", "graph_generator.py", "-d ../../logs/ -s", eps, ">", \
                        "../../setting_plots/eps{}_max_iter{}_stop_tol{}_diag_perturb{}".format(eps, max_iter, stop_tol, diag_perturb)])

if __name__ == "__main__":
    # for i in range(10):
    #     call(["echo", str(i)])
    # call(["ls"])
    main()