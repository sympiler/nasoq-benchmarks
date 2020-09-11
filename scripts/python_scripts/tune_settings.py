from subprocess import call
import sys
import os

def tune(eps, diag_perturb, max_iter, stop_tol, build_folder, dataset):
    """
    draw the plotting based on the tuned settings
    """

    path = "{}_setting_plots/eps{}_max_iter{}_stop_tol{}_diag_perturb{}".format(dataset, eps, max_iter, stop_tol, diag_perturb)

    if os.path.exists(path) and len(os.listdir(path)):
        return

    if not os.path.exists(path):
        os.makedirs(path)

    os.system("rm -f *.png")

    call(["echo", "Running NASOQ-Fixed ..."])
    os.system("bash scripts/NASOQ_bench.sh {}/nasoq/NASOQ-BIN".format(build_folder) + " {} ".format(dataset) + " " + eps + \
        " -p {} -r {} -t {}".format(diag_perturb, max_iter, stop_tol) + " > " + "logs/nasoq-fixed-e{}.csv".format(eps))

    call(["echo", "Running NASOQ-Tuned ..."])
    call(["bash","scripts/NASOQ_bench.sh", "{}/nasoq/NASOQ-BIN".format(build_folder), "{}".format(dataset), eps, \
        "-p {} -r {} -t {} -v tuned".format(diag_perturb, max_iter, stop_tol), ">", "logs/nasoq-tuned-e{}.csv".format(eps)])

    call(["echo", "Running customized NASOQ ..."])
    call(["bash", "scripts/NASOQ_bench.sh", "{}/nasoq/NASOQ-BIN".format(build_folder), "{}".format(dataset), eps, \
        "-p {} -r {} -t {} -v predet".format(diag_perturb, max_iter, stop_tol), ">", "logs/nasoq-custom-e{}.csv".format(eps)])

    os.chdir("scripts/python_scripts")
    call(["python", "graph_generator.py", "-d",  "../../logs/", "-s", eps])
    os.system("rm -f *.txt")
    # os.system("rm -f ../../logs/*.csv")
    os.system("mv *.png ../../{}_setting_plots/eps{}_max_iter{}_stop_tol{}_diag_perturb{}".format(dataset, eps, max_iter, stop_tol, diag_perturb))
    os.chdir("../..")

def main():
    """
    experiment on the influence of different settings on the behaviour
    of the NASOQ solver
    """
    if len(sys.argv) != 3:
        call(["echo", "<dataset>", "<build folder>"])
        sys.exit(1)
    
    dataset, build_folder = sys.argv[1], sys.argv[2]

    if not os.path.exists(dataset + "_setting_plots"):
        # call(["rm", "-rf", "setting_plots"])
        os.makedirs(dataset + "_setting_plots")

    suggested_eps = list(map(str, [-3, -6]))
    suggested_max_iter = list(map(str, [0, 5, 10]))
    suggested_stop_tol = list(map(str, [-13, -15, -17]))
    suggested_diag_perturb = list(map(str, [-6, -9, -12]))

    # os.chdir("scripts/python_scripts")

    for eps in suggested_eps:

        max_iter = suggested_max_iter[0]
        stop_tol = suggested_stop_tol[0]
        for diag_perturb in suggested_diag_perturb:
            tune(eps, diag_perturb, max_iter, stop_tol, build_folder, dataset)
        
        max_iter = suggested_max_iter[0]
        diag_perturb = suggested_diag_perturb[0]
        for stop_tol in suggested_stop_tol:
            tune(eps, diag_perturb, max_iter, stop_tol, build_folder, dataset)
        
        stop_tol = suggested_stop_tol[0]
        diag_perturb = suggested_diag_perturb[0]
        for max_iter in suggested_max_iter:
            tune(eps, diag_perturb, max_iter, stop_tol, build_folder, dataset)

if __name__ == "__main__":
    # for i in range(10):
    #     call(["echo", str(i)])
    # call(["ls"])
    main()