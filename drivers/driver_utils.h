//
// Created by kazem on 6/7/20.
//

#ifndef NASOQ_BENCHMARKS_DRIVER_UTILS_H
#define NASOQ_BENCHMARKS_DRIVER_UTILS_H
namespace nasoq_bench{

 bool parse_args(int argc, char **argv, std::map<std::string, std::string> &qp_args) {
  const char *const short_opts = "i:o:l:d:v:p:r:e:t:h";
  const option long_opts[] = {
    {"input",      required_argument, nullptr, 'i'},
    {"output",      required_argument, nullptr, 'o'},
    {"log",      required_argument, nullptr, 'l'},
    {"header",      required_argument, nullptr, 'd'},
    {"variant",      required_argument, nullptr, 'v'},
    {"perturb",    required_argument, nullptr, 'p'},
    {"refinement", required_argument, nullptr, 'r'},
    {"epsilon",    required_argument, nullptr, 'e'},
    {"toli",       required_argument, nullptr, 't'},
    {"help",       no_argument,       nullptr, 'h'},
    {nullptr,      no_argument,       nullptr, 0}
  };
  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  while (true) {
   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

   if (-1 == opt)
    break;
   switch (opt) {
    case 'i':
     qp_args.insert(std::pair<std::string, std::string>("input", optarg));
     break;
    case 'o':
     qp_args.insert(std::pair<std::string, std::string>("output", optarg));
     break;
    case 'l':
     qp_args.insert(std::pair<std::string, std::string>("log", optarg));
     break;
    case 'd':
     qp_args.insert(std::pair<std::string, std::string>("header", optarg));
     break;
    case 'v':
     qp_args.insert(std::pair<std::string, std::string>("variant", optarg));
     break;
    case 'p':
     qp_args.insert(std::pair<std::string, std::string>("perturbation", optarg));
     break;
    case 'r':
     qp_args.insert(std::pair<std::string, std::string>("iterations", optarg));
     break;
    case 'e':
     qp_args.insert(std::pair<std::string, std::string>("epsilon", optarg));
     break;
    case 't':
     qp_args.insert(std::pair<std::string, std::string>("tolerance", optarg));
     break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
     print_help();
     break;
   }
  }
  if(qp_args.find("input") == qp_args.end()){
   print_help();
   return false;
  }
  return true;
 }


 void print_header(){
  std::cout<<"Tool Name,Problem Name,Hessian dim,Hessian NNZ,# of Eq Constraints,"
             "Eq Constraint NNZ,# of Ineq Const,Ineq Constraint NNZ,# of Threads,"
             "eps_abs,Outer GMRES Iter,Inner GMRES Iter,GMRES Tol,Diagonal Pert,"
             "Status,# of Iterations,Time (s),Active-set Size,Constraint Satisfaction Inf,"
             "Residual Lagrangian inf,Primal Obj,Dual Obj,Obj Value,Non-negativity Inf,Complementarity Inf,"
             "Problem Type,Problem Class,\n";
 }
}

#endif //NASOQ_BENCHMARKS_DRIVER_UTILS_H
