//
// Created by kazem on 10/15/19.
//

#include <iostream>
#include <chrono>
#include <omp.h>
#include "Gurobi.h"
#include <mp_format_converter.h>
#include "Eigen_utils.h"
#include "qp_utils.h"
#include "driver_utils.h"

using namespace Eigen;



/*
 * writing a vector to a file.
 */
int write_vector_mosek(std::string fName, int n,
                       const VectorXd& vec_vals,
                       int prec=30){
 if(n <= 0 )
  return 0;
 double value;
 std::ofstream out_file;
 out_file.precision(prec);
 out_file.open (fName);
 for(int i = 0; i<n; i++){//writing from file row by row
  value=vec_vals[i];
  out_file<<value;
  if(i!=n-1)
   out_file<<"\n";
 }
 out_file.close();
 return 1;
}

void export_to_file(std::string p_name, int m, int n,
                    const VectorXd& dual_vars,
                    const VectorXd& primal_vars,
                    int num_iter, double primal_obj, double tot){
 std::string primal_name = p_name + "_mosek_primal.txt";
 write_vector_mosek(primal_name,n,primal_vars);
 std::string dual_name = p_name + "_mosek_dual.txt";
 write_vector_mosek(dual_name,m,dual_vars);
 VectorXd info(5);
 info[0] = num_iter;
 info[1] = primal_obj;
 info[2] = 0;
 info[3] = primal_obj;
 info[4] = tot;
 std::string other_name = p_name + "_mosek_info.txt";
 write_vector_mosek(other_name,5,info);
}


bool gurobiSolve(SparseMatrix<double>& MDK, const VectorXd& b,
                 SparseMatrix<double>& Aeq, const VectorXd& beq,
                 SparseMatrix<double>& Aineq, const VectorXd& bineq,
                 VectorXd& v, VectorXd& w,VectorXd& ww, double &elapsed_time,
                 const double eps, double &obj, int &num_iter,
                 int num_thread);

int main(int argc, char **argv) {
 double mosek_time=0;

 std::map<std::string, std::string> qp_args;
 if (!nasoq_bench::parse_args(argc, argv, qp_args))
  return -1;

 /// New settings if provided
 std::string input_qp_path = qp_args["input"];
 double reg_diag = pow(10,-9);
 double zero_threshold = reg_diag;
 double eps = 1e-3;
 int inner_iter = 2;
 int outer_iter = 2;
 double stop_tol = 1e-15;
 int solver_mode = 0;
 bool print_header = false;
 if (qp_args.find("variant") != qp_args.end()) {
  std::string osqp_mode = qp_args["variant"];
  if (osqp_mode == "polished") {
   solver_mode = 1;
  }
 }
 if(qp_args.find("epsilon") != qp_args.end())
  eps = pow(10, std::stoi(qp_args["epsilon"]) );
 if(qp_args.find("header") != qp_args.end())
  print_header = true;


 auto *QPFC = new format::QPFormatConverter();
 if(!QPFC->load_smp(input_qp_path))
  return -1;
 QPFC->smp_to_ie();
 std::string  p_name = QPFC->smp_->desc_struct_.name_;


 int n = QPFC->num_var();
 int  me = QPFC->ief_->get_num_eqc();
 int mie = QPFC->ief_->get_num_ineqc();
 int  P_nnz = QPFC->ief_->H->nnz;
 int A_nnz = QPFC->ief_->A ? QPFC->ief_->A->nnz : 0;

 int is_dual = 1;

 int status = 0;
 int *colA, *rowA, *colH, *rowH;
 double *valA, *valH, *q;
 double *b_ineq;
 double *optimal_primal, *optimal_dual;
 size_t sizeH, nnzA, ad1, ad2, nnzH;
 double *init_x = NULL;
 int *init_as = NULL;
 double objective = 0;
 double *x, *u, *w;
 int num_iter=0, num_thread=omp_get_max_threads();

/*
 SparseMatrix<double> H_eigen(QPFC->H->ncol,QPFC->H->ncol);
if(qp_type == 2){
 convert_mtx_to_eigen(QPFC->H_full->ncol,QPFC->H_full->ncol,
                      QPFC->H_full->p,QPFC->H_full->i,QPFC->H_full->x,
                      H_eigen);
}else if(qp_type == 1){ // for type 1, the hessian should be passed in full
 convert_mtx_to_eigen(QPFC->H->ncol,QPFC->H->ncol,
                      QPFC->H->p,QPFC->H->i,QPFC->H->x,
                      H_eigen);
}
*/


 SparseMatrix<double> H_eigen(n,n);
 convert_mtx_to_eigen(n,n,QPFC->ief_->H_general->p,QPFC->ief_->H_general->i,
                      QPFC->ief_->H_general->x,
                      H_eigen);
 SparseMatrix<double> B_eq_eig(me,n);
 if(me>0)
  convert_mtx_to_eigen(me,n,
                       QPFC->ief_->A->p,QPFC->ief_->A->i,QPFC->ief_->A->x,
                       B_eq_eig);

 SparseMatrix<double> A_eig(mie,n);
 if(mie>0)
  convert_mtx_to_eigen(mie,n,
                       QPFC->ief_->C->p,QPFC->ief_->C->i,QPFC->ief_->C->x,
                       A_eig);


 VectorXd q_eig(n);
 for (int i = 0; i < n; ++i) {
  q_eig[i] = QPFC->ief_->q->a[i];
 }

 VectorXd b_eq_eig(me);
 if(me>0){
  for (int i = 0; i < me; ++i) {
   b_eq_eig[i] = QPFC->ief_->b->a[i];
  }
 }

 VectorXd b_eig(mie);
 if(mie>0){
  for (int i = 0; i < mie; ++i) {
   b_eig[i] = QPFC->ief_->d->a[i];
  }
 }

 x = new double[n]();
 if(me>0)
  u = new double[me]();
 else
  u = NULL;
 if(mie>0)
  w = new double[mie]();
 else
  w = NULL;
 VectorXd primal(n);
 VectorXd dual_eq(me);
 VectorXd dual_ineq(mie);



 bool solve_state = gurobiSolve(H_eigen,q_eig,B_eq_eig,b_eq_eig,A_eig,b_eig,
                               primal,dual_eq,dual_ineq,mosek_time,eps,objective,
                               num_iter,num_thread);
/* if(solve_state){
  std::cout<<"Solved in: "<<mosek_time<<";";
 }else{
  std::cout<<"not solved;";
 }*/
 for (int i = 0; i < n; ++i) {
  x[i] = primal[i];
 }

 for (int i = 0; i < me; ++i) {
  u[i] = -dual_eq[i];
 }

 for (int i = 0; i < mie; ++i) {
  w[i] = -dual_ineq[i];
 }

 //TODO: remove this
 auto H_old = nasoq_bench::new_to_old(QPFC->ief_->H);
 auto A_old = nasoq_bench::new_to_old(QPFC->ief_->A);
 auto AT_old = nasoq_bench::new_to_old(QPFC->ief_->AT);
 auto C_old = nasoq_bench::new_to_old(QPFC->ief_->C);
 auto CT_old = nasoq_bench::new_to_old(QPFC->ief_->CT);
 double *q_o = QPFC->ief_->q ? QPFC->ief_->q->a : NULLPNTR;
 double *b_o = QPFC->ief_->b ? QPFC->ief_->b->a : NULLPNTR;
 double *d_o = QPFC->ief_->d ? QPFC->ief_->d->a : NULLPNTR;
 double obj_ql = nasoq::compute_primal_obj(x,H_old, q_o);
 double cs_norm = nasoq::constraint_sat_norm(C_old,A_old,d_o,
                                             b_o,x);
 double nn_norm = nasoq::non_negativity_norm(C_old,w);
 double comp_norm = nasoq::complementarity_norm(C_old,A_old,d_o,b_o,x,w);

 //std::cout<<"dsddddddd\n";
 //print_vec("Bineq_dual",0,QPFC->A_ineq->nrow,w);
// std::cout<<"dsddddddd\n";
 double lag_norm = nasoq::lagrangian_residual_norm(H_old,C_old,CT_old,
                                                   A_old,AT_old,q_o,x,
                                                   w,u);

 if(print_header)
  nasoq_bench::print_header();

 std::cout<<"Gurobi,";
 QPFC->print_log();
 std::cout<<num_thread<<","<<eps<<",N/A,N/A,N/A,N/A,"<<
          solve_state<<
          ","<<num_iter<<","<<mosek_time
          <<",N/A,";
 std::cout<<cs_norm<<","<<lag_norm<<","<<obj_ql<<",N/A,N/A,";
 std::cout<<nn_norm<<","<<comp_norm<<",";
 std::cout<<QPFC->smp_->desc_struct_.application_<<",";
 std::cout<<QPFC->smp_->desc_struct_.category_<<",";
/* std::cout<<"Primal: \n";
 for (int j = 0; j < 20; ++j) {
  std::cout<<primal[j]<<";";
 }
 std::cout<<"\n";*/


 delete []x;
 if(me>0)
  delete []u;
 if(me>0)
  delete []w;
 delete QPFC;
 return 1;
}


bool gurobiSolve(SparseMatrix<double>& MDK, const VectorXd& b,
                 SparseMatrix<double>& Aeq, const VectorXd& beq,
                 SparseMatrix<double>& Aineq, const VectorXd& bineq,
                 VectorXd& v,VectorXd& w,VectorXd& ww, double &elapsed_time,
                 const double eps, double &obj, int &num_iter,
                 int num_thread)
{
 std::chrono::time_point <std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 GurobiSparse qp(b.size(), beq.size(), bineq.size());
 qp.displayOutput(false);

 SparseMatrix<double> Sb(b.sparseView());
 SparseVector<double> Sbeq(beq.sparseView());
 SparseVector<double> Sbineq(bineq.sparseView());

 MDK.makeCompressed();
 Sb.makeCompressed();
 Aeq.makeCompressed();
 Aineq.makeCompressed();

 VectorXd XL, XU;
 double inf = std::numeric_limits<double>::infinity();
 XL.setConstant(b.size(), -inf);
 XU.setConstant(b.size(), inf);
 qp.optimalityTolerance(eps);
 qp.feasibilityTolerance(eps);
 qp.threads(num_thread);
 start = std::chrono::system_clock::now();
 bool success = qp.solve(MDK, Sb,
                         Aeq, Sbeq,
                         Aineq, Sbineq,
                         XL, XU);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 elapsed_time = elapsed_seconds.count();
 v = qp.result();
 w = qp.dualEqResult();
 ww = qp.dualIneqResult();
 num_iter = qp.iter();

 return success;
}
