//
// Created by Kazem on 3/1/19.
//

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <cfloat>
#include <types.h>
#include <chrono>
#include <cmath>
#include <mkl_service.h>
#include "stdio.h"
#include <mp_format_converter.h>
#include <osqp.h>
#include <Util.h>
#include "driver_utils.h"

#define MAX_DBL 1e20
/*
 * writing a vector to a file.
 */
int write_vector_osqp(std::string fName, c_int n, c_float* vec_vals,
                      c_int prec=30){
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

void export_to_file(std::string p_name, c_int m, c_int n,
                   c_float *dual_vars, c_float* primal_vars,
                   c_int num_iter, c_float primal_obj, c_float tot){
 std::string primal_name = p_name + "_osqp_primal.txt";
 write_vector_osqp(primal_name,n,primal_vars);
 std::string dual_name = p_name + "_osqp_dual.txt";
 write_vector_osqp(dual_name,m,dual_vars);
 c_float info[5]={};
 info[0] = num_iter;
 info[1] = primal_obj;
 info[2] = 0;
 info[3] = primal_obj;
 info[4] = tot;
 std::string other_name = p_name + "_osqp_info.txt";
 write_vector_osqp(other_name,5,info);
}

c_float absolute(c_float a){
 double tmp = a>=0 ? a : -a;
 return tmp;
}

double norm_dense
  (
    /* ---- input ---- */
    c_int nrow,
    c_int ncol,
    c_float *X,	/* matrix to compute the norm of */
    c_int norm		/* type of norm: 0: inf. norm, 1: 1-norm, 2: 2-norm */
  )
{
 if (norm < 0 || norm > 2 || (norm == 2 && ncol > 1)){
  std::cout<<"NORM error";
  return -1;
 }
 c_int d = nrow;
 c_float xnorm = 0, s = 0 ;
 /* infinity-norm = max row sum, using stride-1 access of X */
 if (norm == 0) {
  /* infinity-norm = max row sum, using stride-d access of X */
  for (c_int i = 0; i < nrow; i++) {
   s = 0;
   for (c_int j = 0; j < ncol; j++) {
    s += absolute(X[i + j * d]);
   }
   if ((!std::isfinite(s) || s > xnorm) && std::isfinite(xnorm) ) {
    xnorm = s;
   }
  }
 } else if (norm == 1) {
  /* 1-norm = max column sum */
  for (c_int j = 0; j < ncol; j++) {
   s = 0;
   for (c_int i = 0; i < nrow; i++) {
    s += absolute(X[i + j * d]);
   }
   if ((!std::isfinite(s) || s > xnorm) && std::isfinite(xnorm)) {
    xnorm = s;
   }
  }
 } else {
  /* 2-norm = sqrt (sum (X.^2)) */
  for (c_int i = 0; i < nrow; i++) {
   c_float x = X[i];
   xnorm += x * x;
  }
  xnorm = sqrt(xnorm);
 }
 return (xnorm) ;
}

c_int spmv_csc_small_osqp (c_int m, c_int n, c_int *Ap, c_int *Ai, c_float *Ax,
                    c_float *x, c_float *y)
{
 c_int p, j;
 std::fill_n(y,m,0);
 if (!Ap || !x || !y) return (0) ;       /* check inputs */
 for (j = 0 ; j < n ; j++) {
  for (p = Ap [j] ; p < Ap [j+1] ; p++) {
   y [Ai [p]] += Ax [p] * x [j] ;
  }
 }
 return (1) ;
}

c_float compute_non_negativity(c_int m,  c_int n, c_int *Hp,
  c_int *Hi,c_float *Hx, c_float *primal_vars, c_float *y,
  c_float *l, c_float *u, c_float *Ax){
 double nn1, nn2;
 c_float *axly = new c_float[m]();
 c_float *axuy = new c_float[m]();
 c_float *tmp = new c_float[m];
 spmv_csc_small_osqp(m,n,Hp, Hi, Hx, primal_vars,tmp);
 for (int i = 0; i < m; ++i) {
  if(l[i] != u[i] || !(l[i] <= -MAX_DBL && u[i] >= MAX_DBL) ){
   double yn = std::min(y[i],0.0);
   double yp = std::max(y[i],0.0);
   if(l[i] > -MAX_DBL){
    axly[i] = (tmp[i]-l[i]) * yn;
   }
   if(u[i] < MAX_DBL){
    axuy[i] = (tmp[i]-u[i]) * yp;
   }
  }
 }

 nn1 = norm_dense(m, 1, axly,0);
 nn2 = norm_dense(m, 1, axuy,0);
 c_float  non_negativity_infn = std::max(nn1, nn2);
 delete []tmp;
 delete []axly;
 delete []axuy;
 return non_negativity_infn;
}



/*
 *
 */
double compute_slackness_osqp(c_int m, c_int n, c_int *Hp, c_int *Hi,
  c_float *Hx,c_float *primal_vars, c_float *y, c_float *l, c_float *u){
 double slack, slack_l=0, slack_u=0;
 c_float *tmp = new c_float[m];
 c_float *axly = new c_float[m]();
 c_float *axuy = new c_float[m]();
 spmv_csc_small_osqp(m,n,Hp, Hi, Hx, primal_vars,tmp);
 for (int i = 0; i < m; ++i) {
  if(l[i] != u[i]){
    double yn = y[i] < 0 ? y[i] : 0; //min(y,0)
    axly[i] = std::min(std::abs(tmp[i]-l[i]), -yn);
    double yp = y[i] < 0 ? 0 : y[i]; //max(y,0)
    axuy[i] = std::min(std::abs(-tmp[i]+u[i]), yp);
  }
  /*else{
   axuy[i] = std::min(std::abs(tmp[i]-u[i]), y[i]);
   //std::cout<<std::min(std::abs(tmp[i]-l[i]), y[i])<<";";
  // std::cout<< tmp[i]-l[i]<<";";
  }*/
 }
 slack_u = norm_dense(m, 1, axuy,0);
 slack_l = norm_dense(m, 1, axly,0);
 slack = std::max(slack_l,slack_u);
 delete []tmp;
 delete []axly;
 delete []axuy;
 return slack;
}

double dot(c_int n, c_float *a, c_float *b){
 double result = 0.0;
 for (int i = 0; i < n; ++i) {
  result += (a[i]*b[i]);
 }
 return result;
}
/*
 * Computes primal objective
 */
double compute_primal_obj(c_int n, c_int *Hp, c_int *Hi, c_float *Hx,
                          c_float *primal_vars, c_float *q){
 double primal_obj=0;
 c_float *tmp = new c_float[n];
 spmv_csc_small_osqp(n,n,Hp, Hi, Hx, primal_vars,tmp);
 /*CSC *HTT = ptranspose(H,2,NULL,NULL,0,status);
 primal_obj = quad_form(HTT,primal_vars);*/
 primal_obj = 0.5 * dot(n,tmp,primal_vars);
 primal_obj += dot(n,q,primal_vars);
 return primal_obj;
}

/*
* Computes min(b-Bx,0)
*/
c_float constraint_sat_norm_osqp(c_int ad1, c_int ad2,
                           c_int nnzA,
                           c_int* &colA, c_int* &rowA, c_float * &valA,
                           c_float *b, c_float *primal_vars){
 c_float constraint_sat_norm=0;
 c_float *tmp = new c_float[ad2];
 spmv_csc_small_osqp(ad1,ad2,colA, rowA, valA, primal_vars, tmp);
 for (int i = 0; i < ad1; ++i) {
  c_float diff= tmp[i]-b[i];
  //tmp[i] = std::min(tmp[i]-b[i],0.0);
  //tmp[i] = std::max(b[i] - tmp[i], 0.0); //Equal to top
  tmp[i] = diff < 0.0 ? 0.0 : diff;
 }
 constraint_sat_norm = norm_dense(ad1,1,tmp,0);
 delete []tmp;
 return constraint_sat_norm;
}


double error_L2(c_int n, c_float *expected, c_float *x){
 c_float x_error = 0.0;
 c_float exp_sum = 0.0;
 for (int i = 0; i < n; i++ ){
  x_error += (expected[i] - x[i]) * ( expected[i] - x[i]) ;
 }
 x_error = sqrt ( x_error );
 return x_error;
}

/*
 *
 * /home/kazem/Dropbox/UFDB/contact/sparse_version/recontactqps_P.mtx /home/kazem/Dropbox/UFDB/contact/sparse_version/recontactqps_A.mtx   /home/kazem/Dropbox/UFDB/contact/sparse_version/recontactqps_q none /home/kazem/Dropbox/UFDB/contact/sparse_version/recontactqps_l /home/kazem/Dropbox/UFDB/contact/sparse_version/recontactqps_cvx
 * /home/kazem/UFDB/graphics/8/P9 /home/kazem/UFDB/graphics/8/A9 /home/kazem/UFDB/graphics/8/q9 none  /home/kazem/UFDB/graphics/8/l9 /home/kazem/UFDB/graphics/8/primal9
 */

int main(int argc, char **argv) {
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;

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

 // Load problem data
 int is_polish = 0;
 c_float *P_x, *Pu_x;
 c_int   P_nnz = 0, Pu_nnz=0;
 c_int   *P_i, *Pu_i;
 c_int   *P_p, *Pu_p;
 c_float *q;
 c_float *A_x;
 c_int   A_nnz;
 c_int   *A_i;
 c_int   *A_p;
 c_float *l;
 c_float *u;
 c_int n;
 c_int m;
 std::string p_name;

 is_polish = solver_mode;

 //std::string f6 = argv[7];

 auto *QPFC = new format::QPFormatConverter();
 if(!QPFC->load_smp(input_qp_path))
  return -1;
 QPFC->smp_to_bounded();
 p_name = QPFC->smp_->desc_struct_.name_;
 n = QPFC->num_var();
 m = QPFC->bf_->A ? QPFC->bf_->A->m : 0;
 P_nnz = QPFC->bf_->H->nnz;
 A_nnz = QPFC->bf_->A ? QPFC->bf_->A->nnz : 0;
 // To keep the default setting of OSQP

 auto int_to_cint =[](int *t, int n){
  auto *tc = new c_int[n]();
  for (int i = 0; i < n; ++i) {
   tc[i] = t[i];
  }
  return tc;
 };
// Just to avoid double free
 auto double_to_cfloat =[](double *t, int n){
  auto *tc = new c_float[n]();
  for (int i = 0; i < n; ++i) {
   tc[i] = t[i];
  }
  return tc;
 };

 q = double_to_cfloat(QPFC->bf_->q->a, n);
 l = double_to_cfloat(QPFC->bf_->l->a, m);
 u = double_to_cfloat(QPFC->bf_->u->a, m);

 auto *HT = sym_lib::transpose_general(QPFC->bf_->H);
 //auto *Hf = sym_lib::make_full(QPFC->bf_->H);
 //auto *HT = sym_lib::make_half(Hf->n, Hf->p, Hf->i, Hf->x,false);
 P_p = int_to_cint(HT->p,n+1);
 P_i = int_to_cint(HT->i,P_nnz);
 P_x = double_to_cfloat(HT->x,P_nnz);
 //sym_lib::print_csc(n,n,HT->p, HT->i, HT->x);


 A_p = int_to_cint(QPFC->bf_->A->p,n+1);
 A_i = int_to_cint(QPFC->bf_->A->i,A_nnz);
 A_x = double_to_cfloat(QPFC->bf_->A->x,A_nnz);


 c_float *x_exp;
// if (!read_vector_osqp(f6,n,x_exp))
//  return -1;
 // Problem settings
 OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

 c_int exitflag = 0;
 // Structures
 OSQPWorkspace *work; // Workspace
 OSQPData *data;      // OSQPData

 // Populate data
 data    = (OSQPData *)c_malloc(sizeof(OSQPData));
 data->n = n;
 data->m = m;
 data->P = csc_matrix(data->n, data->n, QPFC->bf_->H->nnz, P_x, P_i, P_p);
 data->q = q;
 data->A = csc_matrix(data->m, data->n, QPFC->bf_->A->nnz, A_x, A_i, A_p);
 data->l = l;
 data->u = u;


 // Define Solver settings as default
 osqp_set_default_settings(settings);
 settings->linsys_solver = MKL_PARDISO_SOLVER;
 settings->eps_abs=eps;
 settings->eps_rel=eps;
 settings->max_iter=40000000;
 settings->eps_prim_inf = eps;
 settings->eps_dual_inf = eps;
 settings->polish = is_polish;
 settings->verbose = 0;
 settings->time_limit = 2000.0;
 settings->scaling = 0;
/* settings->polish = 1;
 settings->verbose = 1;
 settings->eps_abs=1e-20;
 settings->eps_rel=1e-3;
 //settings->adaptive_rho = 1;
 settings->scaled_termination = 1;*/

 start = std::chrono::system_clock::now();
 // Setup workspace
 exitflag = osqp_setup(&work, data, settings);
 // Solve Problem
 osqp_solve(work);
 end = std::chrono::system_clock::now();

 elapsed_seconds = end-start;
 double durationSym=elapsed_seconds.count();
 c_float tot = durationSym;
 //c_float abs_nrm2 = error_L2(n,x_exp,work->solution->x);
/* std::cout<<"\n Abs primal Error 2: "<<abs_nrm2<<"\n Time(sec): "<<durationSym<<"\n";
 std::cout<<"constraint sat norm: "
          <<constraint_sat_norm_osqp(data->m,data->n,A_nnz,A_p,A_i,A_x,u,work->solution->x)<<"\n";
 std::cout<<"prime sat norm: "
          <<work->info->pri_res<<"\n";
 std::cout<<"dual residual norm: "
          <<work->info->dua_res<<" \n";
 std::cout<<"Obj: "<<compute_primal_obj(data->n,P_p,P_i,P_x,work->solution->x,q);*/
 /*export_to_file(f3,data->m,data->n,work->solution->y,work->solution->x,
                work->info->iter,work->info->obj_val,tot);*/
 double  nonnegativity, comp_slack1, comp_slack2;
 nonnegativity = compute_non_negativity(data->A->m,
   data->A->n,data->A->p,data->A->i,
   data->A->x,work->solution->x,work->solution->y,
   data->l, data->u, work->Ax);

 if(print_header)
  nasoq_bench::print_header();

 comp_slack2 = compute_slackness_osqp(data->A->m,data->A->n,
   data->A->p,data->A->i,data->A->x,work->solution->x,
   work->solution->y,data->l,data->u);
 if(settings->polish == 1)
  std::cout<<"OSQP-polished,"<<p_name<<",";
 else
  std::cout<<"OSQP,"<<p_name<<",";
 std::cout<<data->n<<","<<data->P->nzmax<<",";
 std::cout<<"N/A,N/A,"<<data->m<<","<<data->A->nzmax<<",";
 std::cout<<work->linsys_solver->nthreads<<",";
 std::cout<<settings->eps_abs<<",N/A,N/A,N/A,N/A,";
 std::cout<<work->info->status<<","<<work->info->iter<<","<< durationSym<<",N/A,";
 std::cout<<work->info->pri_res<<","<<work->info->dua_res<<",";
 std::cout<<work->info->obj_val<<","<<"N/A"<<",N/A,";
 std::cout<<nonnegativity<<","<<comp_slack2<<","<<QPFC->smp_->desc_struct_.application_<<",";
 std::cout<<QPFC->smp_->desc_struct_.category_<<",";

 // Clean workspace
 osqp_cleanup(work);
 c_free(data->A);
 c_free(data->P);
 c_free(data);
 c_free(settings);

 delete HT;
 delete QPFC;


 return exitflag;
}
