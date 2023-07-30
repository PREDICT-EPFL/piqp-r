#include <Rcpp.h>
#include <RcppEigen.h>
#include <piqp/piqp.hpp>

void update_settings(piqp::Settings<double> &s, Rcpp::List rs) {
    s.rho_init = Rcpp::as<double>(rs["rho_init"]);
    s.delta_init = Rcpp::as<double>(rs["delta_init"]);
    s.eps_abs = Rcpp::as<double>(rs["eps_abs"]);
    s.eps_rel = Rcpp::as<double>(rs["eps_rel"]);
    s.check_duality_gap = Rcpp::as<bool>(rs["check_duality_gap"]);
    s.eps_duality_gap_abs = Rcpp::as<double>(rs["eps_duality_gap_abs"]);
    s.eps_duality_gap_rel = Rcpp::as<double>(rs["eps_duality_gap_rel"]);
    s.reg_lower_limit = Rcpp::as<double>(rs["reg_lower_limit"]);
    s.reg_finetune_lower_limit = Rcpp::as<double>(rs["reg_finetune_lower_limit"]);
    s.reg_finetune_primal_update_threshold = Rcpp::as<int>(rs["reg_finetune_primal_update_threshold"]);
    s.reg_finetune_dual_update_threshold = Rcpp::as<int>(rs["reg_finetune_dual_update_threshold"]);
    s.max_iter = Rcpp::as<int>(rs["max_iter"]);
    s.max_factor_retires = Rcpp::as<int>(rs["max_factor_retires"]);
    s.preconditioner_scale_cost = Rcpp::as<bool>(rs["preconditioner_scale_cost"]);
    s.preconditioner_iter = Rcpp::as<int>(rs["preconditioner_iter"]);
    s.tau = Rcpp::as<double>(rs["tau"]);
    s.iterative_refinement_always_enabled = Rcpp::as<bool>(rs["iterative_refinement_always_enabled"]);
    s.iterative_refinement_eps_abs = Rcpp::as<double>(rs["iterative_refinement_eps_abs"]);
    s.iterative_refinement_eps_rel = Rcpp::as<double>(rs["iterative_refinement_eps_rel"]);
    s.iterative_refinement_max_iter = Rcpp::as<int>(rs["iterative_refinement_max_iter"]);
    s.iterative_refinement_min_improvement_rate = Rcpp::as<double>(rs["iterative_refinement_min_improvement_rate"]);
    s.iterative_refinement_static_regularization_eps = Rcpp::as<double>(rs["iterative_refinement_static_regularization_eps"]);
    s.iterative_refinement_static_regularization_rel = Rcpp::as<double>(rs["iterative_refinement_static_regularization_rel"]);
    s.verbose = Rcpp::as<bool>(rs["verbose"]);
    s.compute_timings = Rcpp::as<bool>(rs["compute_timings"]);
}

//' Solve using the piqp dense matrix interface
//'
//' @param P a symmetric positive definite matrix
//' @param c the objective cost
//' @param A the coefficient matrix of primal constraints
//' @param b the right hand side vector of the primal constraints
//' @param G the matrix of inequality constraints
//' @param h the right side vector of the matrix inequality constraints
//' @param x_lb the lower bounds for the solution
//' @param x_ub the upper bounds for the solution
//' @param rs a list of control settings
//' @return named list
//'   - `x` (primal solution)
//'   - `y` (dual solution of equality constraints)
//'   - `z` (dual solution of inequality constraints)
//'   - `z_lb` (dual solution of lower bound box constraints)
//'   - `z_ub` (dual solution of upper bound box constraints)
//'   - `obj` (primal objective value)
//'   - `run_time` (total runtime, if so specified in settings)
// [[Rcpp::export]]
SEXP piqp_dense_solve (Eigen::Map<Eigen::MatrixXd> P,
		       Eigen::Map<Eigen::VectorXd> c,
		       Eigen::Map<Eigen::MatrixXd> A,
		       Eigen::Map<Eigen::VectorXd> b,
		       Eigen::Map<Eigen::MatrixXd> G,
		       Eigen::Map<Eigen::VectorXd> h,
		       Eigen::Map<Eigen::VectorXd> x_lb,
		       Eigen::Map<Eigen::VectorXd> x_ub,
		       Rcpp::List rs) {
    piqp::DenseSolver<double> solver;
    update_settings(solver.settings(), rs);
    solver.setup(P, c, A, b, G, h, x_lb, x_ub);
    piqp::Status status = solver.solve();
    piqp::Result<double> result = solver.result();
    // double run_time = 0.0;

    return Rcpp::List::create(
			      Rcpp::_["status"] = (int) status,			      
			      Rcpp::_["x"] = result.x,
			      Rcpp::_["y"] = result.y,
			      Rcpp::_["z"] = result.z,
			      Rcpp::_["z_lb"] = result.z_lb,
			      Rcpp::_["z_ub"] = result.z_ub,
			      Rcpp::_["pobj"] = result.info.primal_obj,
			      Rcpp::_["dobj"] = result.info.dual_obj,
			      Rcpp::_["dgap"] = result.info.duality_gap,
			      Rcpp::_["dgap_rel"] = result.info.duality_gap_rel,
			      Rcpp::_["run_time"] = solver.settings().compute_timings? result.info.run_time: 0.0);
}

//' Solve using the piqp sparse matrix interface
//'
//' @param P a symmetric positive definite matrix
//' @param c the objective cost
//' @param A the coefficient matrix of primal constraints
//' @param b the right hand side vector of the primal constraints
//' @param G the matrix of inequality constraints
//' @param h the right side vector of the matrix inequality constraints
//' @param x_lb the lower bounds for the solution
//' @param x_ub the upper bounds for the solution
//' @param rs a list of control settings
//' @return named list
//'   - `x` (primal solution)
//'   - `y` (dual solution of equality constraints)
//'   - `z` (dual solution of inequality constraints)
//'   - `z_lb` (dual solution of lower bound box constraints)
//'   - `z_ub` (dual solution of upper bound box constraints)
//'   - `obj` (primal objective value)
//'   - `run_time` (total runtime, if so specified in settings)
// [[Rcpp::export]]
Rcpp::List piqp_sparse_solve (Eigen::Map<Eigen::SparseMatrix<double>> P,
			      Eigen::Map<Eigen::VectorXd> c,
			      Eigen::Map<Eigen::SparseMatrix<double>> A,
			      Eigen::Map<Eigen::VectorXd> b,
			      Eigen::Map<Eigen::SparseMatrix<double>> G,
			      Eigen::Map<Eigen::VectorXd> h,
			      Eigen::Map<Eigen::VectorXd> x_lb,
			      Eigen::Map<Eigen::VectorXd> x_ub,
			      Rcpp::List rs) {
    piqp::SparseSolver<double> solver;
    update_settings(solver.settings(), rs);
    solver.setup(P, c, A, b, G, h, x_lb, x_ub);
    piqp::Status status = solver.solve();
    piqp::Result<double> result = solver.result();
    return Rcpp::List::create(
			      Rcpp::_["status"] = (int) status,			      
			      Rcpp::_["x"] = result.x,
			      Rcpp::_["y"] = result.y,
			      Rcpp::_["z"] = result.z,
			      Rcpp::_["z_lb"] = result.z_lb,
			      Rcpp::_["z_ub"] = result.z_ub,
			      Rcpp::_["pobj"] = result.info.primal_obj,
			      Rcpp::_["dobj"] = result.info.dual_obj,
			      Rcpp::_["dgap"] = result.info.duality_gap,
			      Rcpp::_["dgap_rel"] = result.info.duality_gap_rel,
			      Rcpp::_["run_time"] = solver.settings().compute_timings? result.info.run_time: 0.0);
}
