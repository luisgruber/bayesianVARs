#include <Rcpp.h>
#include <lp_lib.h>

void setverbose(lprec* lp_ptr, int verbose){
  // Retrieve the callable function
  void (*set_verbose_ptr)(lprec*, int) = NULL;
  if(!set_verbose_ptr) set_verbose_ptr = (void (*)(lprec*, int)) R_GetCCallable("lpSolveAPI", "set_verbose");
  // Call the function
  set_verbose_ptr(lp_ptr, verbose);
}

void setmaxim(lprec* lp_ptr){
  void (*set_maxim_ptr)(lprec*) = NULL;
  if(!set_maxim_ptr) set_maxim_ptr = (void (*)(lprec*)) R_GetCCallable("lpSolveAPI", "set_maxim");
  set_maxim_ptr(lp_ptr);
}

void setobj(lprec* lp_ptr, int colnr, LPSREAL value){
  void (*set_obj_ptr)(lprec*, int, LPSREAL) = NULL;
  if(!set_obj_ptr) set_obj_ptr = (void (*)(lprec*, int, LPSREAL)) R_GetCCallable("lpSolveAPI", "set_obj");
  set_obj_ptr(lp_ptr, colnr, value);
}

void setbounds(lprec* lp_ptr, int colnr, LPSREAL lower, LPSREAL upper){
  MYBOOL (*set_bounds_ptr)(lprec*, int, LPSREAL, LPSREAL) = NULL;
  if(!set_bounds_ptr) set_bounds_ptr = (MYBOOL (*)(lprec*, int, LPSREAL, LPSREAL)) R_GetCCallable("lpSolveAPI", "set_bounds");
  set_bounds_ptr(lp_ptr, colnr, lower, upper);
}

void setadd_rowmode(lprec* lp_ptr, MYBOOL turnon){
  MYBOOL (*set_add_rowmode_ptr)(lprec*, MYBOOL) = NULL;
  if(!set_add_rowmode_ptr) set_add_rowmode_ptr = (MYBOOL (*)(lprec*, MYBOOL)) R_GetCCallable("lpSolveAPI", "set_add_rowmode");
  set_add_rowmode_ptr(lp_ptr, turnon);
}

MYBOOL addconstraintex(lprec* lp_ptr, int count, LPSREAL* row, int* colno, int constr_type, LPSREAL rh){
  MYBOOL (*add_constraintex_ptr)(lprec*, int, LPSREAL*, int*, int, LPSREAL) = NULL;
  if(!add_constraintex_ptr) add_constraintex_ptr = (MYBOOL (*)(lprec*, int, LPSREAL*, int*, int, LPSREAL)) R_GetCCallable("lpSolveAPI", "add_constraintex");
  return(add_constraintex_ptr(lp_ptr, count, row, colno, constr_type, rh));
}

int lpsolve(lprec* lp_ptr){
  int (*lp_solve_ptr)(lprec*) = NULL;
  if(!lp_solve_ptr) lp_solve_ptr = (int (*)(lprec*)) R_GetCCallable("lpSolveAPI", "solve");
  return(lp_solve_ptr(lp_ptr));
}

void getvariables(lprec* lp_ptr, LPSREAL* var){
  MYBOOL (*get_variables_ptr)(lprec*, LPSREAL*) = NULL;
  if(!get_variables_ptr) get_variables_ptr = (MYBOOL (*)(lprec*, LPSREAL*)) R_GetCCallable("lpSolveAPI", "get_variables");
  get_variables_ptr(lp_ptr, var);
}

void printlp(lprec* lp_ptr){
  void (*print_lp_ptr)(lprec*) = NULL;
  if(!print_lp_ptr) print_lp_ptr = (void (*)(lprec*)) R_GetCCallable("lpSolveAPI", "print_lp");
  print_lp_ptr(lp_ptr);
}

void deletelp(lprec* lp_ptr){
  void (*delete_lp_ptr)(lprec*) = NULL;
  if(!delete_lp_ptr) delete_lp_ptr = (void (*)(lprec*)) R_GetCCallable("lpSolveAPI", "delete_lp");
  delete_lp_ptr(lp_ptr);
}

void setobj_fnex(lprec* lp_ptr, int count, LPSREAL* row, int* colno){
  MYBOOL (*set_obj_fnex_ptr)(lprec*, int, LPSREAL*, int*) = NULL;
  if(!set_obj_fnex_ptr) set_obj_fnex_ptr = (MYBOOL (*)(lprec*, int, LPSREAL*, int*)) R_GetCCallable("lpSolveAPI", "set_obj_fnex");
  set_obj_fnex_ptr(lp_ptr, count, row, colno);
}

void resizelp(lprec* lp_ptr, int rows, int columns){
  MYBOOL (*resize_lp_ptr)(lprec*, int, int) = NULL;
  if(!resize_lp_ptr) resize_lp_ptr = (MYBOOL (*)(lprec*, int, int)) R_GetCCallable("lpSolveAPI", "resize_lp");
  resize_lp_ptr(lp_ptr, rows, columns);
}

void setbinary(lprec* lp_ptr, int colnr, MYBOOL must_be_bin){
  MYBOOL (*set_binary_ptr)(lprec*, int, MYBOOL) = NULL;
  if(!set_binary_ptr) set_binary_ptr = (MYBOOL (*)(lprec*, int, MYBOOL)) R_GetCCallable("lpSolveAPI", "set_binary");
  set_binary_ptr(lp_ptr, colnr, must_be_bin);
}

void setcol_name(lprec* lp_ptr, int colnr, char* name_cstr){
  MYBOOL (*set_col_name_ptr)(lprec*, int, char*) = NULL;
  if(!set_col_name_ptr) set_col_name_ptr = (MYBOOL (*)(lprec*, int, char*)) R_GetCCallable("lpSolveAPI", "set_col_name");
  set_col_name_ptr(lp_ptr, colnr, name_cstr);
}

int getNrows(lprec* lp_ptr){
  int (*get_Nrows_ptr)(lprec*) = NULL;
  if(!get_Nrows_ptr) get_Nrows_ptr = (int (*)(lprec*)) R_GetCCallable("lpSolveAPI", "get_Nrows");
  return(get_Nrows_ptr(lp_ptr));
}

int getNcolumns(lprec* lp_ptr){
  int (*get_Ncolumns_ptr)(lprec*) = NULL;
  if(!get_Ncolumns_ptr) get_Ncolumns_ptr = (int (*)(lprec*)) R_GetCCallable("lpSolveAPI", "get_Ncolumns");
  return(get_Ncolumns_ptr(lp_ptr));
}
