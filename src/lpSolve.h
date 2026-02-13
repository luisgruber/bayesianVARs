#ifndef LPSOLVE_H
#define LPSOLVE_H

#include <Rcpp.h>
#include <lp_lib.h>

void setverbose(lprec* lp_ptr, int verbose);

void setmaxim(lprec* lp_ptr);

void setobj(lprec* lp_ptr, int colnr, LPSREAL value);

void setbounds(lprec* lp_ptr, int colnr, LPSREAL lower, LPSREAL upper);

void setadd_rowmode(lprec* lp_ptr, MYBOOL turnon);

MYBOOL addconstraintex(lprec* lp_ptr, int count, LPSREAL* row, int* colno, int constr_type, LPSREAL rh);

int lpsolve(lprec* lp_ptr);

void getvariables(lprec* lp_ptr, LPSREAL* var);

void printlp(lprec* lp_ptr);

void deletelp(lprec* lp_ptr);

void setobj_fnex(lprec* lp_ptr, int count, LPSREAL* row, int* colno);

void resizelp(lprec* lp_ptr, int rows, int columns);

void setbinary(lprec* lp_ptr, int colnr, MYBOOL must_be_bin);

void setcol_name(lprec* lp_ptr, int colnr, char* name_cstr);

int getNrows(lprec* lp_ptr);

int getNcolumns(lprec* lp_ptr);

#endif
