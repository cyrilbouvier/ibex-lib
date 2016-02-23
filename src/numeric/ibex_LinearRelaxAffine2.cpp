//============================================================================
//                                  I B E X                                   
// File        : ibex_CtcART.cpp
// Author      : Ignacio Araya, Bertrand Neveu , Gilles Trombettoni
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 1, 2012
// Last Update : Nov 15, 2012
//============================================================================

#include "ibex_LinearRelaxAffine2.h"

namespace ibex {

// the constructor
LinearRelaxAffine2::LinearRelaxAffine2(const System& sys1) :
				LinearRelax(sys1), sys(sys1) {

}

LinearRelaxAffine2::~LinearRelaxAffine2() {

}

bool LinearRelaxAffine2::goal_linearization(const IntervalVector& box, LinearSolver& lp_solver) {
	// Linearization of the objective function by AF2
	Affine2 af2;
	
	sys.goal->eval_affine2(box,af2);

	if (af2.is_empty()) {
		return false;
	}
	try {
		if (af2.size() == sys.nb_var) { // if the affine2 form is valid
			// convert the epsilon variables to the original box
			double tmp=0;
			for (int i =0; i <sys.nb_var; i++) {
				tmp = box[i].rad();
				if (tmp==0) { // sensible case to avoid rowconst[i]=NaN
					if (af2.val(i)==0)
						lp_solver.setVarObj(i, 0);
					else {
						return false; // sensible case to avoid
					}
				} else {
					lp_solver.setVarObj(i, af2.val(i) / tmp);
				}
			}
		}
		else {
			return false;
		}
		return true;

	} catch (LPException&) {
		return false;
	}
}


int LinearRelaxAffine2::inlinearization(const IntervalVector& box, LinearSolver& lp_solver) {
	// TODO a verifier et finir

	Affine2 af2;

	int cont=0;
	Interval ev(0), center(0), err(0);
	Vector rowconst(sys.nb_var);

	// Create the linear relaxation of each constraint
	for (int ctr = 0; ctr < sys.nb_ctr; ctr++) {
		af2 = 0.0;
		CmpOp op = sys.ctrs[ctr].op;
		ev = sys.ctrs[ctr].f.eval_affine2(box, af2);
		//std::cout <<ev<<":::"<< af2<<"  "<<af2.size()<<"  " <<sys.nb_var<< std::endl;

		if (af2.size() == sys.nb_var) { // if the affine2 form is valid
			bool b_abort=false;
			// convert the epsilon variables to the original box
			double tmp=0;
			center =0;
			err =0;
			for (int i =0;(!b_abort) &&(i <sys.nb_var); i++) {
				tmp = box[i].rad();
				if (tmp==0) { // sensible case to avoid rowconst[i]=NaN
					if (af2.val(i)==0)
						rowconst[i]=0;
					else {
						b_abort =true;
					}
				} else {
					rowconst[i] =af2.val(i) / tmp;
					center += rowconst[i]*box[i].mid();
					err += fabs(rowconst[i])*  pow(2,-50);
				}
			}
			if (!b_abort) {
				switch (op) {
				case LEQ:
				case LT: {
					if (0.0 < ev.ub()) {
						try {// TODO TO CHECK
							lp_solver.addConstraint(rowconst, LEQ,	(-(af2.err()+err) - (af2.center()-center)).lb());
							cont++;
						} catch (LPException&) { }
					}
					break;
				}
				case GEQ:
				case GT: {
					if (ev.lb() < 0.0) {
						try {// TODO TO CHECK
							lp_solver.addConstraint(rowconst, GEQ,	((af2.err()+err) - (af2.center()-center)).ub());
							cont++;
						} catch (LPException&) { }
					}
					break;
				}
				case EQ: {
					not_implemented("LinearRelaxAffine2::inlinearization not implemented for equality constrains");
				}
				default:
					break;
				}
			}
		}

	}

	return -1;
}


/*********generation of the linearized system*********/
int LinearRelaxAffine2::linearization(const IntervalVector& box, LinearSolver& lp_solver) {

	Affine2 af2;
	Vector rowconst(sys.nb_var);
	Interval ev(0.0);
	Interval center(0.0);
	Interval err(0.0);
	CmpOp op;
	int cont = 0;

	// Create the linear relaxation of each constraint
	for (int ctr = 0; ctr < sys.nb_ctr; ctr++) {

		af2 = 0.0;
		op = sys.ctrs[ctr].op;

		ev = sys.ctrs[ctr].f.eval_affine2(box, af2);

		if (ev.is_empty()) {
			af2.set_empty();
		}
		//std::cout <<ev<<":::"<< af2<<"  "<<af2.size()<<"  " <<sys.nb_var<< std::endl;

		if (af2.size() == sys.nb_var) { // if the affine2 form is valid
			bool b_abort=false;
			// convert the epsilon variables to the original box
			double tmp=0;
			center =0;
			err =0;
			for (int i =0;(!b_abort) &&(i <sys.nb_var); i++) {
				tmp = box[i].rad();
				if (tmp==0) { // sensible case to avoid rowconst[i]=NaN
					if (af2.val(i)==0)
						rowconst[i]=0;
					else {
						b_abort =true;
					}
				} else {
					rowconst[i] =af2.val(i) / tmp;
					center += rowconst[i]*box[i].mid();
					err += fabs(rowconst[i])*  pow(2,-50);
				}
			}
			if (!b_abort) {
				switch (op) {
				case LT:
					if (ev.lb() == 0.0) return -1;
				case LEQ:
					if (0.0 < ev.lb()) return -1;
					else if (0.0 < ev.ub()) {
						try {
							lp_solver.addConstraint(rowconst, LEQ,	((af2.err()+err) - (af2.center()-center)).ub());
							cont++;
						} catch (LPException&) { }
					}
					break;
				case GT:
					if (ev.ub() == 0.0) return -1;
				case GEQ:
					if (ev.ub() < 0.0) return -1;
					else if (ev.lb() < 0.0) {
						try {
							lp_solver.addConstraint(rowconst, GEQ,	(-(af2.err()+err) - (af2.center()-center)).lb());
							cont++;
						} catch (LPException&) { }
					}
					break;
				case EQ:
					if (!ev.contains(0.0)) return -1;
					else {
						if (ev.diam()>2*lp_solver.getEpsilon()) {
							try {
								lp_solver.addConstraint(rowconst, GEQ,	(-(af2.err()+err) - (af2.center()-center)).lb());
								cont++;
								lp_solver.addConstraint(rowconst, LEQ,	((af2.err()+err) - (af2.center()-center)).ub());
								cont++;
							} catch (LPException&) { }
						}
					}
					break;
				}
			}
		}

	}
	return cont;

}

//void CtcART::convert_back(IntervalVector & box, IntervalVector & epsilon) {
//
//	for (int i = 0; i < box.size(); i++) {
//		box[i] &= box[i].mid() + (box[i].rad() * epsilon[i]);
//	}
//}

}
