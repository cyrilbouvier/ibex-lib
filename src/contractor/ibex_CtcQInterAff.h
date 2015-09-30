//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor specialized by adding a qintersection in the direction of the gradient of the function of the measurement
//  
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 30, 2012
// Last Update : Apr 30, 2012
//============================================================================



#ifndef __IBEX_CTC_Q_INTERAFF_H__
#define __IBEX_CTC_Q_INTERAFF_H__

#include "ibex_CtcQInter.h"

#include "ibex_Function.h"

#include <vector>


using namespace std;

namespace ibex {
class CtcQInterAff : virtual public CtcQInter {

 public:

  CtcQInterAff(int n, const Array<Ctc>& ctc_list, int q,  qintermethod meth=QINTERPROJ,  int kfun=1);
  CtcQInterAff(int n, const Array<Ctc>& ctc_list, int pmax, int q,  qintermethod meth=QINTERPROJ,  int kfun=1);
 CtcQInterAff(int n, const Array<Ctc>& ctc_list, int q, Function*** funlist,qintermethod meth=QINTERPROJ,  int kfun=1 );

   ~CtcQInterAff();


	/**
	 * List of functions 
	 */
	Function*** function_list;
	void contract(IntervalVector& box);

 protected: 
	int pmax;
	double *** alpha;
	double ** valmean;
	double ** err;
        void init_interaff();
	Interval* interaf;
	virtual double slope_compute(int iter, int j, int i ,IntervalVector& box, AffineLin& af);
	vector<double> aff_dir; 
	virtual double err_compute( int iter, int k, IntervalVector& box,AffineLin& af);
	virtual double	valmean_compute(int iter, int i, IntervalVector& box, AffineLin& af);
	void affine_projection (IntervalVector& box, int& p);
	virtual void compute_affine_evaluation( int i, int iter,  AffineLin& af, Interval& af2);
 };


}
#endif // __IBEX_CTC_Q_INTERAFF_H__
