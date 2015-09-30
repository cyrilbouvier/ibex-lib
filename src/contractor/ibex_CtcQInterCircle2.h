//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor specialized for circle detection 
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Feb 4, 2015
// Last Update : Feb 4, 2015
//============================================================================

#ifndef __IBEX_CTC_Q_INTERCIRCLE2_H__
#define __IBEX_CTC_Q_INTERCIRCLE2_H__

#include "ibex_CtcQInterAff.h"


namespace ibex {

class CtcQInterCircle2 : virtual public CtcQInter {
public:
	/**
	 * \brief q-intersection for circle estimation
	 *
	 */

  CtcQInterCircle2(int n, const Array<Ctc>& ctc_list,double*** measure, double epseq, int q, qintermethod QINTERPROJ, int K=1);
  void point_contract(IntervalVector& box,int iter); 
	double *** measure;
	double epseq;
	double compute_err_iter(IntervalVector & box, int iter);
	int activepoints_count(IntervalVector& box);        
	int activepoints_contract_count(IntervalVector& box);        

 protected : 
	  void fwd(IntervalVector & box, int iter);
	  void fwdbwd(IntervalVector & box, int iter);
	
	  Interval eval_ctc(IntervalVector & box, int iter , int k);
	Interval eval_dist(IntervalVector & box, int iter);

};


 class CtcQInterAffCircle2 : virtual public CtcQInter, public CtcQInterCircle2, public CtcQInterAff {

   public :
 CtcQInterAffCircle2(int n, const Array<Ctc>& ctc_list,double*** measure, double epseq, int q,  qintermethod QINTERPROJ, int K=1);

   protected : 
   double slope_compute(int iter, int j, int k, IntervalVector& box, AffineLin& af);
	 double err_compute( int iter, int k, IntervalVector& box,AffineLin& af);
	 double	valmean_compute(int iter, int k, IntervalVector& box, AffineLin& af);
	 void compute_affine_evaluation (int i , int iter, AffineLin& af, Interval & af2);

 };

} // end namespace ibex
#endif // __IBEX_CTC_Q_INTERCIRCLE2_H__
