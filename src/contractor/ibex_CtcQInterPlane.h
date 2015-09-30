//============================================================================

//                                  I B E X                                   
// File        : Q-intersection contractor specialized for plane detection 
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Feb 4, 2015
// Last Update : Feb 4, 2015
//============================================================================

#ifndef __IBEX_CTC_Q_INTERPLANE_H__
#define __IBEX_CTC_Q_INTERPLANE_H__

#include "ibex_CtcQInterAff.h"
#include "ibex_System.h"
#include "ibex_CtcPolytopeHull.h"
#include "ibex_LinearRelaxXTaylor.h"
namespace ibex {

class CtcQInterPlane : virtual public CtcQInter {
public:
	/**
	 * \brief q-intersection on a list of contractors.
	 *
	 * The list itself is not kept by reference.
	 */

  CtcQInterPlane(int n, const Array<Ctc>& ctc_list,double*** linfun, 
		   double epseq, int q,  qintermethod QINTERPROJ, int K=1);
  double compute_err_iter(IntervalVector & box, int iter);
  void point_contract(IntervalVector& box,int iter); 
  double *** linfun;
  double epseq;

  //  int points_count(IntervalVector& box);        
  //  int activepoints_count(IntervalVector& box);        
  //  int activepoints_contract_count(IntervalVector& box);        



 protected :
        void fwdbwd(IntervalVector & box, int iter);
        void fwd(IntervalVector & box, int iter);
        Interval eval_ctc(IntervalVector & box, int iter, int k);
	Interval eval_dist(IntervalVector & box, int iter, int k);

};

 class CtcQInterAffPlane :
   virtual public CtcQInter, 
   public CtcQInterPlane, public CtcQInterAff {

   public :

   CtcQInterAffPlane(int n, const Array<Ctc>& ctc_list,double*** linfun, 
		     double epseq, int q,  qintermethod QINTERPROJ, int K=1);

 protected :
   double valmean_compute(int iter, int i, IntervalVector& box, AffineLin& af);
   double slope_compute(int iter, int j, int i , IntervalVector& box,AffineLin& af);
   double err_compute( int iter, int k, IntervalVector& box,AffineLin& af);
   void compute_affine_evaluation( int i, int iter,  AffineLin& af , Interval& af2);
};

} // end namespace ibex
#endif // __IBEX_CTC_Q_INTERPLANE_H__
