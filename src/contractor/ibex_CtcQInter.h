//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 30, 2012
// Last Update : Apr 30, 2012
//============================================================================

#ifndef __IBEX_CTC_Q_INTER_H__
#define __IBEX_CTC_Q_INTER_H__

#include "ibex_Ctc.h"


#include "ibex_IntervalMatrix.h"
#include "ibex_QInter.h"
#include <vector>
#include <list>


using namespace std;

namespace ibex {

/**
 * \ingroup contractor
 * \brief Q-intersection contractor.
 *
 */
class CtcQInter : public Ctc {
public:
	/**
	 * \brief q-intersection on a list of contractors.
	 *
	 * The list itself is not kept by reference.
	 */
  
  CtcQInter(int n, const Array<Ctc>& ctc_list, int q, qintermethod meth=QINTERPROJ, int kfun=1);
  ~CtcQInter();
	
	/**
	 * \brief Contract the box.
	 */
        virtual void contract(IntervalVector& box);
        virtual void point_contract(IntervalVector& box,int iter); 
	int ctc_contract(IntervalVector& box);

	int midbox_activepoints_number(IntervalVector& box);
	int midbox_activepoints_contract_count(IntervalVector& box);
	//	int activepoints_number(IntervalVector& box);
        virtual        int activepoints_count(IntervalVector& box);
        virtual int activepoints_contract_count(IntervalVector& box);
	double compute_err_sol(IntervalVector & box);
	virtual double compute_err_iter(IntervalVector & box, int iter);
         void updatepoints();
	void updateinterpoints(IntervalVector& box);
	/**
	 * List of contractors
	 */
	Array<Ctc> ctc_list;

	
	/**
	 * The number of contractors we have to intersect the
	 * result.
	 */

	int q;
	int qmax;
        
        qintermethod meth; // The Qinter method 
	list<int>* points;  // the list of current compatible measurements : initialized to all measures and managed by solver when the QInter constraint is used in a parameter estimation solver.
	
        bool points_to_delete;  //to manage shared pointer  points

       virtual  IntervalVector qinter_contract(int & p,int n0) ;
protected:

       int kfun;
        void init();
	IntervalMatrix* boxes; // store boxes for each contraction



};

 

} // end namespace ibex
#endif // __IBEX_CTC_Q_INTER_H__
