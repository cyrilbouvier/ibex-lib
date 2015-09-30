//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 30, 2012
// Last Update : Feb 4, 2015
//============================================================================

using namespace std;
#include "ibex_CtcQInter.h"
#include "ibex_QInter.h"

#include<vector>
#include<algorithm>

namespace ibex {

  CtcQInter::CtcQInter(int n, const Array<Ctc>& ctc_list, int q, qintermethod meth, int kfun) : Ctc(n),ctc_list(ctc_list),  q(q), meth(meth), kfun(kfun){
    init();
  }
  
  

  void CtcQInter::init ()
  {
    points= new list<int> (); 
    boxes= new IntervalMatrix(ctc_list.size(),nb_var);
    for(int i=0;i<ctc_list.size();i++) {points->push_back(i);}
    _side_effects=false;
    points_to_delete=true;
  }

  CtcQInter::~CtcQInter()
   {
     //    cout << " points to delete " << points_to_delete << endl;
     if (points_to_delete) {points->clear(); delete points;}
     delete boxes;
   }

  // this function is called by QInterSolver when a solution is found
  // points are the validpoints
  // the box is the box around the point with the valid measurements
  // to better discriminate the solutions with the same number of valid points , a sum of squares is used
  double CtcQInter::compute_err_sol(IntervalVector & box){
    double err=0;
    list<int>::iterator iter= points->begin();
    while (iter != points->end()){err += std::pow(compute_err_iter(box,*iter),2); iter++;}
    //    cout << " total err " << err<< endl;
    return std::sqrt(err);
  }

  // Individual error : to be implemented for specific problems
  // for the moment implemented in CtcQInterPlane, CtcQInterCircle et CtcQInterCircle2
  double CtcQInter::compute_err_iter(IntervalVector & box, int iter)
  {return 0;}



   

  /* erase from the points list the boxes  becoming empty after intersecting them with box   */
  void CtcQInter::updateinterpoints(IntervalVector & box)
  {
    list<int>::iterator iter= points->begin();
    while (iter != points->end()){
      (*boxes)[*iter]&=box;
      if ((*boxes)[*iter].is_empty())
	{iter=points->erase(iter);}
      else
	iter++;
    }
  }

  /* erase the empty boxes from the points list  */
 void CtcQInter::updatepoints()
  {
    list<int>::iterator iter= points->begin();
    while (iter != points->end()){
      if ((*boxes)[*iter].is_empty())
	{iter=points->erase(iter);}
      else
	iter++;
    }
  }

  /* set the measures to box, contract them ; does not update the points list ; return the number of active measures in box*/
 int CtcQInter::activepoints_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;

    while (iter != points->end()){
      IntervalVector box1=box;
      point_contract(box1,*iter);
      if (!(box1.is_empty())) p++;
      iter++;
    }
    return p;
 }

 /* set the measures to box, contract them ; update the points list ; return the number of active measures in box*/
 int CtcQInter::activepoints_contract_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;

    while (iter != points->end()){
      IntervalVector box1=box;
      point_contract(box1,*iter); 
      if (box1.is_empty())
	iter=points->erase(iter);
      else {
	p++; iter++;
      }
    }
    
    return p;
 }

 



  // returns the number of measurements valid at the middle of the box.
  int CtcQInter::midbox_activepoints_number(IntervalVector& box){
    IntervalVector mid (box.mid());
    return (activepoints_count(mid));
  }


// returns the number of measurements valid at the middle of the box and removes the invalid ones from points
// This function should only be called when a solution is found.
 int CtcQInter::midbox_activepoints_contract_count(IntervalVector& box){
    IntervalVector mid (box.mid());
    return (activepoints_contract_count(mid));
  }


  

 

  // call of the qinter algorithm ; n0 is the first variable to treated by the algorithm 
  
  IntervalVector CtcQInter::qinter_contract(int& p,   int n0) {
    int qproj=points->size();
    IntervalVector box= (*boxes)[0];
    switch (meth)

      {case QINTERPROJ : box=  qinter_projf(*boxes,q,qproj, p,  points,n0); break;
      case QINTERCORE : return  qinter_coref(*boxes,q,p, points,n0);
      case QINTERFULL : return  qinter2(*boxes,q,p, points);
      case QINTERGRID : return qinter(*boxes,q,p,points);
      default : ibex_error("Qinter contract : impossible case");
      }
    
    if (qproj < qmax && side_effects() ) {
      //if (qmax1 < p) cout << "qmax reduction " << p << " " << qmax1; 
      qmax=qproj;}  // during the projection the value of qmax has been decreased
    
    return box;
  }

  void CtcQInter::contract(IntervalVector& box) {
// p is the number of non empty contract results
//    cout << " debut contract " << box << " " << side_effects() << endl;
    int p=ctc_contract(box);

    int p0=p;
    //  cout << " apres contract " << " p " << p << "  " << box << endl;
    if (p<q)  {box.set_empty();return;}
 
    box=qinter_contract(p, 0);
    if (box.is_empty()) return;
    if (p<q) {box.set_empty();return;}
  // updating points only in case of direct contraction call (not in  cid) using  side_effects
  
    //     cout << " side effects " << _side_effects << endl;
    
    if (p<p0 &&  _side_effects) 
      //   updateinterpoints(box);  // mise à jour trop chère 
    updatepoints();
    
    //    ctc_contract(box,1);      trop cher

    
  }

  // contract box with the iterth  measurement : uses the contractor of the measurement  (virtual function , can be redefined if a more simple computation exists
 void  CtcQInter::point_contract(IntervalVector& box, int iter)
  {  //cout << "iter " << iter << " box " << box  << endl;
    ctc_list[iter].contract(box);}



  // contracts all the measurements with box , removes the non compatible from points and returns the number of compatible measurements
  // if ind=0 , the initial box for the measurement is the box (the current box) 
  // if ind=1,  it is intersected with box (in a fixpoint operator)

  /*
int  CtcQInter::ctc_contract(IntervalVector& box, int ind)
{int p=0;  // the number of compatible measurements (no empty box after contraction)
  list<int>::iterator iter = points->begin() ;
  //  cout << " nb points  " << points->size() << endl;
  while (iter != points->end())

    { 

      try {

	if (ind ==0)(*boxes)[*iter]=box;
	else  (*boxes)[*iter]&=box;

	(*boxes)[*iter]&=box;
       	point_contract((*boxes)[*iter],*iter);
	p++;
	iter++;
      }
      catch(EmptyBoxException&) {
	
	assert((*boxes)[*iter].is_empty());
	if (side_effects)    // use of side_effects for not updating points in case of 3bcid
	  iter=points->erase(iter);
	else 
	  iter++;
      }

    }
  return p;
}
  */
  // returns the number of compatible measurements (no empty box after contraction) : the measurements with empty boxes are removed from points.
int  CtcQInter::ctc_contract(IntervalVector& box)
{int p=0; 
  list<int>::iterator iter = points->begin() ;
  //  cout << " nb points  " << points->size() << endl;
  while (iter != points->end())

    { 

      

	(*boxes)[*iter]=box;
	//cout << *iter <<  " box " << box;
       	point_contract((*boxes)[*iter],*iter);  // contraction of the measurement *iter
	//	if (!((*boxes)[*iter][0] == box[0] )) cout << *iter << " " << (*boxes)[*iter][0] <<  " " << box[0] << endl;
	//	cout << 	(*boxes)[*iter] << endl;
	if (! (*boxes)[*iter].is_empty()) {
	    p++;
	    iter++;
	  }
	else
      
	  {
	
	//	cout << " empty box " << endl;
	    if (side_effects())    // use of side_effects for not updating points in case of 3bcid
	     iter=points->erase(iter);
	   else 
	     iter++;
	  }
    }
  return p;
}

 

} // end namespace ibex










