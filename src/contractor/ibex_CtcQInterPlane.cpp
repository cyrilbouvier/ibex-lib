//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor specialized for plane detection 
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Feb 6, 2015
// Last Update : Feb 6, 2015
//============================================================================

using namespace std;
#include "ibex_CtcQInterPlane.h"
#include "ibex_Affine2.h"

// can be used with 2 modelizations that have a linear constraint : 
// the classic one (ax+by+cz+d=eps) with a^2+b^2+c^2=0   nb_var=4
// and with the modelization   (1-b-c)x + by + cz + d  explained in ictai paper     nb_var=3

namespace ibex {

 
  /*
Interval CtcQInterPlane::fwd(IntervalVector & box, int iter)
  {
    Interval eps (-epseq,epseq);
    Interval evald= eps+linfun[iter][0][0];
    
    for (int j=0;j<nb_var-1;j++)
      evald+=linfun[iter][j+1][0]*box[j];
    return evald;
  }
  
  */
 // the evaluation of the iterth equation with uncertainties
  // to check if the iterth measurement is valid (0 must be in the returned interval)

  
  Interval CtcQInterPlane::eval_ctc(IntervalVector & box, int iter, int k){
    Interval eps (-epseq,epseq);
    return (eval_dist(box,iter,k) + eps);
  }

  
  // 
  Interval CtcQInterPlane::eval_dist(IntervalVector & box, int iter, int k)
  {
    Interval evald= linfun[iter][0][k];
    for (int j=0;j<nb_var-1;j++)
      evald+=linfun[iter][j+1][k]*box[j];
    return evald;
  }



  // the individual error at the middle of box for the membership equation of the iterth point to be on the
  // plane 
  // property : it is less than epseq  :  only implemented for kfun=1 
  double CtcQInterPlane::compute_err_iter(IntervalVector & box, int iter){
    IntervalVector mid (box.mid());
    Interval res=  eval_dist(mid,iter,0) - box[nb_var-1].mid();
    //    cout << " " << iter << " " << res << endl;
    return res.mag();
  }


  /* fwd  evaluation of the linear constraint  box[n]= linfun[iter][1][k]* box[0] + linfun[iter][2][k]* box[1]  + ...  
for n=2,3,4 */
  //  void CtcQInterPlane::fwd(IntervalVector & box, int iter)

  // for Plane constraint, simple forward evaluation of the constraint (fwdbwd too expensive)  
  void CtcQInterPlane::point_contract(IntervalVector & box, int iter)
  {
      for (int k=0; k<kfun; k++)
	{ 
	  box[nb_var-1]&=eval_ctc(box, iter, k);
	  if  (box[nb_var-1].is_empty())
	    {box.set_empty();return;}
	}
  }
  

/* fwdbwd for the linear constraint  box[n]= linfun[iter][1][k]* box[0] + linfun[iter][2][k]* box[1]  + ...  
for n=2,3,4 */

 void CtcQInterPlane::fwdbwd(IntervalVector & box, int iter)
  { //cout << " contract planeprojf " << endl;
    int index=1;
    //    IntervalVector old_box(n);
    Interval eps (-epseq,epseq);
    while (index) {
      //old_box=box;

      for (int k=0; k<kfun; k++)
	{  


	  //	Interval evald= eps+linfun[iter][0][k];

	  //	for (int j=0;j<nb_var-1;j++)
	  //	  evald+=linfun[iter][j+1][k]*box[j];
	  box[nb_var-1]&=eval_ctc(box,iter,k);
	if  (box[nb_var-1].is_empty())
	  {box.set_empty();return;}

        if (nb_var==2)
	  {box[0]&=(eps +box[nb_var-1]- linfun[iter][0][k])/linfun[iter][1][k];
	    if  (box[0].is_empty())
	      {box.set_empty();return;}
	  }
	else if (nb_var==3)
	  { 
	    box[0]&=(eps + box[nb_var-1]-  linfun[iter][0][k] - linfun[iter][2][k]*box[1])/linfun[iter][1][k];
	    if  (box[0].is_empty())
	  {box.set_empty();return;}

	    box[1]&=(eps + box[nb_var-1]-  linfun[iter][0][k] - linfun[iter][1][k]*box[0])/linfun[iter][2][k];
	    if  (box[1].is_empty())
	      {box.set_empty();return;}
	  }
	else if (nb_var==4)
	  {box[0]&=(eps +box[nb_var-1] -linfun[iter][0][k]- linfun[iter][2][k]*box[1] -linfun[iter][3][k]*box[2] )/linfun[iter][1][k];
	    if  (box[0].is_empty())
	      {box.set_empty();return;}

	    box[1]&=(eps+ box[nb_var-1] -linfun[iter][0][k]- linfun[iter][1][k]*box[0] - linfun[iter][3][k]*box[2])/linfun[iter][2][k];
	    if  (box[1].is_empty())
	      {box.set_empty();return;}
	  
	    box[2]&=(eps +box[nb_var-1] -linfun[iter][0][k]- linfun[iter][1][k]*box[0] - linfun[iter][2][k]*box[1])/linfun[iter][3][k];
	    if  (box[2].is_empty())
	      {box.set_empty();return;}
	  }
	}
      index=0;
      /*
      if (kfun==1) 
	index=0;
      else if (old_box.rel_distance(box) < 0.1)
	index=0;
      */
      //      cout << " box " << box << endl;
    }
  }
 
 



  CtcQInterPlane::CtcQInterPlane(int n, const Array<Ctc>& ctc_list, double*** linfun, 
				 double epseq, int q, qintermethod meth, int K ) : 
    CtcQInter(n,ctc_list,q,meth,K),linfun(linfun),
       epseq(epseq) {
  }
  /*
  int CtcQInterPlane::points_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;
    int ll= points->size();
    while (iter != points->end()){
      if (!(fwd (box,*iter)&box[nb_var-1]).is_empty()) p++;
      if (p + ll-*iter < q) break; // counting p only if it can be greater than q
      iter++;

    }
    
    return p;
 }
  */


  /*

// no need to redefine activepoints_count and activepoints_contract_count
// 


int CtcQInterPlane::activepoints_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;
    //    int ll= points->size();
    while (iter != points->end()){
      if (!(eval_ctc (box,*iter)&box[nb_var-1]).is_empty()) p++;
      //      if (p + ll-*iter < q) break; // counting p only if it can be greater than q
      iter++;

    }
    
    return p;
 }
 

int CtcQInterPlane::activepoints_contract_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;
    while (iter != points->end()){
      if (!(eval_ctc (box,*iter)&box[nb_var-1]).is_empty()){
	p++;
	iter++;}
      else
	iter=points->erase(iter);
    }
    
    return p;
 }
  */
 

/* for linear constraint, CtcFwdBWd made by hand */

/*
 void  CtcQInterPlane::point_contract(IntervalVector& box, int iter) 
 {
       fwdbwd(box,iter);
       //  fwd(box,iter);   // more efficient ->  fix point useless on CtcQInter constraint
 }
*/

  /* computations for the affine projection : see  CtcQInterAff.cpp */


CtcQInterAffPlane::CtcQInterAffPlane(int n, const Array<Ctc>& ctc_list, double*** linfun, 
				 double epseq, int q,  qintermethod meth, int K ) : 
  CtcQInter(n,ctc_list,q,meth,K),
  CtcQInterPlane (n,ctc_list,linfun,epseq,q,meth,K),
  CtcQInterAff(n,ctc_list,q,meth,K)
        {
  }

 double CtcQInterAffPlane::err_compute( int iter, int k, IntervalVector& box,AffineLin& af)
  {return epseq;}

  double CtcQInterAffPlane::valmean_compute(int iter, int k, IntervalVector& box, AffineLin& af)
  { double valmean = linfun[iter][0][k];
    for (int j =0; j< nb_var-1; j++)
      valmean+=linfun[iter][j+1][k]*box[j].mid();
    return valmean;}
 



  double  CtcQInterAffPlane::slope_compute(int iter, int j , int k , IntervalVector& box,AffineLin& af2)
  {    return linfun[iter][j+1][k];
  }

  
  void CtcQInterAffPlane::compute_affine_evaluation( int i, int iter,  AffineLin& af, Interval& af2) {
   ; }
  
}

