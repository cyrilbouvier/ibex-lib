//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor specialized for circle detection (uncertainty on the equation)
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Feb 6, 2015
// Last Update : Feb 6, 2015
//============================================================================

using namespace std;
#include "ibex_CtcQInterCircle2.h"
#include "ibex_Affine2.h"



namespace ibex {

  /* fwdbwd for the circle constraint  box[n]= (linfun[iter][0][k]- box[0])^2  + (linfun[iter][0][k]- box[1]) ^2
for n=2,3,4 */

Interval CtcQInterCircle2::eval_dist(IntervalVector & box, int iter)
{Interval evald =sqr (measure[iter][0][0] - box[0]) + sqr (measure[iter][1][0] - box[1]);
  return evald;
}

  Interval CtcQInterCircle2::eval_ctc(IntervalVector & box, int iter, int k) {
Interval eps (-epseq,epseq);
 Interval evald= eps + sqr (measure[iter][0][k] - box[0]) + sqr (measure[iter][1][k] - box[1]);
 return evald;
  }
  
  int CtcQInterCircle2::activepoints_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;
    while (iter != points->end()){
      IntervalVector box1=box;
      for (int k=0; k<kfun; k++)
	box1[nb_var-1]&= eval_ctc (box1,*iter,k);
      if (!box1[nb_var-1].is_empty()) p++;
      iter++;
    }
    return p;
 }


int CtcQInterCircle2::activepoints_contract_count(IntervalVector& box){
    int p=0;
    list<int>::iterator iter = points->begin() ;
    while (iter != points->end()){
      IntervalVector box1=box;
      for (int k=0; k<kfun; k++)
	box1[nb_var-1]&= eval_ctc (box1,*iter,k);
      if (!box1[nb_var-1].is_empty()) {p++; iter++;}
      else 
	iter=points->erase(iter);
    }
    return p;
}

  

double CtcQInterCircle2::compute_err_iter(IntervalVector & box, int iter){
    Interval eps (-epseq,epseq);
    IntervalVector mid (box.mid());
    Interval res=  eval_dist(mid,iter) - box[nb_var-1].mid();
    //    cout << " " << iter << " " << eval_dist(mid,iter) << " " <<box[nb_var-1].mid() << " res " <<   res <<  " " << res.mag() << endl;
    return res.mag();
  }

  void CtcQInterCircle2::fwd(IntervalVector & box, int iter)
  { 
    for (int k=0; k<kfun; k++)
      {  Interval eps (-epseq,epseq);
	Interval evald= eps+sqr (measure[iter][0][k] - box[0]) + sqr (measure[iter][1][k] - box[1]);

	box[nb_var-1]&=evald;
	if  (box[nb_var-1].is_empty())
	  {box.set_empty(); return;}
      }
  }


  void CtcQInterCircle2::fwdbwd(IntervalVector & box, int iter)
  { 
    for (int k=0; k<kfun; k++)
      {  Interval eps (-epseq,epseq);
	Interval evald= eps+sqr (measure[iter][0][k] - box[0]) + sqr (measure[iter][1][k] - box[1]);

	box[nb_var-1]&=evald;
	if  (box[nb_var-1].is_empty())
	  {box.set_empty();return;}

       	Interval yb= eps +box[nb_var-1] - sqr (measure[iter][1][k] - box[1]);
	Interval yb1 = sqrt(yb);
	Interval yb2 = (box[0] & (measure[iter][0][k] + yb1)) |( box[0] & (measure[iter][0][k] - yb1));
	box[0]=yb2;

	if  (box[0].is_empty())
	  {box.set_empty();return;}
	Interval xa= eps + box[nb_var-1] - sqr (measure[iter][0][k] - box[0]);
	Interval xa1 = sqrt(xa);
        Interval xa2 = (box[1] & (measure[iter][1][k] + xa1)) | ( box[1] & (measure[iter][1][k] - xa1));
	box[1]=xa2;

	if  (box[1].is_empty())
	  {box.set_empty();return;}
	  
      }
  }

 CtcQInterCircle2::CtcQInterCircle2(int n, const Array<Ctc>& ctc_list, double*** measure, double epseq, int q, qintermethod meth, int K ) : CtcQInter(n,ctc_list,q,meth,K),measure(measure),epseq(epseq) {;}

/* for circle constraint, CtcFwdBWd made by hand */
 void  CtcQInterCircle2::point_contract(IntervalVector& box, int iter) 
 { fwdbwd(box,iter);
   //fwd(box,iter);
}

  double CtcQInterAffCircle2::err_compute(int iter, int k, IntervalVector& box, AffineLin& af)
  {

    double err1=0;
    for (int j=0; j<nb_var-1; j++){
      err1 += ((std::pow(box[j].ub(),2) + std::pow(box[j].lb(),2))/2 - std::pow(box[j].mid(),2))/2;
    }
    return epseq + err1 ;
  }

  double CtcQInterAffCircle2::valmean_compute(int iter, int k, IntervalVector& box, AffineLin& af)
  { double valmean = 0;
    for (int j =0; j< nb_var-1; j++)
      {
	double err1 = ((std::pow(box[j].ub(),2) + std::pow(box[j].lb(),2))/2 - std::pow(box[j].mid(),2))/2;
	valmean+=std::pow(box[j].mid(),2) + std::pow(measure[iter][j][k],2)+ err1 -2 *measure[iter][j][k]* box[j].mid() ;

      }
    return valmean;
  }
 
  double CtcQInterAffCircle2::slope_compute(int iter, int j, int k, IntervalVector& box, AffineLin& af)
  {  return (
	     (std::pow(box[j].ub()-measure[iter][j][k],2) - std::pow(box[j].lb()-measure[iter][j][k],2)) / box[j].diam() 
	     );
}

  void CtcQInterAffCircle2::compute_affine_evaluation (int i , int iter, AffineLin& af, Interval & af2)
  {;}


 CtcQInterAffCircle2::CtcQInterAffCircle2(int n, const Array<Ctc>& ctc_list, double*** measure, double epseq, int q,  qintermethod meth, int K ) :
   CtcQInter(n,ctc_list,q,meth,K),
   CtcQInterCircle2(n,ctc_list,measure,epseq, q,meth,K),
   CtcQInterAff(n,ctc_list,q,meth,K) {;}


}
