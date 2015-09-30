//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor  with projection in the average gradient direction of the measurement equation
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 30, 2012
// Last Update : Feb 4, 2015
//============================================================================


// code ne marche plus en version générique : bug provenant de l'arithmétique affine ec cours de developpement Affine2
// marche avec l'ancienne version AffineLin


using namespace std;
#include "ibex_CtcQInterAff.h"
#include "ibex_QInter.h"  
#include "ibex_Affine2.h"
#include "ibex_AffineProjectionException.h"
#include<vector>
#include<list>
#include<algorithm>

namespace ibex {
  CtcQInterAff::CtcQInterAff(int n, const Array<Ctc>& ctc_list, int q, qintermethod meth, int K) : CtcQInter(n,ctc_list,q,meth,K) {
      for(int i=0;i<nb_var-1;i++) aff_dir.push_back(0);
      pmax=ctc_list.size();
      init_interaff();
    
  }

  CtcQInterAff::CtcQInterAff(int n, const Array<Ctc>& ctc_list,int pmax, int q,  qintermethod meth, int K) : CtcQInter(n,ctc_list,q,meth,K),pmax(pmax){
      for(int i=0;i<nb_var-1;i++) aff_dir.push_back(0);
      init_interaff();
    
  }

  CtcQInterAff::CtcQInterAff(int n, const Array<Ctc>& ctc_list, int q, Function*** funlist, qintermethod meth, int K) :  CtcQInter(n,ctc_list,q,meth,K), function_list(funlist) {

    for(int i=0;i<nb_var-1;i++) aff_dir.push_back(0);
    pmax=ctc_list.size(); 
    init_interaff();
  }
  
  void CtcQInterAff::init_interaff(){
    interaf=new Interval[pmax];
    alpha= new double**[pmax];
    for (int i=0; i<pmax; i++)
      {alpha[i] = new double*[nb_var];
	for (int j=0; j<nb_var; j++)
	  alpha[i][j]= new double[kfun];
      }
    valmean= new double*[pmax];
    for (int i=0; i<pmax; i++)
      valmean[i]=new double[kfun];
    err= new double*[pmax];
    for (int i=0; i<pmax; i++)
      err[i]=new double[kfun];

  }
  
  CtcQInterAff::~CtcQInterAff(){
      

    delete [] interaf;
    for (int i=0; i<pmax; i++)
      {for (int j=0; j<nb_var; j++)
	  delete [] alpha[i][j];
	delete [] alpha[i];
      }
    delete [] alpha;

    for (int i=0; i<pmax; i++)
      delete [] err[i];
    delete [] err;
    for (int i=0; i<pmax; i++)
      delete [] valmean[i];
    delete [] valmean;

      
  }
  

  // code ne marchant plus af.val(j+1) est faux
double CtcQInterAff::slope_compute(int iter, int j , int i, IntervalVector& box,AffineLin& af)
{  //cout << iter << " " << j << "  " << af.val(j+1) << endl;
  return 2*af.val(j+1)/((*boxes)[iter][j].diam());}
  double CtcQInterAff::valmean_compute(int iter, int i, IntervalVector& box, AffineLin& af)
  {return af.val(0);}
  double CtcQInterAff::err_compute( int iter, int k, IntervalVector& box,AffineLin& af)
  {return af.err();}
 
 // compute the boxes corresponding to a projection in a direction computed as the average (or median)  of the coefficients of the affine 
  // formula  giving the last variable of box in function of the  others
  // This formula is computed by an affine evaluation .
  // This is also useful by adding a nth+1 variable in case of an implicit constraint f(x1,...,xn) =0

  // code limite à p=50000 ??

  void CtcQInterAff::compute_affine_evaluation( int i, int iter,  AffineLin& af, Interval& af2) {
    af2= function_list[i][iter]->eval_affine2((*boxes)[iter],af);
  }


  void  CtcQInterAff::affine_projection(IntervalVector& box, int& pp)  {
 
    for (int j=0; j<nb_var-1; j++) aff_dir[j]=0;

    int p=ctc_list.size();
    int pmax=points->size();
    
    
    int k=0;
    AffineLin af;
    //    vector<double> affineslope[n-1];   // variante médiane
    list<int>::iterator iter = points->begin() ; 
    int p1=0;
    Interval af2 (0.0,0.0);
    while (iter != points->end()){
      //      (*boxes)[*iter] &=box;  // ne sert à rien.
      if ((*boxes)[*iter].is_empty() ) {iter++;p1++;continue;}
      int sortie_ok=1;
      for (int i=0; i<kfun;i++){
        
	//	cout  << *iter << " " << (*boxes)[*iter] << endl;
       
	compute_affine_evaluation( i, *iter, af, af2);
	  //	  cout << " afval " << af.val(0) << " " << af.val(1) << " " << af.val(2) << endl;
	  //          boxes[*iter][n-1]&=af2;  //inutile 
	
	  if (
	      af2.is_empty()
	      //	      boxes[*iter][n-1].is_empty()//inutile 
	      )
	    {//iter=points->erase(iter);
              if (!((*boxes)[*iter].is_empty()))
		{(*boxes)[*iter].set_empty();pp--;}
	      iter++;sortie_ok=0;p1++;
	      
              break;}  // evaluation affine vide exemple sqrt nb négatif
	
      
	  for (int j =0; j< nb_var-1; j++) {
	    if ((*boxes)[*iter][j].diam() < 1.e-20) throw AffineProjectionException();
	  //	  cout << "af.val"  << af.val(j+1) << endl;
	    alpha[p1][j][i]= slope_compute(*iter,j,i,(*boxes)[*iter],af);
	    //	    cout << " slope "  << *iter << " " << j << "  " <<  alpha[p1][j][i] <<  "  box " << (*boxes)[*iter] << endl;
	  //	    affineslope[j].push_back(alpha[*iter][j][i]);    // variante médiane
	  aff_dir[j] = (aff_dir[j]*k + alpha[p1][j][i])/(k+1);  // calcul incrémental de la pente moyenne
	  }
	valmean[p1][i]=valmean_compute(*iter,i, (*boxes)[*iter],af);
	err[p1][i]=err_compute(*iter,i, (*boxes)[*iter],af);
	k++;
      }
      if (sortie_ok) {iter++;p1++;}
    }
    

    // essai : utiliser pente médiane au lieu de pente moyenne : trop cher à cause du sort
    /*
    for (int j =0; j< n-1; j++)
      {sort(affineslope[j].begin(),affineslope[j].end());
	aff_dir[j]= affineslope[j][affineslope[j].size()/2];
      }
    */

    // calcul de l'intervalle dans la direction à projeter  aff_dir
    /*
    for (int j=0; j<n-1 ; j++)
      cout << j << " affdir  " << aff_dir[j] << endl;
    cout << "p " << pp << endl;
    */
    iter = points->begin() ; p1=0;
    while (iter != points->end()){
      if ((*boxes)[*iter].is_empty() ) {iter++;p1++;continue;}
      for (int i=0; i<kfun;i++){
	double errj=0;double valmeanj=0;
	  for (int j =0; j< nb_var-1; j++){
	    double alpha1= alpha[p1][j][i]-aff_dir[j];
	    errj +=fabs(alpha1)*(*boxes)[*iter][j].diam();
	    valmeanj +=aff_dir[j] *  (*boxes)[*iter][j].mid();
	    //	    cout << j << " " << "affdir " <<aff_dir[j] << " " << " midbox " << boxes[*iter][j].mid() << endl;
	  }
	  double valmean1= valmean[p1][i] - valmeanj;
	  //          cout << "  valmeani  " << valmean[*iter][i] << " valmeanj " << valmeanj << endl;
	  double err1=2*err[p1][i]+errj;
	  //	  cout << " valmean1 " << valmean1 << " err1 " << err1 << endl;
	  if (i==0)
	    interaf[p1]= Interval(valmean1-(err1/2), valmean1+(err1/2));
	  else
	    interaf[p1]&= Interval(valmean1-(err1/2), valmean1+(err1/2));
	  //	  cout << *iter << "  " << interaf[p1] << endl;
	  //	  cout << (*boxes)[*iter][0].lb() << " " << (*boxes)[*iter][0].ub() <<  " " << valmean1 << " " << alpha[*iter] <<  "  "  << err1 << endl;


      }
      if (interaf[p1].is_empty()){ // cas se produisant uniquement en cas de mesure multiple (kfun >1)
	//iter=points->erase(iter);
	if (!( (*boxes)[*iter].is_empty())) {(*boxes)[*iter].set_empty();pp--;}
	iter++;p1++;}
      else
	{(*boxes)[*iter][nb_var-1]=interaf[p1]; 
	  //	  cout << boxes[*iter][n-1]  << endl;
	  iter++;p1++;
	}
    }
    /*
    
    */
  }

void CtcQInterAff::contract(IntervalVector& box) {
  //  cout << " ctc aff " <<  " " << side_effects()  << endl;
  CtcQInter::contract(box);

  if (box.is_empty()) return;
  int p=points->size();

  int  p0=p;

  // la projection affine n'est appelée que par la contraction directe (pas dans 3bcid)
  // 
  if (
        side_effects ()
      //  && p < 2*q
      ) {

    //     cout << " avant proj affine " << " p " << p <<  "  box " << box << endl;
    try {
    affine_projection (box,p);

    // cout << " apres proj affine " << " p " << p <<  "  box " << box << endl;
    //    p0=p;
    //    int res2=1;


    // while(res2){
    // res2=0;
    // p updated by previous call to qinter_projf
    //    cout << " iter proj affine " << p << endl;

    int n2=nb_var-1;
    IntervalVector box1 = qinter_contract(p,n2);
    //        cout << " apres 2me proj" << " p " << p << " " << box << "   " << box1 << endl;
      if (box1.is_empty()) {box.set_empty();return;}
      for (int i=n2; i<nb_var;i++) {if (box1[i].is_empty()) {box.set_empty();return;}}
      if (p<q) {box.set_empty();return;}
    //    cout << box1[0].lb() << " " << box1[0].ub() << "  " << box1[1].mid() <<  " "<<  0 << " " << box1[1].diam() << endl;
      Interval res1;
      Interval resi(0,0);
      // the n-1 first directions
      for (int i=n2; i<nb_var-1;i++)
	box[i]&=box1[i];
      // projecting the result of the new nth direction (the mean gradient) on the old nth direction 
      for (int i=0;i<nb_var-1;i++)
	resi+=aff_dir[i]*box[i];
      res1=box1[nb_var-1] + resi;
      //      IntervalVector initbox(n);
      //      initbox=box;
      //    cout << box[0].lb() << " " << box[0].ub() << "  " << res1.mid() <<  " "<<  0 << " " << res1.diam() << endl;
      //    cout << " box avant " << box << " reduc " << res1 <<  " ";
      
      box[nb_var-1]&=res1;
      for (int i=n2; i<nb_var;i++) {
	if (box[i].is_empty())  { box.set_empty();return;}
      }
      
      //    cout << box << " " << box1 << " " << box2 << endl;
      //      if (initbox.rel_distance(box) > 0.1) res2=1;
      
      // updating points only in case of direct contraction call (not in  cid) using  side_effects
      

      if ( p<p0 &&   side_effects() ) 
	//	updateinterpoints(box1);  // inutile ?
	updatepoints();

      //      cout << " fin contract proj affine " << endl;
    }
    catch (AffineProjectionException err) {;}
      //      if (points->size()<q) {box.set_empty();}
  }
      
      
  

}


}
