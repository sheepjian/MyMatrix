#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Matrix_data.h"
#include "j_algorithm.cpp"
using namespace std;

template <typename T>
void output(Matrix<T> A,int acu=2){
	for(int i = 0; i< A.R(); i++){
		for(int j = 0; j< A.C(); j++){
			cout<<setw(acu+2)<<setprecision(acu)<<setiosflags(ios::left)<<A.M(i,j)<<" ";	
		}
		cout<<endl;	
	}
}


template <typename T>
Matrix<T> init_file(int r, int c, char *file_name){
	Matrix<T> A(r,c);
	T tmp;
	fstream file;
	file.open(file_name,ios::in);
	for(int i = 0; i< A.R(); i++){
		for(int j = 0; j< A.C(); j++){
			file>>tmp;
			A.set_M(i,j,tmp);
		}
	}
	file.close();
	return A;
}


template <typename T>
void transpose(Matrix<T> &A){
    T tmp;
	for(int i = 0; i< A.R(); i++){
		for(int j = i+1; j< A.C(); j++){
			tmp = A.M(i,j);
			A.set_M(i,j, A.M(j,i));
			A.set_M(j,i,tmp);
		}
	}
}


#ifndef max      
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define TINY 1.0e-8

template <typename T>
int LU_decmp(Matrix<T> &A,vector<int> &H, vector<int> &V){
    int i,j,k;

	int rev = 0;
	
	int n = A.R();

    int pi,pj;
    int tmp1,tmp2;
	T max;
    T ftmp,dtmp;
    for (k=0;k<n;k++){
        pi=-1,pj=-1,max=0;
        //find pivot in submatrix a(k:n,k:n)
        for (i=k;i<n;i++) {
            for (j=k;j<n;j++) {
				if (fabs(A.M(i,j))>max){
                    max = fabs(A.M(i,j));
                    pi=i;
                    pj=j;
                }
            }
        }
        //Swap Row
		tmp1 = H.at(k);
		tmp2 = H.at(pi);
		H.at(k)=tmp2;
		H.at(pi)=tmp1;
		if(pi!=k)
		rev++;

        for (j=0;j<n;j++){
            ftmp=A.M(k,j);
			dtmp=A.M(pi,j);

			A.set_M(k,j,dtmp);
			A.set_M(pi,j,ftmp);
        }
        //Swap Col
		tmp1 = V.at(k);
		tmp2 = V.at(pj);
		V.at(k)=tmp2;
		V.at(pj)=tmp1;
		if(pj!=k)
		rev++;

        for (i=0;i<n;i++){
            ftmp=A.M(i,k);
			dtmp=A.M(i,pj);

			A.set_M(i,k,dtmp);
			A.set_M(i,pj,ftmp);
        }
        //END PIVOT

        //check pivot size and decompose
        if ((fabs(A.M(k,k))>TINY)){
            for (i=k+1;i<n;i++){
                //Column normalisation
                ftmp=A.M(i,k)/A.M(k,k);
				A.set_M(i,k,ftmp);
                for (j=k+1;j<n;j++){
                    //a(ik)*a(kj) subtracted from lower right submatrix elements
                   dtmp =  A.M(i,j)-(ftmp*A.M(k,j));
				   A.set_M(i,j,dtmp);
                }
            }
        }
	}

	return rev;
}


template <typename T>
T det(Matrix<T> A){
	// only square matrix has determinant

     Matrix<T> tmp = A;
	 int dim = tmp.R();
	 int rev;

	 vector<int> v_t ;
	 vector<int> h_t;

	 for(int i=0; i<dim; i++){
		 v_t.push_back(i);
		 h_t.push_back(i);
	 }

	 rev=LU_decmp(tmp,h_t,v_t);

	 T prod = 1.0;
	 for(int i=0; i < dim; i++){
		 prod*= tmp.M(i,i);
	 }

	 int pre;
	 if(rev%2==0){
	      pre = 1;
	 }else{
	      pre = -1;
	 }
	 return prod*pre;
}


template <typename T>
Matrix<T> inv(Matrix<T> A){
	int dim = A.C();

	Matrix<T> inver(dim,dim);

	Matrix<T> tmp(dim-1,dim-1);

	T det_m;

	int i,j,p=0,q=0,pre;

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(p=0;p<i;p++){
			    for(q=0;q<j;q++){
					tmp.set_M(p,q,A.M(p,q));
			    }
			}

			for(p=i+1;p<dim;p++){
			    for(q=0;q<j;q++){
					tmp.set_M(p-1,q,A.M(p,q));
			    }
			}
		    
			for(p=0;p<i;p++){	
				for(q=j+1;q<dim;q++){		
					tmp.set_M(p,q-1,A.M(p,q));
				}
			}

			for(p=i+1;p<dim;p++){
			    for(q=j+1;q<dim;q++){
					tmp.set_M(p-1,q-1,A.M(p,q));
			    }
			}

			det_m = det<T>(tmp);

			 if((i+j)%2==0){
				  pre = 1;
			 }else{
				  pre = -1;
			 }
			 det_m*=pre;
			 inver.set_M(i,j,det_m);
	    }
	}

	T det_A =  det<T>(A);
	T t_inv;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			t_inv = inver.M(i,j);
			inver.set_M(i,j,t_inv/det_A);
		}
	}
	return inver;
}
