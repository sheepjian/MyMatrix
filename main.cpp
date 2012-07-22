#include "Matrix_data.h"
#include "h_func.cpp"
#include <iostream>

using namespace std;


void main(){

      Matrix<float> A = init_file<float>(5,5,"matrix.txt");
	  output<float>(A,3);

	 // transpose<float>(A);
	 // output<float>(A,7);

	  cout<<endl;
	  transpose<float>(A);

	  cout<<endl<<"the det of A"<<det<float>(A)<<endl;
	
	  Matrix<float> B(5,5);

	  B = inv<float>(A);

	  output<float>(B,3);
	  system("pause");

}