#ifndef  M_H
#define M_H

#include <vector>
using namespace std;

template <typename T>
struct Matrix{
private:
     vector<T> data;
	 int col;
	 int row;

public:
	 Matrix(int new_col=0, int new_row=0){
	    col = new_col;
		row = new_row;

		data.resize(col*row);

	}


	T M(int x,int y){
		return data.at(x*row+y);
	}

	void set_M(int x,int y,T val){
	     data.at(x*row+y) = val;
	}


	int R(){
	    return row;
 	}

	int C(){
	    return col;
	}

	
	 ~Matrix() { data.clear();}
};

#endif