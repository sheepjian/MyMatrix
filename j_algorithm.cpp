#include <vector>
using namespace std;

template <typename T>
int MergeSort(vector<T> &data,int first, int last)
{
	int change = 0;

	if( first < last){
	int mid =(first+last)/2;

	change += MergeSort(data,first,mid);

	change += MergeSort(data,mid+1,last);

	change += Merge(data,first,mid,last);
	return change;
	} else {
	   return 0;
	}	
}


template <typename T>
int Merge(vector<T> data,int first, int mid, int last){
	
    int change =0;
	int i=0;

	vector<T> tmp(data);
	
	int b_A=first, end_A = mid, b_B = mid+1, end_B = last;

	while(b_A<=end_A && b_B <= end_B){
		if(data.at(b_A)<data.at(b_B)){
			data.at(i+first) = tmp.at(b_A);
			b_A++;
			i++;
		} else{
		    data.at(i+first) = tmp.at(b_B);
			b_B++;
			i++;
			change = end_A - b_A+1;
		}
	}

	while(b_A<=end_A){
	    data.at(i+first) = tmp.at(b_A);
	    b_A++;
		i++;
	}

	while(b_B<=end_B){
	    data.at(i+first) = tmp.at(b_B);
	    b_B++;
		i++;
	}
	return change;
}
