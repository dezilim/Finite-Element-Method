#include "Elem.hpp"
#include "Node.hpp"
#include <vector>
#include <iostream>

using namespace std;

Elem::Elem(int elem_index){
    //node = NULL;
    this->elem_index_ = elem_index;
	//node = NULL;
	//elemNodes.resize(0);

}
void Elem::PrintNodes(){
    vector<Node>::iterator itn = elemNodes.begin();
	for (; itn != elemNodes.end(); itn++)
	{
		//cout << "Elem: printing nodes"  << "...\n" ;
		cout << *itn << endl;
	}

}



// // wanna print each 
// ostream& operator<<(ostream& os, const Elem& E) {
//     os << "--------------------------------------------------\n";

//     vector<Node>::iterator itn = E.elemNodes.begin();
// 	for (; itn != E.elemNodes.end(); itn++)
// 	{
// 		os << "Node in element"  << "..." ;
// 		os << *itn << endl;
// 	}

//     return os;	

	
// }