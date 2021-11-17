#ifndef Elem_hpp
#define Elem_hpp

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include <vector>
#include "Matrix.hpp"
#include "Node.hpp"
using namespace std;


class Elem{
    public: 
        // want to be able to access nodes in an elemenet in no particular order...
        vector<Node> elemNodes;
        int elem_index_;
        // pointer of nodes in this element. 
        //Node* node;

        //simple constructor 
        Elem(int elem_index);

        void PrintNodes();
        

};

// ostream& operator<<(ostream& os, const Elem& E);


#endif 