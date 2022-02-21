#ifndef Grid_hpp
#define Grid_hpp

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include <vector>
#include "Matrix.hpp"
#include "Node.hpp"
#include "Elem.hpp"
using namespace std;


class Grid{
    public: 
        // want to be able to access nodes in an elemenet in no particular order...
        vector<Elem> gridElems; 
        vector<Node> gridNodes;
        int gridSize;
        // pointer of nodes in this element. 
        //Node* node;

        //simple constructor 
        Grid(int Size);

        //void PrintNodes();
        

};

// ostream& operator<<(ostream& os, const Elem& E);


#endif 