#ifndef Node_hpp
#define Node_hpp

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include "Matrix.hpp"
using namespace std;


class Node{
    public: 
        Matrix position_;
        int elem_index_;

        //simple constructor 
        Node(Matrix& position, int elem_index);

        friend ostream& operator<<(ostream& os, const Node& N);

};


#endif 