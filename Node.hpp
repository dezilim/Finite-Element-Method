#ifndef Node_hpp
#define Node_hpp

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include <Eigen/Dense>
#include <array>
#include <math.h>

using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;


class Node{
    public: 
        Vector3d position_;
        int elem_index_;

        //simple constructor 
        Node(Vector3d& position, int elem_index);

        friend ostream& operator<<(ostream& os, const Node& N);

};


#endif 