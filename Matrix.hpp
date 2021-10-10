#ifndef Matrix_hpp
#define Matrix_hpp

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include "Vector.hpp"
using namespace std;


class Matrix {
    public: 
        int rows_, cols_;
        double** val;

        // create matrix with constant value, if on value specified, default is 0 matrix 
        Matrix(int rows_=2, int cols_=2, double a=0.0);
        Matrix(const Matrix& Q);


        // += and -= only if matrices are same dim. 
        Matrix& operator+=(const Matrix& Q);
        Matrix& operator-=(const Matrix& Q);

        // multiplication by matrix must have compatible dimensions
        Matrix& operator*=(const Matrix& Q);
        Matrix& operator*=(double a);

        

        
        Matrix transpose();
        Matrix fillSym();
        void alloc();
        


};
// A*B  = C ; a*A ; A*a
Matrix operator*(const Matrix&, const Matrix&);  // CHECK COMPATIBILITY
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);



// CHECK SAME DIM 
Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);

ostream& operator<<(ostream& os, const Matrix& M);
bool checkCompatible(const Matrix&, const Matrix&); // check if compatible for multiplication
bool checkSameDim(const Matrix&, const Matrix&); // check if both matrices have exactly same dimensions

Matrix identity(int size);
Matrix colVect(int size, int i); // create col vect with 1 at the ith location. Check 0 <= i < size
Matrix rowVect(int size, int j); // create col vect with 1 at the ith location. Check 0 <= i < size
Matrix joinLR(const Matrix&, const Matrix&);
Matrix joinUD(const Matrix&, const Matrix&);



#endif 