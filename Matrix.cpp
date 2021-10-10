#include "Matrix.hpp"
// functions from the class
Matrix::Matrix(int rows_, int cols_, double a) 
{
	this->rows_ = rows_;
	this->cols_ = cols_;
    alloc();
    for (int i = 0; i < rows_; ++i) 
    {
        for (int j = 0; j < cols_; ++j) 
        {
            val[i][j] = a;
        }
    }
}

Matrix::Matrix(const Matrix& Q)
{
    this->rows_ = Q.rows_;
    this->cols_ = Q.cols_;
    alloc();
    for (int i = 0; i < rows_; ++i) 
    {
        for (int j = 0; j < cols_; ++j) 
        {
            val[i][j] = Q.val[i][j];
        }
    }
} 

//  IN CLASS OPERATOR OVERLOADS 
// ----------------------------

Matrix& Matrix::operator+=(const Matrix& Q)
{
    if (!checkSameDim(*this, Q)){
        cout << "!! Dimensions not equal for += " << endl;
        exit(-1); 
    }
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            val[i][j] += Q.val[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator-=(const Matrix& Q)
{
    if (!checkSameDim(*this, Q)){
        cout << "!! Dimensions not equal for += " << endl;
        exit(-1); 
    }
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            val[i][j] -= Q.val[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(double a)
{
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            val[i][j] *= a;
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& Q)
{
    if (!checkCompatible(*this, Q)){
        cout << "!! Dimensions incompatible for *= " << endl;
        exit(-1); 
    }
    Matrix TMP(rows_, Q.cols_);
    for (int i = 0; i < TMP.rows_; ++i) {
        for (int j = 0; j < TMP.cols_; ++j) {
            for (int k = 0; k < cols_; ++k) {
                TMP.val[i][j] += (val[i][k] * Q.val[k][j]);
            }
        }
    }
    return (*this = TMP);
}

Matrix Matrix::transpose()
{
    // create temporary matrix 
    Matrix TMP(cols_, rows_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            TMP.val[j][i] = val[i][j]; // val refers to the reference matrix this function acts on
        }
    }
    return TMP;
}

Matrix Matrix:: fillSym()
{
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < i; j++) {
            val[i][j] = val[j][i];
        }
    }
    return *this;
}
void Matrix::alloc()
{
    val = new double*[rows_];
    for (int i = 0; i < rows_; i ++)
    {
        val[i] = new double[cols_];
    }
}

//-------------- fuctions ouside the class --------------------
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//-------------------------------------------------------------
ostream& operator<<(ostream& os, const Matrix& M) {
	os << "\n";
	for(int i=0; i<M.rows_; i++){
        for(int j=0; j<M.cols_; j++){
            os << M.val[i][j] << " ";     
        }
        os << "\n" ;
    } 
	return os;
}

bool checkCompatible(const Matrix& P, const Matrix& Q){
    if(P.cols_ == Q.rows_){
        return (true);
    }
    return (false);
}

bool checkSameDim(const Matrix& P, const Matrix& Q){
    if((P.rows_ == Q.rows_) && (P.cols_ == Q.cols_)){
        return (true);
    }
    return (false);
}

Matrix operator* (const Matrix& P, const Matrix& Q) {
    if (!checkCompatible(P, Q)){
        cout << "!! Dimensions incompatible for *= " << endl;
        exit(-1); 
    }
	Matrix TMP(P);
	return (TMP*= Q);
}

Matrix operator*(const Matrix& P, double a){
    Matrix TMP(P);
    return (TMP*=a);
}
Matrix operator*(double a, const Matrix& P){
    Matrix TMP(P);
    return (TMP*=a);
}



Matrix operator+(const Matrix& P, const Matrix& Q)
{
    if (!checkSameDim(P, Q)){
        cout << "!! Dimensions not equal for += " << endl;
        exit(-1); 
    }
    Matrix TMP(P);
    return (TMP += Q);
}

Matrix operator-(const Matrix& P, const Matrix& Q)
{
    if (!checkSameDim(P, Q)){
        cout << "!! Dimensions not equal for += " << endl;
        exit(-1); 
    }
    Matrix TMP(P);
    return (TMP -= Q);
}

Matrix identity(int size)
{
    Matrix TMP(size, size);
    for (int i = 0; i < TMP.rows_; i++) {
        TMP.val[i][i] = 1;
    }
    return TMP;
}

Matrix colVect(int size, int i)
{
    Matrix TMP(size, 1);
    for (int k = 0; k < TMP.rows_; k++) {
        TMP.val[k][0] = 0;
    }
    TMP.val[i-1][0] = 1;
    return TMP;
}


Matrix rowVect(int size, int j)
{
    Matrix TMP(1, size);
    for (int k = 0; k < TMP.cols_; k++) {
        TMP.val[0][k] = 0;
    }
    TMP.val[0][j-1] = 1;
    return TMP;
}

Matrix joinLR(const Matrix& P, const Matrix& Q)
{
    if (P.rows_ != Q.rows_ ){
        cout << "!! Rows not equal for combining LR " << endl;
        exit(-1); 
    }
    Matrix TMP(P.rows_, P.cols_+Q.cols_);
    for (int i = 0; i < P.rows_; i++) {
        for (int j = 0; j < P.cols_ + Q.cols_; j++){
            if (j < P.cols_) {
                TMP.val[i][j] = P.val[i][j];
            }
            else {
                TMP.val[i][j] = Q.val[i][j-P.cols_];
            }
            
        }
    }
    return TMP;

}


Matrix joinUD(const Matrix& P, const Matrix& Q)
{
    if (P.cols_ != Q.cols_ ){
        cout << "!! Columns not equal for combining UD " << endl;
        exit(-1); 
    }
    Matrix TMP(P.rows_ + Q.rows_, P.cols_);
    for (int j = 0; j < P.cols_; j++) {
        for (int i = 0; i < P.rows_ + Q.rows_; i++){
            if (i < P.rows_) {
                TMP.val[i][j] = P.val[i][j];
            }
            else {
                TMP.val[i][j] = Q.val[i- P.rows_][j];
            }
            
        }
    }
    return TMP;

}
