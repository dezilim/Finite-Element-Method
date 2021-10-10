//#include <ctime>
#include "Vector.hpp"
#include "Matrix.hpp"
// #include "Particule.hpp"
// #include "Boite.hpp"

#include <iostream>	// pour pouvoir utiliser des objets ostreram
using namespace std;

int main( int argc, char * argv[] ) {
	//initialisation du générateur aleatoire
	//srand (time(0));
	
	/*    Proper test of matrix class   */
	Matrix M1(4,2,3); Matrix M1copy(M1); Matrix I4 = identity(4);
	Matrix M2(4,2,5); Matrix M3(2,4,4); 
	cout << "M1(4,2,3): " << M1 << endl;
	cout << "M1copy(M1):" << M1copy << endl;
	cout << "I4:" << I4 << endl;
	cout << "M2(4,2,5): " << M2 << endl;
	cout << "M3(2,4,4): " << M3 << endl;
	cout << "----- Test Matrix overload operations  -----" << endl;
	M1 += M2;
	cout << "M1 +=  M2:" << M1 << endl;
	cout << "M1 *= 3:" << (M1 *= 3) << endl;
	cout << "M1 *= M3:" << (M1 *= M3) << endl; 

	cout << "----- Test Matrix fillSym and transpose  -----" << endl;
	I4.val[1][0] = 2; I4.val[2][3] = -3; 
	cout << "new I4:" << 2*I4 << endl; 
	Matrix I4T = 2*I4.transpose();
	cout << "I4.transpose():" << I4T << endl;
	cout << "I4.fillSym():" << I4T.fillSym() << endl;

	cout << "----- Test Matrix external functions  -----" << endl;
	cout << "Mcopy + M2:" << M1copy + M2 << endl;
	cout << "M2 * M3:" << M2 * M3 << endl;
	cout << "3*Mcopy:" << 3*M1copy << endl;
	cout << "Mcopy*3:" << M1copy*3 << endl;

	// vectors of dim 4 
	Matrix E_1 = colVect(4,1); Matrix E_1T = rowVect(4,1);
	cout << "Vector E_1:" << E_1 << endl;
	cout << "Vector E_1T:" << E_1T << endl;
	cout << "Vector E_1T:" << E_1.transpose() << endl;

	// vector matrix multiplication 
	cout << "E_1T*M1copy:" << E_1T*M1copy << endl;
	cout << "M1copyT*E_1:" << (M1copy.transpose())*E_1 << endl;

	// join vector with matrix 
	cout << "joinLR(E_1,M1copy):"<< joinLR(E_1, M1copy) << endl;

	

	// ------------------------------ END OF TESTS ------------------------------ 
	/*    create Me   */
	Matrix I3 = identity(3); Matrix I6 = identity(6);
	Matrix Mright = joinUD(joinLR(I3,I3),joinLR(I3,I3));
	Matrix Mleft = joinUD(joinLR(2*I3,I3),joinLR(I3,2*I3));
	cout << "Matrix test, Mright :" << Mright << endl;
	cout << "Matrix test, Mleft = 2*I6*Mright:" << Mleft << endl;
	Matrix Me = joinUD(joinLR(Mleft,Mright), joinLR(Mright,Mleft));
	cout << "Elementary matrix Me:" << Me << endl;


	/*   create Ke   */
	double E = 110; double  v = 0.34; double G = 46;
	double c11 = E*(1-v)/((1-2*v)*(1+v)); 
	double c12 = E*v/((1-2*v)*(1+v));
	double F = (c11 - c12)/2;
	cout << "c11 : " << c11 << endl;
	cout << "c12 : " << c12 << endl;
	cout << "G : " << F << endl;

	Matrix C = c11*I6;
	cout << "Matrix test, c11*I6 :" << C << endl;
	for(int i = 3; i < 6; i++){
		C.val[i][i] = G; 
	} 
	C.val[0][1] = c12; C.val[0][2] = c12; C.val[1][1] = c12;
	C.fillSym();
	cout << "Matrix test, C :" << C << endl;




	

	/*   Test of matrix and vectors   */
	// Matrix E2_1 = colVect(2,1); Matrix E5_3 = rowVect(3, 5);

	





	// cout << "Matrix test, W = colVect(2,1):" << W << endl;

	



	// //test de la classe Vector
	// Vector O(3), S(3,1.0);
	// Vector P = Base(3, 0), Q = Base(3, 1), R = Base(3, 2);
	// cout << "O=" << O << endl;
	// cout << "P=" << P << endl;
	// cout << "Q=" << Q << endl;
	// cout << "R=" << R << endl;
	// Vector T = O+P+2.*Q-R/2.;
	// cout << "O+P+2.*Q-R/2.=" << T << endl;

	// //test class Matrix 
	// Matrix M(4,2,3); Matrix N(M); 
	// Matrix W = colVect(2,1);


	// Matrix I4 = identity(4);Matrix I2 = identity(2);
	// cout << "Matrix test, M(4,2,3): " << M << endl;
	// cout << "Matrix test, I(4):" << I4 << endl;
	// cout << "Matrix test, N(M):" << N << endl;
	// cout << "Matrix test, W = colVect(2,1):" << W << endl;
	// cout << "Matrix test, MW :" << M*W << endl;

	// Matrix J1 = joinLR(M,I4);Matrix J2 = joinUD(M,I2);
	// cout << "Matrix test, J1 = [M I]:" << J1 << endl;
	// cout << "Matrix test, J2 = [M; I]:" << J2 << endl;
	// J2.val[1][1] = -5;
	// cout << "Matrix test, J2 = [M; I]:" << J2 << endl;



	// // create Ke 
	// double E = 110; double  v = 0.34; double G = 46;
	// double c11 = E*(1-v)/((1-2*v)*(1+v)); 
	// double c12 = E*v/((1-2*v)*(1+v));
	// double F = (c11 - c12)/2;
	// cout << "c11 : " << c11 << endl;
	// cout << "c12 : " << c12 << endl;
	// cout << "G : " << F << endl;

	// Matrix C = c11*I6;
	// cout << "Matrix test, c11*I6 :" << C << endl;
	// for(int i = 3; i < 6; i++){
	// 	C.val[i][i] = G; 
	// } 
	// C.val[0][1] = c12; C.val[0][2] = c12; C.val[1][1] = c12;
	// C.fillSym();
	// cout << "Matrix test, C :" << C << endl;


	// Matrix MUP = joinLR(M1,M2); 
	// cout << "Elementary matrix MUP:" << MUP << endl;
	// cout << "Elementary matrix MUP^T:" << MUP.transpose() << endl;





	



	// test de la classe Particule
	//Particule Soleil(3), Mercure(P), Venus(Q), Terre(R), Mars(S);
	/*cout << Soleil << Mercure << Venus << Terre << endl;
	for(int i=0; i<10; i++) Terre.saute_mouton(1.0, 1e-4);
	cout << Soleil << Mercure << Venus << Terre << endl;*/
	
	// test de la classe Boite
	// Boite B(3);
	// B.Ajouter(Soleil);
	// B.Ajouter(Mercure);
	// B.Ajouter(Venus);
	// B.Ajouter(Terre);
	// B.Ajouter(Mars);
	// cout << B << endl;
	// B.Update();
	// cout << B << endl;
	
	// B.Appliquer(B);
	// cout << B << endl;
	// B.saute_mouton(1.0);
	// cout << B << endl;
	// Reconstruire
	
	//B.Ajouter(Mercure);
	//cout << B << endl;
	/*//test de la classe Numeros
	Numeros ns(1,2,3);
	cout << "Numeros ns(1,2,3) -> " << ns << endl;
	cout << "ns(2)=" << ns(2) << endl;
	cout<<"ns.shift(2) -> "<<ns.shift(2)<<endl;
	
	//test de la classe Maillage
	Maillage mail(4,2);
	cout<<"Maillage du carre unite"<<endl;
	mail.affiche();
	
	cout<<"Copie du maillage du carre unite"<<endl;
	Maillage mail2(mail);
	mail2.affiche();
	
	cout<<"Maillage du rectangle [2,4]x[10,20]"<<endl;
	Maillage mailr(2.0,4.0,10.0,20.0,4,5);
	mailr.affiche();
	
	cout<<"Transformation affine du maillage du carre unite"<<endl;
	vector<double> A(4,0); A[0]=1; A[3]=1;
	vector<double> t(2,0); t[0]=1; t[1]=1;
	mail2.tf_affine(A,t); //transformation affine de maillage
	mail2.affiche();*/
	
	/*Maillage mail3(2,2);	//maillage carré unité
	t[0]=-1;A[0]=1;
	mail3.tf_affine(A,t); //transformation affine de maillage
	mail3.affiche();
	cout<<"Concatenation de maillage"<<endl;
	Maillage mail4=mail2+mail3; //concaténation de maillages
	mail4.affiche();
	mail4.saveToFile("mail4.dat");
	*/
	return 0;
}
