//#include <ctime>
#include "Matrix.hpp"
#include "Node.hpp"
#include "Elem.hpp"
#include "Grid.hpp"
// #include "Particule.hpp"
// #include "Boite.hpp"

#include <fstream>
#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include <Eigen/Dense>
#include <Eigen/SparseLU> 
#include <Eigen/Sparse> 

// #include <unsupported/Eigen/IterativeSolvers>

#include <Eigen/IterativeLinearSolvers>

#include <vector>
#include <array>
#include <math.h>
#include <string>


using Eigen::MatrixXd;

using namespace Eigen;
using namespace std;

// declare functions 
vector<Elem> generateGrid(int N);
Grid makeGrid(int N);

int main( int argc, char * argv[] ) {

	// string node_fileName =  "nodeList.txt";
	// string soln_fileName =  "solution.txt";	
	// string K_fileName =  "K.txt";
	// string F_fileName =  "F.txt";

	ofstream nodeListTxt, elemListTxt, solnTxt, K_Txt, F_Txt;
	nodeListTxt.open("nodeList.txt");
	elemListTxt.open("elemList.txt");
	solnTxt.open("solution.txt");
	K_Txt.open("K.txt");
	F_Txt.open("F.txt");

	// ------------------------------ GENERATE GRID AND CHECK ------------------------------ 


	// NOTE: Change the number in makeGrid and also in the nb_nodes 
	Grid mainGrid = makeGrid(3);
	int N = 3;
	int nb_nodes = pow((N+1),3);
	int nb_nodes_bound = pow((N+1),2);

	cout << "Grid of size " << N << " with "<< nb_nodes << " nodes generated.\n" << endl;


	cout << "Running through global node list in mainGrid to check: ...\n" << endl;
	int count = 1;
	vector<Node>::iterator itglobal_node = mainGrid.gridNodes.begin();
	//cout << (*itglobal) << endl;
	for (; itglobal_node != mainGrid.gridNodes.end();itglobal_node++)
	{	
		// print node in myNodes
		// cout << "node number: " << count << endl;
		// cout << (*itglobal) << endl;
		nodeListTxt << (*itglobal_node) ;
		//nodeListTxt << "\n";

		count++;
	}

	// Run through global element list to save which nodes are in each elem
	vector<Elem>::iterator itglobal_elem = mainGrid.gridElems.begin();
	//cout << (*itglobal) << endl;
	for (; itglobal_elem != mainGrid.gridElems.end();itglobal_elem++)
	{	
		vector<Node>::iterator nodeInElem = (*itglobal_elem).elemNodes.begin();
		for (; nodeInElem != (*itglobal_elem).elemNodes.end();nodeInElem++){

			elemListTxt << (*nodeInElem);
			// elemListTxt << "\n";
		}
		elemListTxt << "NEXT\n";
	}





	// cout << "Running through elements list in mainGrid to check: ...\n" << endl;
	// vector<Elem>::iterator itelem = mainGrid.gridElems.begin();
	// //cout << (*itglobal) << endl;
	// for (; itelem != mainGrid.gridElems.end();itelem++)
	// {	
	// 	//cout << "this element " << endl;
	// 	// print node in myNodes
	// 	//cout << (*itglobal) << endl;
	// 	//(*itelem).PrintNodes();
	// }

	// ------------------------------ CREATE LOCAL MATRICES ------------------------------ 
    /*    create Me   */



    MatrixXd I3 = MatrixXd::Identity(3,3); 
    MatrixXd I6 = MatrixXd::Identity(6,6); 
    MatrixXd I12 = MatrixXd::Identity(12,12); 
    cout << 2*I3 << endl;

    /*  I3I3
    1 0 0 1 0 0
    0 1 0 0 1 0
    0 0 1 0 0 1
    */
    MatrixXd I3I3(I3.rows(), I3.cols()+I3.cols());
    I3I3 << I3, I3; 
    MatrixXd I3I3_I3I3(I3I3.rows()+I3I3.rows(), I3I3.cols()); // <-- D(A.rows() + B.rows(), ...)
    I3I3_I3I3 << I3I3, I3I3; // <-- syntax is the same for vertical and horizontal concatenation 

    cout << I3I3_I3I3 << endl;

    // create the 12 by 12 matrix Me but without the 2 yet 
    MatrixXd Me_half(I3I3_I3I3.rows(), I3I3_I3I3.cols()+I3I3_I3I3.cols());
    Me_half << I3I3_I3I3, I3I3_I3I3; 
    MatrixXd Me_full(Me_half.rows()+Me_half.rows(), Me_half.cols());
    Me_full << Me_half, Me_half; 
    //cout << Me_full << endl;

    MatrixXd Me = (Me_full + I12)/20;
    cout <<  "Elementary matrix Me:\n" << Me << endl;

	// =================================
    /*     create Ke = Ve B^t C B     */
	// =================================

	//double E = 110; double  v = 0.34; double G = 46;
	double E = 200; double  v = 0.3; double G = 46;
	double c11 = E*(1-v)/((1-2*v)*(1+v)); 
	double c12 = E*v/((1-2*v)*(1+v));
	double c13 = E/(2*(1+v));
	// //double F = (c11 - c12)/2;
	cout << "c11 : " << c11 << endl;
	cout << "c12 : " << c12 << endl;
	// //cout << "G : " << F << endl;
    MatrixXd C = c11*I6; 

	// Matrix C = c11*I6;
	// cout << "Matrix test, c11*I6 :\n" << C << endl;
	for(int i = 3; i < 6; i++){
		C(i,i) = c13; 
	} 
	C(0,1) = c12; C(0,2) = c12; C(1,2) = c12;
    C(1,0) = c12; C(2,0) = c12; C(2,1) = c12;
	// C.fillSym();
	cout << "Matrix test, C :\n" << C << endl;   


	// creation of matrix B = LN 
	// ========================
	MatrixXd B1up = -1*I3; MatrixXd B1down = -1*MatrixXd::Ones(3,3);
	// B1up(1,1) = 1;
	B1down(1,1) = 0;B1down(0,2) = 0;B1down(2,0) = 0;

	cout << "B1up :\n" << B1up << endl;
	cout << "B1down :\n" << B1down << endl;


	MatrixXd Bleft(2*B1up.rows(), B1up.cols());
	Bleft << B1up, B1down ;
	cout << "Bleft :\n" << Bleft << endl;

	//MatrixXd Bright = MatrixXd::Zero(6,9);
	MatrixXd Bright = MatrixXd::Zero(6,9);
	Bright(0,0) = -1; 
	Bright(1,4) = 1;
	Bright(2,8) = -1;
	Bright(3,1) = -1;
	Bright(3,3) = 1;
	Bright(4,2) = -1;
	Bright(4,6) = -1;
	Bright(5,5) = 1;
	Bright(5,7) = -1;

	cout << "Bright :\n" << Bright << endl; 

	// this is B for reference tetrahedron. 
	MatrixXd B(6,12);
    B << Bleft,Bright; 
	cout << "B = LN :\n" << B << endl; 






	// another way to make B by chadrupatla
	// =====================================
	// MatrixXd Aref = MatrixXd::Identity(3,3);

	// MatrixXd B_byA_ref = MatrixXd::Zero(6,12); 
	// for(int i = 0; i< 3; i++){
	// 	cout << "A(0,i): \n" << Aref(0,i) << endl;
	// 	B_byA_ref(0,0 + 3*i) = Aref(0,i); 
	// 	B_byA_ref(1,1 + 3*i) = Aref(1,i);
	// 	B_byA_ref(2,2 + 3*i) = Aref(2,i);
	// 	B_byA_ref(3,1 + 3*i) = Aref(2,i); B_byA_ref(3,2 + 3*i) = Aref(1,i); 
	// 	B_byA_ref(4,0 + 3*i) = Aref(2,i); B_byA_ref(4,2 + 3*i) = Aref(0,i);
	// 	B_byA_ref(5,0 + 3*i) = Aref(1,i); B_byA_ref(5,1 + 3*i) = Aref(0,i);
	// }
	// double A1r = Aref(0,0) + Aref(0,1) + Aref(0,2); 
	// double A2r = Aref(1,0) + Aref(1,1) + Aref(1,2); 
	// double A3r = Aref(2,0) + Aref(2,1) + Aref(2,2); 
	// B_byA_ref(0,9) = -A1r; B_byA_ref(1,10) = -A2r; B_byA_ref(2,11) = -A3r; 
	// B_byA_ref(3,10) = -A3r; B_byA_ref(3,11) = -A2r; 
	// B_byA_ref(4,9) = -A3r; B_byA_ref(4,11) = -A1r; 
	// B_byA_ref(5,9) = -A2r; B_byA_ref(5,10) = -A1r; 

	// cout << "B_byA reference : \n" << B_byA_ref << endl;



	// MatrixXd Ke_byA = (B_byA_ref.transpose())*C*B_byA_ref/6;
	// cout << "Ke_byA :\n" << Ke_byA << endl; 

	// =====================================


	// CEHCKING WTH THE GENERAL FORM 9.50 of AFEM ch 9 colorado
	// ========================================================
	// MatrixXd Ke_check = MatrixXd::Zero(12,3); 
	// double E_tilde = E/(12*(1-2*v)*(1+v)); double v_tilde = 1-2*v; double v_hat = 1- v;
	// MatrixXd Ke_34L = -v_tilde*MatrixXd::Identity(3,4); MatrixXd Ke_35R = -v_tilde*MatrixXd::Identity(3,5); MatrixXd Ke_99 = v_tilde*MatrixXd::Identity(9,9);
	// cout << "v, vhat, vtilde: \n" << v << "  " << v_hat << "  " << v_tilde << "  \n" <<endl;
	// Ke_34L(0,0) = -2*v_hat;
	// Ke_34L(0,1) = -v_tilde;Ke_34L(0,2) = -v_tilde;Ke_34L(0,3) = -v_tilde;
	// Ke_34L(1,3) = -v_tilde;
	// Ke_34L(1,0) = -2*v;Ke_34L(2,0) = -2*v;

	// Ke_35R(0,0) = -2*v;
	// Ke_35R(0,2) = -v_tilde;
	// Ke_35R(2,3) = -v_tilde;Ke_35R(2,1) = -v_tilde;
	// Ke_35R(1,3) = -v_tilde;
	// Ke_35R(1,0) = -2*v_hat;Ke_35R(2,0) = -2*v;
	// Ke_35R(0,4) = -2*v;Ke_35R(1,4) = -2*v;Ke_35R(2,4) = -2*v_hat;


	// cout << "Ke_34L:\n" << Ke_34L << endl;
	// cout << "Ke_35R:\n" << Ke_35R << endl;

	// Ke_99(0,0) = 2*v_hat; Ke_99(0,4) = 2*v;  Ke_99(0,8) = 2*v; 
	// Ke_99(4,4) = 2*v_hat; Ke_99(4,0) = 2*v;  Ke_99(4,8) = 2*v; 
	// Ke_99(8,8) = 2*v_hat; Ke_99(8,0) = 2*v;  Ke_99(8,4) = 2*v; 
	// Ke_99(1,3) = v_tilde; Ke_99(3,1) = v_tilde; Ke_99(2,6) = v_tilde; Ke_99(6,2) = v_tilde; Ke_99(7,5) = v_tilde;  Ke_99(5,7) = v_tilde;  
 	// cout << "Ke_99:\n" << Ke_99 << endl;

	

	// // first 3 rows
	// Ke_check(0,0) = 4 - 6*v;Ke_check(0,1) = 1;Ke_check(0,2) = 1;
	// Ke_check(1,1) = 4 - 6*v;Ke_check(1,0) = 1;Ke_check(1,2) = 1;
	// Ke_check(2,2) = 4 - 6*v;Ke_check(2,0) = 1;Ke_check(2,1) = 1;

	// Ke_check(3,0) = -2*v_hat;Ke_check(4,0) = -v_tilde;Ke_check(5,0) = -v_tilde;Ke_check(6,0) = -v_tilde;
	// Ke_check(3,1) = -2*v;Ke_check(4,1) = -v_tilde;Ke_check(6,1) = -v_tilde;
	// Ke_check(3,2) = -2*v;Ke_check(5,2) = -v_tilde;

	// Ke_check(7,0) = -2*v;Ke_check(9,0) = -v_tilde;
	// Ke_check(7,1) = -2*v_hat;Ke_check(8,1) = -v_tilde;Ke_check(10,1) = -v_tilde;
	// Ke_check(7,2) = -2*v;Ke_check(8,2) = -v_tilde;Ke_check(9,2) = -v_tilde;Ke_check(10,2) = -v_tilde;

	// Ke_check(11,0) = -2*v;
	// Ke_check(11,1) = -2*v;
	// Ke_check(11,2) = -2*v_hat;
	// cout << "KeCHECK:\n" << Ke_check << endl;


	// MatrixXd Ke_check_tot(12,12); MatrixXd Ke_right_top(3,9); MatrixXd Ke_right_tot(12,9);
	// Ke_right_top << Ke_34L, Ke_35R; 
	// Ke_right_tot << Ke_right_top, Ke_99;
	// Ke_check_tot << Ke_check, Ke_right_tot;

	// cout << "Ke_check_TOTAL:\n" << E_tilde*Ke_check_tot <<endl;
	// cout << "=================\n\n" <<endl;


	// ========================================================



	MatrixXd Ke = (B.transpose())*C*B/6;
	cout << "Ke :\n" << Ke << endl; 

	Vector3d Fe1 = Vector3d::Zero();
	Fe1(2) = -9.81/10;
	//cout << "Fe1 :\n" << Fe1 << endl;
	VectorXd Fe = VectorXd::Zero(12);
	Fe << Fe1, Fe1, Fe1, Fe1;

	cout << "Fe :\n" << Fe << endl;

 



    int nb_fixedNodes = 4;
    // // ------------------------------ CREATE GLOBAL KE ------------------------------ 
    int elem_index = 0;
	// pseudo code 
	// for each element, get the coordinates of Xi and then calculate the Jacobian then get the element matrice Ke-global 
    // for each element of type Elem, it has a vector list of Nodes that are arranged (indexed). Can access from 0 to 4
    vector<Elem>::iterator itn = mainGrid.gridElems.begin();
	for (; itn != mainGrid.gridElems.end(); itn++){
		// cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
		// cout << "--------------------------------------- " << endl;
        
		//(*itn).PrintNodes();
        // this is of type Node

        // cout << (*itn).elemNodes[0] << endl;
        // cout << (*itn).elemNodes[1] << endl;
        // cout << (*itn).elemNodes[2] << endl;

		MatrixXd Ke_global; MatrixXd Ke_globalBB; VectorXd Fe_global; VectorXd D_e = VectorXd::Zero(12); MatrixXd D = MatrixXd::Zero(3,4);
		MatrixXd Ke_global_inv; MatrixXd Fe_global_inv;
		MatrixXd nodeSelec(3,4);
		nodeSelec << (*itn).elemNodes[0].position_, (*itn).elemNodes[1].position_,  (*itn).elemNodes[2].position_,  (*itn).elemNodes[3].position_;
		//cout << "nodeSelec: \n" << nodeSelec << endl;





        // Initialise Jacobian matrix 
		//Matrix Je(3,3);
        // assign to X0 the coordinates of the first node of the element
        Vector3d X0((*itn).elemNodes[0].position_.transpose());
        //cout << "X0:\n" << X0 << endl;
		Vector3d X1((*itn).elemNodes[1].position_.transpose());
		//cout << "X1:\n" << X1 << endl;
		Vector3d X2((*itn).elemNodes[2].position_.transpose());
		//cout << "X2:\n" << X2 << endl;
		Vector3d X3((*itn).elemNodes[3].position_.transpose());
		//cout << "X3:\n" << X3 << "\n" << endl;

		D_e << X0, X1, X2, X3;
		D << X0, X1, X2, X3;
		// cout << "De:\n" << D_e << "\n" << endl;		
		//cout << "D:\n" << D << "\n" << endl;		

		Vector3d diffX1X0 = X1-X0;
		Vector3d diffX2X0 = X2-X0;
		Vector3d diffX3X0 = X3-X0;	


		// cout << "X1-X0:\n" << diffX1X0 << endl;
		// cout << "X2-X0:\n" << diffX2X0 << endl;
		// cout << "X3-X0:\n" << diffX3X0 << endl;





		// // ========== Try to make Ke_global bty explicitly calculating B. 

		MatrixXd B_global = MatrixXd::Zero(6,12); 
		for(int i = 0; i < 4; i++){
			// calculate b_i, c_i, d_i. Do this by creating the matrix to take det first. 
			MatrixXd B_TMP = MatrixXd::Ones(3,3);
			MatrixXd C_TMP = MatrixXd::Ones(3,3);
			MatrixXd D_TMP = MatrixXd::Ones(3,3);

			// cout << "checking (i+1) %4: \n" << (i+1)%4 << endl;
			// first arg of node selec is to say either x, y or z. second arg is to choose the subscript. 
			// b_0 = -det (1 y1 z1)
			//             1 y2 z2)
			//             1 y3 z3)
			B_TMP(0,1) = nodeSelec(1,(i+1)%4); B_TMP(0,2) = nodeSelec(2,(i+1)%4); 
			B_TMP(1,1) = nodeSelec(1,(i+2)%4); B_TMP(1,2) = nodeSelec(2,(i+2)%4); 
			B_TMP(2,1) = nodeSelec(1,(i+3)%4); B_TMP(2,2) = nodeSelec(2,(i+3)%4); 
			// cout << "B_TMP: \n" << B_TMP << endl;
			
			C_TMP(0,0) = nodeSelec(0,(i+1)%4); C_TMP(0,2) = nodeSelec(2,(i+1)%4); 
			C_TMP(1,0) = nodeSelec(0,(i+2)%4); C_TMP(1,2) = nodeSelec(2,(i+2)%4); 
			C_TMP(2,0) = nodeSelec(0,(i+3)%4); C_TMP(2,2) = nodeSelec(2,(i+3)%4); 
			// cout << "C_TMP: \n" << C_TMP << endl;
						
			D_TMP(0,0) = nodeSelec(0,(i+1)%4); D_TMP(0,1) = nodeSelec(1,(i+1)%4); 
			D_TMP(1,0) = nodeSelec(0,(i+2)%4); D_TMP(1,1) = nodeSelec(1,(i+2)%4); 
			D_TMP(2,0) = nodeSelec(0,(i+3)%4); D_TMP(2,1) = nodeSelec(1,(i+3)%4); 
			// cout << "D_TMP: \n" << D_TMP << endl;


			double b_i =  -1*B_TMP.determinant();
			double c_i =  -1*C_TMP.determinant();
			double d_i =  -1*D_TMP.determinant();
			// cout << "b_i, c_i, d_i:\n" << b_i << " " << c_i << " " << d_i << endl;

			// fil in the matrix B for this element. 
			B_global(0,0 + 3*i) = b_i; 
			B_global(1,1 + 3*i) = c_i;
			B_global(2,2 + 3*i) = d_i;
			B_global(3,0 + 3*i) = c_i; B(3,1 + 3*i) = b_i; 
			B_global(4,1 + 3*i) = d_i; B(4,2 + 3*i) = c_i;
			B_global(5,0 + 3*i) = d_i; B(5,2 + 3*i) = b_i;

		}
		//cout << "B_global: \n" << B_global << endl;










		// ================================================================
		// Creating B_byA (chandrupatla)
		// ================================================================


		MatrixXd Je(3,3); 
		Je << diffX1X0, diffX2X0, diffX3X0;

		// MatrixXd Zero3 = MatrixXd::Zero(3,3); Vector3d ZeroCol3 = Vector3d::Zero(3);
		// MatrixXd A(3,3);
		// MatrixXd JeT = Je.transpose();
		// A = JeT.inverse();
		// MatrixXd Tv = Je.inverse();
		// // cout << "Je inverse = Tv: \n" << Tv << endl;
		// //cout << "Je : \n" << Je << endl;
		// //cout << "A : \n" << Je << endl;
		// // cout << "JeT : \n" << JeT << endl;




		// MatrixXd B_byA = MatrixXd::Zero(6,12); 
		// for(int i = 0; i< 3; i++){
		// 	// cout << "A(0,i): \n" << A(0,i) << endl;
		// 	B_byA(0,0 + 3*i) = A(0,i); 
		// 	B_byA(1,1 + 3*i) = A(1,i);
		// 	B_byA(2,2 + 3*i) = A(2,i);
		// 	B_byA(3,1 + 3*i) = A(2,i); B_byA(3,2 + 3*i) = A(1,i); 
		// 	B_byA(4,0 + 3*i) = A(2,i); B_byA(4,2 + 3*i) = A(0,i);
		// 	B_byA(5,0 + 3*i) = A(1,i); B_byA(5,1 + 3*i) = A(0,i);
		// }
		// double A1 = A(0,0) + A(0,1) + A(0,2); 
		// double A2 = A(1,0) + A(1,1) + A(1,2); 
		// double A3 = A(2,0) + A(2,1) + A(2,2); 
		// B_byA(0,9) = -A1; B_byA(1,10) = -A2; B_byA(2,11) = -A3; 
		// B_byA(3,10) = -A3; B_byA(3,11) = -A2; 
		// B_byA(4,9) = -A3; B_byA(4,11) = -A1; 
		// B_byA(5,9) = -A2; B_byA(5,10) = -A1; 

		//cout << "B_byA : \n" << B_byA << endl;
		
		// MatrixXd T0(3,12); MatrixXd T1(3,12); MatrixXd T2(3,12); MatrixXd T3(3,12);
		// MatrixXd W0(3,12); MatrixXd W1(3,12); MatrixXd W2(3,12); MatrixXd W3(3,12);

		// T0 << Zero3, Zero3, Zero3, Zero3; 
    	// T1 << -Tv, Tv, Zero3, Zero3;  
		// T2 << -Tv, Zero3, Tv, Zero3; 
		// T3 << -Tv, Zero3, Zero3, Tv; 

		// W0 << Zero3, Zero3, Zero3, Zero3; 
    	// W1 << -Je, Je, Zero3, Zero3;  
		// W2 << -Je, Zero3, Je, Zero3; 
		// W3 << -Je, Zero3, Zero3, Je; 
		
		// MatrixXd T(12,12); MatrixXd W(12,12);
		// T << T0, T1, T2, T3;
		// W << W0, W1, W2, W3;

		// cout << "T :\n" << T << endl; 
		// cout << "W :\n" << W << endl; 
		// cout << "T.inv :\n" << T.inverse() << endl; 
		// cout << "T.transpose :\n" << T.transpose() << endl; 







		// Using the brute force method of calculating B
		//Ke_globalBB = (Je.determinant())* (B_global.transpose())*C*B_global/6;
		Ke_globalBB = (B_global.transpose())*C*B_global/6;
		// Ke_global = (B_byA.transpose())*C*B_byA/6;

		// Using the transformation method
		// Ke_global = T.transpose()*Ke*T;
		// Fe_global = T.transpose()*Fe;
		// Ke_global_inv = T.inverse()*Ke*T;
		// Fe_global_inv = T.inverse()*Fe;
		// Fe_global = Fe + W*Fe;
		Fe_global = Fe;
		//cout << "W*Fe: \n" << W*Fe << endl;

		// assign Ke and Fe for the element
		(*itn).Ke_global = Ke_globalBB;
		(*itn).Fe_global = Fe_global;

		cout << " ------ Summary for element ------ "<< elem_index + 1 << endl;
		// cout << "de = T De:\n" << T*D_e << "\n" << endl;		
		// cout << "Jacobian: \n" << Je << endl;
		// cout << "Det (Je) : \n" << Je.determinant() << endl;
		// cout << "Ke global by A: \n" << Ke_global << endl;

		cout << "Ke global BB : \n" << Ke_globalBB << endl;
		// cout << "Fe : \n" << Fe_global << endl;
		// cout << "Ke global INV: \n" << Ke_global_inv << endl;
		// cout << "Fe INV: \n" << Fe_global_inv << endl;


        elem_index++;

	}
	// //SparseMatrix<double> K(3*nb_nodes, 3*nb_nodes);
	MatrixXd K = MatrixXd::Zero(3*nb_nodes, 3*nb_nodes);
	VectorXd F = VectorXd::Zero(3*nb_nodes);

	// Assemble F 
	// for force in z - direction, 3*globalNodeIndex_j+2. !!!!!!!!!!!!!!!!!!!

	//cout << "Initializing K and F...\n "<< K << "\n" << F << endl;


	// Matrix K(3*nb_nodes, 3*nb_nodes, 0); Matrix F(3*nb_nodes, 1, 0);
	// cout << "Initialising global stiffness matrix K... " << K << endl;
	// // run through each node vertex and check in each element the node exists there


	cout << "RUNNING THROUGH NODES FOR FILLING MATRIX K\n" << endl;
	int globalNodeIndex_i = 0;
	vector<Node>::iterator itnode = mainGrid.gridNodes.begin();
	for (; itnode != mainGrid.gridNodes.end();itnode++)
	{	
		// print node in myNodes
		//cout << "node number: " << count << endl;
		// this is
		globalNodeIndex_i = itnode - mainGrid.gridNodes.begin();
		// cout << "Looking at a current node with global index " << globalNodeIndex_i << endl;  
		// cout << (*itnode) << endl;
		// cout << "========================================== " << endl;
		// cout << "========================================== " << endl;

		
		//given this current node, we run through the elements. 
		vector<Elem>::iterator itelem = mainGrid.gridElems.begin();
		for (; itelem != mainGrid.gridElems.end();itelem++)
		{	
			// first we need to check if this node (*itnode) is in this element. 
			//---------------------------------------
			
			// if itnode.pos == 
			// run through nodes of this element and check if itnode is inside
			bool inElem = false;
			int pos_inElem_i = 0;
			vector<Node>::iterator itelemnode = (*itelem).elemNodes.begin();
			for (; itelemnode != (*itelem).elemNodes.end();itelemnode++)
			{	
				// print node in myNodes
				//cout << "node number: " << count << endl;
				// this is
				// GLOBAL NODE i IS IN THE ELEMENT.   
				// ===============================
				// ===============================
				if((*itelemnode).position_ == (*itnode).position_){
					inElem = true; 
					pos_inElem_i = itelemnode - (*itelem).elemNodes.begin();
					// cout << "========================================== " << endl;
					// cout << (*itnode) << "is found in the element at position " << pos_inElem_i << endl;
					// (*itelem).PrintNodes();
					// cout << " ==== Run through other nodes in this elem... ==== " << endl;
					
					// now run through all the other nodes of this element and get the local and global index (j)
					int pos_inElem_j = 0; int globalNodeIndex_j = 0;
					vector<Node>::iterator itelemnode2 = (*itelem).elemNodes.begin();
					for(; itelemnode2 != (*itelem).elemNodes.end();itelemnode2++)
					{
						// I HAVE A FEEELING THE LINE AFTER THIS IS ACTL NOT TAKING INTO ACCOUNT THE SELF ADDIDIOTN
						//=================================================== !!!!!!!!!!!!!
						//if (!((*itelemnode2).position_ == (*itnode).position_)){
						//=================================================== !!!!!!!!!!!!!

						// cout << "j in this element: " << (*itelemnode2).position_ << endl;
						
						// get local and global positions of this node. 
						pos_inElem_j = itelemnode2 - (*itelem).elemNodes.begin();
						// get global node position of itelemnode2
						vector<Node>::iterator itnode2 = mainGrid.gridNodes.begin();
						for (; itnode2 != mainGrid.gridNodes.end();itnode2++)
						{
							if((*itnode2).position_ == (*itelemnode2).position_){
								globalNodeIndex_j = itnode2 - mainGrid.gridNodes.begin();
							}
						}
						// cout << "global i:" << globalNodeIndex_i <<endl;
						// cout << "global j:" << globalNodeIndex_j <<endl;

						// cout << "local j :" << pos_inElem_j << endl;
						// cout << " ------ end of resume for j ------ " << endl;


						// little loop because we need to fill in 3 by 3 matrix 
						// // for z = 0, 3*globalNodeIndex_j+2. !!!!!!!!!!!!!!!!!!!
						// F(3*globalNodeIndex_j+2) -= 2210*9.81/24;
						//F(3*globalNodeIndex_j+2) -= 9.81;
						// F(3*globalNodeIndex_j+2) += (*itelem).Fe_global(2);
						

						for (int i_local = 0; i_local < 3; i_local++){
							for (int j_local = 0; j_local < 3; j_local++){
								// i and j for the Ke global should be the positions found for the nodes
								// for example you can have i = 5, but the global element matrix has only Kij for i,j = 1,2,3,4
								// and each Kij is a 3 by 3 matrix
								
								//K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) = K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) +  Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
								K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) += (*itelem).Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
								//+= Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
								F(3*globalNodeIndex_i+j_local) += (*itelem).Fe_global(3*pos_inElem_j + j_local);
								//cout << "updating K : ..............." << K << endl;
								

							}
						}
						// cout << "updating K : ...............\n" << K << endl;


						// }
					}

				}
			}			
		}
	}

	
	cout << "============= Summary =============" << endl; 
	// cout << "Number of nodes: " << nb_nodes << endl; 
	// cout << "Global stiffness matrix: \n" << K << endl;
	// cout << "Global force vector: \n" << F << endl;




	// run through each node, check if z= 0. If yes, modify that node index row. For example, if node i = 6 is a node of z= 0, 
	// then we need u_i = 0, so we need Kii = 1, Fi = 0, Kij = 0 if i != j
	cout << "Taking into account dirichlet conditions ..." << endl;
	vector<Node>::iterator itDir = mainGrid.gridNodes.begin();
	int i_Dir;
	for (; itDir != mainGrid.gridNodes.end();itDir++)
	{	
		// cout << "--------------" << endl;
		// cout << (*itDir) << endl;

		// check if the z component of this node's position is 0.
		// cout << (*itDir).position_(2) << endl;
		if ((*itDir).position_(2) == 0) {
			i_Dir = itDir - mainGrid.gridNodes.begin();
			// cout << "i_Dir = " << i_Dir << endl;
			// cout << "This node is at the boundary, need to modify the Kij andFi" << endl;

			for (int i_local = 0; i_local <3; i_local++ ){
				// K.row(3*i_Dir+2) *= 0; K.col(3*i_Dir+2) *= 0;
				// K(3*i_Dir+2,3*i_Dir+2) =  1;
				// F(3*i_Dir+2) = 0;				
				K.row(3*i_Dir+i_local) *= 0; K.col(3*i_Dir+i_local) *= 0;
				K(3*i_Dir+i_local,3*i_Dir+i_local) =  1;
				F(3*i_Dir+i_local) = 0;
			}
		}

	}
	cout << "============= Summary After Dirichlet conditions =============" << endl; 
	// cout << "Number of nodes: " << nb_nodes << endl; 
	// cout << "Global stiffness matrix: \n" << K << endl;
	// cout << "Global force vector: \n" << F << endl;

	// MatrixXd K_modified = MatrixXd::Zero(3*nb_nodes-3*nb_nodes_bound, 3*nb_nodes-3*nb_nodes_bound);
	// VectorXd F_modified = VectorXd::Zero(3*nb_nodes-3*nb_nodes_bound);

	for(int i = 0; i < 3*nb_nodes; i++){
		// if the row does not contain 0 0 0 0 1 0 0 0 ... , copy to K_modified
		F_Txt << F(i);
		for(int j = 0; j < 3*nb_nodes; j++){
			
			K_Txt<< K(i,j) << " " ;

			// if((i > 3*nb_nodes_bound-1) && (j > 3*nb_nodes_bound-1)){
			// 	K_modified(i-3*nb_nodes_bound, j-3*nb_nodes_bound) = K(i,j);
			// 	F_modified(i-3*nb_nodes_bound) = F(i);


			// }
		}
		K_Txt <<"\n" ;
		F_Txt <<"\n" ;
	}


	// Saving Kmod and Fmod instead 
	// for(int i = 0; i < 3*nb_nodes-3*nb_nodes_bound; i++){
	// 	F_Txt << F_modified(i);
	// 	for(int j = 0; j < 3*nb_nodes-3*nb_nodes_bound; j++){
	// 		K_Txt<< K_modified(i,j) << " " ;
	// 	}
	// 	K_Txt <<"\n" ;
	// 	F_Txt <<"\n" ;
	// }

	// cout << "K Modified: \n" << K_modified << endl;
	// cout << "F Modified: \n" << F_modified << endl;




	// cout << "Copying matrix K to testM " << endl;
	// // just testing out the matrix class from Eigen 
	// // Eigen::MatrixXf is the default general matrix one but for sparse matrices it is better to use Eigen::SparseMatrix

	// int Ksize = (int) nb_nodes*3;
	// Eigen::SparseMatrix<double> testK(Ksize, Ksize);
	// for(int i = 0; i< nb_nodes*3; i++){
	// 	for(int j = 0; j< nb_nodes*3; j++){
	// 		testK.coeffRef(i,j) = K.val[i][j];
	// 		//cout << testK.coeffRef(i,j) << " " ;
	// 	}
	// 	//cout << "\n" << endl;
	// }

	// cout << "Copying matrix K to testM_matrix " << endl;

	// Eigen::MatrixXf testK_matrix(Ksize, Ksize);
	// for(int i = 0; i< nb_nodes*3; i++){
	// 	for(int j = 0; j< nb_nodes*3; j++){
	// 		testK_matrix(i,j) = K.val[i][j];
	// 		//cout << testK_matrix(i,j) << " " ;
	// 	}
	// 	//cout << "\n" << endl;
	// }

	// cout << "Copying matrix F to testF " << endl;

	// Eigen::VectorXf testF(Ksize);
	// for(int i = 0; i< nb_nodes*3; i++){
	// 	testF(i) = F.val[i][0];
	// 	//cout << testF(i) << " " ;
	// }
	


	//============================================================
	//========= SOLVING USING HOUSEHOLDER QR APPROX (WORKS) ======
	//Eigen::VectorXd x_sol_inv = K_modified.inverse()*F_modified;
	//cout << "Kinv: \n" << K_modified.inverse();


	//Eigen::VectorXd x_sol_inv = K_modified.inverse()*F_modified;
	// cout << "Kmodinv: \n" << K_modified.inverse();
	// cout << "===============END=============" << endl;;

	// cout << "Kmodinv*Fmod: \n" << x_sol_inv;
	// cout << "===============END=============" << endl;;


	// Eigen::VectorXd x_sol_full = K.colPivHouseholderQr().solve(F);
    // //cout << "The full solution is:\n" << x_sol_full << endl;
	// for(int i = 0; i < 3*nb_nodes; i++){
	// 	solnTxt<< x_sol_full(i) ;
	// 	solnTxt << "\n";
	// }


    // // SOLVING BY PARTITION
	// Eigen::VectorXd x_sol = K_modified.colPivHouseholderQr().solve(F_modified);
	// Eigen::VectorXd x_sol = (K_modified.inverse())*F_modified;
    // cout << "The solution is:\n" << x_sol << endl;

	// for(int i = 0; i < 3*nb_nodes; i++){
	// 	if(i < 3*nb_nodes_bound){
	// 		solnTxt<< 0;
	// 	}
	// 	else{
	// 		solnTxt<< x_sol(i-3*nb_nodes_bound);
	// 	}
	// 	solnTxt << "\n";
	// }



	// cout << "x_sol = :\n" << x_sol << endl;
	// cout << "x_sol_inv = :\n" << x_sol_inv << endl;
	// cout << "Check solution, KU = :\n" << K*x_sol << endl;











	nodeListTxt.close();
	elemListTxt.close();
	solnTxt.close();
	K_Txt.close();
	F_Txt.close();
	return 0;
}






// double get_b(Vector3d X0, Vector3d& X1, Vector3d& X2, Vector3d& X3, int i){
// 	MatrixXd TMP = MatrixXd::Identity(3,3);
// }


Grid makeGrid(int N){
	int Size = N;
	Grid grid(Size);
    //vector<Elem> elements_grid_list;
	//vector<Node> nodes_grid_list;
	for(int level = 0; level < N; level++){
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				// Going through a cubic element
				// cout << "i, j = " << i << ", " << j << endl;
				int elem_index = i*N + (j+1) + N*level;
				// cout << "*******************************" << endl;

				// cout << "*** Elem_cube_index = "<< elem_index << " *** " <<  endl;
				// cout << "*******************************" << endl;

				// this displays the left side of cubic element on grid. We should create 6 tetrahedral elements each 
				// containing respective nodes. 

				vector<Elem> elements_cube_list;
				vector<Node> nodes_cube_list;

				// create nodes and simultaneously push them into the elements 
				int modo = 1;
				int i_loc = i-1; int j_loc = j; int z_loc = 0;

				// for each node in cubic element
				
				for(int k = 0; k<8; k++){

					// test 
					// cout << "k%4: " << k%4 << endl;
					// cout << "k/4: " << k/4 << endl;
					z_loc = k/4; 
					// cout << "=======" << endl;

					modo += k%2 ;
					j_loc += pow(-1, modo)*(k%2);
					i_loc += pow(-1, modo)*((k-1)%2);
					// cout << "Node "<< k +1 << ": " << j_loc << ", " <<  i_loc << ", "  << z_loc << endl;

					//cout << "Node "<< k +1<<  "[modo, x, y, z]: "<<modo << ", " << j_loc << ", " <<  i_loc << ", "  << z_loc << endl;
					Vector3d pos_loc( (1.0)*j_loc,  (1.0)*i_loc, (1.0)*z_loc + level);
					cout << "CHECKING VECTOR3d pos_loc: \n" << pos_loc << endl; 
					//Vect(3,0); pos_loc.val[0][0] = (1.0)*j_loc; pos_loc.val[0][1] = (1.0)*i_loc; pos_loc.val[0][2] = (1.0)*z_loc;
					
					// !!!!!!!!!!!!!!!!!!!!!!!!!!
					// tried to make a workaround by creating two nodes, one to push into nodes cube list and another to put into gridNodes. 
					Node node_loc(pos_loc, elem_index); Node node_glob(pos_loc, elem_index);
					// cout << "Pushing local node into nodes_cube_list\n" << endl;
					nodes_cube_list.push_back(node_loc);
					// pushing all nodes generated for a cube into this list. In fact this overlaps.
					// cout << "Pushing global node into gridNodes_list\n" << endl;
					// we want to push it into gridNodes only if we cannot find the element in the existing list. 
					bool found_global_node = false; 
					vector<Node>::iterator itcheck = grid.gridNodes.begin();
					//cout << (*itglobal) << endl;
					for (; itcheck != grid.gridNodes.end();itcheck++)
					{	
						// print node in myNodes
						//cout << "node number: " << count << endl;
						// cout << (*itcheck) << endl;
						if ((*itcheck).position_ == pos_loc)
						{
							// cout << "node already registered in global node list" << endl;
							found_global_node = true;
						}
					}
					if (!found_global_node) {
						grid.gridNodes.push_back(node_glob);
					}
					// if(std::find(grid.gridNodes.begin(), grid.gridNodes.end(), node_glob) != grid.gridNodes.end()){
					// 	cout<< "Node globa found in existing list";
					// 	//grid.gridNodes.push_back(node_glob);
					// } 
					// else{
					// 	grid.gridNodes.push_back(node_glob);
					// }
					

				}
				// cout << "=======" << endl;
				// create 6 local elements and stuff them into list. For each element we do need to say which nodes they have
				for(int l = 0; l<6 ; l++){
					Elem elem_loc((elem_index-1)*6 + l);
					// every tetrahedral element has node 1 in it
					elem_loc.elemNodes.push_back(nodes_cube_list[0]);
					if (l == 0){
						elem_loc.elemNodes.push_back(nodes_cube_list[1]);
						elem_loc.elemNodes.push_back(nodes_cube_list[2]);
					}
					if (l == 1){
						elem_loc.elemNodes.push_back(nodes_cube_list[5]);
						elem_loc.elemNodes.push_back(nodes_cube_list[1]);
					}
					if (l == 2){
						elem_loc.elemNodes.push_back(nodes_cube_list[4]);
						elem_loc.elemNodes.push_back(nodes_cube_list[5]);
					}
					if (l == 3){
						elem_loc.elemNodes.push_back(nodes_cube_list[7]);
						elem_loc.elemNodes.push_back(nodes_cube_list[4]);
					}
					if (l == 4){
						elem_loc.elemNodes.push_back(nodes_cube_list[3]);
						elem_loc.elemNodes.push_back(nodes_cube_list[7]);
					}
					if (l == 5){
						elem_loc.elemNodes.push_back(nodes_cube_list[2]);
						elem_loc.elemNodes.push_back(nodes_cube_list[3]);
					}
					elem_loc.elemNodes.push_back(nodes_cube_list[6]);

					// put the local element into the elements cube list and element grid list. 
					// element cube list is only for when in this loop going through each cube
					elements_cube_list.push_back(elem_loc); 
					grid.gridElems.push_back(elem_loc); 

				}

				// check process for a single cubic element (created 8 nodes and 6 elements)
				// vector<Elem>::iterator itn = elements_cube_list.begin();
				// for (; itn != elements_cube_list.end(); itn++)
				// {
					








				// 	cout << "Run through element in a cube"  << "..." <<endl;
				// 	//cout << "\t" << endl;
				// 	(*itn).PrintNodes();
				// 	//cout << *itn << endl;
				// }


			}
		}
	}

	
    return grid;
}