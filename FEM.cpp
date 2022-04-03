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
#include <chrono>

#include <stdio.h>
#include <set>

#include <C:\Users\lin-d\AppData\Local\Programs\Python\Python310\include\Python.h>



using Eigen::MatrixXd;

using namespace Eigen;
using namespace std;

// declare functions 
vector<Elem> generateGrid(int N);
void removeCubes(int* cubesToRemove, int size, vector<Elem>& vectElem);
void createVerticalHole(int* hole_index, int nb_holes, int N, vector<Elem>& vectElem);
void addCube(set<int>& cubesToAdd, int size, vector<Elem>& vectRef, vector<Elem>& vectElem);
void simulateFEM(vector<Node>& thisGridNodes,vector<Elem>& thisGridElems, int nb_nodes, ofstream& K_txt, ofstream& F_txt);
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
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	/* code you want to measure */


	// NOTE: Change the number in makeGrid and also in the nb_nodes 
	Grid mainGrid = makeGrid(8);
	int N = 8;
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
	// vector<Elem>::iterator itglobal_elem = mainGrid.gridElems.begin();
	// //cout << (*itglobal) << endl;
	// for (; itglobal_elem != mainGrid.gridElems.end();itglobal_elem++)
	// {	
	// 	vector<Node>::iterator nodeInElem = (*itglobal_elem).elemNodes.begin();
	// 	for (; nodeInElem != (*itglobal_elem).elemNodes.end();nodeInElem++){

	// 		elemListTxt << (*nodeInElem);
	// 		// elemListTxt << "\n";
	// 	}
	// 	elemListTxt << "NEXT\n";
	// }




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
	// cout << "c11 : " << c11 << endl;
	// cout << "c12 : " << c12 << endl;
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
	// cout << "Matrix test, C :\n" << C << endl;   


	// creation of matrix B = LN 
	// ========================
	MatrixXd B1up = -1*I3; MatrixXd B1down = -1*MatrixXd::Ones(3,3);
	// B1up(1,1) = 1;
	B1down(1,1) = 0;B1down(0,2) = 0;B1down(2,0) = 0;

	// cout << "B1up :\n" << B1up << endl;
	// cout << "B1down :\n" << B1down << endl;


	MatrixXd Bleft(2*B1up.rows(), B1up.cols());
	Bleft << B1up, B1down ;
	// cout << "Bleft :\n" << Bleft << endl;

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

	// cout << "Bright :\n" << Bright << endl; 

	// this is B for reference tetrahedron. 
	MatrixXd B(6,12);
    B << Bleft,Bright; 
	// cout << "B = LN :\n" << B << endl; 


	MatrixXd Ke = (B.transpose())*C*B/6;
	// cout << "Ke :\n" << Ke << endl; 

	Vector3d Fe1 = Vector3d::Zero();
	Fe1(2) = -9.81/5;
	//cout << "Fe1 :\n" << Fe1 << endl;
	VectorXd Fe = VectorXd::Zero(12);
	Fe << Fe1, Fe1, Fe1, Fe1;

	// cout << "Fe :\n" << Fe << endl;

    int nb_fixedNodes = 4;

	// can try running through every thing and get rid of unphysical solution 


    // // ------------------------------ CREATE GLOBAL KE FOR EACH ELEMENT ------------------------------ 
    int elem_index = 0;
	// pseudo code 
	// for each element, get the coordinates of Xi and then calculate the Jacobian then get the element matrice Ke-global 
    // for each element of type Elem, it has a vector list of Nodes that are arranged (indexed). Can access from 0 to 4
    vector<Elem>::iterator itn = mainGrid.gridElems.begin();
	for (; itn != mainGrid.gridElems.end(); itn++){
		cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
		// cout << "--------------------------------------- " << endl;
        
		//(*itn).PrintNodes();
        // this is of type Node

        // cout << (*itn).elemNodes[0] << endl;
        // cout << (*itn).elemNodes[1] << endl;
        // cout << (*itn).elemNodes[2] << endl;

		MatrixXd Ke_global; MatrixXd Ke_globalBB; VectorXd Fe_global; VectorXd D_e = VectorXd::Zero(12); MatrixXd D = MatrixXd::Zero(3,4);
		MatrixXd nodeSelec(3,4);
		nodeSelec << (*itn).elemNodes[0].position_, (*itn).elemNodes[1].position_,  (*itn).elemNodes[2].position_,  (*itn).elemNodes[3].position_;
		//cout << "nodeSelec: \n" << nodeSelec << endl;




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


		Ke_globalBB = (B_global.transpose())*C*B_global/6;

		// Fe_global = Fe;
		//cout << "W*Fe: \n" << W*Fe << endl;

		// assign Ke and Fe for the element
		(*itn).Ke_global = Ke_globalBB;

		// cout << " ------ Summary for element ------ "<< elem_index + 1 << endl;
		// cout << "Ke global BB : \n" << Ke_globalBB << endl;
		// cout << "Fe : \n" << Fe_global << endl;
        elem_index++;

	}


	//cout << "Initializing K and F...\n "<< K << "\n" << F << endl;

	// --------------------------------------
	// --------------------------------------
	// checking original elements
	vector<Elem>::iterator init = mainGrid.gridElems.begin();
	for (; init != mainGrid.gridElems.end(); init++){
		cout << "Original element index: \n" << (*init).elem_index_ << endl;
		// (*itn).PrintNodes();
		// cout << "###################\n" << endl;
	}

	// FULL GRID 
	// simulateFEM(mainGrid.gridNodes, mainGrid.gridElems, nb_nodes, K_Txt, F_Txt);

	// --------------------------------------
	// TESTING REMOVAL OF CUBES FROM A LIST 
	// --------------------------------------

	// int cubesToRemove[3] = {1,4,5};
	// removeCubes(cubesToRemove, 3, mainGrid.gridElems);


	// vector<Elem>::iterator fin = mainGrid.gridElems.begin();
	// for (; fin != mainGrid.gridElems.end(); fin++){
	// 	cout << "Final element index: \n" << (*fin).elem_index_ << endl ;
	// 	// (*itn).PrintNodes();
	// 	// cout << "###################\n" << endl;
	// }
	// --------------------------------------
	// CHECK REMOVAL OF CUBES FROM A LIST 
	// --------------------------------------

	// // run through each node vertex and check in each element the node exists there
	// for now i try to run thorugh removing all cubes from the top layer. So there is n^2 cubes to remove
	// in each iteration it removes one cube.


	// ======================================
	// --------------------------------------
	// 1. REMOVING CUBES ONE BY ONE
	// --------------------------------------
	// ======================================


	// for (int iteCube = 0; iteCube < N*N*N - N*N; iteCube++){
	// 	cout << "========= REMOVING CUBE!!!! =========\n" ;

	// 	mainGrid.gridElems.erase(mainGrid.gridElems.end()-6,mainGrid.gridElems.end());

	// 	// ============================
	// 	// check after cube removal
	// 	// ============================

	// 	// vector<Elem>::iterator itn = mainGrid.gridElems.begin();
	// 	// for (; itn != mainGrid.gridElems.end(); itn++){
	// 	// 	cout << "This element index is: \n" << (*itn).elem_index_ << endl ;
	// 	// 	// cout << "###################\n" << endl;
	// 	// }
	// 	simulateFEM(mainGrid.gridNodes, mainGrid.gridElems, nb_nodes, K_Txt, F_Txt);	
	// }
	// --------------------------------------
	// ======================================



	// ======================================
	// --------------------------------------
	// 2. REMOVING LAYER ONE BY ONE
	// --------------------------------------
	// ======================================
	
	
	// for (int iteCube = 0; iteCube < N-1; iteCube++){
	// 	cout << "========= REMOVING CUBE!!!! =========\n" ;

	// 	mainGrid.gridElems.erase(mainGrid.gridElems.end()-6*N*N,mainGrid.gridElems.end());

	// 	// ============================
	// 	// check after cube removal
	// 	// ============================

	// 	// vector<Elem>::iterator itn = mainGrid.gridElems.begin();
	// 	// for (; itn != mainGrid.gridElems.end(); itn++){
	// 	// 	cout << "This element index is: \n" << (*itn).elem_index_ << endl ;
	// 	// 	// cout << "###################\n" << endl;
	// 	// }
	// 	simulateFEM(mainGrid.gridNodes, mainGrid.gridElems, nb_nodes, K_Txt, F_Txt);	
	// }
	// --------------------------------------
	// ======================================




	// ======================================
	// --------------------------------------
	//
	//
	// 3. REMOVING A VERTICAL AND HORIZONTAL HOLE 
	//
	//
	// --------------------------------------
	// ======================================
	// set<int> finalGridIndices;
	// // Vertical hole, given a given first layer index. 
	// // this is for N = 5
	// int hole_indices[4] = {8,9,13,14};
	// // int hole_indices[4] = {1,2,6,7};
	// // int hole_indices[4] = {8,9,13,14};
	// int cubesToRemove[20] = {8, 33, 58, 83, 108, 9, 34, 59, 84, 109, 13, 38, 63, 88, 113, 14, 39, 64, 89, 114};
	// createVerticalHole(hole_indices, 4, N, mainGrid.gridElems);
	// // removeCubes(cubesToRemove, 20, mainGrid.gridElems);
	// simulateFEM(mainGrid.gridNodes, mainGrid.gridElems, nb_nodes, K_Txt, F_Txt);
	// vector<Elem>::iterator fin = mainGrid.gridElems.begin();
	// for (; fin != mainGrid.gridElems.end(); fin++){
	// 	// cout << "Final element index: \n" << (*fin).elem_index_ << endl ;
	// 	// (*itn).PrintNodes();
	// 	// cout << "###################\n" << endl;
	// 	finalGridIndices.insert((*fin).elem_index_);
	// }

	// set<int>::iterator itr;
	// cout << "Final grid elements: ====== " << endl;
	// // Displaying set elements
	// // for (itr = toAdd.begin(); itr != toAdd.end(); itr++) 
	// for (itr = finalGridIndices.begin(); itr != finalGridIndices.end(); itr++) 
	// {
	// 	cout << *itr << " ";
	// 	elemListTxt << *itr << endl;
	// }
	// elemListTxt << "===\n";
	// --------------------------------------
	// ======================================
	





	// ======================================
	// --------------------------------------
	//
	//
	// 4. CREATE STICKS GIVEN INDICES . I try to make a standalone hole with thickness
	//
	//
	// --------------------------------------
	// ======================================
	// ======================================
	set<int> finalGridIndices;
	set<int> toAdd;
	vector<Elem> thisGridElems;	
	int nb_base = 12;
	// this little one works with N = 5
	// int base_indices[nb_base] = {7,8,9,12,14,17,18,19};

	// this little one works with N = 8, 28 pieces
	// int base_indices[nb_base] = {11,12,13,14,18,19,20,21,22,23,26,27,30,31,34,35,38,39,42,43,44,45,46,47,51,52,53,54};
	// for N = 8, 12 pieces of hole meat
	int base_indices[nb_base] = {19,20,21,22, 27,30,35,38, 43,44,45,46};
	int hole_indices[4] = {28,29,36,37};
	for(int i = 0; i < nb_base; i++){
		for(int j = 0; j < N; j++){
			cout << "Adding making the hole meat. ";
			toAdd.insert(base_indices[i]+j*N*N);
		}
	}



	// check 
	set<int>::iterator itr;
	cout << "====== To add based on base indices ====== " << endl;
	// Displaying set elements 
	// for (itr = toAdd.begin(); itr != toAdd.end(); itr++) 
	for (itr = toAdd.begin(); itr != toAdd.end(); itr++) 
	{
		cout << *itr << " ";
		elemListTxt << *itr << endl;
	}

	addCube(toAdd,12,mainGrid.gridElems, thisGridElems);
	simulateFEM(mainGrid.gridNodes, thisGridElems, nb_nodes, K_Txt, F_Txt);

	// ======================================
	// BELOW THIS IS A FAILED EXPLORATION 
	






	// cout << "\n====== Now pass the toAdd into the method ====== " << endl;
	// // this is just the normal "Hole meat". works 
	// // addCube(toAdd, 28*N, mainGrid.gridElems, thisGridElems);
	// // simulateFEM(mainGrid.gridNodes, thisGridElems, nb_nodes, K_Txt, F_Txt);

	// // after making the hole thingy, i would like to randomly generate stuff around it! 
	// // to do that, we just have to add extra stuff into the toAdd set. 
	// vector<Elem> currGridElems;
	// for (int iteCube = 0; iteCube < 10; iteCube++){
	// 	cout << "===== TEST CHOOSING BY RANDOM ===== " << endl;
	// 	// currGridElems.insert(currGridElems.end(),thisGridElems.begin(), thisGridElems.end());
	// 	int loop_counter = 0;

	// 	// 10 is the extra we wanna add. 28 is making the hole meat. N is the hole.
	// 	int hole_meat = nb_base*N;
	// 	int hole = 4*N;
	// 	// int nb_cubes = N*N*N-3;
	// 	// int nb_cubes = N*N*N;
	// 	// set <int> extraToAdd;
	// 	// IT IS IN this list that you include in the roots so dont include the hole! 
	// 	int randomToAdd[hole_meat + 20] = {1};
	// 	// store data of if the cubes are already taken or nah. 
	// 	int takenCubes[N*N*N+1] = {0};

	// 	// cout  << "nb_cubes: " << nb_cubes << endl;

	// 	// we have to say that the hole thingy is all taken liao. 
	// 	// accounting for the meat making the hole
	// 	int i = 0;
	// 	for (itr = toAdd.begin(); itr != toAdd.end(); itr++) 
	// 	{
	// 		// cout << "i = " <<  i << "; (*itr) = " << *itr << " ";
	// 		cout << *itr << ", ";
	// 		takenCubes[*itr] = 1;
	// 		randomToAdd[i] = *itr;
	// 		i+= 1;
	// 	}
	// 	cout << "total number making up hole meat: " << i << endl;
	// 	// accounting for the hole - we dont add the hole to toAdd.Only need to keep track of it being taken.
	// 	for(int j = 0; j < N; j++){
	// 		for (int k = 0; k < 4; k++){
	// 			cout << "Accounting for hole: " << hole_indices[k]+j*N*N <<  endl;
	// 			takenCubes[hole_indices[k]+j*N*N] =1;
	// 			// randomToAdd[i] = hole_indices[k]+j*N*N;
	// 			// i += 1;
	// 		}
	// 		// takenCubes[28+j*N*N] =1;
	// 		// takenCubes[29+j*N*N] =1;
	// 		// takenCubes[36+j*N*N] =1;
	// 		// takenCubes[37+j*N*N] =1;
	// 		// randomToAdd[i] = 28+j*N*N;
	// 		// randomToAdd[i] = 29+j*N*N;
	// 		// randomToAdd[i] = 36+j*N*N;
	// 		// randomToAdd[i] = 37+j*N*N;
	// 	}

	// 	cout << "Accounted for the previous cube making the hole " << endl;
	// 	// select first root from any of the taken ones. At this point, theres i things taken. 
	// 	// int root = rand()%(N*N) + 1;
	// 	int root, random_index, original_root;
	// 	// int original_root = root; 
	// 	// randomToAdd[0] = root;
	// 	// cout << "First root: " << root << endl;
	// 	for(int i = hole_meat; i < hole_meat + 20; i++){
	// 		cout << "Entering loop to choose next cubes, i =  " << i << " ..................... " << endl;
	// 		// generate the next cube. first, pick a new root out of all the already picked cubes. 
	// 		// if (i == 0){random_index = 0;} 
	// 		// else{
	// 		// 	random_index = rand()%i;
	// 		// 	// if the random index chosen for i != 0 is such that the randomtoadd is 0, then choose another number. 
	// 		// 	while(randomToAdd[random_index] == 0){
	// 		// 		random_index = rand()%i;
	// 		// 	}
	// 		// }
	// 		random_index = rand()%i;
	// 		cout << "Random index: ====" << random_index << endl;
	// 		// if the random index chosen for i != 0 is such that the randomtoadd is 0, then choose another number. 
	// 		loop_counter = 0;
	// 		while(randomToAdd[random_index] == 0 or randomToAdd[random_index]%(N*N) == 28 or randomToAdd[random_index]%(N*N) == 29 or randomToAdd[random_index]%(N*N) == 36 or randomToAdd[random_index]%(N*N) == 37){
	// 			cout << "Stuck in loop:( " << endl;
	// 			if(loop_counter = 5){
	// 				random_index = 0;
	// 				break;
	// 			}
	// 			if(loop_counter > 3){
	// 				random_index = i-1; 
	// 			}
	// 			else{
	// 				random_index = (random_index+1)%i ;
	// 			}
				
	// 			loop_counter += 1;
	// 		}
	// 		cout << "Final random index: ====" << random_index << endl;
	// 		root = randomToAdd[random_index];
	// 		// temporary placeholder. 
	// 		original_root = root; 
	// 		cout << "================================" << endl;
	// 		cout << "    Selected for next root: " << root << endl;
	// 		cout << "================================" << endl;
	// 		cout << "root%N: ==== " << root%N << endl;
	// 		// the root that is selected is obviously taken alr. 
	// 		//now depending on where root is, we choose the next cube by chance. also need to check if the cube is same as before a not.
	// 		if (root%N == 0){
	// 			cout << "root%N: ==== " << root%N << endl;
	// 			cout << "Root is at right side" << endl;

	// 			// top right corner
	// 			if ((root-1)%(N*N) < N){
	// 				cout << "Root is at right side, top" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "((root-1) - (root-1)%(N*N))/(N*N) == 0" << endl;
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,D or Higher
	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
				
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "((root-1) - (root-1)%(N*N))/(N*N) == N-1" << endl;
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,D or Lower

	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can L,D,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;


	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}

	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 

	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
					
	// 			}
	// 			// bottom right corner
	// 			else if ((root-1)%(N*N) >= N*(N-1)){
	// 				cout << "Root is at right side, bottom" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,U or Higher
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;

	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if(choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,U or Lower
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can L,U,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}

	// 			else{
	// 				cout << "Root is at right side, center" << endl;
	// 				// right and not a corner 
	// 				// choose up or down or left 
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,U,D or Higher
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,U, D or Lower
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter >= 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can L,U, D or Lower or higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%5; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 5){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter >= 0){ 
	// 							choice = (choice + 1)%5;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 4){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}

	// 		}	
	// 		// LEFT SIDE (SYMMETRIC TO RIGHT I JUST REPLACED ALL R WITH L )
	// 		else if (root%N == 1){
	// 			cout << "Root is at left side" << endl;

	// 			// top left corner
	// 			if ((root-1)%(N*N) < N){
	// 				cout << "Root is at left side, top" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "((root-1) - (root-1)%(N*N))/9 == 0" << endl;
	// 					cout << "Root is at first floor" << endl;
	// 					// choose R,D or Higher
	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
				
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "((root-1) - (root-1)%(N*N))/9 == N-1" << endl;
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose R,D or Lower

	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can R,D,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;


	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}

	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 

	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
					
	// 			}
	// 			// bottom left corner
	// 			else if ((root-1)%(N*N) >= N*(N-1)){
	// 				cout << "Root is at left side, bottom" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,U or Higher
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;

	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if(choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,U or Lower
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can R,U,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}

	// 			else{
	// 				cout << "Root is at left side, center" << endl;
	// 				// left and not a corner 
	// 				// choose up or down or right 
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose R,U,D or Higher
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose R,U, D or Lower
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can R,U, D or Lower or higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%5; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 5){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%5;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 4){
	// 							root += N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}
	// 		}

	// 		// TOP AND BOTTOM CENTERS
	// 		// ======================
	// 		else if((root-1)%(N*N) < N){
	// 			cout << "Root is at top, center" << endl;
	// 			// top with no corners
	// 			// choose left rigt or down
	// 			if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 				cout << "Root is at first floor" << endl;
	// 				// choose L,R,D or Higher
    //                 int choice = rand()%4; // select 0 or 1 
	// 				cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                 // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root += N;
	// 					}
	// 					if (choice == 3){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);

	// 			}
	// 			else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 				cout << "Root is at highest floor" << endl;
	// 				// choose L,R,D or Lower
    //                 int choice = rand()%4; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root += N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}
    //                     loop_counter += 1; 

    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 			// otherwise, means can L,R,D or Lower or higher
	// 			else{ 
	// 				cout << "Not top or first layer" << endl;
    //                 int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}
	// 					if (choice == 2){
	// 						root += N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}						
	// 					if (choice == 4){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 		}
	// 		// BOTTOM
	// 		else if((root-1)%(N*N) >= N*(N-1)){
	// 			cout << "Root is at bottom, center" << endl;
	// 			// bottom with no corners
	// 			// choose left rigt or up
	// 			if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 				cout << "Root is at first floor" << endl;
	// 				// choose L,R,U or Higher
	// 				int choice = rand()%4; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root -= N;
	// 					}
	// 					if (choice == 3){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 			else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 				cout << "Root is at highest floor" << endl;
	// 				// choose L,R,U or Lower
    //                 int choice = rand()%4; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root -= N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 			// otherwise, means can L,R,U or Lower or higher
	// 			else{ 
	// 				cout << "Not top or first layer" << endl;
	// 				int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}
	// 					if (choice == 2){
	// 						root -= N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}						
	// 					if (choice == 4){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
	// 			}
	// 		}
	// 		else{
	// 			cout << "Root is at center" << endl;
				
	// 			// choose L,R,U,D
	// 			if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 				cout << "Root is at first floor" << endl;
	// 				// choose L,R,U,D or Higher
	// 				int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}	
	// 					if (choice == 2){
	// 						root -= N;
	// 					}					
	// 					if (choice == 3){
	// 						root += N;
	// 					}
	// 					if (choice == 4){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
    //                 // 
	// 			}
	// 			else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 				cout << "Root is at highest floor" << endl;
	// 				// choose L,R,U,D or Lower
	// 				int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root -= N;
	// 					}					
	// 					if (choice == 3){
	// 						root += N;
	// 					}
	// 					if (choice == 4){
	// 						root -= N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
    //                 // 
	// 			}
	// 			// otherwise, means can L,R,U,D,Lower,Higher
	// 			else{ 
	// 				cout << "Not top or first layer" << endl;
	// 				int choice = rand()%6; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 6){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%6;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}
	// 					if (choice == 2){
	// 						root -= N;
	// 					}					
	// 					if (choice == 3){
	// 						root += N;
	// 					}
	// 					if (choice == 4){
	// 						root -= N*N;
	// 					}						
	// 					if (choice == 5){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 		}


	// 		takenCubes[root] = 1;
	// 		cout << "-----------------------" << endl;
	// 		cout << "Next cube selected is: " << root << endl;
	// 		cout << "-----------------------" << endl;

	// 		// the one that has zero alr in the spot means it was unsuccessful. 
	// 		// if (randomToAdd[i]!= 0){randomToAdd[i] = root;}
	// 		// if(randomToAdd[i] == 0){
	// 		// 	cout << "!!! ==== This index got a repeated cube." << endl;
	// 		// }
	// 		// else{
	// 		// 	cout << "!!! ==== Unique cube." << endl;
				
	// 		// }
	// 		int checkRoot = root%64;
	// 		if(checkRoot == 28 or checkRoot == 29 or checkRoot == 36 or checkRoot == 37){
	// 			cout << "!!! ==== This index is part of a hole" << endl;
	// 		}
	// 		else{
	// 			randomToAdd[i] = root;
	// 			toAdd.insert(root);
	// 		}
			
	// 		// LR UD = 1,2,3,4. LR is +- 1, UD is +-N
	// 		// at right side. 
			
	// 	}
	// 	cout << "--------- Summary of extra cubes to add for this iteration ---------" << endl;
	// 	for(int i = hole_meat; i < hole_meat + 20; i++){
	// 		// cout << "RandomToAdd[" << i << "]: " << randomToAdd[i] << endl;  
	// 		cout <<  randomToAdd[i] << ", ";  
	// 	}

	// 	set<int>::iterator itr;
	// 	cout << "\nSet toAdd: ====== " << endl;
	// 	// Displaying set elements
	// 	int toAddSize = 0;
	// 	for (itr = toAdd.begin(); itr != toAdd.end(); itr++) 
	// 	{
	// 		cout << *itr << " ";
	// 		toAddSize += 1;
	//      //elemListTxt << *itr << endl;
	// 	}
	// 	cout << "\nnb_cubes, toAddSize = " << hole_meat + 20 << ", " << toAddSize << endl;
	// 	addCube(toAdd, toAddSize, mainGrid.gridElems, currGridElems);
	// 	// interested to check what is the currGridElems. 
	// 	vector<Elem>::iterator fin = currGridElems.begin();
	// 	for (; fin != currGridElems.end(); fin += 6){
	// 		cout << "Final element index: \n" << (*fin).elem_index_ << endl ;
	// 		// (*itn).PrintNodes();
	// 		// cout << "###################\n" << endl;
	// 		finalGridIndices.insert((*fin).elem_index_);
	// 	}

	// 	simulateFEM(mainGrid.gridNodes,currGridElems, nb_nodes, K_Txt, F_Txt);	
	// 	cout << "Done with the simulation" << endl;
	//  	// elemListTxt << "===\n";
	// 	currGridElems.clear();
	// }







	// addCube(toAdd, nb_base*N, mainGrid.gridElems, thisGridElems);
	// simulateFEM(mainGrid.gridNodes, thisGridElems, nb_nodes, K_Txt, F_Txt);
	// elemListTxt << "===\n" << endl;		






	// for (int iteCube = 0; iteCube < N-1; iteCube++){
	// 	cout << "========= REMOVING CUBE!!!! =========\n" ;

	// 	mainGrid.gridElems.erase(mainGrid.gridElems.end()-6*N*N,mainGrid.gridElems.end());

	// 	// ============================
	// 	// check after cube removal
	// 	// ============================

	// 	// vector<Elem>::iterator itn = mainGrid.gridElems.begin();
	// 	// for (; itn != mainGrid.gridElems.end(); itn++){
	// 	// 	cout << "This element index is: \n" << (*itn).elem_index_ << endl ;
	// 	// 	// cout << "###################\n" << endl;
	// 	// }
	// 	simulateFEM(mainGrid.gridNodes, mainGrid.gridElems, nb_nodes, K_Txt, F_Txt);	
	// }



	// ======================================
	// --------------------------------------
	// 3. CHOOSING A SET OF CUBES BY RANDOM
	// --------------------------------------
	// ======================================

	// works ok, picks a different root each time. PROBLEM WITH THE PYTHON CODE NOT SAVING THE FILES. !! 17 march
	// if ((root-1) - (root-1)%(N*N))/9 = 0, means is first layer cannot go down. If = N-1, means last layer cannot go up. 
	// implementing the algo of choose root each time 

	// vector<Elem> currGridElems;	
	// for (int iteCube = 0; iteCube < 10; iteCube++){
	// 	cout << "===== TEST CHOOSING BY RANDOM ===== " << endl;
	// 	int loop_counter = 0;
	// 	// pick a root from the first floor. 
	// 	// how many cubes will be adding for the object
	// 	// put extra cos keep getting like duplicates 
	// 	int nb_cubes = N*N*N-5;
	// 	// int nb_cubes = N*N*N-3;
	// 	// int nb_cubes = N*N*N;
	// 	set <int> toAdd;
	// 	int randomToAdd[nb_cubes] = {1};
	// 	// store data of if the cubes are already taken or nah. 
	// 	int takenCubes[N*N*N] = {0};
	// 	int random_index = 0;
	// 	int again = 0;

	// 	// select first root from first floor. 
	// 	int root = rand()%(N*N) + 1;
	// 	int original_root = root; 
	// 	toAdd.insert(root);
	// 	randomToAdd[0] = root;
	// 	takenCubes[root] = 1;
	// 	cout << "First root: " << root << endl;
	// 	for(int i = 1; i < nb_cubes; i++){
	// 		again = 0;
	// 		cout << "Entering loop to choose next cubes, i =  " << i << " ..................... " << endl;
	// 		// generate the next cube. first, pick a new root out of all the already picked cubes. 
	// 		// if (i == 0){random_index = 0;} 
	// 		// else{
	// 		// 	random_index = rand()%i;
	// 		// 	// if the random index chosen for i != 0 is such that the randomtoadd is 0, then choose another number. 
	// 		// 	while(randomToAdd[random_index] == 0){
	// 		// 		random_index = rand()%i;
	// 		// 	}
	// 		// }
	// 		random_index = rand()%i;
	// 		cout << "Random index: ====" << random_index << endl;
	// 		// if the random index chosen for i != 0 is such that the randomtoadd is 0, then choose another number. 
	// 		loop_counter = 0;
	// 		while(randomToAdd[random_index] == 0){
	// 			cout << "Stuck in loop:( " << endl;
	// 			if(loop_counter = 5){
	// 				random_index = 0;
	// 				break;
	// 			}
	// 			if(loop_counter > 3){
	// 				random_index = i-1; 
	// 			}
	// 			else{
	// 				random_index = (random_index+1)%i;
	// 			}
				
	// 			loop_counter += 1;
	// 		}
	// 		cout << "Final random index: ====" << random_index << endl;
	// 		root = randomToAdd[random_index];
	// 		// temporary placeholder. 
	// 		original_root = root; 
	// 		cout << "================================" << endl;
	// 		cout << "    Selected for next root: " << root << endl;
	// 		cout << "================================" << endl;
	// 		cout << "root%N: ==== " << root%N << endl;
	// 		// the root that is selected is obviously taken alr. 
	// 		//now depending on where root is, we choose the next cube by chance. also need to check if the cube is same as before a not.
	// 		if (root%N == 0){
	// 			cout << "root%N: ==== " << root%N << endl;
	// 			cout << "Root is at right side" << endl;

	// 			// top right corner
	// 			if ((root-1)%(N*N) < N){
	// 				cout << "Root is at right side, top" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "((root-1) - (root-1)%(N*N))/(N*N) == 0" << endl;
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,D or Higher
	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
				
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "((root-1) - (root-1)%(N*N))/(N*N) == N-1" << endl;
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,D or Lower

	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can L,D,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;


	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}

	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 

	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
					
	// 			}
	// 			// bottom right corner
	// 			else if ((root-1)%(N*N) >= N*(N-1)){
	// 				cout << "Root is at right side, bottom" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,U or Higher
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;

	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if(choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,U or Lower
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can L,U,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}

	// 			else{
	// 				cout << "Root is at right side, center" << endl;
	// 				// right and not a corner 
	// 				// choose up or down or left 
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,U,D or Higher
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,U, D or Lower
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter >= 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can L,U, D or Lower or higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%5; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 5){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter >= 0){ 
	// 							choice = (choice + 1)%5;
	// 						}
	// 						if (choice == 0){
	// 							root -= 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 4){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}

	// 		}	
	// 		// LEFT SIDE (SYMMETRIC TO RIGHT I JUST REPLACED ALL R WITH L )
	// 		else if (root%N == 1){
	// 			cout << "Root is at left side" << endl;

	// 			// top left corner
	// 			if ((root-1)%(N*N) < N){
	// 				cout << "Root is at left side, top" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "((root-1) - (root-1)%(N*N))/9 == 0" << endl;
	// 					cout << "Root is at first floor" << endl;
	// 					// choose R,D or Higher
	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
				
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "((root-1) - (root-1)%(N*N))/9 == N-1" << endl;
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose R,D or Lower

	// 					int choice = rand()%3; 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}

	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can R,D,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;


	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}

	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 

	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
					
	// 			}
	// 			// bottom left corner
	// 			else if ((root-1)%(N*N) >= N*(N-1)){
	// 				cout << "Root is at left side, bottom" << endl;
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose L,U or Higher
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;

	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if(choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root += N*N;
	// 						}
	// 						loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose L,U or Lower
	// 					int choice = rand()%3; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 3){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%3;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can R,U,Lower or Higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root -= N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}

	// 			else{
	// 				cout << "Root is at left side, center" << endl;
	// 				// left and not a corner 
	// 				// choose up or down or right 
	// 				if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 					cout << "Root is at first floor" << endl;
	// 					// choose R,U,D or Higher
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root += N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 					cout << "Root is at highest floor" << endl;
	// 					// choose R,U, D or Lower
	// 					int choice = rand()%4; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 4){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%4;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}						
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 				// otherwise, means can R,U, D or Lower or higher
	// 				else{ 
	// 					cout << "Not top or first layer" << endl;
	// 					int choice = rand()%5; // select 0 or 1 
	// 					cout << "choice = " << choice << endl;
	// 					loop_counter = 0;
	// 					do{
	// 						// EXIT LOOP 
	// 						if (loop_counter >= 5){
	// 							cout << "Hit max choice permutation." << endl;
	// 							randomToAdd[i] = 0;
	// 							break;
	// 						}

	// 						root = original_root; 
	// 						if (loop_counter > 0){ 
	// 							choice = (choice + 1)%5;
	// 						}
	// 						if (choice == 0){
	// 							root += 1;
	// 						} 
	// 						if (choice == 1){
	// 							root += N;
	// 						}
	// 						if (choice == 2){
	// 							root -= N;
	// 						}
	// 						if (choice == 3){
	// 							root -= N*N;
	// 						}						
	// 						if (choice == 4){
	// 							root += N*N;
	// 						}
    //                         loop_counter += 1; 
	// 					}
	// 					while (takenCubes[root] == 1);
						
	// 				}
	// 			}
	// 		}

	// 		// TOP AND BOTTOM CENTERS
	// 		// ======================
	// 		else if((root-1)%(N*N) < N){
	// 			cout << "Root is at top, center" << endl;
	// 			// top with no corners
	// 			// choose left rigt or down
	// 			if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 				cout << "Root is at first floor" << endl;
	// 				// choose L,R,D or Higher
    //                 int choice = rand()%4; // select 0 or 1 
	// 				cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                 // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root += N;
	// 					}
	// 					if (choice == 3){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);

	// 			}
	// 			else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 				cout << "Root is at highest floor" << endl;
	// 				// choose L,R,D or Lower
    //                 int choice = rand()%4; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root += N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}
    //                     loop_counter += 1; 

    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 			// otherwise, means can L,R,D or Lower or higher
	// 			else{ 
	// 				cout << "Not top or first layer" << endl;
    //                 int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}
	// 					if (choice == 2){
	// 						root += N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}						
	// 					if (choice == 4){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 		}
	// 		// BOTTOM
	// 		else if((root-1)%(N*N) >= N*(N-1)){
	// 			cout << "Root is at bottom, center" << endl;
	// 			// bottom with no corners
	// 			// choose left rigt or up
	// 			if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 				cout << "Root is at first floor" << endl;
	// 				// choose L,R,U or Higher
	// 				int choice = rand()%4; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root -= N;
	// 					}
	// 					if (choice == 3){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 			else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 				cout << "Root is at highest floor" << endl;
	// 				// choose L,R,U or Lower
    //                 int choice = rand()%4; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 4){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%4;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root -= N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 			// otherwise, means can L,R,U or Lower or higher
	// 			else{ 
	// 				cout << "Not top or first layer" << endl;
	// 				int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}
	// 					if (choice == 2){
	// 						root -= N;
	// 					}
	// 					if (choice == 3){
	// 						root -= N*N;
	// 					}						
	// 					if (choice == 4){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
	// 			}
	// 		}
	// 		else{
	// 			cout << "Root is at center" << endl;
				
	// 			// choose L,R,U,D
	// 			if (((root-1) - (root-1)%(N*N))/(N*N) == 0){
	// 				cout << "Root is at first floor" << endl;
	// 				// choose L,R,U,D or Higher
	// 				int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}	
	// 					if (choice == 2){
	// 						root -= N;
	// 					}					
	// 					if (choice == 3){
	// 						root += N;
	// 					}
	// 					if (choice == 4){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
    //                 // 
	// 			}
	// 			else if (((root-1) - (root-1)%(N*N))/(N*N) == N-1){
	// 				cout << "Root is at highest floor" << endl;
	// 				// choose L,R,U,D or Lower
	// 				int choice = rand()%5; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 5){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%5;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}						
	// 					if (choice == 2){
	// 						root -= N;
	// 					}					
	// 					if (choice == 3){
	// 						root += N;
	// 					}
	// 					if (choice == 4){
	// 						root -= N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
    //                 // 
	// 			}
	// 			// otherwise, means can L,R,U,D,Lower,Higher
	// 			else{ 
	// 				cout << "Not top or first layer" << endl;
	// 				int choice = rand()%6; // select 0 or 1 
    //                 cout << "choice = " << choice << endl;
	// 				loop_counter = 0;
    //                 do{
    //                     // EXIT LOOP 
    //                     if (loop_counter >= 6){
	// 						cout << "Hit max choice permutation." << endl;
    //                         randomToAdd[i] = 0;
    //                         break;
    //                     }

    //                     root = original_root; 
    //                     if (loop_counter > 0){ 
    //                         choice = (choice + 1)%6;
    //                     }
    //                     if (choice == 0){
	// 						root -= 1;
	// 					} 
	// 					if (choice == 1){
	// 						root += 1;
	// 					}
	// 					if (choice == 2){
	// 						root -= N;
	// 					}					
	// 					if (choice == 3){
	// 						root += N;
	// 					}
	// 					if (choice == 4){
	// 						root -= N*N;
	// 					}						
	// 					if (choice == 5){
	// 						root += N*N;
	// 					}
    //                     loop_counter += 1; 
    //                 }
	// 				while (takenCubes[root] == 1);
                    
	// 			}
	// 		}


	// 		takenCubes[root] = 1;
	// 		cout << "-----------------------" << endl;
	// 		cout << "Next cube selected is: " << root << endl;
	// 		cout << "-----------------------" << endl;

	// 		// the one that has zero alr in the spot means it was unsuccessful. 
	// 		// if (randomToAdd[i]!= 0){randomToAdd[i] = root;}
	// 		if(randomToAdd[i] == 0){
	// 			cout << "!!! ==== This index got a repeated cube." << endl;
	// 		}
	// 		else{
	// 			cout << "!!! ==== Unique cube." << endl;
	// 		}
	// 		randomToAdd[i] = root;
	// 		toAdd.insert(root);
	// 		// LR UD = 1,2,3,4. LR is +- 1, UD is +-N
	// 		// at right side. 
			
	// 	}
	// 	cout << "--------- Summary of cubes to add for this iteration ---------" << endl;
	// 	for(int i = 0; i < nb_cubes; i++){
	// 		// cout << "RandomToAdd[" << i << "]: " << randomToAdd[i] << endl;  
	// 		cout <<  randomToAdd[i] << ", ";  
	// 	}

	// 	set<int>::iterator itr;
	// 	cout << "Set toAdd: ====== " << endl;
	// 	// Displaying set elements
	// 	for (itr = toAdd.begin(); itr != toAdd.end(); itr++) 
	// 	{
	// 		cout << *itr << " ";
	//      elemListTxt << *itr << endl;
	// 	}
	// 	addCube(toAdd, nb_cubes, mainGrid.gridElems, currGridElems);

	// 	// currGridElems.insert(currGridElems.end(), mainGrid.gridElems.begin()+(randomToAdd[i]-1)*6,mainGrid.gridElems.begin()+(randomToAdd[i]-1)*6+6);

	// 	// vector<Elem>::iterator checkCurr = currGridElems.begin();
	// 	// for (; checkCurr != currGridElems.end(); checkCurr++){
	// 	// 	cout << "Curr element index: \n" << (*checkCurr).elem_index_ << endl ;
	// 	// }
	// 	simulateFEM(mainGrid.gridNodes,currGridElems, nb_nodes, K_Txt, F_Txt);	
	// 	cout << "Done with the simulation" << endl;
	//  	elemListTxt << "===\n";
	// 	currGridElems.clear();
	// }
	





	
	
	// cout << "============= Summary =============" << endl; 
	// cout << "Number of nodes: " << nb_nodes << endl; 
	// cout << "Global stiffness matrix: \n" << K << endl;
	// cout << "Global force vector: \n" << F << endl;


	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	cout << "Time difference (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0  << endl;
	// cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << endl;
	// cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << endl;	// time taken to generate and store the mesh and perform the simulation
	// cout << "============= Summary After Dirichlet conditions =============" << endl; 


	// cout << "Number of nodes: " << nb_nodes << endl; 
	// cout << "Global stiffness matrix: \n" << K << endl;
	// cout << "Global force vector: \n" << F << endl;

	// MatrixXd K_modified = MatrixXd::Zero(3*nb_nodes-3*nb_nodes_bound, 3*nb_nodes-3*nb_nodes_bound);
	// VectorXd F_modified = VectorXd::Zero(3*nb_nodes-3*nb_nodes_bound);

	


	nodeListTxt.close();
	elemListTxt.close();
	solnTxt.close();
	K_Txt.close();
	F_Txt.close();


	// =====================
	// RUN PYTHON FILE 
	// =====================
    // FILE* file;
    // int argc;
    // wchar_t* argv[1];
	// argv[0] = "solveFEM.py";


    // argc = 3;
    // argv[1] = "-m";
    // argv[2] = "/tmp/targets.list";

    // Py_SetProgramName("solveFEM.py");
    // 


    // FILE* file = fopen("solveFEM.py","r");

    // // Py_SetProgramName(argv[0]);



	// Py_Initialize();
    // PyRun_SimpleFile(file, "solveFEM.py");
    // Py_Finalize();





	return 0;
}



// we need a method that can remove an array of indexed cubes from the mainGridElem.
// CUBES TO REMOVED NEEDS TO BE SORTED!!!!!
void removeCubes(int* cubesToRemove, int size, vector<Elem>& vectElem){
	int count = 0;
	cout << "SIZE IN REMOVE CUBES: " << size << endl;
	for (int i = 0; i < size; i++){
		// the values in cubesToRemove is indexed starting from 1 
		cout << "Removing i = " << i << " cube: " << cubesToRemove[i] << endl;
		// cout << (vectElem.begin()+(cubesToRemove[i]-1)*6 - count*6).elem_index << endl;
		// cout << *(vectElem.begin()+(cubesToRemove[i]-1)*6 + 6 - count*6).elem_index << endl;
		vectElem.erase(vectElem.begin()+(cubesToRemove[i]-1)*6 - count*6,vectElem.begin()+(cubesToRemove[i]-1)*6 + 6 - count*6);
		// CEHCKING 
		vector<Elem>::iterator init = vectElem.begin();
		cout << "Removed cube for iteration i = " << i << endl;
		for (; init != vectElem.end(); init+= 6){
			cout << (*init).elem_index_ << ", " ;
			// (*itn).PrintNodes();
			// cout << "###################\n" << endl;
		}
		count += 1;
		cout << "\n" << endl;
	}
}

void createVerticalHole(int* hole_indices, int nb_holes, int N, vector<Elem>& vectElem){
	int cubes[N*nb_holes];
	cout << "Creating list of cubes for hole" << endl;
	for(int j = 0; j < N; j++){
		// for each hole index
		for(int i = 0; i < nb_holes; i++){
			cubes[i + j*nb_holes] = hole_indices[i] + N*N*j;
			// cout << cubes[i] ;
		}
	}
	// CHECKING
	cout << "Checking the list for the cubes" << endl;
	for(int i = 0; i < N*nb_holes; i++){
		cout << cubes[i] << ", " ;
	}
	
	cout << "\nCreated hole, calling removeCubes" << endl;
	removeCubes(cubes, N*nb_holes, vectElem);
}

void addCube(set<int>& cubesToAdd, int size, vector<Elem>& vectRef, vector<Elem>& vectElem){

	set<int>::iterator itr;
	// Displaying set elements
	for (itr = cubesToAdd.begin(); itr != cubesToAdd.end(); itr++) 
	{
		// cout << *itr << " ";
		if((*itr) == 0) {continue;}
		vectElem.insert(vectElem.end(), vectRef.begin()+((*itr)-1)*6,vectRef.begin()+((*itr)-1)*6+6);
	}
	// for(int i = 0; i < size; i++){
	// 	// this accounts for the case when there are no other good options. 
	// 	if(cubesToAdd[i] == 0) {continue;}
	// 	vectElem.insert(vectElem.end(), vectRef.begin()+(cubesToAdd[i]-1)*6,vectRef.begin()+(cubesToAdd[i]-1)*6+6);
	// }
}

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
				int elem_index = i*N + (j+1) + N*N*level;
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
					// cout << "CHECKING VECTOR3d pos_loc: \n" << pos_loc << endl; 
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
					elem_loc.elem_index_ = elem_index;

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




void simulateFEM(vector<Node>& thisGridNodes,vector<Elem>& thisGridElems, int nb_nodes, ofstream& K_Txt, ofstream& F_Txt){
	MatrixXd K = MatrixXd::Zero(3*nb_nodes, 3*nb_nodes);
	VectorXd F = VectorXd::Zero(3*nb_nodes);

	cout << "SIMULATING\n" << endl;
	//cout << "SIMULATING\n RUNNING THROUGH NODES FOR FILLING MATRIX K\n" << endl;
	int globalNodeIndex_i = 0;
	vector<Node>::iterator itnode = thisGridNodes.begin();

	// CALCULATIONS 
	for (; itnode !=  thisGridNodes.end();itnode++)
	{	
		globalNodeIndex_i = itnode - thisGridNodes.begin();
		// cout << "Looking at a current node with global index " << globalNodeIndex_i << endl;  
		// cout << (*itnode) << endl;
		// cout << "========================================== " << endl;
		// cout << "========================================== " << endl;


		//given this current node, we run through the elements. 
		vector<Elem>::iterator itelem = thisGridElems.begin();
		for (; itelem !=  thisGridElems.end();itelem++)
		{	
			// first we need to check if this node (*itnode) is in this element. 
			//---------------------------------------
			// run through nodes of this element and check if itnode is inside
			bool inElem = false;
			int pos_inElem_i = 0;
			vector<Node>::iterator itelemnode = (*itelem).elemNodes.begin();
			for (; itelemnode != (*itelem).elemNodes.end();itelemnode++)
			{	

				// GLOBAL NODE i IS IN THE ELEMENT.   
				// ===============================
				// ===============================
				if((*itelemnode).position_ == (*itnode).position_){
					inElem = true; 
					pos_inElem_i = itelemnode - (*itelem).elemNodes.begin();
					
					// now run through all the other nodes of this element and get the local and global index (j)
					int pos_inElem_j = 0; int globalNodeIndex_j = 0;
					vector<Node>::iterator itelemnode2 = (*itelem).elemNodes.begin();
					for(; itelemnode2 != (*itelem).elemNodes.end();itelemnode2++)
					{

						// cout << "j in this element: " << (*itelemnode2).position_ << endl;
						// get local and global positions of this node. 
						pos_inElem_j = itelemnode2 - (*itelem).elemNodes.begin();
						// get global node position of itelemnode2
						vector<Node>::iterator itnode2 =  thisGridNodes.begin();
						for (; itnode2 !=  thisGridNodes.end();itnode2++)
						{
							if((*itnode2).position_ == (*itelemnode2).position_){
								globalNodeIndex_j = itnode2 -  thisGridNodes.begin();
							}
						}
						
						// F(3*globalNodeIndex_i+2) -= 9.81/5;
						F(3*globalNodeIndex_i+2) -= 9.81;
						for (int i_local = 0; i_local < 3; i_local++){
							for (int j_local = 0; j_local < 3; j_local++){

								K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) += (*itelem).Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
								
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
	// WRITE INTO FILE 
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

	K_Txt <<"===\n" ;
	F_Txt <<"===\n" ;

}





