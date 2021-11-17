//#include <ctime>
#include "Matrix.hpp"
// #include "Particule.hpp"
// #include "Boite.hpp"

#include <iostream>	// pour pouvoir utiliser des objets ostreram
using namespace std;

int main( int argc, char * argv[] ) {
	//initialisation du générateur aleatoire
	//srand (time(0));
	
	// /*    Proper test of matrix class   */
	// Matrix M1(4,2,3); Matrix M1copy(M1); Matrix I4 = identity(4);
	// Matrix M2(4,2,5); Matrix M3(2,4,4); 
	// cout << "M1(4,2,3): " << M1 << endl;
	// cout << "M1copy(M1):" << M1copy << endl;
	// cout << "I4:" << I4 << endl;
	// cout << "M2(4,2,5): " << M2 << endl;
	// cout << "M3(2,4,4): " << M3 << endl;
	// cout << "----- Test Matrix overload operations  -----" << endl;
	// M1 += M2;
	// cout << "M1 +=  M2:" << M1 << endl;
	// cout << "M1 *= 3:" << (M1 *= 3) << endl;
	// cout << "M1 *= M3:" << (M1 *= M3) << endl; 

	// cout << "----- Test Matrix fillSym and transpose  -----" << endl;
	// I4.val[1][0] = 2; I4.val[2][3] = -3; 
	// cout << "new I4:" << 2*I4 << endl; 
	// Matrix I4T = 2*I4.transpose();
	// cout << "I4.transpose():" << I4T << endl;
	// cout << "I4.fillSym():" << I4T.fillSym() << endl;

	// cout << "----- Test Matrix external functions  -----" << endl;
	// cout << "Mcopy + M2:" << M1copy + M2 << endl;
	// cout << "M2 * M3:" << M2 * M3 << endl;
	// cout << "3*Mcopy:" << 3*M1copy << endl;
	// cout << "Mcopy*3:" << M1copy*3 << endl;

	// // vectors of dim 4 
	// Matrix E_1 = colVect(4,1); Matrix E_1T = rowVect(4,1);
	// cout << "Vector E_1:" << E_1 << endl;
	// cout << "Vector E_1T:" << E_1T << endl;
	// cout << "Vector E_1T:" << E_1.transpose() << endl;

	// // vector matrix multiplication 
	// cout << "E_1T*M1copy:" << E_1T*M1copy << endl;
	// cout << "M1copyT*E_1:" << (M1copy.transpose())*E_1 << endl;

	// // join vector with matrix 
	// cout << "joinLR(E_1,M1copy):"<< joinLR(E_1, M1copy) << endl;

	

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
	//double F = (c11 - c12)/2;
	cout << "c11 : " << c11 << endl;
	cout << "c12 : " << c12 << endl;
	//cout << "G : " << F << endl;

	Matrix C = c11*I6;
	cout << "Matrix test, c11*I6 :" << C << endl;
	for(int i = 3; i < 6; i++){
		C.val[i][i] = G; 
	} 
	C.val[0][1] = c12; C.val[0][2] = c12; C.val[1][1] = c12;
	C.fillSym();
	cout << "Matrix test, C :" << C << endl;

	// creation of matrix B = LN 
	Matrix e_11 = elemSelec(3,1,1); Matrix e_22 = elemSelec(3,2,2);Matrix e_33 = elemSelec(3,3,3);
	Matrix Low1(3,3,-1); Matrix Low2(3,3,0);Matrix Low3(3,3,0);Matrix Low4(3,3,0);
	Low1 += identity(3);
	cout << "Low1: " << Low1 << endl;
	Low2.val[1][2] = 1; Low3.val[0][2] = 1; Low4.val[0][1] = 1;

	Low2.fillSym();Low3.fillSym();Low4.fillSym();
	Matrix LowTotal = joinLR(joinLR(Low1, Low2), joinLR(Low3,Low4));
	Matrix High = joinLR(joinLR(-1*I3,e_11),joinLR(e_22,e_33));
	Matrix B = joinUD(High,LowTotal);
	cout << "B = LN :" << B << endl;

	Matrix Ke = (B.transpose())*C*B;
	cout << "Ke :" << Ke << endl;

	double coords[12][3] = {{-2,-2,-2}, {2,-2,-2}, {2,2,-2}, {-2,2,-2}, {-2,-2,2}, {2,-2,2}, {2,2,2} , {-2,2,2}, {-2,-2,6}, {2,-2,6}, {2,2,6} , {-2,2,6}};
	int elems[12][4] = {{1,2,3,7},{1,6,2,7}, {1,5,6,7}, {1,8,5,7}, {1,4,8,7}, {1,3,4,7},{5,6,7,11}, {5,10,6,11}, {5,9,10,11}, {5,12,9,11}, {5,8,12,11}, {5,7,8,11}};
	int nb_vertices = 12;
	// this may expand later, need a way to detect which nodes are fixed
	int nb_fixedNodes = 4;

	int fixedNodes[4] = {1,2,3,4};

	Matrix Ke_global;
	// pseudo code 
	// for each element, get the coordinates of Xi and then calculate the Jacobian then get the element matrice Ke-global 
	for(int elem_index = 0; elem_index < 12; elem_index++){
		cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
		cout << "--------------------------------------- " << endl;
		int ref_vertex = elems[elem_index][0];
		// Jacobian matrix 
		Matrix Je(3,3);
		// col vect (x0, y0, z0)^T
		Matrix X0 = colVect(3,1);
		X0.val[0][0] = coords[ref_vertex-1][0]; X0.val[1][0] = coords[ref_vertex-1][1]; X0.val[2][0] = coords[ref_vertex-1][2];
		cout << "Initialised Jacobian:" << Je << endl;
		cout << "Reference vector X0:" << X0 <<endl;


		// print the vertices used for each element 
		cout << "Relevant vertices: " << endl;
		for (int vertex_index = 1; vertex_index < 4; vertex_index++){
			Matrix Xi = colVect(3,1);
			int curr_vertex = elems[elem_index][vertex_index];
			//cout << curr_vertex << endl;
			Xi.val[0][0] = coords[curr_vertex-1][0]; Xi.val[1][0] = coords[curr_vertex-1][1]; Xi.val[2][0] = coords[curr_vertex-1][2];

			cout << "Coordinates of current vertex "<< curr_vertex << " : --- "<< endl;
			cout << Xi << endl;
			Matrix diffXiX0 = Xi-X0;
			cout << "Xi - X0:" << diffXiX0 << endl;
			for(int i = 0; i < 3; i++){
				Je.val[i][vertex_index-1] = diffXiX0.val[i][0];
			}
			cout << "Jacobian: " << Je << endl;
		}
		cout << "=================================" << endl;
		cout << "Summary for element "<< elem_index + 1 << endl;
		cout << "=================================" << endl;
		cout << "Jacobian Je: " << Je << endl; 
		cout << "Det (Je) : " << Je.getDet3() << endl;
		Ke_global = (Je.getDet3())*Ke;
		//cout << "Ke global :" <<Ke_global << endl;

	}
	cout << "Ke global :" << Ke_global  << endl;

	// For now, each element has the same Ke global sa they are all quite standard shapes. But in the next step we have to efficiently store this data...


	Matrix K(3*nb_vertices, 3*nb_vertices, 0); Matrix F(3*nb_vertices, 1, 0);
	cout << "Initialising global stiffness matrix K... " << K << endl;
	// run through each node vertex and check in each element the node exists there
	for (int i = 0; i < nb_vertices; i++){
		cout << "------- Node index i :  " << i + 1 << " ------- " << endl;
		cout << "-------------------------------------- " << endl;
		for(int elem_index = 0; elem_index < 12; elem_index++){
		 		cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
		 		cout << "---------------------------------------- " << endl;
				// create pointer to find in this element if vertex i is inside
				int *find_i = std::find(std::begin(elems[elem_index]), std::end(elems[elem_index]), i+1);

				// we will only check for which vertex j is inside if we already  know that i is in the element.
				// this way  we do nont have to check all j for all elements.
				// if we find an element that conntains node i, we need to run through all its nodes (of this element) labelled j and fill in the Kij
				// because if node i is in element elem, then all nodes of this element is a neighour of node i
				if (find_i != std::end(elems[elem_index])) {
					int i_pos = std::distance(elems[elem_index], find_i);
					cout << "Node i = " << i+1 << " is found in element " << elem_index+1 << " at position " << i_pos << endl;
					for (int vertex_index = 0; vertex_index < 4; vertex_index++){
						int j = elems[elem_index][vertex_index]-1;
						int *find_j = std::find(std::begin(elems[elem_index]), std::end(elems[elem_index]), j+1);
						if (find_j != std::end(elems[elem_index])) {
							int j_pos = std::distance(elems[elem_index], find_j);
							cout << "\t Node j = " << j+1 << "at position "<< j_pos << endl;
							// modify part of the global stiffness matrix K
							// ex . i = 0, j = 2
							// for now not so sure about assembling the force vect, i think its similar to assemblinng the K matrix. 
							// dk what exact value but I know that the first entry is a constant. etc fb = (g,0,0) that makes force in the x direction
							F.val[3*j+2][0] += 9.81*16/6;
							for (int i_local = 0; i_local < 3; i_local++){
								for (int j_local = 0; j_local < 3; j_local++){
									// i and j for the Ke global should be the positions found for the nodes
									// for example you can have i = 5, but the global element matrix has only Kij for i,j = 1,2,3,4
									// and each Kij is a 3 by 3 matrix
								  
									K.val[3*i + i_local][3*j + j_local] += Ke_global.val[3*i_pos + i_local][3*j_pos + j_local];
									
									//cout << "updating K : ..............." << K << endl;

								}
							}
						// ex i = 2, j = 5 
						//Kij += 
						// j is the actual value but i is the computer convention value 
						}
					}
				}

			}


	}
	cout << "============= Summary =============" << endl; 
	cout << "Number of elements: " << 12 << endl; 
	cout << "Global stiffness matrix: " << K << endl;
	cout << "Global force vector: " << F << endl;



		// for(int elem_index = 0; elem_index < 12; elem_index++){
		// 	cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
		// 	cout << "--------------------------------------- " << endl;
			

		// 	for (int vertex_index = 1; vertex_index < 4; vertex_index++){
		// 		int curr_vertex = elems[elem_index][vertex_index];
		// 		if (vertex == )

	// }
	// apply dirichlet conditions
	cout << "Taking into account dirichlet conditions ..." << endl;
	for (int vertex_index = 0; vertex_index < nb_fixedNodes; vertex_index++){
		int i = fixedNodes[vertex_index]-1;

		for (int i_local = 0; i_local <3; i_local++ ){
			K.changeAllInRow(3*i+i_local,0);
			K.changeAllInCol(3*i+i_local,0);
			K.val[3*i+i_local][3*i+i_local] =  1;
			F.val[3*i+i_local][0] = 0;
		}

		// cout << "Global stiffness matrix: " << K << endl;


	}
	cout << "============= Summary (after dirichlet conditions) =============" << endl; 
	cout << "Number of elements: " << 12 << endl; 
	cout << "Global stiffness matrix: " << K << endl;
	cout << "Global force vector: " << F << endl;

	// Solve KU = F
	Matrix linProb = joinLR(K,-1*F);
	Matrix linProbElim = linProb.gaussianElim();
	cout << "linProb.gaussianElim.sol :" <<linProbElim.gaussianSol() << endl;

	/*
	// test gaussian elim
	Matrix gTest(3,3,1);
	gTest.val[1][0] = 2;gTest.val[1][1] = 3; gTest.val[1][2] = 7;
	gTest.val[2][1] = 3; gTest.val[2][2] = -2;
	cout << "gTest: " << gTest << endl;
	Matrix fTest(3,1,1);
	fTest.val[0][0] = 3; fTest.val[1][0] = 0; fTest.val[2][0] = 17;
	Matrix linProb = joinLR(gTest, fTest);

	
	cout << "fTest " << fTest << endl;

	cout << "gTest.gaussianElim :" << gTest.gaussianElim() << endl;
	cout << "linProb :" <<linProb << endl;
	cout << "linProb.gaussianElim :" <<linProb.gaussianElim() << endl;
	Matrix linProbElim = linProb.gaussianElim();
	cout << "linProb.gaussianElim.sol :" <<linProbElim.gaussianSol() << endl;
	*/

	




	return 0;
}
