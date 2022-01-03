//#include <ctime>
#include "Matrix.hpp"
#include "Node.hpp"
#include "Elem.hpp"
#include "Grid.hpp"
// #include "Particule.hpp"
// #include "Boite.hpp"

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include <Eigen/Dense>
#include<Eigen/SparseLU> 
#include <vector>
#include <array>
#include <math.h>

using Eigen::MatrixXd;
using namespace std;

// declare functions 
vector<Elem> generateGrid(int N);
Grid makeGrid(int N);

int main( int argc, char * argv[] ) {

	// ------------------------------ GENERATE GRID AND CHECK ------------------------------ 



	Grid mainGrid = makeGrid(3);
	cout << "Grid generated.\n" << endl;
	int nb_nodes = (3+1)*(3+1)*2;

	cout << "Running through global node list in mainGrid to check: ...\n" << endl;
	int count = 1;
	vector<Node>::iterator itglobal = mainGrid.gridNodes.begin();
	//cout << (*itglobal) << endl;
	for (; itglobal != mainGrid.gridNodes.end();itglobal++)
	{	
		// print node in myNodes
		cout << "node number: " << count << endl;
		cout << (*itglobal) << endl;
		count++;
	}

	// cout << "Running through elements list in mainGrid to check: ...\n" << endl;
	// vector<Elem>::iterator itelem = mainGrid.gridElems.begin();
	// //cout << (*itglobal) << endl;
	// for (; itelem != mainGrid.gridElems.end();itelem++)
	// {	
	// 	cout << "this element " << endl;
	// 	// print node in myNodes
	// 	//cout << (*itglobal) << endl;
	// 	(*itelem).PrintNodes();
	// }

	// ------------------------------ CREATE LOCAL MATRICES ------------------------------ 
    /*    create Me   */
	Matrix I3 = identity(3); Matrix I6 = identity(6);
	Matrix Mright = joinUD(joinLR(I3,I3),joinLR(I3,I3));
	Matrix Mleft = joinUD(joinLR(2*I3,I3),joinLR(I3,2*I3));
	cout << "Matrix test, Mright :\n" << Mright << endl;
	cout << "Matrix test, Mleft = 2*I6*Mright:\n" << Mleft << endl;
	Matrix Me = joinUD(joinLR(Mleft,Mright), joinLR(Mright,Mleft));
	cout << "Elementary matrix Me:\n" << Me << endl;


	/*   create Ke   */
	double E = 110; double  v = 0.34; double G = 46;
	double c11 = E*(1-v)/((1-2*v)*(1+v)); 
	double c12 = E*v/((1-2*v)*(1+v));
	//double F = (c11 - c12)/2;
	cout << "c11 : " << c11 << endl;
	cout << "c12 : " << c12 << endl;
	//cout << "G : " << F << endl;

	Matrix C = c11*I6;
	cout << "Matrix test, c11*I6 :\n" << C << endl;
	for(int i = 3; i < 6; i++){
		C.val[i][i] = G; 
	} 
	C.val[0][1] = c12; C.val[0][2] = c12; C.val[1][1] = c12;
	C.fillSym();
	cout << "Matrix test, C :\n" << C << endl;

	// creation of matrix B = LN 
	Matrix e_11 = elemSelec(3,1,1); Matrix e_22 = elemSelec(3,2,2);Matrix e_33 = elemSelec(3,3,3);
	Matrix Low1(3,3,-1); Matrix Low2(3,3,0);Matrix Low3(3,3,0);Matrix Low4(3,3,0);
	Low1 += identity(3);
	cout << "Low1: \n" << Low1 << endl;
	Low2.val[1][2] = 1; Low3.val[0][2] = 1; Low4.val[0][1] = 1;

	Low2.fillSym();Low3.fillSym();Low4.fillSym();
	Matrix LowTotal = joinLR(joinLR(Low1, Low2), joinLR(Low3,Low4));
	Matrix High = joinLR(joinLR(-1*I3,e_11),joinLR(e_22,e_33));
	Matrix B = joinUD(High,LowTotal);
	cout << "B = LN :\n" << B << endl;

	Matrix Ke = (B.transpose())*C*B;
	cout << "Ke :\n" << Ke << endl;


    int nb_fixedNodes = 4;
    // ------------------------------ CREATE GLOBAL KE ------------------------------ 
    Matrix Ke_global;
    int elem_index = 0;
	// pseudo code 
	// for each element, get the coordinates of Xi and then calculate the Jacobian then get the element matrice Ke-global 
    // for each element of type Elem, it has a vector list of Nodes that are arranged (indexed). Can access from 0 to 4
    vector<Elem>::iterator itn = mainGrid.gridElems.begin();
	for (; itn != mainGrid.gridElems.end(); itn++){
		cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
		cout << "--------------------------------------- " << endl;
        //(*itn).PrintNodes();
        // this is of type Node
        cout << (*itn).elemNodes[0] << endl;
        cout << (*itn).elemNodes[1] << endl;
        cout << (*itn).elemNodes[2] << endl;
        // Initialise Jacobian matrix 
		//Matrix Je(3,3);
        // assign to X0 the coordinates of the first node of the element
        Matrix X0((*itn).elemNodes[0].position_.transpose());
        cout << "X0:" << X0 << endl;
		Matrix X1((*itn).elemNodes[1].position_.transpose());
		cout << "X1:" << X1 << endl;
		Matrix X2((*itn).elemNodes[2].position_.transpose());
		cout << "X2:" << X2 << endl;
		Matrix X3((*itn).elemNodes[3].position_.transpose());
		cout << "X3:" << X3 << "\n" << endl;

		Matrix diffX1X0(X1-X0);Matrix diffX2X0(X2-X0);Matrix diffX3X0(X3-X0);
		cout << "X1-X0:" << diffX1X0 << endl;
		cout << "X2-X0:" << diffX2X0 << endl;
		cout << "X3-X0:" << diffX3X0 << endl;

		Matrix Je = joinLR(diffX1X0,joinLR(diffX2X0,diffX3X0));
		Ke_global = (Je.getDet3())*Ke;
		cout << "Summary for element "<< elem_index + 1 << endl;
		cout << "Jacobian: " << Je << "\n" << endl;
		cout << "Det (Je) : " << Je.getDet3() << endl;
		cout << "Ke global :" << Ke_global  << endl;

        elem_index++;

	}
	
	Matrix K(3*nb_nodes, 3*nb_nodes, 0); Matrix F(3*nb_nodes, 1, 0);
	cout << "Initialising global stiffness matrix K... " << K << endl;
	// run through each node vertex and check in each element the node exists there


	cout << "RUNNING THROUGH NODES FOR FILLING MATRIX K\n" << endl;
	int globalNodeIndex_i = 0;
	vector<Node>::iterator itnode = mainGrid.gridNodes.begin();
	for (; itnode != mainGrid.gridNodes.end();itnode++)
	{	
		// print node in myNodes
		//cout << "node number: " << count << endl;
		// this is
		globalNodeIndex_i = itnode - mainGrid.gridNodes.begin();
		cout << "Looking at a current node with global index " << globalNodeIndex_i << endl;  
		cout << (*itnode) << endl;
		cout << "========================================== " << endl;
		cout << "========================================== " << endl;

		
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
					cout << "========================================== " << endl;
					cout << (*itnode) << "is found in the element at position " << pos_inElem_i << endl;
					(*itelem).PrintNodes();
					cout << " ==== Run through other nodes in this elem... ==== " << endl;
					// now run through all the other nodes of this element and get the local and global index (j)
					int pos_inElem_j = 0; int globalNodeIndex_j = 0;
					vector<Node>::iterator itelemnode2 = (*itelem).elemNodes.begin();
					for(; itelemnode2 != (*itelem).elemNodes.end();itelemnode2++)
					{
						if (!((*itelemnode2).position_ == (*itnode).position_)){
							cout << "j in this element: " << (*itelemnode2).position_ << endl;
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
							cout << "global j:" << globalNodeIndex_j <<endl;
							cout << "local j :" << pos_inElem_j << endl;
							cout << " ------ end of resume for j ------ " << endl;
							// little loop because we need to fill in 3 by 3 matrix 
							F.val[3*globalNodeIndex_j+2][0] += 9.81*16/6;
							for (int i_local = 0; i_local < 3; i_local++){
								for (int j_local = 0; j_local < 3; j_local++){
									// i and j for the Ke global should be the positions found for the nodes
									// for example you can have i = 5, but the global element matrix has only Kij for i,j = 1,2,3,4
									// and each Kij is a 3 by 3 matrix
								  
									K.val[3*globalNodeIndex_i + i_local][3*globalNodeIndex_j + j_local] += Ke_global.val[3*pos_inElem_i + i_local][3*pos_inElem_j + j_local];
									
									//cout << "updating K : ..............." << K << endl;

								}
							}


						}
					}

				}
			}
			// if (*itnode) is indeed in this element, we have to 
			

		}
	}
	
	
	cout << "============= Summary =============" << endl; 
	cout << "Number of nodes: " << nb_nodes << endl; 
	cout << "Global stiffness matrix: " << K << endl;
	cout << "Global force vector: " << F << endl;

	// run through each node, check if z= 0. If yes, modify that node index row. For example, if node i = 6 is a node of z= 0, 
	// then we need u_i = 0, so we need Kii = 1, Fi = 0, Kij = 0 if i != j
	cout << "Taking into account dirichlet conditions ..." << endl;
	vector<Node>::iterator itDir = mainGrid.gridNodes.begin();
	int i_Dir;
	for (; itDir != mainGrid.gridNodes.end();itDir++)
	{	
		cout << "--------------" << endl;
		cout << (*itDir) << endl;
		// check if the z component of this node's position is 0.
		cout << (*itDir).position_.val[0][2] << endl;
		if ((*itDir).position_.val[0][2] == 0) {
			i_Dir = itDir - mainGrid.gridNodes.begin();
			cout << "i_Dir = " << i_Dir << endl;
			cout << "This node is at the boundary, need to modify the Kij andFi" << endl;

			for (int i_local = 0; i_local <3; i_local++ ){
				K.changeAllInRow(3*i_Dir+i_local,0);
				K.changeAllInCol(3*i_Dir+i_local,0);
				K.val[3*i_Dir+i_local][3*i_Dir+i_local] =  1;
				F.val[3*i_Dir+i_local][0] = 0;
			}
		}

	}
	cout << "============= Summary After Dirichlet conditions =============" << endl; 
	cout << "Number of nodes: " << nb_nodes << endl; 
	cout << "Global stiffness matrix: " << K << endl;
	cout << "Global force vector: " << F << endl;

	// Solve KU = F
	// Matrix linProb = joinLR(K,-1*F);
	// cout << "linProb: "<< linProb << "\n\n" << endl;


	cout << "Copying matrix K to testM " << endl;
	// just testing out the matrix class from Eigen 
	// Eigen::MatrixXf is the default general matrix one but for sparse matrices it is better to use Eigen::SparseMatrix

	int Ksize = (int) nb_nodes*3;
	Eigen::SparseMatrix<double> testK(Ksize, Ksize);
	for(int i = 0; i< nb_nodes*3; i++){
		for(int j = 0; j< nb_nodes*3; j++){
			testK.coeffRef(i,j) = K.val[i][j];
			//cout << testK.coeffRef(i,j) << " " ;
		}
		//cout << "\n" << endl;
	}

	cout << "Copying matrix K to testM_matrix " << endl;

	Eigen::MatrixXf testK_matrix(Ksize, Ksize);
	for(int i = 0; i< nb_nodes*3; i++){
		for(int j = 0; j< nb_nodes*3; j++){
			testK_matrix(i,j) = K.val[i][j];
			//cout << testK_matrix(i,j) << " " ;
		}
		//cout << "\n" << endl;
	}

	cout << "Copying matrix F to testF " << endl;

	Eigen::VectorXf testF(Ksize);
	for(int i = 0; i< nb_nodes*3; i++){
		testF(i) = F.val[i][0];
		//cout << testF(i) << " " ;
	}
	
	Eigen::VectorXf x_sol = testK_matrix.colPivHouseholderQr().solve(testF);
    cout << "The solution is:\n" << x_sol << endl;

	// // try to solve using SparseLU factorisation
	//    =================================
	// Eigen::SparseLU<Eigen::SparseMatrix<double> > solverA;
	// testK.makeCompressed();
	// solverA.analyzePattern(testK);
	// solverA.factorize(testK);

	// if(solverA.info()!=Eigen::Success) {
	// 	std::cout << "Oh: Very bad" <<"\n";
	// }
	// else{
	// 	std::cout<<"okay computed"<<"\n";
	// }
	// Eigen::VectorXd solnew = solverA.Eigen::SparseLU<Eigen::SparseMatrix<double> >::solve(testF);

	
	// linProb is a matrix, no problem. But for some reason the gaussian elim is wrong, i think at some points they divide by zero 
	// Matrix linProbElim = linProb.gaussianElim();
	// cout << "linProb.gaussianElim() :" <<linProbElim<< "\n\n" << endl;
	// cout << "linProb.gaussianElim.sol :" <<linProbElim.gaussianSol() << endl;

	return 0;
}

//
vector<Elem> generateGrid(int N){
    vector<Elem> elements_grid_list;
	// vector<Node> nodes_grid_list;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
			// Going through a cubic element
            cout << "i, j = " << i << ", " << j << endl;
            int elem_index = i*N + (j+1);
			cout << "*******************************" << endl;

            cout << "*** Elem_cube_index = "<< elem_index << " *** " <<  endl;
			cout << "*******************************" << endl;
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
				cout << "=======" << endl;

				modo += k%2 ;
				j_loc += pow(-1, modo)*(k%2);
				i_loc += pow(-1, modo)*((k-1)%2);
				cout << "Node "<< k +1 << ": " << j_loc << ", " <<  i_loc << ", "  << z_loc << endl;
				//cout << "Node "<< k +1<<  "[modo, x, y, z]: "<<modo << ", " << j_loc << ", " <<  i_loc << ", "  << z_loc << endl;
				Matrix pos_loc = rowVect(3,0); pos_loc.val[0][0] = (1.0)*j_loc; pos_loc.val[0][1] = (1.0)*i_loc; pos_loc.val[0][2] = (1.0)*z_loc;
				
				Node node_loc(pos_loc, elem_index);
				nodes_cube_list.push_back(node_loc);
				//nodes_grid_list.push_back(node_loc);

			}
			cout << "=======" << endl;
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
				elements_grid_list.push_back(elem_loc); 

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
    return elements_grid_list;
}







Grid makeGrid(int N){
	int Size = N;
	Grid grid(Size);
    //vector<Elem> elements_grid_list;
	//vector<Node> nodes_grid_list;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
			// Going through a cubic element
            cout << "i, j = " << i << ", " << j << endl;
            int elem_index = i*N + (j+1);
			cout << "*******************************" << endl;

            cout << "*** Elem_cube_index = "<< elem_index << " *** " <<  endl;
			cout << "*******************************" << endl;
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
				cout << "=======" << endl;

				modo += k%2 ;
				j_loc += pow(-1, modo)*(k%2);
				i_loc += pow(-1, modo)*((k-1)%2);
				cout << "Node "<< k +1 << ": " << j_loc << ", " <<  i_loc << ", "  << z_loc << endl;
				//cout << "Node "<< k +1<<  "[modo, x, y, z]: "<<modo << ", " << j_loc << ", " <<  i_loc << ", "  << z_loc << endl;
				Matrix pos_loc = rowVect(3,0); pos_loc.val[0][0] = (1.0)*j_loc; pos_loc.val[0][1] = (1.0)*i_loc; pos_loc.val[0][2] = (1.0)*z_loc;
				
				// !!!!!!!!!!!!!!!!!!!!!!!!!!
				// tried to make a workaround by creating two nodes, one to push into nodes cube list and another to put into gridNodes. 
				Node node_loc(pos_loc, elem_index); Node node_glob(pos_loc, elem_index);
				cout << "Pushing local node into nodes_cube_list\n" << endl;
				nodes_cube_list.push_back(node_loc);
				// pushing all nodes generated for a cube into this list. In fact this overlaps.
				cout << "Pushing global node into gridNodes_list\n" << endl;
				// we want to push it into gridNodes only if we cannot find the element in the existing list. 
				bool found_global_node = false; 
				vector<Node>::iterator itcheck = grid.gridNodes.begin();
				//cout << (*itglobal) << endl;
				for (; itcheck != grid.gridNodes.end();itcheck++)
				{	
					// print node in myNodes
					//cout << "node number: " << count << endl;
					cout << (*itcheck) << endl;
					if ((*itcheck).position_ == pos_loc)
					{
						cout << "node already registered in global node list" << endl;
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
			cout << "=======" << endl;
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
	
    return grid;
}