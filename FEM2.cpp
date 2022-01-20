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

	string node_fileName =  "nodeList.txt";
	string soln_fileName =  "solution.txt";
	ofstream nodeListTxt, solnTxt;
	nodeListTxt.open("nodeList.txt");
	solnTxt.open("solution.txt");


	// ------------------------------ GENERATE GRID AND CHECK ------------------------------ 


	// NOTE: Change the number in makeGrid and also in the nb_nodes 
	Grid mainGrid = makeGrid(2);
	int N = 2;
	int nb_nodes = pow((N+1),3);

	cout << "Grid of size " << N << " with "<< nb_nodes << " nodes generated.\n" << endl;


	cout << "Running through global node list in mainGrid to check: ...\n" << endl;
	int count = 1;
	vector<Node>::iterator itglobal = mainGrid.gridNodes.begin();
	//cout << (*itglobal) << endl;
	for (; itglobal != mainGrid.gridNodes.end();itglobal++)
	{	
		// print node in myNodes
		// cout << "node number: " << count << endl;
		// cout << (*itglobal) << endl;
		nodeListTxt << (*itglobal) ;
		nodeListTxt << "\n";

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

    /*     create Ke = Ve B^t C B     */
	double E = 110; double  v = 0.34; double G = 46;
	double c11 = E*(1-v)/((1-2*v)*(1+v)); 
	double c12 = E*v/((1-2*v)*(1+v));
	// //double F = (c11 - c12)/2;
	cout << "c11 : " << c11 << endl;
	cout << "c12 : " << c12 << endl;
	// //cout << "G : " << F << endl;
    MatrixXd C = c11*I6; 

	// Matrix C = c11*I6;
	// cout << "Matrix test, c11*I6 :\n" << C << endl;
	for(int i = 3; i < 6; i++){
		C(i,i) = G; 
	} 
	C(0,1) = c12; C(0,2) = c12; C(1,2) = c12;
    C(1,0) = c12; C(2,0) = c12; C(2,1) = c12;
	// C.fillSym();
	cout << "Matrix test, C :\n" << C << endl;   


	// creation of matrix B = LN 
	// ========================
	MatrixXd B1up = -1*I3;
	cout << "B1down :\n" << B1up << endl;

	// MatrixXd B1down = MatrixXd::Identity(3,3);
	// for (int i = 0; i< 3; i++){
	// 	for(int j = 0; j < 3; j++){
	// 		B1down(i,j) -= 1;
	// 	}
	// }
	// cout << "B1down :\n" << B1down << endl;

	MatrixXd Bleft(B1up.rows() + B1up.rows(), B1up.cols());
	Bleft << B1up, B1up ;
	Bleft(3,1) = -1;
	Bleft(5,0) = -1;
	cout << "Bleft :\n" << Bleft << endl;

	MatrixXd Bright = MatrixXd::Zero(6,9);
	Bright(0,0) = 1; 
	Bright(1,4) = 1;
	Bright(2,8) = 1;
	Bright(3,1) = 1;
	Bright(3,3) = 1;
	Bright(4,7) = 1;
	Bright(5,2) = 1;
	Bright(5,6) = 1;

	cout << "Bright :\n" << Bright << endl; 

	MatrixXd B(Bleft.rows(), Bleft.cols()+Bright.cols());
    B << Bleft, Bright; 
	cout << "B = LN :\n" << B << endl; 

	MatrixXd Ke = (B.transpose())*C*B/6;
	cout << "Ke :\n" << Ke << endl; 


 



    int nb_fixedNodes = 4;
    // // ------------------------------ CREATE GLOBAL KE ------------------------------ 
    MatrixXd Ke_global;
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


        // Initialise Jacobian matrix 
		//Matrix Je(3,3);
        // assign to X0 the coordinates of the first node of the element
        Vector3d X0((*itn).elemNodes[0].position_.transpose());
        // cout << "X0:" << X0 << endl;
		Vector3d X1((*itn).elemNodes[1].position_.transpose());
		// cout << "X1:" << X1 << endl;
		Vector3d X2((*itn).elemNodes[2].position_.transpose());
		// cout << "X2:" << X2 << endl;
		Vector3d X3((*itn).elemNodes[3].position_.transpose());
		// cout << "X3:" << X3 << "\n" << endl;

		Vector3d diffX1X0 = X1-X0;
		Vector3d diffX2X0 = X2-X0;
		Vector3d diffX3X0 = X3-X0;
		// cout << "X1-X0:\n" << diffX1X0 << endl;
		// cout << "X2-X0:\n" << diffX2X0 << endl;
		// cout << "X3-X0:\n" << diffX3X0 << endl;

		MatrixXd Je(3,3);
		Je << diffX1X0, diffX2X0, diffX3X0;
		Ke_global = (Je.determinant())*Ke;
		cout << "Summary for element "<< elem_index + 1 << endl;
		cout << "Jacobian: \n" << Je << endl;
		cout << "Det (Je) : \n" << Je.determinant() << endl;
		cout << "Ke global : \n" << Ke_global << endl;

        elem_index++;

	}
	//SparseMatrix<double> K(3*nb_nodes, 3*nb_nodes);
	MatrixXd K = MatrixXd::Zero(3*nb_nodes, 3*nb_nodes);
	VectorXd F = VectorXd::Zero(3*nb_nodes);

	// Assemble F 
	// for force in z - direction, 3*globalNodeIndex_j+2. !!!!!!!!!!!!!!!!!!!
	for (int i_loc = 0; i_loc < nb_nodes; i_loc++){
		// exaggerating to check results 
		//F(3*i_loc+2) -= 9.81;
		//F(3*i_loc+2) -= 9.81/24;
	}
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
					//(*itelem).PrintNodes();
					// cout << " ==== Run through other nodes in this elem... ==== " << endl;
					
					// now run through all the other nodes of this element and get the local and global index (j)
					int pos_inElem_j = 0; int globalNodeIndex_j = 0;
					vector<Node>::iterator itelemnode2 = (*itelem).elemNodes.begin();
					for(; itelemnode2 != (*itelem).elemNodes.end();itelemnode2++)
					{
						if (!((*itelemnode2).position_ == (*itnode).position_)){
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
							// cout << "global j:" << globalNodeIndex_j <<endl;
							// cout << "local j :" << pos_inElem_j << endl;
							// cout << " ------ end of resume for j ------ " << endl;


							// little loop because we need to fill in 3 by 3 matrix 
							// // for z = 0, 3*globalNodeIndex_j+2. !!!!!!!!!!!!!!!!!!!
							F(3*globalNodeIndex_j) -= 2210*9.81/24;
							for (int i_local = 0; i_local < 3; i_local++){
								for (int j_local = 0; j_local < 3; j_local++){
									// i and j for the Ke global should be the positions found for the nodes
									// for example you can have i = 5, but the global element matrix has only Kij for i,j = 1,2,3,4
									// and each Kij is a 3 by 3 matrix
								  
									//K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) = K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) +  Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
									K(3*globalNodeIndex_i + i_local,3*globalNodeIndex_j + j_local) += Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
									//+= Ke_global(3*pos_inElem_i + i_local,3*pos_inElem_j + j_local);
									
									//cout << "updating K : ..............." << K << endl;

								}
							}


						}
					}

				}
			}			
		}
	}

	
	cout << "============= Summary =============" << endl; 
	// cout << "Number of nodes: " << nb_nodes << endl; 
	// cout << "Global stiffness matrix: " << K << endl;
	// cout << "Global force vector: " << F << endl;

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
				K.row(3*i_Dir+2) *= 0; // K.col(3*i_Dir+2) *= 0;
				//K.col(3*i_Dir+i_local) *= 0;
				K(3*i_Dir+2,3*i_Dir+2) =  1;
				F(3*i_Dir+2) = 0;
			}
		}

	}
	cout << "============= Summary After Dirichlet conditions =============" << endl; 
	// cout << "Number of nodes: " << nb_nodes << endl; 
	cout << "Global stiffness matrix: " << K << endl;
	cout << "Global force vector: " << F << endl;









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
	// Eigen::VectorXd x_sol = K.colPivHouseholderQr().solve(F);
    // //cout << "The solution is:\n" << x_sol << endl;
	// for(int i = 0; i < 3*nb_nodes; i++){
	// 	solnTxt<< x_sol(i) ;
	// 	solnTxt << "\n";
	// }
	//cout << "Check solution, KU = :\n" << K*x_sol << endl;





	// ========== SOLVING USING LDLT solver (WORKS but funky) ===========
	// SparseMatrix<double> A(3*nb_nodes,3*nb_nodes);
	// for(int i = 0; i< 3*nb_nodes; i++){
	// 	for(int j = 0; j< 3*nb_nodes; j++){
	// 		if(K(i,j) != 0) A.insert(i,j) = K(i,j);

	// 	}
	// }


	// Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT;
    // solverLDLT.compute(A);
    // Eigen::VectorXd x_sol = solverLDLT.solve(F);


	// //cout << "The solution is:\n" << x_sol << endl;
	// for(int i = 0; i < 3*nb_nodes; i++){
	// 	solnTxt<< x_sol(i) ;
	// 	solnTxt << "\n";
	// }



	// ========== SOLVING USING LU solver 
 	// Eigen::SparseLU<Eigen::SparseMatrix<double>> solverA;
	// solverA.factorize(A);
	// Eigen::VectorXd x_sol = solverA.solve(F);
	// =========

	// ========= SOLVING FOR SPD using Gradient conjugee (iterative algo is reccomended)
	// test by copying NOTE: may be faster to only insert non zero items 

	// VectorXd x_sol(3*nb_nodes);
	// SparseMatrix<double> A(3*nb_nodes,3*nb_nodes);
	// for(int i = 0; i< 3*nb_nodes; i++){
	// 	for(int j = 0; j< 3*nb_nodes; j++){
	// 		if(K(i,j) != 0) A.insert(i,j) = K(i,j);

	// 	}
	// }
	// fill A and b
	// Eigen::ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
	// cg.compute(A);
	// x_sol = cg.solve(F);
	// std::cout << "#iterations:     " << cg.iterations() << std::endl;
	// std::cout << "estimated error: " << cg.error()      << std::endl;
	// // update b, and solve again
	// x_sol = cg.solve(F);

	//cout << "Check solution, KU = :\n" << K*x_sol << endl;

	// ========= SOLVING FOR SPD using GMRES (iterative algo is reccomended)
	// int n = 3*nb_nodes;
	// VectorXd x_sol(n);
	// SparseMatrix<double> A(n,n);
	// for(int i = 0; i< 3*nb_nodes; i++){
	// 	for(int j = 0; j< 3*nb_nodes; j++){
	// 		if(K(i,j) != 0) A.insert(i,j) = K(i,j);
	// 	}
	// }
	// // fill A and b
	// Eigen::GMRES<SparseMatrix<double> > solver(A);
	// x_sol = solver.solve(F);
	// std::cout << "#iterations:     " << solver.iterations() << std::endl;
	// std::cout << "estimated error: " << solver.error()      << std::endl;
	// // update b, and solve again
	// x_sol = solver.solve(F);

	// //cout << "The solution is:\n" << x_sol << endl;
	// for(int i = 0; i < 3*nb_nodes; i++){
	// 	solnTxt<< x_sol(i) ;
	// 	solnTxt << "\n";
	// }



	// ========= SOLVING BY BRUTE FORCE 

	// MatrixXd Kinv = K.inverse();
	// // MatrixXd I = MatrixXd::Identity(3*nb_nodes,3*nb_nodes);
	// // MatrixXd Kinv = K.llt().solve(I);
	// VectorXd x_sol = Kinv*F;
	// cout << "Kinv:\n" << Kinv << endl; 
	// cout << "KinvF = U :\n" << x_sol << endl;
	// for(int i = 0; i < 3*nb_nodes; i++){
	// 	solnTxt<< x_sol(i) ;
	// 	solnTxt << "\n";
	// }

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
	nodeListTxt.close();
	solnTxt.close();
	return 0;
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