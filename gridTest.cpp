//#include <ctime>
#include "Matrix.hpp"
#include "Node.hpp"
#include "Elem.hpp"
// #include "Particule.hpp"
// #include "Boite.hpp"

#include <iostream>	// pour pouvoir utiliser des objets ostreram
#include <vector>
#include <array>

using namespace std;

int main( int argc, char * argv[] ) {

	// ------------------------------ END OF TESTS ------------------------------ 
	


	

	double coords[12][3] = {{-2,-2,-2}, {2,-2,-2}, {2,2,-2}, {-2,2,-2}, {-2,-2,2}, {2,-2,2}, {2,2,2} , {-2,2,2}, {-2,-2,6}, {2,-2,6}, {2,2,6} , {-2,2,6}};
	int elems[12][4] = {{1,2,3,7},{1,6,2,7}, {1,5,6,7}, {1,8,5,7}, {1,4,8,7}, {1,3,4,7},{5,6,7,11}, {5,10,6,11}, {5,9,10,11}, {5,12,9,11}, {5,8,12,11}, {5,7,8,11}};
	int nb_vertices = 12;
	// this may expand later, need a way to detect which nodes are fixed
	int nb_fixedNodes = 4;

	int fixedNodes[4] = {1,2,3,4};
	// vector of pointers : doesnt work
	//vector<array<int, 3> > NODES;

	// manual. works
	vector<Node> myNodes;

	// create an element of index 5
	Elem myElem(5); 
    //vector<int*> NODES;
    int N = 10;
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
			// Going through a cubic element
            cout << "i, j = " << i << ", " << j << endl;
            int elem_index = i*5 + (j+1);
            cout << "Elem_index = "<< elem_index << endl;
			// this displays the left side of cubic element on grid. We should create 6 tetrahedral elements each 
			// containing respective nodes. 

			/* For each cubic element, we need to create 8 nodes. node1 will be the left bottom, z= 0 */
            //initialise a position for this node: from left bottom (x,y) , ACW, then z = 1 
			Matrix node1_ = rowVect(3,0); node1_.val[0][0] = j; node1_.val[0][1] = i; node1_.val[0][2] = 0; 
			Matrix node2_ = rowVect(3,0); node2_.val[0][0] = j+1; node2_.val[0][1] = i; node2_.val[0][2] = 0; 			
			Matrix node3_ = rowVect(3,0); node3_.val[0][0] = j+1; node3_.val[0][1] = i+1; node3_.val[0][2] = 0; 
			Matrix node4_ = rowVect(3,0); node4_.val[0][0] = j; node4_.val[0][1] = i+1; node4_.val[0][2] = 0; 			
			Matrix node5_ = rowVect(3,0); node5_.val[0][0] = j; node5_.val[0][1] = i; node5_.val[0][2] = 1; 
			Matrix node6_ = rowVect(3,0); node6_.val[0][0] = j+1; node6_.val[0][1] = i; node6_.val[0][2] = 1; 			
			Matrix node7_ = rowVect(3,0); node7_.val[0][0] = j+1; node7_.val[0][1] = i+1; node7_.val[0][2] = 1; 
			Matrix node8_ = rowVect(3,0); node8_.val[0][0] = j; node8_.val[0][1] = i+1; node8_.val[0][2] = 1; 

			Node node1(node1_, elem_index);
			Node node2(node2_, elem_index);
			Node node3(node3_, elem_index);
			Node node4(node4_, elem_index);
			Node node5(node5_, elem_index);
			Node node6(node6_, elem_index);
			Node node7(node7_, elem_index);
			Node node8(node8_, elem_index);

			// make tetrahedral elements with these nodes
			//Elem elem1(); 



			// testing for just creating a node at left bottom for each cubic element, z = 0
			Matrix position = rowVect(3,0);
			position.val[0][0] = j; position.val[0][1] = i; position.val[0][2] = 0; 
			Node node(position, elem_index);

            //node[0] = j; node[1] = i; node[2] = 0;
            //myNodes.push_back(node);

			// PROBLEM WITH THIS LINE 
			myElem.elemNodes.push_back(node);


        }
    }

	myElem.PrintNodes();

	// print the nodes from the manual vector of nodes. No problemo
	// vector<Node>::iterator itn = myNodes.begin();
	// for (; itn != myNodes.end(); itn++)
	// {
	// 	cout << "Node added"  << "..." <<endl;
	// 	cout << *itn << endl;
	// }
    
    // for (vector<array<int, 3> >::iterator it = NODES.begin() ; it != NODES.end(); ++it) {

    //     cout << "----"<< endl;
    //     cout << "\t" << *it->0 << ", " <<  *it->1 << ", " << *it->2;
    //     cout << "}; \n";
    // }
    

    
	// Matrix Ke_global;
	// // pseudo code 
	// // for each element, get the coordinates of Xi and then calculate the Jacobian then get the element matrice Ke-global 
	// for(int elem_index = 0; elem_index < 12; elem_index++){
	// 	cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
	// 	cout << "--------------------------------------- " << endl;
	// 	int ref_vertex = elems[elem_index][0];
	// 	// Jacobian matrix 
	// 	Matrix Je(3,3);
	// 	// col vect (x0, y0, z0)^T
	// 	Matrix X0 = colVect(3,1);
	// 	X0.val[0][0] = coords[ref_vertex-1][0]; X0.val[1][0] = coords[ref_vertex-1][1]; X0.val[2][0] = coords[ref_vertex-1][2];
	// 	cout << "Initialised Jacobian:" << Je << endl;
	// 	cout << "Reference vector X0:" << X0 <<endl;


	// 	// print the vertices used for each element 
	// 	cout << "Relevant vertices: " << endl;
	// 	for (int vertex_index = 1; vertex_index < 4; vertex_index++){
	// 		Matrix Xi = colVect(3,1);
	// 		int curr_vertex = elems[elem_index][vertex_index];
	// 		//cout << curr_vertex << endl;
	// 		Xi.val[0][0] = coords[curr_vertex-1][0]; Xi.val[1][0] = coords[curr_vertex-1][1]; Xi.val[2][0] = coords[curr_vertex-1][2];

	// 		cout << "Coordinates of current vertex "<< curr_vertex << " : --- "<< endl;
	// 		cout << Xi << endl;
	// 		Matrix diffXiX0 = Xi-X0;
	// 		cout << "Xi - X0:" << diffXiX0 << endl;
	// 		for(int i = 0; i < 3; i++){
	// 			Je.val[i][vertex_index-1] = diffXiX0.val[i][0];
	// 		}
	// 		cout << "Jacobian: " << Je << endl;
	// 	}
	// 	cout << "=================================" << endl;
	// 	cout << "Summary for element "<< elem_index + 1 << endl;
	// 	cout << "=================================" << endl;
	// 	cout << "Jacobian Je: " << Je << endl; 
	// 	cout << "Det (Je) : " << Je.getDet3() << endl;
	// 	Ke_global = (Je.getDet3())*Ke;
	// 	//cout << "Ke global :" <<Ke_global << endl;

	// }
	// cout << "Ke global :" << Ke_global  << endl;

	// // For now, each element has the same Ke global sa they are all quite standard shapes. But in the next step we have to efficiently store this data...


	// Matrix K(3*nb_vertices, 3*nb_vertices, 0); Matrix F(3*nb_vertices, 1, 0);
	// cout << "Initialising global stiffness matrix K... " << K << endl;
	// // run through each node vertex and check in each element the node exists there
	// for (int i = 0; i < nb_vertices; i++){
	// 	cout << "------- Node index i :  " << i + 1 << " ------- " << endl;
	// 	cout << "-------------------------------------- " << endl;
	// 	for(int elem_index = 0; elem_index < 12; elem_index++){
	// 	 		cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
	// 	 		cout << "---------------------------------------- " << endl;
	// 			// create pointer to find in this element if vertex i is inside
	// 			int *find_i = std::find(std::begin(elems[elem_index]), std::end(elems[elem_index]), i+1);

	// 			// we will only check for which vertex j is inside if we already  know that i is in the element.
	// 			// this way  we do nont have to check all j for all elements.
	// 			// if we find an element that conntains node i, we need to run through all its nodes (of this element) labelled j and fill in the Kij
	// 			// because if node i is in element elem, then all nodes of this element is a neighour of node i
	// 			if (find_i != std::end(elems[elem_index])) {
	// 				int i_pos = std::distance(elems[elem_index], find_i);
	// 				cout << "Node i = " << i+1 << " is found in element " << elem_index+1 << " at position " << i_pos << endl;
	// 				for (int vertex_index = 0; vertex_index < 4; vertex_index++){
	// 					int j = elems[elem_index][vertex_index]-1;
	// 					int *find_j = std::find(std::begin(elems[elem_index]), std::end(elems[elem_index]), j+1);
	// 					if (find_j != std::end(elems[elem_index])) {
	// 						int j_pos = std::distance(elems[elem_index], find_j);
	// 						cout << "\t Node j = " << j+1 << "at position "<< j_pos << endl;
	// 						// modify part of the global stiffness matrix K
	// 						// ex . i = 0, j = 2
	// 						// for now not so sure about assembling the force vect, i think its similar to assemblinng the K matrix. 
	// 						// dk what exact value but I know that the first entry is a constant. etc fb = (g,0,0) that makes force in the x direction
	// 						F.val[3*j+2][0] += 9.81*16/6;
	// 						for (int i_local = 0; i_local < 3; i_local++){
	// 							for (int j_local = 0; j_local < 3; j_local++){
	// 								// i and j for the Ke global should be the positions found for the nodes
	// 								// for example you can have i = 5, but the global element matrix has only Kij for i,j = 1,2,3,4
	// 								// and each Kij is a 3 by 3 matrix
								  
	// 								K.val[3*i + i_local][3*j + j_local] += Ke_global.val[3*i_pos + i_local][3*j_pos + j_local];
									
	// 								//cout << "updating K : ..............." << K << endl;

	// 							}
	// 						}
	// 					// ex i = 2, j = 5 
	// 					//Kij += 
	// 					// j is the actual value but i is the computer convention value 
	// 					}
	// 				}
	// 			}

	// 		}


	// }
	// cout << "============= Summary =============" << endl; 
	// cout << "Number of elements: " << 12 << endl; 
	// cout << "Global stiffness matrix: " << K << endl;
	// cout << "Global force vector: " << F << endl;



	// 	// for(int elem_index = 0; elem_index < 12; elem_index++){
	// 	// 	cout << "------- Going through element " << elem_index + 1 << " ------- " << endl;
	// 	// 	cout << "--------------------------------------- " << endl;
			

	// 	// 	for (int vertex_index = 1; vertex_index < 4; vertex_index++){
	// 	// 		int curr_vertex = elems[elem_index][vertex_index];
	// 	// 		if (vertex == )

	// // }
	// // apply dirichlet conditions
	// cout << "Taking into account dirichlet conditions ..." << endl;
	// for (int vertex_index = 0; vertex_index < nb_fixedNodes; vertex_index++){
	// 	int i = fixedNodes[vertex_index]-1;

	// 	for (int i_local = 0; i_local <3; i_local++ ){
	// 		K.changeAllInRow(3*i+i_local,0);
	// 		K.changeAllInCol(3*i+i_local,0);
	// 		K.val[3*i+i_local][3*i+i_local] =  1;
	// 		F.val[3*i+i_local][0] = 0;
	// 	}

	// 	// cout << "Global stiffness matrix: " << K << endl;


	// }
	// cout << "============= Summary (after dirichlet conditions) =============" << endl; 
	// cout << "Number of elements: " << 12 << endl; 
	// cout << "Global stiffness matrix: " << K << endl;
	// cout << "Global force vector: " << F << endl;

	// // Solve KU = F
	// Matrix linProb = joinLR(K,-1*F);
	// Matrix linProbElim = linProb.gaussianElim();
	// cout << "linProb.gaussianElim.sol :" <<linProbElim.gaussianSol() << endl;

	// /*
	// // test gaussian elim
	// Matrix gTest(3,3,1);
	// gTest.val[1][0] = 2;gTest.val[1][1] = 3; gTest.val[1][2] = 7;
	// gTest.val[2][1] = 3; gTest.val[2][2] = -2;
	// cout << "gTest: " << gTest << endl;
	// Matrix fTest(3,1,1);
	// fTest.val[0][0] = 3; fTest.val[1][0] = 0; fTest.val[2][0] = 17;
	// Matrix linProb = joinLR(gTest, fTest);

	
	// cout << "fTest " << fTest << endl;

	// cout << "gTest.gaussianElim :" << gTest.gaussianElim() << endl;
	// cout << "linProb :" <<linProb << endl;
	// cout << "linProb.gaussianElim :" <<linProb.gaussianElim() << endl;
	// Matrix linProbElim = linProb.gaussianElim();
	// cout << "linProb.gaussianElim.sol :" <<linProbElim.gaussianSol() << endl;
	// */

	

    


	return 0;
}
