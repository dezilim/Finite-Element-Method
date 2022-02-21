#include "Grid.hpp"
#include "Elem.hpp"
#include "Node.hpp"
#include <vector>
#include <iostream>

using namespace std;

Grid::Grid(int Size){
    //node = NULL;
    this->gridSize = Size;
	//node = NULL;
	//elemNodes.resize(0);

}

// void Grid::PrintNodes(){
//     vector<Node>::iterator itn = gridNodes.begin();
// 	for (; itn != gridNodes.end(); itn++)
// 	{
// 		//cout << "Elem: printing nodes"  << "...\n" ;
// 		cout << *itn << endl;
// 	}
// }
