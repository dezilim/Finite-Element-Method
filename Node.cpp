#include "Node.hpp"

Node::Node(Vector3d& position, int elem_index){
    this->position_ = position;
    this->elem_index_ = elem_index;

}


ostream& operator<<(ostream& os, const Node& N) {
	//os << "--node-position--\n" << N.position_ << "\n" << endl;
	os << "--node-position--\n" << N.position_ << endl;
	return os;
}