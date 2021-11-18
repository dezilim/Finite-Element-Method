#include "Node.hpp"

Node::Node(Matrix& position, int elem_index){
    this->position_ = position;
    this->elem_index_ = elem_index;

}


ostream& operator<<(ostream& os, const Node& N) {
	os << "[" << N.position_ << "]" << endl;
	return os;
}