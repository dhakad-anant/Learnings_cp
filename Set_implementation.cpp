#include <bits/stdc++.h>
using namespace std;

/* A node of a BST */
template <typename T>
struct Node {

	T val; /* Node's val */
	Node* left;	/* Left child's Pointer */
	Node* right; /* Right child's Pointer */

public:
	/* To display Inorder traversal of the BST */
	void inorder(Node* root){
		if (root == NULL)return;
		inorder(root->left);
		cout << root->val << " ";
		inorder(root->right);
	}

	/**
		Function to check if BST contains a node with the given val
		
		@param root pointer to the root node
		@param key the val to search
		@return if present -> 1, if not present -> 0
	*/
	int haveNode(Node* root, T key){
		if (root == NULL)return 0;
		int x = root->val == key ? 1 : 0;
		return x | haveNode(root->left, key) | haveNode(root->right, key);
	}

	/**
		Function to insert a node with given val into the BST
		
		@param root pointer to the root node
		@param key the val to insert
		@return root of new bst 
	*/
	Node* insert(Node* root, T key){

		/* When a NULL node is discovered, add the node */
		if (root == NULL){
			Node<T>* newNode = new Node<T>;
			newNode->val = key;
			newNode->left = newNode->right = NULL;
			return newNode;
		}

		/* If val is smaller than the current node, go to the left subtree. */
		if (key < root->val){
			root->left = insert(root->left, key);
			return root;
		}

		/* If val is greater than the current node, traverse the right subtree. */
		else if (key > root->val) {
			root->right = insert(root->right, key);
			return root;
		}
		
		return root;
	}
};

/* Set class (implemented using BST) */
template <typename T>
class Set {

	Node<T>* root; /* pointer to the root the BST tree */
	int size; /* number of elements in the set */

public:
	/* Default constructor initializing the root and size. */
	Set(){
		root = NULL;
		size = 0;
	}

	Set(const Set& s){
		root = s.root;
		size = s.size;
	}

	/**
        function to insert an element to the set.

		@param val the element to add to the set
	*/
	void insert(const T val){
        if (root->haveNode(root, val) == 1)return;

        root = root->insert(root, val);
		size++;
	}

	/**
        This function returns the union of two sets
		
		@param s set to find union with
		@return the union set
	*/
	Set _union(Set& s){
		Set<T> res;

		/* Return 2nd set if 1st set is empty */
		if (root == NULL) return s;
        /* Return 1st set if 2nd set is empty */
		if (s.root == NULL)return *this;

        /* Adding the elements of the 1st set in the resultant set */
		stack<Node<T>*> utilstack;
		utilstack.push(root);

		/* Preorder traversal of the BST */
		while (!utilstack.empty()){
			Node<T>* node;
			node = utilstack.top();
			utilstack.pop();

			// Adding val to the resultant set
			res.insert(node->val);

			if(node->right)utilstack.push(node->right);
			if(node->left)utilstack.push(node->left);
		}

		/* Adding the elements of the 2nd set in the resultant set */
		stack<Node<T>*> utilstack1;
		utilstack1.push(s.root);

		while (!utilstack1.empty()){
			Node<T>* node;
			node = utilstack1.top();
			utilstack1.pop();

            // Adding val to the resultant set
			res.insert(node->val);

			if(node->right)utilstack1.push(node->right);
			if(node->left)utilstack1.push(node->left);
		}

		return res;
	}

	/**
        This function returns the intersection of two sets
		
		@param s the set to find intersection with
		@return the intersection set
	*/
	Set _intersection(Set& s){
		Set<T> res;
		stack<Node<T>*> utilstack;
		utilstack.push(root);

		while (!utilstack.empty()){
			Node<T>* node;
			node = utilstack.top();
			utilstack.pop();
			if (s.contains(node->val)){
				res.insert(node->val);
			}
			if(node->right)utilstack.push(node->right);
			if(node->left)utilstack.push(node->left);
		}
		return res;
	}

	/**
        This function returns the complement of the sets
		
		@param U the universal set
		@return the complement set
	*/
	Set _complement(Set& U){
		return (U - *this);
	}

	/**
		Overloading '-' operator to find the difference of two sets
		
		@param s the set to be subtracted
		@return the Difference set
	*/
	Set operator-(Set& s){
		Set<T> res;
		stack<Node<T>*> utilstack;
		utilstack.push(this->root);

		while (!utilstack.empty()){
			Node<T>* node;
			node = utilstack.top();
			utilstack.pop();
			if (s.contains(node->val) == 0)res.insert(node->val);

			if(node->right)utilstack.push(node->right);
			if(node->left)utilstack.push(node->left);
		}
		return res;
	}

	/**
		Function to check if set are equal
		
		@param s set to check equality with
		@return if equal -> true, if not equal -> false
	*/
	bool operator==(Set& s){
		if (s._size() != size)return false;

		stack<Node<T>*> utilstack;
		utilstack.push(this->root);

		while (!utilstack.empty()){
			Node<T>* node;
			node = utilstack.top();
			utilstack.pop();
			if (!s.contains(node->val))return false;

			if(node->right)utilstack.push(node->right);
			if(node->left)utilstack.push(node->left);
		}
		return true;
	}

	/**
		Function to display the cartesian product of two sets
		
		@param s the set to find product with
	*/
	void printProduct(Set& s){
        int size2 = s._size();
		T *A = _array(), *B = s._array();

		cout << "{ ";
		for(int i = 0; i < size; i++){
			for(int j = 0; j < size2 ; j++){
				cout << "{ " << A[i] << " " << B[j] << " } ";
			}
		}
		cout << "}" << endl;
	}

    /**
		Function to display the power set of the set
	*/
	void printPowerSet(){
		int n = pow(2, size);
		T* A = _array();
		while(n-- > 0){
			cout << "{ ";
			for(int i = 0; i < size; i++){
				if((n & (1 << i)) == 0){
					cout << A[i] << " ";
				}
			}
			cout << "}" << endl;
		}
	}

	/**
		Function to convert the set into an array
		
		@return array of set elements
	*/
	T* _array(){
		T* A = new T[size];
		int i = 0;
		stack<Node<T>*> utilstack;
		utilstack.push(this->root);

		while(!utilstack.empty()){
			Node<T>* node;
			node = utilstack.top();
			utilstack.pop();

			A[i++] = node->val;

			if(node->right)utilstack.push(node->right);
			if(node->left)utilstack.push(node->left);
		}
		return A;
	}

	/**
		Function to check whether the set contains an element
		
		@param key the element to search
		@return if exist -> true, not exist -> false
	*/
	bool contains(T key){
		return root->haveNode(root, key) ? true : false;
	}

	// Function to print the contents of the set
	void printSet(){
		cout << "{ ";
		root->inorder(root);
		cout << "}" << endl;
	}

	/**
        get the current size of the set
		
		@return size of set
	*/
	int _size(){
		return size;
	}
};

int main(){
	Set<int> A;

    /* inserting elements in Set A */
	A.insert(1); A.insert(2); A.insert(3); A.insert(2);

    /* printing the contents of Set A */
	cout << "A = ";
	A.printSet();
	cout << "P(A) = " << '\n';
	A.printPowerSet();

	/* Check if Set A contains some elements */
	cout << "A " << (A.contains(7) ? "contains" : "does not contain") << " 7" << '\n';
	cout << "A " << (A.contains(9) ? "contains" : "does not contain") << " 9" << '\n';
	cout << '\n';

	Set<int> B;

	/* inserting elements in Set B */
	B.insert(1); B.insert(2); B.insert(4);

    /* printing the contents of Set B */
	cout << "B = ";
	B.printSet();
	cout << "P(B) = " << '\n';
	B.printPowerSet();
	cout << '\n';

	Set<int> C;
	C.insert(1); C.insert(2); C.insert(4);

    /* printing the contents of Set C */
	cout << "C = ";
	C.printSet();
	cout << '\n';

    /* F = A - B */
	Set<int> F = A - B;
	cout << "A - B = ";
	F.printSet();
	cout << '\n';

    /* D = A Union B */
	Set<int> D = A._union(B);
	cout << "A union B = ";
	D.printSet();
	cout << '\n';

    /* E = A intersection B */
	Set<int> E = A._intersection(B);
	cout << "A intersection B = ";
	E.printSet();
	cout << '\n';

    /* printing the product */
	cout << "A x B = ";
	A.printProduct(B);
	cout << '\n';

	// Equality tests
	cout << "Equality of Sets:" << '\n';

	cout << "A " << ((A == B) ? "" : "!") << "= B" << '\n';
	cout << "B " << ((B == C) ? "" : "!") << "= C" << '\n';
	cout << "A " << ((A == C) ? "" : "!") << "= C" << '\n';
	cout << '\n';

	Set<int> U;
	U.insert(1); U.insert(2); U.insert(3); U.insert(4);
	U.insert(5); U.insert(6); U.insert(7);

	// Complements of the respective Sets
	Set<int> A1 = A._complement(U);
	Set<int> B1 = B._complement(U);
	Set<int> C1 = C._complement(U);

	cout << "A' = ";
	A1.printSet();
	cout << "B' = ";
	B1.printSet();
	cout << "C' = ";
	C1.printSet();

	return 0;
}
