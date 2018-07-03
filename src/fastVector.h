// fastVector classes
// Author: Dirk Stratmann
// email: dirk.stratmann at upmc.fr


/*
 This file is part of TEF.

    TEF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TEF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TEF.  If not, see <http://www.gnu.org/licenses/>.


(C) 2018 Stratmann D, Pathmanathan JS, Postic G, Rey J, Chomilier J
Contact: dirk.stratmann@sorbonne-universite.fr
*/

#ifndef fastVectorH
#define fastVectorH

#include <vector>
#include <deque>
#include <set>
#include <stdint.h>
#include "MyError.h"

template<typename T> class fastVector;
template<typename T> class fastDeque;

template<typename T> ostream &operator<<( ostream &out, const fastVector<T> &v );
template<typename T> ostream &operator<<( ostream &out, const fastDeque<T> &L );


template<typename T> class fastVector
{
      friend ostream &operator<< <T>( ostream &out, const fastVector<T> &L );
public:
	fastVector() : data(NULL), dataSize(0){
	}
	fastVector(unsigned size_) : data(NULL){
		allocateMem(size_);
	}
	fastVector(unsigned size_, const T& initialObj) : data(NULL){
		allocateMem(size_);
		initialize(initialObj);
	}

	fastVector(const fastVector<T>& obj) : data(NULL){
		allocateMem(obj.dataSize);
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj.data[i];
		}
	}

	fastVector(const vector<T>& obj) : data(NULL){
		allocateMem(obj.size());
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj[i];
		}
	}

	fastVector(const fastDeque<T>& obj) : data(NULL){
		allocateMem(obj.size());
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj[i];
		}
	}

	fastVector<T>& operator=(const fastDeque <T>& obj){
		allocateMem(obj.size());
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj[i];
		}
		return *this;
	}

	fastVector<T>& operator=(const fastVector <T>& obj){
		allocateMem(obj.dataSize);
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj.data[i];
		}
		return *this;
	}

	fastVector<T>& operator=(const vector<T>& obj){
		allocateMem(obj.size());
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj[i];
		}
		return *this;
	}

	fastVector<T>& operator=(const deque<T>& obj){
		allocateMem(obj.size());
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj[i];
		}
		return *this;
	}

	fastVector<T>& operator=(const set<T>& obj){
//		typename set<T>::const_iterator IterType;
		allocateMem(obj.size());
		unsigned ct = 0;
		for (typename set<T>::const_iterator iter = obj.begin(); iter != obj.end(); iter++) {
			data[ct] = *iter;
			ct++;
		}
		return *this;
	}


	void assign(const fastVector<T>& obj){
		if (dataSize != obj.dataSize) {
			throw MyError("ERROR in fastVector::assign");
		}
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = obj.data[i];
		}
	}

	bool operator==(const fastVector<T>& right) const{
		if (dataSize != right.dataSize) {
			return false;
		}
		for (unsigned i=0; i < dataSize; i++) {
			if (!(data[i] == right.data[i])) {
				return false;
			}
		}
		return true;
	}

	bool equalWithoutOrder(const fastVector<T>& right) const{
		if (dataSize != right.dataSize) {
			return false;
		}
		set<T> leftSet;
		set<T> rightSet;
		generateSet(leftSet);
		right.generateSet(rightSet);
		return (leftSet == rightSet);
	}


	void generateSet(set<T>& setObj) const{
		for (unsigned i=0; i < dataSize; i++) {
			setObj.insert(data[i]);
		}
	}

	void getVector(vector<T>& vec){
		vec.clear();
		vec.resize(dataSize);
		for (int i=0; i < dataSize; i++) {
			vec[i] = data[i];
		}
	}

	void resize(unsigned size_){
		allocateMem(size_);
	}

	void resize(unsigned size_, const T& initialObj){
		allocateMem(size_);
		initialize(initialObj);
	}

	~fastVector(){
		freeMem();
	}

	void freeMem(){
		if (data != NULL) {
			delete [] data;
			data = NULL;
		}
		dataSize = 0;
	}

	void allocateMem(unsigned size_){
		freeMem();
		dataSize = size_;
		data = new T[dataSize];
	}

	void initialize(const T& initialObj){
		for (unsigned i=0; i < dataSize; i++) {
			data[i] = initialObj;
		}
	}


	T& operator[](unsigned index){
		return data[index];
	}

	const T& operator[](unsigned index) const{
		return data[index];
	}

	T* getPointer(unsigned index){
        return &(data[index]);
	}

	unsigned size() const{
		return dataSize;
	}
	bool empty() const{
		if (dataSize == 0) {
			return true;
		} else{
			return false;
		}
	}

	int find(const T& obj) const{
		for (int i=0; i < dataSize; i++) {
			if (data[i] == obj) {
				return i;
			}
		}
		return -1;
	}

protected:
	T* data;
	unsigned dataSize;
};


template<typename T> ostream &operator<<( ostream &out, const fastVector<T> &v ) {
  out << "[";
  for (unsigned k=0; k< v.size(); k++) {
    out << v[k] << ' ';
  }
  out << "]";
  return out;
}


class bitVector : public fastVector<uint32_t>
{
public:
	bitVector() : fastVector<uint32_t>() { }
	bitVector(unsigned bitsize) : fastVector<uint32_t>(){
		unsigned size = bitsize / 32;
		if (bitsize % 32 > 0) {
			size++;
		}
		allocateMem(size);
		initialize(0);
	}
	bitVector(unsigned bitsize, const uint32_t& initialObj) : fastVector<uint32_t>(){
		unsigned size = bitsize / 32;
		if (bitsize % 32 > 0) {
			size++;
		}
		allocateMem(size);
		initialize(initialObj);
	}

	void setBit(unsigned bitNo){ // set the bit with bitNo to 1.
//		unsigned word = bitNo >> 5; // bitNo / 32
//		unsigned pos = bitNo & 31; // bitNo % 32
		data[bitNo >> 5] |= (1 << (bitNo & 31));
	}
	void clearBit(unsigned bitNo){ // set the bit with bitNo to 0.
//		unsigned word = bitNo >> 5; // bitNo / 32
//		unsigned pos = bitNo & 31; // bitNo % 32
		data[bitNo >> 5] &= ~(1 << (bitNo & 31));
	}
	unsigned getBit(unsigned bitNo){ // get the bit at bitNo
		return (data[bitNo >> 5] & (1 << (bitNo & 31)));
	}
	unsigned getBitSum(){
		unsigned bitSum = 0;
		unsigned sizeOfVec = size();
		unsigned bitIndex = 1;
		for (unsigned i=0; i < sizeOfVec; i++) {
			unsigned vecElement = data[i];
			for (unsigned j=0; j < 32; j++) {
				if(vecElement & bitIndex){
					bitSum++;
				}
				bitIndex = bitIndex << 1;
			}
		}
		return bitSum;
	}

	int ANDcomparison(const bitVector& second) const{
//		if (second.dataSize != dataSize) {
//			return -1;
//		}
		uint32_t value;
		unsigned bitIndex = 1;
		unsigned bitSum = 0;
		int bitID = -1;
		int resID = -1;
		for (unsigned i=0; i < dataSize; i++) {
			value = data[i] & second.data[i];
			if (value != 0){
				if (resID == -1) {
					bitIndex = 1;
					bitSum = 0;
					bitID = -1;
					for (unsigned j=0; j < 32; j++) {
						if(value & bitIndex){
							bitSum++;
							if (bitSum > 1) {
								bitID = -1;
								return -3; // more than one possibility
							} else{
								bitID = j;
							}
						}
						bitIndex = bitIndex << 1;
					}
					if (bitID != -1) {
						resID = (i << 5) + bitID; // i << 5 == i * 32
					}
				} else{
					return -3; // more than one possibility
				}
			}
		}
		if (resID == -1) {
			return -2; // no possibility found
		}
		return resID; // one single possiblity found
	}
private:

};


#define fastDeque_SIZE 5  // initial Size of a fastDeque .

template<typename T> class fastDeque
{
      friend ostream &operator<< <T>( ostream &out, const fastDeque<T> &L );
public:
	fastDeque() : numIncrement(fastDeque_SIZE), fastVec(fastDeque_SIZE), numElements(0) {
    	fastVecSize = fastVec.size();
	}
	fastDeque(unsigned numIncr) : numIncrement(numIncr), fastVec(numIncr), numElements(0) {
    	fastVecSize = fastVec.size();
	}
	fastDeque(unsigned numIncr, const T& initialObj) : numIncrement(numIncr), fastVec(numIncr, initialObj), numElements(0) {
    	fastVecSize = fastVec.size();
	}

	fastDeque<T>& operator=(const fastVector <T>& obj){
		fastVec = obj;
		numElements = fastVec.size();
		return *this;
	}

	fastDeque<T>& operator=(const fastDeque <T>& obj){
		fastVec = obj.fastVec;
		fastVecSize = fastVec.size();
		numElements = obj.numElements;
		numIncrement = obj.numIncrement;
		return *this;
	}

	const T& operator[] (unsigned offset) const {
		if (offset < numElements) {
			return fastVec[offset];
		} else{
			throw MyError("ERROR in fastDeque::operator[]");
		}
	}

	T& operator[] (unsigned offset) {
		if (offset < numElements) {
			return fastVec[offset];
		} else{
			throw MyError("ERROR in fastDeque::operator[]");
		}
	}

	T* getPointer(unsigned index){
        return fastVec.getPointer(index);
	}

	const T& back() const{
		return fastVec[numElements - 1];
	}
	T& back(){
		return fastVec[numElements - 1];
	}

	bool operator==(const fastDeque<T>& right) const{
		if (numElements != right.numElements) {
			return false;
		}
		for (unsigned i=0; i < numElements; i++) {
			if (!(fastVec[i] == right.fastVec[i])) {
				return false;
			}
		}
		return true;
	}

	void getVector(vector<T>& vec){
		vec.clear();
		vec.resize(numElements);
		for (int i=0; i < numElements; i++) {
			vec[i] = fastVec[i];
		}
	}


	unsigned size() const{
		return numElements;
	}
	bool empty() const{
		if (numElements == 0) {
			return true;
		} else{
			return false;
		}
	}

	void push_back(const T& obj){
		if (numElements < fastVecSize) {
			fastVec[numElements++] = obj;
		} else{ // need a larger array :
			fastVector<T> dummyVec = fastVec;
			fastVec.resize(fastVec.size() + numIncrement);
			for (unsigned i=0; i < numElements; i++) {
				fastVec[i] = dummyVec[i];
			}
			fastVec[numElements++] = obj;
			fastVecSize = fastVec.size();
		}
//		numElements++;
	}


	void insert(unsigned pos, const T& obj){
		if (pos > numElements) {
			throw MyError("ERROR in fastDeque::insert: pos > numElements");
		} else if (pos == numElements) {
			push_back(obj);
		}
		if (numElements < fastVecSize) {
			for (int i=numElements; i > pos; i--) {
				fastVec[i] = fastVec[i-1];
			}
			fastVec[pos] = obj;
		} else{
			fastVector<T> dummyVec = fastVec;
			fastVec.resize(fastVec.size() + numIncrement);
			fastVecSize = fastVec.size();
			for (int i=numElements; i > pos; i--) {
				fastVec[i] = dummyVec[i-1];
			}
			fastVec[pos] = obj;
			if (pos > 0) {
				for (int i=pos-1; i >= 0; i--) {
					fastVec[i] = dummyVec[i];
				}
			}
		}
		numElements++;
	}


	T pop_back(){
		if (numElements == 0) {
			throw MyError("ERROR in fastDeque::pop_back: numElements == 0");
		}
		numElements--;
		return fastVec[numElements];
	}

	void cutElement(unsigned offset){
		for (unsigned i=offset; i < numElements - 1; i++) {
            fastVec[i] = fastVec[i+1];
		}
		numElements--;
	}

	void clear(){ // the objects are not destroyed by clear() .
		numElements = 0;
	}

	int find(const T& obj) const{
		for (int i=0; i < numElements; i++) {
			if (fastVec[i] == obj) {
				return i;
			}
		}
		return -1;
	}

protected:
	unsigned numIncrement;
	fastVector<T> fastVec;
	unsigned numElements;
	unsigned fastVecSize;
};

template<typename T> ostream &operator<<( ostream &out, const fastDeque<T> &d ) {
  out << "[";
  for (unsigned k=0; k< d.size(); k++) {
    out << d[k] << ' ';
  }
  out << "]";
  return out;
}


template<typename T> class fastDoubleList
{
public:
	struct node
	  { T item; node* previous; node* next;
		node(const T& x, node* prev_, node* next_)
		  { item = x; previous = prev_, next = next_; }
	  };
	typedef node *link;

	fastDoubleList(){
		firstLink = NULL;
		lastLink = NULL;
		numElements = 0;
	}

	~fastDoubleList(){
		while (numElements > 0){
			deleteNode(lastLink);
		}
	}
/*	fastList(unsigned initialSize){
		construct(initialSize);
	}*/

	void insert(link where, const T& item){ // add before the position where
		link newEntry = new node(item, where->previous, where); // allocate Memory for new Node
		if (newEntry->previous != NULL) {
			newEntry->previous->next = newEntry;
		}
		where->previous = newEntry;
		numElements++;
	}

	void push_back(const T& item){
		link newEntry = new node(item, lastLink, NULL);
		if (lastLink != NULL) {
			lastLink->next = newEntry;
		}
		lastLink = newEntry;
		numElements++;
	}

	void deleteNode(link x){
		if (x != NULL) {
			if (x->previous != NULL) {
				x->previous->next = x->next;
			} else{
                firstLink = x->next;
            }
			if (x->next != NULL) {
				x->next->previous = x->previous;
			} else{
				lastLink = x->previous;
            }
			delete x;
			numElements--;
		}
	}

	link begin()
		{ return firstLink; }
	link end()
		{ return lastLink; }
	T& back()
		{ return lastLink->item; }
	const T& back() const
		{ return lastLink->item; }

	unsigned size() const{
		return numElements;
    }

	T& item(link x)
	  { return x->item; }
	link next(link x)
	  { return x->next; }

protected:

/*	void construct(unsigned N)
	  {
		freelist = new node[N+1];
		for (int i = 0; i < N; i++)
		  freelist[i].next = &freelist[i+1];
		freelist[N].next = 0;
	  }
	link newNode(int i)
	  { link x = remove(freelist);
		x->item = i; x->next = x;
		return x;
	  }
	void deleteNode(link x)
	  { insert(freelist, x); }
	void insert(link x, link t)
	  { t->next = x->next; x->next = t; }
	link remove(link x)
	  { link t = x->next; x->next = t->next; return t; }*/

	link firstLink;
	link lastLink;
	unsigned numElements;
};


// singly-linked list:
template<typename T> class fastList
{
public:
	struct node
	  { T item; node* next;
		node(const T& x, node* next_)
		  { item = x; next = next_; }
	  };
	typedef node *link;

	fastList(){
		firstLink = NULL;
	}

	~fastList(){
		while (pop_front()){
			// ...
		}
	}
/*	fastList(unsigned initialSize){
		construct(initialSize);
	}*/

	void insert(link where, const T& item){ // add after the position where
		link newEntry = new node(item, where->next); // allocate Memory for new Node
		where->next = newEntry;
	}

	bool deleteEntry(link afterWhich){ // delete the entry after the link to afterWhich
									// set afterWhich to NULL, if you want to delete the first entry
		if (afterWhich == NULL) {
            return pop_front();
		}
		link deletePnt = afterWhich->next;
		afterWhich->next = deletePnt->next;
		delete deletePnt;
		return true;
	}

	void push_front(const T& item){
		link newEntry = new node(item, firstLink);
		firstLink = newEntry;
	}

	bool pop_front(){
		if (firstLink != NULL) {
        	link nextLink = firstLink->next;
			delete firstLink;
			firstLink = nextLink;
			return true;
		} else{
            return false;
        }
	}

	link begin()
		{ return firstLink; }

	T& firstItem(){
		if (firstLink == NULL) {
			throw MyError("ERROR in fastList::firstItem");
		}
		return firstLink->item;
	}
	const T& firstItem() const{
		if (firstLink == NULL) {
			throw MyError("ERROR in fastList::firstItem");
		}
		return firstLink->item;
	}

	bool empty() const{
		if (firstLink == NULL) {
            return true;
		} else{
            return false;
		}
	}

	unsigned size() const{
		unsigned numElements = 0;
		if (firstLink != NULL) {
        	numElements++;
			link nextLink = firstLink->next;
			while (nextLink != NULL){
				numElements++;
				nextLink = nextLink->next;
			}
		}
		return numElements;
	}

	T& item(link x)
	  { return x->item; }
	link next(link x)
	  { return x->next; }

protected:

/*	void construct(unsigned N)
	  {
		freelist = new node[N+1];
		for (int i = 0; i < N; i++)
		  freelist[i].next = &freelist[i+1];
		freelist[N].next = 0;
	  }
	link newNode(int i)
	  { link x = remove(freelist);
		x->item = i; x->next = x;
		return x;
	  }
	void deleteNode(link x)
	  { insert(freelist, x); }
	void insert(link x, link t)
	  { t->next = x->next; x->next = t; }
	link remove(link x)
	  { link t = x->next; x->next = t->next; return t; }*/

	link firstLink;
};



template<typename T> class linesMatrix
{
public:
	linesMatrix() : data(NULL) {
		dataSize = 0;
		lineNum = 0;
	}
	linesMatrix(const linesMatrix<T>& obj) : data(NULL){
    	(*this) = obj;
	}

	linesMatrix(const vector<unsigned>& lineLengths_, const T& initialObj)
	: lineLengths(lineLengths_), data(NULL)
	{
		allocateMem();
		initialize(initialObj);
	}

	~linesMatrix(){
		freeMem();
	}

	linesMatrix<T>& operator =(const linesMatrix<T>& obj){
		lineLengths = obj.lineLengths;
		allocateMem();
		for (unsigned i=0; i < lineNum; i++) {
			for (unsigned j=0; j < lineLengths[i]; j++) {
				data[i][j] = obj.data[i][j];
			}
		}
		return *this;
	}

	void freeMem(){
		if (data != NULL) {
			for (unsigned i=0; i < lineNum; i++) {
				if (data[i] != NULL) {
//					try{
						delete [] data[i];
//					} catch(...){
//						cout << "ERROR in linesMatrix::freeMem()" << endl;
//						throw;
//					}
					data[i] = NULL;
				}
			}
			delete [] data;
			data = NULL;
		}
		dataSize = 0;
		lineNum = 0;
	}

	void allocateMem(){
		freeMem();
		data = new T*[lineLengths.size()];
		for (int i=0; i < lineLengths.size(); i++) {
			data[i] = new T[lineLengths[i]];
			dataSize += lineLengths[i];
			lineNum++;
		}
	}

	void initialize(const T& initialObj){
		for (unsigned i=0; i < lineNum; i++) {
			unsigned lineLen = lineLengths[i];
			T* linePnt = data[i];
			for (unsigned j=0; j < lineLen; j++) {
				linePnt[j] = initialObj;
			}
		}
	}

	void resize(const vector<unsigned>& lineLengths_){
		lineLengths = lineLengths_;
		allocateMem();
	}

	void resize(const vector<unsigned>& lineLengths_, const T& initialObj){
		lineLengths = lineLengths_;
		allocateMem();
		initialize(initialObj);
	}

	void resize(unsigned rows, unsigned cols, const T& initialObj){  // creates a classical Matrix with equal line lengths .
		lineLengths = vector<unsigned>(rows, cols);
		allocateMem();
		initialize(initialObj);
	}

	const T& operator()(unsigned row, unsigned col) const
	{
		#ifdef DEBUG_LINESMATRIX
		verifyIndex(row, col);
		#endif
		return data[row][col];
	}
	T& operator()(unsigned row, unsigned col)
	{
		#ifdef DEBUG_LINESMATRIX
		verifyIndex(row, col);
		#endif
		return data[row][col];
	}

//	T* operator[](unsigned row)   // could lead to confusion if a xxx[a][b]  is not replaced
//	{									// and T is a pointer to something
//    	return getLinePnt(row);
//	}

	T* getLinePnt(unsigned row){
		#ifdef DEBUG_LINESMATRIX
		verifyIndex(row, 0);
		#endif
		return data[row];
	}

	bool setRow(unsigned row, const fastVector<T>& lineVec){
		if (row >= lineNum) {
			return false;
		}
		if (lineVec.size() != lineLengths[row]) {
			return false;
		}
		for (unsigned i=0; i < lineVec.size(); i++) {
        	data[row][i] = lineVec[i];
		}
	}

	bool verifyIndex(unsigned row, unsigned col) const{
		if (row >= lineNum) {
			printf("ERROR"); // to be replaces by a throw
			return false;
		} else if(col >= lineLengths[row]){
			printf("ERROR");
			return false;
		}
		return true;
	}


	unsigned size() const{
		return lineNum;
	}
	int size(unsigned row) const{
		#ifdef DEBUG_LINESMATRIX
		verifyIndex(row, 0);
		#endif
		return lineLengths[row];
	}

protected:
	T** data;
	unsigned dataSize;
	unsigned lineNum;
	fastVector<unsigned> lineLengths;
};


// the binaryMatrix must only have 0 or 1 entries
/*void BinaryMatrix2LinesMatrix(const matrix<int>& binaryMatrix, linesMatrix<unsigned>& linesMat)
{
	vector<unsigned> lineLengths(binaryMatrix.columns());
	for (unsigned i=0; i < binaryMatrix.columns(); i++) {
		lineLengths[i] = binaryMatrix.sum_column(i);
	}
	linesMatrix<unsigned> dummy(lineLengths, 0);
	linesMat = dummy;

	for (unsigned i=0; i < binaryMatrix.columns(); i++) {
		unsigned lineCt = 0;
		for (unsigned j=0; j < binaryMatrix.rows(); j++) {
			if (binaryMatrix(j, i) == 1) {
				if (lineCt >= lineLengths[i]) {
					throw MyError("ERROR in BinaryMatrix2LinesMatrix()");
				}
				linesMat(i, lineCt) = j;
				lineCt++;
			} else if (binaryMatrix(j, i) != 0) {
				throw MyError("ERROR in BinaryMatrix2LinesMatrix()");
			}
		}
		if (lineCt != lineLengths[i]) {
           	throw MyError("ERROR in BinaryMatrix2LinesMatrix()");
		}
	}
}*/

#endif
