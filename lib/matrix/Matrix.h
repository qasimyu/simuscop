// ***************************************************************************
// Matrix.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <cmath>

using namespace std;

template <class numtype>
class Matrix {
	private:
		int ROWS, COLS;
		static double ZERO_FINAL;
		numtype* m_matrix;
	public:
		Matrix();
		Matrix(int rows, int cols);
		Matrix(int rows, int cols, numtype value);
		Matrix(const Matrix<numtype> &mat);
		~Matrix();
		void clear();
		void Print() const;
		int getROWS() const {return ROWS;}
		int getCOLS() const {return COLS;}
		numtype* getEntrance() const {return m_matrix;}
		
		Matrix<numtype> transpose() const;
		Matrix<numtype> inverse() const;
		numtype determinant() const;
		numtype sum() const;
		Matrix<numtype> sumRows() const;
		Matrix<numtype> sumCols() const;
		Matrix<numtype> Row(int row) const;
		Matrix<numtype> Col(int col) const;
		Matrix<numtype> Rows(vector<int> index) const;
		Matrix<numtype> Cols(vector<int> index) const;
		void setRow(int row, Matrix<numtype> &mat);
		void setRow(int row, numtype value);
		void setCol(int col, Matrix<numtype> &mat);
		void setCol(int col, numtype value);
		void setCols(vector<int> colIndxs, Matrix<numtype> &mat);
		Matrix<numtype> repeat(int M, int N) const;
		Matrix<numtype> logValue(double base) const;
		Matrix<numtype> reshape(int rows, int cols) const;
		void resize(int rows, int cols, bool reset);
		void normalize(bool direction);
		Matrix<numtype> cumsum() const;
		Matrix<numtype> concat(Matrix<numtype> &mat, int direction) const;
		
		numtype get(int row, int col) const;
		void set(int row, int col, numtype value);
		void set(numtype value);
		
		void apply(numtype (*func)(numtype));
		
		inline Matrix<numtype> dotProduct(const Matrix<numtype> &mat) const;
		inline Matrix<numtype> dotDivide(const Matrix<numtype> &mat) const;
		
		inline Matrix<numtype> operator+(const Matrix<numtype> &mat) const;
		inline Matrix<numtype> operator+(numtype a) const;
		inline void operator+=(Matrix<numtype> &mat);
		inline Matrix<numtype> operator-(const Matrix<numtype> &mat) const;
		inline Matrix<numtype> operator-(numtype a) const;
		inline void operator-=(Matrix<numtype> &mat);
		inline Matrix<numtype> operator*(const Matrix<numtype> &mat) const;
		inline Matrix<numtype> operator*(numtype a) const;
		inline void operator*=(Matrix<numtype> &mat);
		inline Matrix<numtype> operator/(const Matrix<numtype> &mat) const;
		inline Matrix<numtype> operator/(numtype a) const;
		inline void operator/=(Matrix<numtype> &mat);
		
		inline void operator=(const Matrix<numtype> &mat);
};

template <class numtype>
double Matrix<numtype>::ZERO_FINAL = 2.2204e-16;

template <class numtype>
Matrix<numtype>::Matrix() {
	ROWS = 0;
	COLS = 0;
	m_matrix = NULL;
	//ZERO_FINAL = 2.2204e-16;
	//cerr << "default construction function used" << endl;
}

template <class numtype>
Matrix<numtype>::Matrix(int rows,int cols) : ROWS(rows), COLS(cols) {
	assert(rows >= 0);
	assert(cols >= 0);
	if(ROWS*COLS > 0) {
		m_matrix = new numtype[ROWS*COLS];
	}
	else {
		m_matrix = NULL;
	}
	
}

template <class numtype>
Matrix<numtype>::Matrix(int rows, int cols, numtype value) : ROWS(rows), COLS(cols) {
	assert(rows >= 0);
	assert(cols >= 0);
	if(ROWS*COLS > 0) {
		m_matrix = new numtype[ROWS*COLS];
		int n = ROWS*COLS;
		for(int i = 0; i < n; i++) {
			m_matrix[i] = value;
		}
	}
	else {
		m_matrix = NULL;
	}
	
}

template <class numtype>
Matrix<numtype>::Matrix(const Matrix<numtype> &mat) {
	ROWS = mat.ROWS;
	COLS = mat.COLS;
	if(ROWS*COLS > 0) {
		m_matrix = new numtype[ROWS*COLS];
		int n = ROWS*COLS;
		for(int i = 0; i < n; i++) {
			m_matrix[i] = mat.m_matrix[i];
		}
	}
	else {
		m_matrix = NULL;
	}
	//ZERO_FINAL = 2.2204e-16;
}

template <class numtype>
Matrix<numtype>::~Matrix() {
    clear();
}

template <class numtype>
void Matrix<numtype>::clear() {
	if(m_matrix == NULL) {
		return;	
	}
    delete[] m_matrix;
	m_matrix = NULL;
	ROWS = COLS = 0;
}

template <class numtype>
void Matrix<numtype>::Print() const {
	if(m_matrix == NULL) {
		return;
	}
	for(int i = 0; i < ROWS;  i++) {
		for(int j = 0; j < COLS;  j++) {
			cerr << m_matrix[i*COLS+j] << '\t';
		}
		cerr << endl;
	}
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::transpose() const {
	Matrix<numtype> ret(COLS, ROWS);
	for(int i = 0; i < COLS;  i++) {
		for(int j = 0; j < ROWS;  j++) {
			ret.m_matrix[i*ROWS+j] = m_matrix[j*COLS+i];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::inverse() const {
	assert(ROWS == COLS);
	numtype det = determinant();
	assert(det != 0);
	
	Matrix<numtype> ret(ROWS, COLS);
	
	int i, j, k, l, m, n;
	for(i = 0; i < ROWS; i++) {
		for(j = 0; j < COLS; j++) {
			Matrix<numtype> a(ROWS-1, COLS-1);
			for(m = -1, k = 0; k < ROWS-1; k++) {
				m++;
				if(m == i) {
					m++;
				}
				for(n = -1, l = 0; l < COLS-1; l++) {
					n++;
					if(n == j) {
						n++;
					}
					a.m_matrix[k*(COLS-1)+l] = m_matrix[m*COLS+n];
				}
			}
			numtype temp = a.determinant();
			if((i+j)%2 == 0) {
				ret.m_matrix[j*COLS+i] = (numtype) temp/det;
			}
			else {
				ret.m_matrix[j*COLS+i] = (numtype) -temp/det;
			}
			if(fabs(ret.m_matrix[j*COLS+i]) < ZERO_FINAL) {
				ret.m_matrix[j*COLS+i] = 0;
			}
		}
	}
	
	return ret;
}

template <class numtype>
numtype Matrix<numtype>::determinant() const {
	assert(ROWS == COLS);
	Matrix<numtype> cur_mat(*this);
	int i, j, k, m;
	int flag = 0, switchcount = 0;
	numtype a, ret = 1;
	for(i = 0; i < ROWS-1; i++) {
		j = i+1;
		if(cur_mat.m_matrix[i*COLS+i] == 0) {
			while(cur_mat.m_matrix[j*COLS+i] == 0) {
				j++;
				if(j == ROWS) {
					flag = 1;
					break;
				}
			}
			if(flag == 1) {
				continue;
			}
			switchcount++;
			for(k = 0; k < COLS; k++) {
				numtype temp = cur_mat.m_matrix[i*COLS+k];
				cur_mat.m_matrix[i*COLS+k] = cur_mat.m_matrix[j*COLS+k];
				cur_mat.m_matrix[j*COLS+k] = temp;
			}
		}
		for(j = i+1; j < ROWS; j++) {
			if(cur_mat.m_matrix[j*COLS+i] == 0) {
				continue;
			}
			a = cur_mat.m_matrix[j*COLS+i]/cur_mat.m_matrix[i*COLS+i];
			for(k = 0; k < COLS; k++) {
				cur_mat.m_matrix[j*COLS+k] -= a*cur_mat.m_matrix[i*COLS+k];
			}
		}
	}
	for(i = 0; i < ROWS; i++) {
		ret *= cur_mat.m_matrix[i*COLS+i];
	}
	if(switchcount%2) {
		ret = -ret;
	}
	return ret;
}

/*
template <class numtype>
numtype Matrix<numtype>::determinant() const {
	assert(ROWS == COLS);
	int i, j, k, m;
	numtype a, ret = 0;
	int *temp = new int[COLS];
	for(i = 0; i < COLS; i++) {
		temp[i] = i;
	}
	do {
		a = 1;
		for(i = 0; i < ROWS; i++) {
			a *= m_matrix[i][temp[i]];
		}
		m = 0;
		for(j = 0; j < COLS-1; j++) {
			for(k = j+1; k < COLS; k++) {
				if(temp[j] > temp[k]) {
					m++;
				}
			}
		}
		if(m%2 == 0) {
			ret += a;
		}
		else {
			ret -= a;
		}
	}while(next_permutation(temp,temp+COLS));
	
	delete[] temp;
	return ret;
}
*/

template <class numtype>
numtype Matrix<numtype>::sum() const {
	numtype ret = 0;
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret += m_matrix[i];
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::sumRows() const {
	Matrix<numtype> ret(1, COLS, 0);
	for(int i = 0; i < ROWS; i++) {
		for(int j = 0; j < COLS; j++) {
			ret.m_matrix[j] += m_matrix[i*COLS+j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::sumCols() const {
	Matrix<numtype> ret(ROWS, 1, 0);
	for(int i = 0; i < ROWS; i++) {
		for(int j = 0; j < COLS; j++) {
			ret.m_matrix[i] += m_matrix[i*COLS+j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Row(int row) const {
	assert(row >= 0 && row < ROWS);
	Matrix<numtype> ret(1, COLS);
	for(int j = 0; j < COLS; j++) {
		ret.m_matrix[j] = m_matrix[row*COLS+j];
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Rows(vector<int> index) const {
	Matrix<numtype> ret(index.size(), COLS);
	for(int i = 0; i < index.size(); i++) {
		assert(index[i] >= 0 && index[i] < ROWS);
		for(int j = 0; j < COLS; j++) {
			ret.m_matrix[i*COLS+j] = m_matrix[index[i]*COLS+j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Col(int col) const {
	assert(col >= 0 && col < COLS);
	Matrix<numtype> ret(ROWS, 1);
	for(int i = 0; i < ROWS; i++) {
		ret.m_matrix[i] = m_matrix[i*COLS+col];
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Cols(vector<int> index) const {
	Matrix<numtype> ret(ROWS, index.size(), false);
	for(int j = 0; j < index.size(); j++) {
		assert(index[j] >= 0 && index[j] < COLS);
		for(int i = 0; i < ROWS; i++) {
			ret.m_matrix[i*index.size()+j] = m_matrix[i*COLS+index[j]];
		}
	}
	return ret;
}

template <class numtype>
void Matrix<numtype>::setRow(int row, Matrix<numtype> &mat) {
	assert(row >= 0 && row < ROWS);
	assert(mat.ROWS == 1 && COLS == mat.COLS);
	for(int j = 0; j < COLS; j++) {
		m_matrix[row*COLS+j] = mat.m_matrix[j];
	}
}

template <class numtype>
void Matrix<numtype>::setRow(int row, numtype value) {
	assert(row >= 0 && row < ROWS);
	for(int j = 0; j < COLS; j++) {
		m_matrix[row*COLS+j] = value;
	}
}

template <class numtype>
void Matrix<numtype>::setCol(int col, Matrix<numtype> &mat) {
	assert(col >= 0 && col < COLS);
	assert(mat.COLS == 1 && ROWS == mat.ROWS);
	for(int i = 0; i < ROWS; i++) {
		m_matrix[i*COLS+col] = mat.m_matrix[i];
	}
}

template <class numtype>
void Matrix<numtype>::setCol(int col, numtype value) {
	assert(col >= 0 && col < COLS);
	for(int i = 0; i < ROWS; i++) {
		m_matrix[i*COLS+col] = value;
	}
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::repeat(int M, int N) const {
	Matrix<numtype> ret(M*ROWS, N*COLS);
	for(int m = 0; m < M; m++) {
		for(int n = 0; n < N; n++) {
			for(int i = 0; i < ROWS; i++) {
				for(int j = 0; j < COLS; j++) {
					ret.m_matrix[(i+m*ROWS)*N*COLS+j+n*COLS] = m_matrix[i*COLS+j];
				}
			}
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::reshape(int rows, int cols) const {
	assert(rows > 0 || cols > 0);
	if(rows <= 0)
		rows = ROWS*COLS/cols;
	if(cols <= 0)
		cols = ROWS*COLS/rows;
	
	assert(rows*cols == ROWS*COLS);
	Matrix<numtype> ret(rows, cols);
	int m = 0, n = 0;
	for(int j = 0; j < cols; j++) {
		for(int i = 0; i < rows; i++) {
			ret.m_matrix[i*cols+j] = m_matrix[m*COLS+n];
			m++;
			if(m == ROWS) {
				m = 0;
				n++;
			}
		}
	}
	return ret;
}

template <class numtype>
void Matrix<numtype>::resize(int rows, int cols, bool reset) {
	assert(rows >= 0 && cols >= 0);
	
	if(rows*cols == 0) {
		delete[] m_matrix;
		m_matrix = NULL;
		ROWS = COLS = 0;
		return;
	}
	
	if(ROWS != rows || COLS != cols) {
		numtype *m_matrix_n = new numtype[rows*cols];
		if(reset) {
			for(int i = 0; i < rows;  i++) {
				for(int j = 0; j < cols; j++) {
					m_matrix_n[i*cols+j] = (i < ROWS && j < COLS)? m_matrix[i*COLS+j]:0;
				}
			}
		}
		ROWS = rows;
		COLS = cols;
		delete[] m_matrix;
		m_matrix = m_matrix_n;
	}
}

template <class numtype>
void Matrix<numtype>::normalize(bool direction) {
	if(m_matrix == NULL) {
		return;
	}
	if(direction == 1) {
		Matrix<numtype> temp = sumRows();
		for(int i = 0; i < ROWS; i++) {
			for(int j = 0; j < COLS; j++) {
				m_matrix[i*COLS+j] /= (ZERO_FINAL + temp.m_matrix[j]);
			}
		}
	}
	else {
		Matrix<numtype> temp = sumCols();
		for(int i = 0; i < ROWS; i++) {
			for(int j = 0; j < COLS; j++) {
				m_matrix[i*COLS+j] /= (ZERO_FINAL + temp.m_matrix[i]);
			}
		}
	}
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::cumsum() const {
	if(m_matrix == NULL) {
		return *this;
	}
	Matrix<numtype> ret(ROWS, COLS, 0);
	for(int i = 0; i < ROWS; i++) {
		for(int j = 0; j < COLS; j++) {
			if(j > 0) {
				ret.m_matrix[i*COLS+j] = m_matrix[i*COLS+j] + ret.m_matrix[i*COLS+j-1];
			}
			else {
				ret.m_matrix[i*COLS+j] = m_matrix[i*COLS+j];
			}
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::concat(Matrix<numtype> &mat, int direction) const {
	int rows, cols;
	int i, j;
	if(m_matrix == NULL) {
		Matrix<numtype> ret = mat;
		return ret;
	}
	if(direction == 1) {
		assert(COLS == mat.COLS);
		rows = ROWS+mat.ROWS;
		cols = COLS;
		Matrix<numtype> ret(rows, cols);
		for(j = 0; j < cols; j++) {
			for(i = 0; i < ROWS; i++) {
				ret.m_matrix[i*cols+j] = m_matrix[i*COLS+j];
			}
			for(; i < rows; i++) {
				ret.m_matrix[i*cols+j] = mat.m_matrix[(i-ROWS)*mat.COLS+j];
			}
		}
		return ret;
	}
	else {
		assert(ROWS == mat.ROWS);
		rows = ROWS;
		cols = COLS+mat.COLS;
		Matrix<numtype> ret(rows, cols);
		for(i = 0; i < rows; i++) {
			for(j = 0; j < COLS; j++) {
				ret.m_matrix[i*cols+j] = m_matrix[i*COLS+j];
			}
			for(; j < cols; j++) {
				ret.m_matrix[i*cols+j] = mat.m_matrix[i*mat.COLS+j-COLS];
			}
		}
		return ret;
	}
}

template <class numtype>
numtype Matrix<numtype>::get(int row, int col) const {
	assert(row >= 0 && row < ROWS);
	assert(col >= 0 && col < COLS);
	return m_matrix[row*COLS+col];
}

template <class numtype>
void Matrix<numtype>::set(int row, int col, numtype value) {
	assert(row >= 0 && row < ROWS);
	assert(col >= 0 && col < COLS);
	m_matrix[row*COLS+col] = value;
}

template <class numtype>
void Matrix<numtype>::set(numtype value) {
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		m_matrix[i] = value;
	}
}

template <class numtype>
void Matrix<numtype>::apply(numtype (*func)(numtype)) {
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		m_matrix[i] = (*func)(m_matrix[i]);
	}
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::dotProduct(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]*mat.m_matrix[i];
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::dotDivide(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]/(ZERO_FINAL+mat.m_matrix[i]);
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator+(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]+mat.m_matrix[i];
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator+(numtype a) const {
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]+a;
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator+=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		m_matrix[i] += mat.m_matrix[i];
	}
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator-(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]-mat.m_matrix[i];
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator-(numtype a) const {
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]-a;
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator-=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		m_matrix[i] -= mat.m_matrix[i];
	}
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator*(const Matrix<numtype> &mat) const {
	assert(COLS == mat.ROWS);
	Matrix<numtype> ret(ROWS, mat.COLS, 0);
	for(int i = 0; i < ROWS; i++) {
		for(int j = 0; j < mat.COLS; j++) {
			for(int k = 0; k < COLS; k++) {
				ret.m_matrix[i*mat.COLS+j] += m_matrix[i*COLS+k]*mat.m_matrix[k*mat.COLS+j];
			}
			if(fabs(ret.m_matrix[i*mat.COLS+j]) < ZERO_FINAL) {
				ret.m_matrix[i*mat.COLS+j] = 0;
			}
		}
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator*(numtype a) const {
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]*a;
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator*=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		m_matrix[i] *= mat.m_matrix[i];
	}
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator/(const Matrix<numtype> &mat) const {
	Matrix<numtype> temp = mat.inverse();
	assert(COLS == temp.ROWS);
	Matrix<numtype> cur_mat(*this);
	return cur_mat*temp;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator/(numtype a) const {
	assert(a != 0);
	Matrix<numtype> ret(ROWS, COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		ret.m_matrix[i] = m_matrix[i]/a;
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator/=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	int n = ROWS*COLS;
	for(int i = 0; i < n; i++) {
		m_matrix[i] /= (ZERO_FINAL+mat.m_matrix[i]);
	}
}

template <class numtype>
inline void Matrix<numtype>::operator=(const Matrix<numtype> &mat) {
	int i, j;
	if(ROWS != mat.ROWS || COLS != mat.COLS) {
		delete[] m_matrix;
		ROWS = mat.ROWS;
		COLS = mat.COLS;
		m_matrix = new numtype[ROWS*COLS];
	}
	int n = ROWS*COLS;
	for(i = 0; i < n; i++) {
		m_matrix[i] = mat.m_matrix[i];
	}
	
}


#endif

