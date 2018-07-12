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
		numtype** m_matrix;
	public:
		Matrix();
		Matrix(int rows, int cols, bool reset);
		Matrix(const Matrix<numtype> &mat);
		~Matrix();
		void clear();
		void Print() const;
		int getROWS() const {return ROWS;}
		int getCOLS() const {return COLS;}
		numtype** getEntrance() const {return m_matrix;}
		
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
		void setCol(int col, Matrix<numtype> &mat);
		Matrix<numtype> repeat(int M, int N) const;
		Matrix<numtype> logValue(double base) const;
		Matrix<numtype> reshape(int rows, int cols) const;
		void resize(int rows, int cols, bool reset);
		void normalize(bool direction);
		Matrix<numtype> cumsum() const;
		Matrix<numtype> concat(Matrix<numtype> &mat, int direction) const;
		
		numtype get(int row, int col) const;
		void set(int row, int col, numtype value);
		
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
Matrix<numtype>::Matrix(int rows,int cols, bool reset) : ROWS(rows), COLS(cols) {
	assert(rows >= 0);
	assert(cols >= 0);
	if(ROWS > 0) {
		m_matrix = new numtype*[ROWS];
		for(size_t i = 0; i < ROWS;  i++) {
			m_matrix[i] = new numtype[COLS];
			if(reset) {
				memset(m_matrix[i], 0, COLS*sizeof(numtype));
			}
		}
	}
	else {
		m_matrix = NULL;
	}
	//ZERO_FINAL = 2.2204e-16;
	
}

template <class numtype>
Matrix<numtype>::Matrix(const Matrix<numtype> &mat) : ROWS(mat.ROWS), COLS(mat.COLS){
	assert(mat.ROWS >= 0);
	assert(mat.COLS >= 0);
	if(ROWS > 0) {
		m_matrix = new numtype*[ROWS];
		for(size_t i = 0; i < ROWS;  i++) {
			m_matrix[i] = new numtype[COLS];
			for(size_t j = 0; j < COLS; j++) {
				m_matrix[i][j] = mat.m_matrix[i][j];
			}
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
	for(size_t i = 0; i < ROWS; i++) {
		delete[] m_matrix[i];
		m_matrix[i] = NULL;
	}
    	delete[] m_matrix;
	m_matrix = NULL;
	ROWS = COLS = 0;
	m_matrix = NULL;
}

template <class numtype>
void Matrix<numtype>::Print() const {
    for(size_t i = 0; i < ROWS;  i++) {
		for(size_t j = 0; j < COLS;  j++) {
			cerr << m_matrix[i][j] << '\t';
		}
		cerr << endl;
	}
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::transpose() const {
	Matrix<numtype> ret(COLS, ROWS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[j][i] = m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::inverse() const {
	assert(ROWS == COLS);
	numtype det = determinant();
	assert(det != 0);
	
	Matrix<numtype> ret(ROWS, COLS, false);
	
	size_t i, j, k, l, m, n;
	for(i = 0; i < ROWS; i++) {
		for(j = 0; j < COLS; j++) {
			Matrix<numtype> a(ROWS-1, COLS-1, false);
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
					a.m_matrix[k][l] = m_matrix[m][n];
				}
			}
			numtype temp = a.determinant();
			if((i+j)%2 == 0) {
				ret.m_matrix[j][i] = (numtype) temp/det;
			}
			else {
				ret.m_matrix[j][i] = (numtype) -temp/det;
			}
			if(fabs(ret.m_matrix[j][i]) < ZERO_FINAL) {
				ret.m_matrix[j][i] = 0;
			}
		}
	}
	
	return ret;
}

template <class numtype>
numtype Matrix<numtype>::determinant() const {
	assert(ROWS == COLS);
	Matrix<numtype> cur_mat(*this);
	size_t i, j, k, m;
	size_t flag = 0, switchcount = 0;
	numtype a, ret = 1;
	for(i = 0; i < ROWS-1; i++) {
		j = i+1;
		if(cur_mat.m_matrix[i][i] == 0) {
			while(cur_mat.m_matrix[j][i] == 0) {
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
			numtype *temp = cur_mat.m_matrix[i];
			cur_mat.m_matrix[i] = cur_mat.m_matrix[j];
			cur_mat.m_matrix[j] = temp;
		}
		for(j = i+1; j < ROWS; j++) {
			if(cur_mat.m_matrix[j][i] == 0) {
				continue;
			}
			a = cur_mat.m_matrix[j][i]/cur_mat.m_matrix[i][i];
			for(k = 0; k < COLS; k++) {
				cur_mat.m_matrix[j][k] -= a*cur_mat.m_matrix[i][k];
			}
		}
	}
	for(i = 0; i < ROWS; i++) {
		ret *= cur_mat.m_matrix[i][i];
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
	size_t i, j, k, m;
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
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret += m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::sumRows() const {
	Matrix<numtype> ret(1, COLS, true);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[0][j] += m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::sumCols() const {
	Matrix<numtype> ret(ROWS, 1, true);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][0] += m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Row(int row) const {
	assert(row >= 0 && row < ROWS);
	Matrix<numtype> ret(1, COLS, false);
	for(size_t j = 0; j < COLS; j++) {
		ret.m_matrix[0][j] = m_matrix[row][j];
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Rows(vector<int> index) const {
	Matrix<numtype> ret(index.size(), COLS, false);
	for(size_t i = 0; i < index.size(); i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[index[i]][j];
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Col(int col) const {
	assert(col >= 0 && col < COLS);
	Matrix<numtype> ret(ROWS, 1, false);
	for(size_t i = 0; i < ROWS; i++) {
		ret.m_matrix[i][0] = m_matrix[i][col];
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::Cols(vector<int> index) const {
	Matrix<numtype> ret(ROWS, index.size(), false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < index.size(); j++) {
			ret.m_matrix[i][j] = m_matrix[i][index[j]];
		}
	}
	return ret;
}

template <class numtype>
void Matrix<numtype>::setRow(int row, Matrix<numtype> &mat) {
	assert(row >= 0 && row < ROWS);
	for(size_t j = 0; j < COLS; j++) {
		m_matrix[row][j] = mat.m_matrix[0][j];
	}
}

template <class numtype>
void Matrix<numtype>::setCol(int col, Matrix<numtype> &mat) {
	assert(col >= 0 && col < COLS);
	for(size_t i = 0; i < ROWS; i++) {
		m_matrix[i][col] = mat.m_matrix[i][0];
	}
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::repeat(int M, int N) const {
	Matrix<numtype> ret(M*ROWS, N*COLS, false);
	for(size_t m = 0; m < M; m++) {
		for(size_t n = 0; n < N; n++) {
			for(size_t i = 0; i < ROWS; i++) {
				for(size_t j = 0; j < COLS; j++) {
					ret.m_matrix[i+m*ROWS][j+n*COLS] = m_matrix[i][j];
				}
			}
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::logValue(double base) const {
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = log(m_matrix[i][j]+ZERO_FINAL)/log(base);
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::reshape(int rows, int cols) const {
	assert(rows >= 0);
	assert(cols >= 0);
	if(rows == 0)
		rows = ROWS*COLS/cols;
	if(cols == 0)
		cols = ROWS*COLS/rows;
	
	assert(rows*cols == ROWS*COLS);
	Matrix<numtype> ret(rows, cols, false);
	size_t m = 0, n = 0;
	for(size_t j = 0; j < cols; j++) {
		for(size_t i = 0; i < rows; i++) {
			ret.m_matrix[i][j] = m_matrix[m][n];
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
	
	numtype** m_matrix_new = new numtype*[rows];
	for(size_t i = 0; i < rows;  i++) {
		m_matrix_new[i] = new numtype[cols];
		if(reset) {
			for(size_t j = 0; j < cols; j++) {
				if(i < ROWS && j < COLS) {
					m_matrix_new[i][j] = m_matrix[i][j];
				}
				else {
					m_matrix_new[i][j] = 0;
				}
			}
		}
	}
	
	for(size_t i = 0; i < ROWS; i++) {
		delete[] m_matrix[i];
	}
	delete[] m_matrix;
	
	ROWS = rows;
	COLS = cols;
	
	m_matrix = m_matrix_new;
}

template <class numtype>
void Matrix<numtype>::normalize(bool direction) {
	if(m_matrix == NULL) {
		return;
	}
	if(direction == 1) {
		Matrix<numtype> temp = sumRows();
		for(size_t i = 0; i < ROWS; i++) {
			for(size_t j = 0; j < COLS; j++) {
				m_matrix[i][j] /= (ZERO_FINAL + temp.m_matrix[0][j]);
			}
		}
	}
	else {
		Matrix<numtype> temp = sumCols();
		for(size_t i = 0; i < ROWS; i++) {
			for(size_t j = 0; j < COLS; j++) {
				m_matrix[i][j] /= (ZERO_FINAL + temp.m_matrix[i][0]);
			}
		}
	}
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::cumsum() const {
	if(m_matrix == NULL) {
		return *this;
	}
	Matrix<numtype> ret(ROWS, COLS, true);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			if(j > 0) {
				ret.m_matrix[i][j] = m_matrix[i][j] + ret.m_matrix[i][j-1];
			}
			else {
				ret.m_matrix[i][j] = m_matrix[i][j];
			}
		}
	}
	return ret;
}

template <class numtype>
Matrix<numtype> Matrix<numtype>::concat(Matrix<numtype> &mat, int direction) const {
	int rows, cols;
	size_t i, j;
	if(m_matrix == NULL) {
		Matrix<numtype> ret = mat;
		return ret;
	}
	if(direction == 1) {
		assert(COLS == mat.COLS);
		rows = ROWS+mat.ROWS;
		cols = COLS;
		Matrix<numtype> ret(rows, cols, false);
		for(j = 0; j < cols; j++) {
			for(i = 0; i < ROWS; i++) {
				ret.m_matrix[i][j] = m_matrix[i][j];
			}
			for(; i < rows; i++) {
				ret.m_matrix[i][j] = mat.m_matrix[i-ROWS][j];
			}
		}
		return ret;
	}
	else {
		assert(ROWS == mat.ROWS);
		rows = ROWS;
		cols = COLS+mat.COLS;
		Matrix<numtype> ret(rows, cols, false);
		for(i = 0; i < rows; i++) {
			for(j = 0; j < COLS; j++) {
				ret.m_matrix[i][j] = m_matrix[i][j];
			}
			for(; j < cols; j++) {
				ret.m_matrix[i][j] = mat.m_matrix[i][j-COLS];
			}
		}
		return ret;
	}
}

template <class numtype>
numtype Matrix<numtype>::get(int row, int col) const {
	assert(row >= 0 && row < ROWS);
	assert(col >= 0 && col < COLS);
	return m_matrix[row][col];
}

template <class numtype>
void Matrix<numtype>::set(int row, int col, numtype value) {
	assert(row >= 0 && row < ROWS);
	assert(col >= 0 && col < COLS);
	m_matrix[row][col] = value;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::dotProduct(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]*mat.m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::dotDivide(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]/(ZERO_FINAL+mat.m_matrix[i][j]);
		}
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator+(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]+mat.m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator+(numtype a) const {
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]+a;
		}
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator+=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			m_matrix[i][j] = m_matrix[i][j]+mat.m_matrix[i][j];
		}
	}
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator-(const Matrix<numtype> &mat) const {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]-mat.m_matrix[i][j];
		}
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator-(numtype a) const {
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]-a;
		}
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator-=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			m_matrix[i][j] = m_matrix[i][j]-mat.m_matrix[i][j];
		}
	}
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator*(const Matrix<numtype> &mat) const {
	assert(COLS == mat.ROWS);
	Matrix<numtype> ret(ROWS, mat.COLS, true);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < mat.COLS; j++) {
			for(size_t k = 0; k < COLS; k++) {
				ret.m_matrix[i][j] += m_matrix[i][k]*mat.m_matrix[k][j];
			}
			if(fabs(ret.m_matrix[i][j]) < ZERO_FINAL) {
				ret.m_matrix[i][j] = 0;
			}
		}
	}
	return ret;
}

template <class numtype>
inline Matrix<numtype> Matrix<numtype>::operator*(numtype a) const {
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]*a;
		}
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator*=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			m_matrix[i][j] = m_matrix[i][j]*mat.m_matrix[i][j];
		}
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
	if(a == 0) {
		a = ZERO_FINAL;
	}
	Matrix<numtype> ret(ROWS, COLS, false);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			ret.m_matrix[i][j] = m_matrix[i][j]/a;
		}
	}
	return ret;
}

template <class numtype>
inline void Matrix<numtype>::operator/=(Matrix<numtype> &mat) {
	assert(ROWS == mat.ROWS);
	assert(COLS == mat.COLS);
	for(size_t i = 0; i < ROWS; i++) {
		for(size_t j = 0; j < COLS; j++) {
			m_matrix[i][j] = m_matrix[i][j]/(ZERO_FINAL+mat.m_matrix[i][j]);
		}
	}
}

template <class numtype>
inline void Matrix<numtype>::operator=(const Matrix<numtype> &mat) {
	size_t i, j;
	for(i = 0; i < ROWS; i++) {
		delete[] m_matrix[i];
	}
	delete[] m_matrix;
	
	ROWS = mat.ROWS;
	COLS = mat.COLS;
	
	//m_matrix = mat.m_matrix;
	
	m_matrix = new numtype*[ROWS];
	for(i = 0; i < ROWS; i++) {
		m_matrix[i] = new numtype[COLS];
		for(j = 0; j < COLS; j++) {
			m_matrix[i][j] = mat.m_matrix[i][j];
		}
	}
	
}


#endif

