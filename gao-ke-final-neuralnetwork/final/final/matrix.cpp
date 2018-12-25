#include "matrix.h"
#include <sstream>
#include <iostream>
#include <random>


//constructor with no argument will create a matrix of 3 x 3 with all entries 0. 
Matrix::Matrix() :ncol(3), nrow(3) {
	data = new double[ncol*nrow];
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = 0;
	}
}
//constructor with one argument, say n, will create an n x n matrix with all entry values 0. 
Matrix::Matrix(int n) :nrow(n), ncol(n) {
	data = new double[ncol*nrow];
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = 0;
	}
};
//constructor with two arguments, say n and m, will create an n x m matrix with all entries 0. 
Matrix::Matrix(int n, int m) :nrow(n), ncol(m) {
	data = new double[ncol*nrow];
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = 0;
	}
};
//constructor with three arguments, say n, m, and v, will create an n x m matrix with all values v. 
Matrix::Matrix(int n, int m, double v) :nrow(n), ncol(m) {
	data = new double[ncol*nrow];
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = v;
	}
};
//constructor with n,m,seq will create an n x m matrix with sequence seq
Matrix::Matrix(int n, int m, double* seq) :nrow(n), ncol(m) {
	data = new double[ncol*nrow];
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = seq[i];
	}
};
Matrix::Matrix(int n, int m, vector<double> vec) :nrow(n), ncol(m) {
	data = new double[ncol*nrow];
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = vec[i];
	}
};
//destructor
Matrix::~Matrix() {
	if (data != nullptr) {
		delete[] data;
		data = nullptr;
	}
};
vector<double> Matrix::MatrixToVector() {
	vector<double> result;
	for (int i = 0; i < ncol*nrow; i++) {
		result.push_back(data[i]);
	}
	return result;
};
//getter
int Matrix::GetNRow() const {
	return nrow;
};
//getter
int Matrix::GetNCol() const {
	return ncol;
};
//getter
double Matrix::GetElement(int n, int m) const {
	return data[ncol*(n - 1) + m - 1];
};
//setter
void Matrix::SetElement(int n, int m, double v) {
	data[ncol*(n - 1) + m - 1] = v;
};
//get column
double* Matrix::GetColumn(int n) const {
	double* result;
	result = new double[nrow];
	for (int i = 0; i < nrow; i++) {
		result[i] = this->GetElement(i + 1, n);
	};
	return result;
};
//get row
double* Matrix::GetRow(int n) const {
	double* result;
	result = new double[ncol];
	for (int i = 0; i < ncol; i++) {
		result[i] = this->GetElement(n, i + 1);
	};
	return result;
};
//overload << operator
ostream& operator<<(ostream& os, const Matrix& inputmatrix) {
	for (int i = 0; i < (inputmatrix.ncol*inputmatrix.nrow); i++) {
		os << inputmatrix.data[i];
		os << "  \t";
		if ((i + 1) % inputmatrix.ncol == 0) {
			os << "\n";
		}
	}
	return os;
};
//overload >> operator
istream& operator>>(istream& is, Matrix& inputmatrix) {
	cout << "Please input nrow, ncol" << endl;
	is >> inputmatrix.nrow >> inputmatrix.ncol;
	return is;
};
Matrix Matrix::PointwiseProduct(const Matrix& inputmatrix) {
	if (ncol == inputmatrix.ncol&&nrow == inputmatrix.nrow) {
		Matrix resultmatrix(inputmatrix);
		for (int i = 0; i < ncol*nrow; i++) {
			resultmatrix.data[i] = inputmatrix.data[i] * data[i];
		}
		return resultmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
//overload assignment operator
Matrix& Matrix::operator= (const Matrix& inputmatrix) {
	if (this == &inputmatrix) {
		return *this;
	}
	nrow = inputmatrix.GetNRow();
	ncol = inputmatrix.GetNCol();
	if (data != nullptr) {
		delete[] data;
	}
	data = new double[nrow*ncol];
	for (int i = 0; i < nrow*ncol; i++) {
		data[i] = inputmatrix.data[i];
	}
	return *this;
};
//overload + operator
Matrix Matrix::operator+ (const Matrix& inputmatrix) {
	if (ncol == inputmatrix.ncol&&nrow == inputmatrix.nrow) {
		Matrix resultmatrix(inputmatrix);
		for (int i = 0; i < ncol*nrow; i++) {
			resultmatrix.data[i] = inputmatrix.data[i] + data[i];
		}
		return resultmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
//overload - operator
Matrix Matrix::operator- (const Matrix& inputmatrix) {
	if (ncol == inputmatrix.ncol&&nrow == inputmatrix.nrow) {
		Matrix resultmatrix(inputmatrix);
		for (int i = 0; i < ncol*nrow; i++) {
			resultmatrix.data[i] = -inputmatrix.data[i] + data[i];
		}
		return resultmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
//overload * operator
Matrix Matrix::operator* (const Matrix& inputmatrix) {
	if (ncol == inputmatrix.nrow) {
		Matrix newmatrix(nrow, inputmatrix.ncol);
		for (int i = 0; i < nrow; i++) {
			double* a;
			a = this->GetRow(i + 1);
			for (int j = 0; j < inputmatrix.ncol; j++) {
				double *b;
				b = inputmatrix.GetColumn(j + 1);
				double sum = 0;
				for (int k = 0; k < ncol; k++) {
					sum = sum + a[k] * b[k];
				}
				newmatrix.SetElement(i + 1, j + 1, sum);
			}
		}
		return newmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
//return number multiply
Matrix Matrix::operator* (const double multiplier) {
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] *= multiplier;
	}
	return *this;
};
//copy constructor
Matrix::Matrix(const Matrix& inputmatrix) {
	nrow = inputmatrix.GetNRow();
	ncol = inputmatrix.GetNCol();
	data = new double[nrow*ncol];
	for (int i = 0; i < nrow*ncol; i++) {
		data[i] = inputmatrix.data[i];
	}
};
//return identity matrix
Matrix Matrix::IdentityMatrix() {
	if (ncol == nrow) {
		Matrix newmatrix(nrow);
		for (int i = 0; i < nrow; i++) {
			newmatrix.SetElement(i + 1, i + 1, 1);
		}
		return newmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
//return transpose matrix
Matrix Matrix::TransposeMatrix() {
	Matrix newmatrix(ncol, nrow);
	for (int i = 0; i < ncol; i++) {
		for (int j = 0; j < nrow; j++) {
			newmatrix.data[j + i * nrow] = *(GetColumn(i + 1) + j);
		}
	}
	return newmatrix;
}
//return adjugate matrix
Matrix Matrix::AdjugateMatrix(int n, int m) {
	if (ncol == nrow) {
		Matrix newmatrix(nrow - 1, ncol - 1);
		int p = 0;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				if (n != i + 1 && m != j + 1) {
					newmatrix.data[p] = GetElement(i + 1, j + 1);
					p++;
				}
			}
		}
		return newmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
}
//return determinant 
double Matrix::DetMatrix() {
	if (nrow == ncol) {
		if (nrow == 1) {
			return GetElement(1, 1);
		}
		else {
			double det = 0;
			Matrix newmatrix;
			for (int i = 0; i < ncol; i++) {
				newmatrix = AdjugateMatrix(1, i + 1);
				det += pow(-1, i + 2)*GetElement(1, i + 1)*newmatrix.DetMatrix();
			}
			return det;
		}
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
//return inverse matrix
Matrix Matrix::InverseMatrix() {
	if (nrow == ncol) {
		double det = DetMatrix();
		if (det == 0) {
			cout << "Error: Det = 0" << endl;
			throw - 1;
		}
		Matrix resultmatrix(ncol);
		Matrix newmatrix;
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				newmatrix = AdjugateMatrix(i + 1, j + 1);
				resultmatrix.SetElement(j + 1, i + 1, pow(-1, i + j + 2)*newmatrix.DetMatrix() / det);// transpose matrix and change +/-
			}
		}
		return resultmatrix;
	}
	else {
		cout << "Matrix Size Error" << endl;
		throw - 1;
	}
};
Matrix Matrix::GetRandomMatrix(int lower_bound, int upper_bound) {
	random_device seed;
	default_random_engine generator(seed());
	uniform_real_distribution<double> distribution(lower_bound, upper_bound);
	for (int i = 0; i < ncol*nrow; i++) {
		data[i] = distribution(generator);
	}
	return *this;
};
Matrix Matrix::EraseLastRow() {
	Matrix newmatrix(GetNRow() - 1, GetNCol());
	for (int i = 0; i < nrow*ncol-ncol; i++) {
		newmatrix.data[i] = data[i];
	}
	return newmatrix;
};