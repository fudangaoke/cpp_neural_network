#pragma once
#include <string>
#include <vector>
using namespace std;
class Matrix {
private:
	int ncol;
	int nrow;
	double* data;
public:
	Matrix();
	Matrix(int n);
	Matrix(int n, int m);
	Matrix(int n, int m, double v);
	Matrix(int n, int m, double* seq);
	Matrix(int n, int m, vector<double> vec);
	~Matrix();
	vector<double> MatrixToVector();
	int GetNRow() const;
	int GetNCol() const;
	double GetElement(int n, int m) const;
	void SetElement(int n, int m, double v);
	double* GetColumn(int n) const;
	double* GetRow(int m)const;
	Matrix& operator= (const Matrix& inputmatrix);
	Matrix operator+ (const Matrix& inputmatrix);
	Matrix operator- (const Matrix& inputmatrix);
	Matrix operator* (const Matrix& inputmatrix);
	Matrix operator* (const double multiplier);
	Matrix(const Matrix& inputmatrix);
	Matrix PointwiseProduct(const Matrix& inputmatrix);
	Matrix IdentityMatrix();
	Matrix TransposeMatrix();
	double DetMatrix();
	Matrix AdjugateMatrix(int n, int m);
	Matrix InverseMatrix();
	friend ostream& operator<<(ostream& os, const Matrix& inputmatrix);
	friend istream& operator>>(istream& is, Matrix& inputmatrix);
	Matrix GetRandomMatrix(int lower_bound, int upper_bound);
	Matrix EraseLastRow();
};