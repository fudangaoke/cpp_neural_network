#include <iostream>
#include <random>
#include "matrix.h"
#include "neuron.h"
#include "layer.h"
#include "neuralnetwork.h"
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

//read CSV file and return data.
vector<vector<double>> ReadCSV(string input_filename) {
	ifstream input_file(input_filename, ios::in);
	string stock_string_line;
	vector<vector<double>> vector_of_doubles;
	while (getline(input_file, stock_string_line)) {
		stringstream ss(stock_string_line);
		string stock_string_unit;
		double stock_double;
		vector<double> vector_double_line;
		while (getline(ss, stock_string_unit, ',')) {
			stringstream ss1(stock_string_unit);
			ss1 >> stock_double;
			vector_double_line.push_back(stock_double);
		}
		vector_of_doubles.push_back(vector_double_line);
	}
	return vector_of_doubles;
}


int main() {

	//set topology of neural network
	vector<int> topology;	
	topology.push_back(4);			//size of input layer = 4
	topology.push_back(4);			//size of hidden layer I = 4
	topology.push_back(3);			//size of hidden layer II = 3
	topology.push_back(2);			//size of hidden layer III = 2
	topology.push_back(1);			//size of output layer = 1

	vector<vector<double>> vector_of_inputs; 
	vector<vector<double>> vector_of_targets;
	//read training file from csv
	vector_of_inputs = ReadCSV("training_input.csv");
	vector_of_targets = ReadCSV("training_target.csv");

	//set parameters for neural network
	int cycle = 0;
	int max_cycle = 1000;
	double error = 1;					//error = 0.5*(output - target)^2
	double derived_error = 1;	//derived error = output - target
	double average_error = 1;	//average of error
	double tolerance_error = 0.0000001;	
	double learning_rate = 1;	//learning rate
	
	//before training
	NeuralNetwork nn(topology, learning_rate);			
	cout << "Before training" << endl<<endl;
	cout << nn << endl;

	//training 
	while (cycle< max_cycle && average_error > tolerance_error) {
		average_error = 0;
		for (int i = 0; i < vector_of_inputs.size(); i++) {
			nn.SetInputValue(vector_of_inputs[i]);
			nn.SetTargetValue(vector_of_targets[i]);
			nn.ForwardPropagation();
			nn.BackPropagation();
			error = nn.GetTotalError();
			average_error += error / vector_of_inputs.size();
		}
		cycle += 1;
		cout << "cycle = " << cycle << "  average error = " << average_error << endl;
	}

	//after training
	cout << endl<<endl;
	cout << "After trainining" << endl << endl;
	nn.ForwardPropagation();
	cout << nn<< endl;

	//output neural network and weight matrices
	ofstream output_file;
	output_file.open("weight_matrices_output.txt", ios::trunc);
	output_file << nn;
	output_file.close();

	//prepare to output prediction results.
	output_file.open("prediction_output.txt", ios::trunc);
	
	//testing
	vector_of_inputs = ReadCSV("test_input.csv");
	vector_of_targets = ReadCSV("test_target.csv");
	cout << "Predict=    \tTarget=    \tPredict - Target =   "<< endl;
	output_file << "Predict=    \tTarget=    \tPredict - Target =   " << endl;
	for (int i = 0; i < vector_of_inputs.size(); i++) {
		nn.SetInputValue(vector_of_inputs[i]);
		nn.SetTargetValue(vector_of_targets[i]);
		nn.ForwardPropagation();
		cout.width(8);
		cout.setf(ios::showpoint);
		cout << setprecision(6);
		cout << nn.GetPredictValue(0) << "    \t" << nn.GetTargetValue(0)<< "    \t" << nn.GetTotalDerivedError()<< endl;
		output_file << nn.GetPredictValue(0) << "    \t" << nn.GetTargetValue(0) << "    \t" << nn.GetTotalDerivedError() << endl;
	}
	output_file.close();
	return 0;
}