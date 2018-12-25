#pragma once
#include <iostream>
#include "matrix.h"
#include "layer.h"
using namespace std;
class NeuralNetwork {
private:
	vector<int> topology;						//The topology of neural network
	vector<Layer> layers;						//The vector of layers
	vector<Matrix> weight_matrices;	//The weight matrices vector.
	int topology_size;							//The number of layers in neural network
	vector<double> input_value;			//The vector of input value of neural network
	vector<double> target_value;		//target value of output layer.
	double total_error;							//total error of neural network
	vector<double> errors;					//the vector of errors for each neuron output in the output layer.
	vector<double> derived_errors;					//the vector of derived errors for each neuron output in the output layer.
	double total_derived_error;
	vector<double> historical_total_errors;	//historical total errors.
	double learning_rate;						     //learning rate
	double bias;
public:
	NeuralNetwork(vector<int> input_topology, double input_learning_rate);		//Build a neural network with certain topology
	void SetInputValue(vector<double> input_value);	//set input value of neural network
	void SetTargetValue(vector<double> target_value);	//set target value of neural network
	//void SetInputTargetValue(vector<vector<double>> input_vector_of_inputs, vector<vector<double>> input_vector_of_targets);
	int GetSizeNeuralNetwork();									//get number of layers.
	vector<Layer> GetLayersNeuralNetwork();			//get vector of layers.
	vector<Matrix> GetWeightMatricesNeuralNetwork();	//get vector of weight matrices
	void ForwardPropagation();												//do forward propagation.
	void BackPropagation();											//do backward propagation
	Matrix GetNeuronValueMatrix(int index);				//get neuron value of layer[index]
	Matrix GetNeuronActivatedValueMatrix(int index);
	Matrix GetNeuronActivatedValueMatrixPlusBias(int index);//get neuron activated value of layer[index]
	Matrix GetNeuronDerivedValueMatrix(int index);//get neuron derived value of layer[index]
	Matrix GetWeightMatrix(int index);					//get weight matrix of weight_matrices[index]
	void SetValueNeuralNetwork(int layer_index, int neuron_index, double neuron_value);	//set value of neuron of layers[layer_index][neuron_index] = neuron value.
	friend ostream& operator<< (ostream& os, NeuralNetwork& nn);		//print neural network
	double GetTotalError();								//get total error
	vector<double> GetErrors();						//get error vector
	vector<double> GetDerivedErrors();
	void SetErrors();											//calculate errors of target and output values.
	void PrintHistoricalError();							//print historical error
	double GetTotalDerivedError();
	double GetPredictValue(int index);
	double GetTargetValue(int index);
};