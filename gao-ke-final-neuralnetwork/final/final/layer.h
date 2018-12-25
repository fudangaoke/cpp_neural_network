#pragma once
#include "neuron.h"
#include "matrix.h"
#include <vector>
#include <iostream>
using namespace std;
class Layer {
private:
	int layer_size;								//The number of neurons in each layer.
	vector<Neuron> neurons;	//Neuron vector in each layer.
public:
	Layer(int input_layer_size);
	void SetValueLayer(int neuron_index, double neuron_value);		//set value for neuron in each layer.
	int GetLayerSize();
	vector<Neuron> GetNeuronLayer();
	Matrix MatrixifyValues();     //make layer value a (1*size) matrix.
	Matrix MatrixifyActivatedValues();     //make layer activated value a (1*size) matrix.
	Matrix MatrixifyActivatedValuesPlusBias();     //make layer activated value a (1*size+1) matrix.
	Matrix MatrixifyDerivedValues();     //make layer derived value a (1*size) matrix.
};