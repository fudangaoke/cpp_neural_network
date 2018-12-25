#include "layer.h"
#include <iostream>
using namespace std;
Layer::Layer(int input_layer_size) {
	layer_size = input_layer_size;
	for (int i = 0; i < layer_size; i++) {
		Neuron input_neuron(0.00);
		neurons.push_back(input_neuron);
	}
};

void Layer::SetValueLayer(int neuron_index, double neuron_value) {
	neurons[neuron_index].SetValueNeuron(neuron_value);
};

Matrix Layer::MatrixifyValues() {
	vector<double> layer_values_vector;
	for (int i = 0; i < layer_size; i++) {
		layer_values_vector.push_back(neurons[i].GetNeuronValue());
	}
	Matrix layer_values_matrix(1, layer_size,layer_values_vector);
	return layer_values_matrix;
};

Matrix Layer::MatrixifyActivatedValues() {
	vector<double> layer_activated_values_vector;
	for (int i = 0; i < layer_size; i++) {
		layer_activated_values_vector.push_back(neurons[i].GetActivatedValue());
	}
	Matrix layer_activated_values_matrix(1, layer_size, layer_activated_values_vector);
	return layer_activated_values_matrix;
};
Matrix Layer::MatrixifyActivatedValuesPlusBias() {
	vector<double> layer_activated_values_vector_plus_bias;
	for (int i = 0; i < layer_size; i++) {
		layer_activated_values_vector_plus_bias.push_back(neurons[i].GetActivatedValue());
	}
	Neuron bias_neuron(0);
	bias_neuron.MakeBiasNeuron();
	layer_activated_values_vector_plus_bias.push_back(bias_neuron.GetActivatedValue());
	Matrix layer_activated_values_matrix_plus_bias(1, layer_size+1, layer_activated_values_vector_plus_bias);		//+1 stands for bias
	return layer_activated_values_matrix_plus_bias;
};

Matrix Layer::MatrixifyDerivedValues() {
	vector<double> layer_derived_values_vector;
	for (int i = 0; i < layer_size; i++) {
		layer_derived_values_vector.push_back(neurons[i].GetDerivedValue());
	}
	Matrix layer_derived_values_matrix(1, layer_size, layer_derived_values_vector);
	return layer_derived_values_matrix;
};
int Layer::GetLayerSize() {
	return layer_size;
};
vector<Neuron> Layer::GetNeuronLayer() {
	return neurons;
};