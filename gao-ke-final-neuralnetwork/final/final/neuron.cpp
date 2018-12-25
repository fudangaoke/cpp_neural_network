#include "neuron.h"
#include <iostream>
using namespace std;
Neuron::Neuron(double input_value) {
	value = input_value;
	Activate();
	Derive();
};
//bipolar sigmoid function
void Neuron::Activate() {
	activatedValue = (2.0 / (1 + exp(-value))) - 1;
};
void Neuron::Derive() {
	derivedValue = 0.5 * (1 + activatedValue)*(1 - activatedValue);
};
double Neuron::GetNeuronValue() {
	return value;
};
double Neuron::GetActivatedValue() {
	return activatedValue;
};
double Neuron::GetDerivedValue() {
	return derivedValue;
};
void Neuron::SetValueNeuron(double input_value) {
	value = input_value;
	Activate();
	Derive();
};
void Neuron::MakeBiasNeuron() {
	value = 0;
	activatedValue = 1;
	derivedValue = 0;
};