#pragma once
class Neuron {
private:
	double value;							//value
	double activatedValue;			//activated value
	double derivedValue;				//derived value
public:
	Neuron(double input_value);
	void Activate();								//bipolar sigmoid function
	void Derive();									//derive bipolar sigmoid function
	double GetNeuronValue();	
	double GetActivatedValue();
	double GetDerivedValue();
	void SetValueNeuron(double input_value);
	void MakeBiasNeuron();

};