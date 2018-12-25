#include "neuralnetwork.h"
#include <random>
NeuralNetwork::NeuralNetwork(vector<int> input_topology, double input_learning_rate) {
	topology = input_topology;
	topology_size = input_topology.size();
	learning_rate = input_learning_rate;
	bias = 1;

	// push layer class into layers vector
	for (int i = 0; i < topology_size; i++) {
		Layer input_layer(topology[i]);
		layers.push_back(input_layer);
	}
	//push weight matrix into weight_matrices vector.
	for (int i = 0; i < topology_size-1; i++) {
		//Matrix input_matrix(topology[i], topology[i + 1]);			//if it is input layer, we don't need bias
		//input_matrix = input_matrix.GetRandomMatrix(0, 1);
		//weight_matrices.push_back(input_matrix);
		if (i != 0) {
			Matrix input_matrix(topology[i]+1, topology[i + 1]);		//topology[i]+1, which plus weights for bias.
			input_matrix = input_matrix.GetRandomMatrix(0, 1);
			weight_matrices.push_back(input_matrix);
		}
		else if(i==0){
			Matrix input_matrix(topology[i], topology[i + 1]);			//if it is input layer, we don't need bias
			input_matrix = input_matrix.GetRandomMatrix(0, 1);
			weight_matrices.push_back(input_matrix);
		}
	}
}

void NeuralNetwork::SetInputValue(vector<double> input_value) {
	this->input_value = input_value;
	if (topology[0] != input_value.size()) {
		cout << "input value size != input layer size" << endl;
	}
	for (int i = 0; i < this->input_value.size(); i++) {
		layers[0].SetValueLayer(i, input_value[i]);
	}
}
void NeuralNetwork::SetTargetValue(vector<double> target_value) {
	this->target_value = target_value;
};
int NeuralNetwork::GetSizeNeuralNetwork() {
	return topology_size;
};
vector<Layer> NeuralNetwork::GetLayersNeuralNetwork() {
	return layers;
};
vector<Matrix> NeuralNetwork::GetWeightMatricesNeuralNetwork() {
	return weight_matrices;
};
Matrix NeuralNetwork::GetNeuronValueMatrix(int index) {
	return layers[index].MatrixifyValues();
};
Matrix NeuralNetwork::GetNeuronActivatedValueMatrix(int index) {
	return layers[index].MatrixifyActivatedValues();
};
Matrix NeuralNetwork::GetNeuronActivatedValueMatrixPlusBias(int index) {
	return layers[index].MatrixifyActivatedValuesPlusBias();
};
Matrix NeuralNetwork::GetNeuronDerivedValueMatrix(int index) {
	return layers[index].MatrixifyDerivedValues();
};
Matrix NeuralNetwork::GetWeightMatrix(int index) {
	return weight_matrices[index];
};
void NeuralNetwork::SetValueNeuralNetwork(int layer_index, int neuron_index, double neuron_value) {
	layers[layer_index].SetValueLayer(neuron_index, neuron_value);
};
void NeuralNetwork::ForwardPropagation() {
	for (int i = 0; i < GetSizeNeuralNetwork()-1; i++) {
		Matrix a(GetNeuronValueMatrix(i));					//if i=0, which means input layer
		if (i != 0) {
			a=GetNeuronActivatedValueMatrixPlusBias(i);			//if i!=0, which means hidden layer
		}
		Matrix b=GetWeightMatrix(i);
		Matrix c(a*b);	
		vector<double> result(c.MatrixToVector());			//make result matrix a vector
		for (int j = 0; j < result.size(); j++) {
			SetValueNeuralNetwork(i + 1, j, result[j]);			//push value to next layer neurons (+bias)
		}
	}
	SetErrors();		//renew errors vector
};

void NeuralNetwork::BackPropagation() {
	int output_layer_index = layers.size() - 1;
	int last_hidden_layer_index = output_layer_index - 1;	

	//part 1: handle output layer to hidden layer.
	Matrix derived_output_to_hidden(layers[output_layer_index].MatrixifyDerivedValues());		//create derived values matrix
	Matrix error_output_to_hidden(1, errors.size(), derived_errors);
	Matrix gradient_output_to_hidden(derived_output_to_hidden.PointwiseProduct(error_output_to_hidden));		// gradient matrix
	//calculate gradient matrix output to hidden
	Matrix delta_output_to_hidden((layers[last_hidden_layer_index].MatrixifyActivatedValuesPlusBias().TransposeMatrix())*gradient_output_to_hidden);
	//refresh weight matrices.
	weight_matrices[last_hidden_layer_index] = weight_matrices[last_hidden_layer_index] - delta_output_to_hidden*learning_rate;
	

	//hidden layer a: left side                hidden layer b: right side
	//gradient value of layer b
	Matrix gradient_hidden(gradient_output_to_hidden);

	//part 2: calcuate gradient matrix from hidden layer to input layer
	for (int i = last_hidden_layer_index; i > 0; i--) {
		//derived hidden value of layer b
		Matrix derived_hidden(layers[i].MatrixifyDerivedValues());
		//weight matrix from a to b
		Matrix weight_matrix_b;
		weight_matrix_b = weight_matrices[i].EraseLastRow();
		//calculate new gradient value of layer b
		gradient_hidden = derived_hidden.PointwiseProduct(gradient_hidden*(weight_matrix_b.TransposeMatrix()));
		Matrix neuron_value_a(layers[i - 1].MatrixifyActivatedValuesPlusBias());
		if (i == 1) {	//which means the layer a is input layer
			neuron_value_a = layers[0].MatrixifyValues();
		}
		Matrix delta_hidden(neuron_value_a.TransposeMatrix()*gradient_hidden);
		weight_matrices[i-1] = weight_matrices[i-1] - delta_hidden*learning_rate;
	}	
};



ostream& operator<< (ostream& os, NeuralNetwork& nn) {
	os << "*************BEGIN**************" << endl;
	for (int i = 0; i < nn.GetSizeNeuralNetwork(); i++) {
		os << "Layer = " << i << endl<<endl;
		if (i == 0) {
			os << "NeuronValue = "<< endl;
			os <<(nn.GetLayersNeuralNetwork())[0].MatrixifyValues();
			os << endl;
			os << "WeightMatrix = " << endl;
			os<<nn.GetWeightMatricesNeuralNetwork()[0];
			os << endl;
			os << "********************************"<<endl<<endl;
		}
		else {
			//os << "NeuronValue = " << endl;
			//os << (nn.GetLayersNeuralNetwork())[i].MatrixifyValues();
			//os << endl;
			os << "NeuronActivatedValue= "<<endl;
			os << (nn.GetLayersNeuralNetwork())[i].MatrixifyActivatedValues();
			os << endl;
			if (i < nn.GetSizeNeuralNetwork() -1) {
				os << "WeightMatrix (nrow + 1 because of bias value ) = " << endl;
				os << nn.GetWeightMatricesNeuralNetwork()[i];
			}
			os << endl;
			os << "********************************" << endl << endl;
		}
	}
	os << "Total Error = " << nn.GetTotalError() << endl<<endl;
	os << "**************END***************" << endl<<endl<<endl;
	return os;
};
double NeuralNetwork::GetTotalError() {
	return total_error;
};
double NeuralNetwork::GetTotalDerivedError() {
	return total_derived_error;
};
vector<double> NeuralNetwork::GetErrors() {
	return errors;
};
vector<double> NeuralNetwork::GetDerivedErrors() {
	return derived_errors;
};
void NeuralNetwork::SetErrors() {
	if (target_value.size() == 0){
		cout << "No target value for neural network!" << endl;
	}
	if (target_value.size() != layers[layers.size()-1].GetLayerSize()) {
		cout << "Target size do not match with output layer size!" << endl;
	}
	total_error = 0.0;
	total_derived_error = 0.0;
	int output_layer_index = layers.size() - 1;
	vector<Neuron> output_neurons = layers[output_layer_index].GetNeuronLayer();
	vector<double> result_errors;
	vector<double> result_derived_errors;
	for (int i = 0; i < target_value.size(); i++) {
		double temp_error = 0.5*(output_neurons[i].GetActivatedValue() - target_value[i])*(output_neurons[i].GetActivatedValue() - target_value[i]);
		double temp_derived_error = (output_neurons[i].GetActivatedValue() - target_value[i]);
		result_errors.push_back(temp_error);
		result_derived_errors.push_back(temp_derived_error);
		total_error += temp_error;
		total_derived_error += temp_derived_error;
	}
	errors = result_errors;
	derived_errors = result_derived_errors;

	total_error = total_error / topology[topology.size() - 1];
	total_derived_error = total_derived_error / topology[topology.size() - 1];
	historical_total_errors.push_back(total_error);
};

void NeuralNetwork::PrintHistoricalError() {
	for (int i = 0; i < historical_total_errors.size(); i++) {
		cout << i+1<<" iteration, total error = "<<historical_total_errors[i] << endl;
	}
};
double NeuralNetwork::GetPredictValue(int index) {
	int output_layer_index = layers.size() - 1;
	return layers[output_layer_index].GetNeuronLayer()[index].GetActivatedValue();
};
double NeuralNetwork::GetTargetValue(int index) {
	return target_value[index];
};