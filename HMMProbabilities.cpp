/*
 * HMMProbablities.cpp
 *
 *	This is the cpp file for the HMMProbabilties object. HMMProbabilities
 *  is a collection of all probabilties needed for a hidden markov model. This
 *  includes initiation, emission and transition probabilties.  There are 
 *	convenience methods for setting and retriving probabilties as well as
 *  the log value of each probabilty.
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "HMMProbabilities.h"
#include <cmath>
#include <sstream>
#include <string>

// Constuctors
// ==============================================
HMMProbabilities::HMMProbabilities() {
}

HMMProbabilities::HMMProbabilities(int numOfStates) {
	numStates = numOfStates;
	createEmissionResidueMap();

	// Initialize all probabilities to zero
	for (int i = 0; i <= numStates; i++) {
		setInitiationProbability(i, 0);
		for (int j = 0; j <= numStates; j++) {
			setTransitionProbability(i, j, 0);
		}
		for (pair<string, int> mapPair : emissionResidueMap) {
			string& residue = mapPair.first;
			setEmissionProbability(i, residue, 0);
		}
	}
}


// Destructor
// =============================================
HMMProbabilities::~HMMProbabilities() {
}

// Public Class Methods
// =============================================

// HMMProbabilities* initialProbabilities()
//  Purpose: 
//		Returns a probabilites object initialzed to the initial
//		probabilites required by genome540 homework #7	
HMMProbabilities* HMMProbabilities::initialProbabilities() {

	HMMProbabilities* probs = new HMMProbabilities(12);

	// initiation probabilties
	probs->setInitiationProbability(6, 1.0);

	// transition probabilities
	probs->setTransitionProbability(1, 3, 1.0);
	probs->setTransitionProbability(2, 3, 1.0);
	probs->setTransitionProbability(3, 4, 1.0);
	probs->setTransitionProbability(4, 2, 0.99);
	probs->setTransitionProbability(4, 5, 0.01);
	probs->setTransitionProbability(5, 6, 1.0);
	probs->setTransitionProbability(6, 1, 0.8);
	probs->setTransitionProbability(6, 6, 0.1);
	probs->setTransitionProbability(6, 11, 0.1);
	probs->setTransitionProbability(7, 6, 1.0);
	probs->setTransitionProbability(8, 9, 1.0);
	probs->setTransitionProbability(9, 10, 1.0);
	probs->setTransitionProbability(10, 7, 0.01);
	probs->setTransitionProbability(10, 8, 0.99);
	probs->setTransitionProbability(11, 9, 1.0);

	// emission probabilities
	// State 1 (Top-strand Start Codon)
	double state1DefaultValue = 1/(double)4;
	probs->setEmissionProbability(1, "ATG", state1DefaultValue);
	probs->setEmissionProbability(1, "CTG", state1DefaultValue);
	probs->setEmissionProbability(1, "GTG", state1DefaultValue);
	probs->setEmissionProbability(1, "TTG", state1DefaultValue);

	// State 2 (1st base of top-strand internal codon)
	double state2DefaultValue = 1/(double)61;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(2, residue, state2DefaultValue);
	}
	probs->setEmissionProbability(2, "TAA", 0);
	probs->setEmissionProbability(2, "TGA", 0);
	probs->setEmissionProbability(2, "TAG", 0);

	// State 3 (2nd Base of top-strand start or internal codon)
	double state3DefaultValue = 1/(double)64;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(3, residue, state3DefaultValue);
	}

	// State 4 (3rd Base of top-strand start or internal codon)
	double state4DefaultValue = 1/(double)64;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(4, residue, state4DefaultValue);
	}

	// State 5 (Top-strand Stop Codon)
	double state5DefaultValue = 1/(double)3;
	probs->setEmissionProbability(5, "TAA", state5DefaultValue);
	probs->setEmissionProbability(5, "TGA", state5DefaultValue);
	probs->setEmissionProbability(5, "TAG", state5DefaultValue);

	// State 6 (intergenic)
	double state6DefaultValue = 1/(double)64;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(6, residue, state6DefaultValue);
	}
	
	// State 7 (Bottom-strand Start Codon)
	double state7DefaultValue = 1/(double)4;
	probs->setEmissionProbability(7, "CAA", state7DefaultValue);
	probs->setEmissionProbability(7, "CAC", state7DefaultValue);
	probs->setEmissionProbability(7, "CAG", state7DefaultValue);
	probs->setEmissionProbability(7, "CAT", state7DefaultValue);

	// State 8 (1st base of bottom-strand internal codon)
	double state8DefaultValue = 1/(double)61;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(8, residue, state8DefaultValue);
	}
	probs->setEmissionProbability(8, "TTA", 0);
	probs->setEmissionProbability(8, "TCA", 0);
	probs->setEmissionProbability(8, "CTA", 0);

	// State 9 (2nd Base of bottom-strand start or internal codon)
	double state9DefaultValue = 1/(double)64;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(9, residue, state9DefaultValue);
	}

	// State 10 (3rd Base of bottom-strand start or internal codon)
	double state10DefaultValue = 1/(double)64;
	for (pair<string, int> mapPair : probs->emissionResidueMap) {
		string& residue = mapPair.first;
		probs->setEmissionProbability(10, residue, state10DefaultValue);
	}

	// State 11 (Bottom-strand Stop Codon)
	double state11DefaultValue = 1/(double)3;
	probs->setEmissionProbability(11, "TTA", state11DefaultValue);
	probs->setEmissionProbability(11, "TCA", state11DefaultValue);
	probs->setEmissionProbability(11, "CTA", state11DefaultValue);

	return probs;

}

// Public Methods
// =============================================

// double emissionProbability(int state, char residue)
//  Purpose: 
//		Returns the emission probability for the state and residue
long double HMMProbabilities::emissionProbability(int state, string residue) {
	return emissionProbabilities.at(state).at(getEmissionResidueIndex(residue));
}

// double initiationProbability(int state)
//  Purpose: 
//		Returns the initiation probability for the state
long double HMMProbabilities::initiationProbability(int state) {
	return initiationProbabilities[state];
}
	
// double transitionProbability(int beginState, int endState)
//  Purpose: 
//		Returns the transition probability for transition from beginState
//		to endState
long double HMMProbabilities::transitionProbability(int beginState, int endState) {
	return transitionProbabilities[beginState][endState];
}

// double logEmissionProbability(int state, char residue)
//  Purpose: 
//		Returns the log of the emission probability for the state and residue
long double HMMProbabilities::logEmissionProbability(int state, string residue) {
	return logEmissionProbabilities.at(state).at(getEmissionResidueIndex(residue));
}

// double logInitiationProbability(int state)
//  Purpose: 
//		Returns the log of the initiation probability for the state
long double HMMProbabilities::logInitiationProbability(int state) {
	return logInitiationProbabilities[state];
}
	
// double logTransitionProbability(int beginState, int endState)
//  Purpose: 
//		Returns the log of the transition probability for transition from
//		beginState to endState
long double HMMProbabilities::logTransitionProbability(int beginState, int endState) {
	return logTransitionProbabilities[beginState][endState];
}

// setEmissionProbability(int state, char residue, double value)
//  Purpose: 
//		Sets the emission probability for the state and residue to value
//	Postconditions:
//		emissionProbabilites - value set for state/residue
//		logEmissionProbabilites - value set for state/residue
void HMMProbabilities::setEmissionProbability(int state, string residue, long double value) {
	emissionProbabilities[state][getEmissionResidueIndex(residue)] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logEmissionProbabilities[state][getEmissionResidueIndex(residue)] = logVal;
}

// setInitiationProbability(int state, double value)
//  Purpose: 
//		Sets the initiation probability for the state to value
//	Postconditions:
//		initiationProbabilites - value set for state
//		logInitiationProbabilites - value set for state
void HMMProbabilities::setInitiationProbability(int state, long double value) {
	initiationProbabilities[state] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logInitiationProbabilities[state] = logVal;
}

// setTransitionProbability(int beginState, int endState, double value)
//  Purpose: 
//		Sets the transition probability from the beginState to the endState
//		to value
//	Postconditions:
//		transitionProbabilites - value set for beginState to endState
//		logTransitionProbabilites - value set for beginState to endState
void HMMProbabilities::setTransitionProbability(int beginState, int endState, long double value) {
	transitionProbabilities[beginState][endState] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logTransitionProbabilities[beginState][endState] = logVal;
}

// string probabilitiesResultsString()
//  Purpose:
//		Returns a string representing the probabilites
//
//		format:
//			<<statesResultsString>>
//			<<initiationProbabilitesResultsString>>
//			<<transmissionProbabilitesResultsString>>
//			...
//			<<emissionProbabilitesResultsString>>
//			...
string HMMProbabilities::probabilitiesResultsString() {
	stringstream ss;

	// Begin Model
	ss << "      <model type=\"hmm\">\n";

	// States
	ss << statesResultsString();

	// Probabiltiies
	ss << intitiationProbabiltiesResultsString();
	for (int i = 1; i < numStates; i++)
		ss << transitionProbablitiesResultsString(i);
	for (int i = 1; i < numStates; i++)
		ss << emissionProbablitiesResultsString(i);
	
	// End Model
	ss << "      </model>\n";

	return ss.str();
}

// string statesResultsString()
//  Purpose:
//		Returns a string representing the states
//
//		format:
//			<result type="states">
//				<<state1>>,<<state2>>,...
//			</result>
string HMMProbabilities::statesResultsString() {
	stringstream ss;

	// Header 
	ss << "        <states>";

	// States
	for (int i = 1; i < numStates; i++) {
		ss << i ;

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</states>\n";

	return ss.str();
}

// string intitiationProbabiltiesResultsString()
//  Purpose:
//		Returns a string representing the initiation probablities
//
//		format:
//			<result type="initiation_probabilites">
//				<<state>>=<<initiation probability>>,
//			</result>
string HMMProbabilities::intitiationProbabiltiesResultsString() {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <initial_state_probabilities>";

	// States
	for (int i = 1; i < numStates; i++) {
		ss
			<< i 
			<< "="
			<< initiationProbability(i);

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</initial_state_probabilities>\n";

	return ss.str();
}

// transitionProbablitiesResultsString(int state)
//  Purpose:
//		Returns a string representing the transition probablities for a state
//
//		format:
//			<result type="transition_probabilites" state="<<state>>">
//				<<to state>>=<<transition probability>>,
//			</result>
string HMMProbabilities::transitionProbablitiesResultsString(int state) {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <transition_probabilities state=\"" << state + 1 << "\">";

	// States
	for (int i = 1; i < numStates; i++) {
		ss
			<< i
			<< "="
			<< transitionProbability(state, i);

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</transition_probabilities>\n";

	return ss.str();
}

// emissionProbablitiesResultsString(int state)
//  Purpose:
//		Returns a string representing the emission probablities for a state
//
//		format:
//			<result type="emission_probabilites" state="<<state>>">
//				<<residue>>=<<emission probability>>,
//			</result>
string HMMProbabilities::emissionProbablitiesResultsString(int state) {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <emission_probabilities state=\"" << state + 1 << "\">";

	// Residues
	for (pair<string, int> mapPair : emissionResidueMap) {
		string& residue = mapPair.first;
		ss << residue << "=" << emissionProbability(state, residue) << ",";
	}

	// Footer
	ss << "</emission_probabilities>\n";

	return ss.str();
}

// map<string, int> createEmissionMap()
//  Purpose: 
//		Creates a map of the index location for a trinucleotide emission
//		in the emission probabilities array
void HMMProbabilities::createEmissionResidueMap() {

	emissionResidueMap["AAA"]= 0;
	emissionResidueMap["AAC"]= 1;
	emissionResidueMap["AAG"]= 2;
	emissionResidueMap["AAT"]= 3;
	emissionResidueMap["ACA"]= 4;
	emissionResidueMap["ACC"]= 5;
	emissionResidueMap["ACG"]= 6;
	emissionResidueMap["ACT"]= 7;
	emissionResidueMap["AGA"]= 8;
	emissionResidueMap["AGC"]= 9;
	emissionResidueMap["AGG"]= 10;
	emissionResidueMap["AGT"]= 11;
	emissionResidueMap["ATA"]= 12;
	emissionResidueMap["ATC"]= 13;
	emissionResidueMap["ATG"]= 14;
	emissionResidueMap["ATT"]= 15;
	emissionResidueMap["CAA"]= 16;
	emissionResidueMap["CAC"]= 17;
	emissionResidueMap["CAG"]= 18;
	emissionResidueMap["CAT"]= 19;
	emissionResidueMap["CCA"]= 20;
	emissionResidueMap["CCC"]= 21;
	emissionResidueMap["CCG"]= 22;
	emissionResidueMap["CCT"]= 23;
	emissionResidueMap["CGA"]= 24;
	emissionResidueMap["CGC"]= 25;
	emissionResidueMap["CGG"]= 26;
	emissionResidueMap["CGT"]= 27;
	emissionResidueMap["CTA"]= 28;
	emissionResidueMap["CTC"]= 29;
	emissionResidueMap["CTG"]= 30;
	emissionResidueMap["CTT"]= 31;
	emissionResidueMap["GAA"]= 32;
	emissionResidueMap["GAC"]= 33;
	emissionResidueMap["GAG"]= 34;
	emissionResidueMap["GAT"]= 35;
	emissionResidueMap["GCA"]= 36;
	emissionResidueMap["GCC"]= 37;
	emissionResidueMap["GCG"]= 38;
	emissionResidueMap["GCT"]= 39;
	emissionResidueMap["GGA"]= 40;
	emissionResidueMap["GGC"]= 41;
	emissionResidueMap["GGG"]= 42;
	emissionResidueMap["GGT"]= 43;
	emissionResidueMap["GTA"]= 44;
	emissionResidueMap["GTC"]= 45;
	emissionResidueMap["GTG"]= 46;
	emissionResidueMap["GTT"]= 47;
	emissionResidueMap["TAA"]= 48;
	emissionResidueMap["TAC"]= 49;
	emissionResidueMap["TAG"]= 50;
	emissionResidueMap["TAT"]= 51;
	emissionResidueMap["TCA"]= 52;
	emissionResidueMap["TCC"]= 53;
	emissionResidueMap["TCG"]= 54;
	emissionResidueMap["TCT"]= 55;
	emissionResidueMap["TGA"]= 56;
	emissionResidueMap["TGC"]= 57;
	emissionResidueMap["TGG"]= 58;
	emissionResidueMap["TGT"]= 59;
	emissionResidueMap["TTA"]= 60;
	emissionResidueMap["TTC"]= 61;
	emissionResidueMap["TTG"]= 62;
	emissionResidueMap["TTT"]= 63;
}

// int getIndex(char residue)
//  Purpose: 
//	  Returns the index in the emission probabilities for the residue
int HMMProbabilities::getEmissionResidueIndex(string residue) {
	return emissionResidueMap.at(residue);
}
