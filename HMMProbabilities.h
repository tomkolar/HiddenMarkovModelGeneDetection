/*
 * HMMProbablities.h
 *
 *	This is the header file for the HMMProbabilties object. HMMProbabilities
 *  is a collection of all probabilties needed for a hidden markov model. This
 *  includes initiation, emission and transition probabilties.  There are 
 *	convenience methods for setting and retriving probabilties as well as
 *  the log value of each probabilty.
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */

#ifndef HMMPROBABILITIES_H
#define HMMPROBABILITIES_H
#include <map>
#include <string>
using namespace std;

class HMMProbabilities
{
public:
	// Constuctors
	// ==============================================
	HMMProbabilities();

	// Destructor
	// =============================================
	~HMMProbabilities();

	// Public Class Methods
	// =============================================

	// HMMProbabilities* initialProbabilities()
	//  Purpose: 
	//		Returns a probabilites object initialzed to the initial
	//		probabilites required by genome540 homework #5	
	static HMMProbabilities* initialProbabilities();

	// HMMProbabilities* testProbabilities()
	//  Purpose: 
	//		Returns a probabilites object initialzed to the probabilites
	//		required by the viterbi toy example
	//			homepages.ulb.ac.be/~dgonze/TEACHING/viterbi.pdf
	static HMMProbabilities* testProbabilities();

	// Public Methods
	// =============================================

	// double emissionProbability(int state, char residue)
	//  Purpose: 
	//		Returns the emission probability for the state and residue
	long  double emissionProbability(int state, char residue);

	// double initiationProbability(int state)
	//  Purpose: 
	//		Returns the initiation probability for the state
	long  double initiationProbability(int state);

	// double transitionProbability(int beginState, int endState)
	//  Purpose: 
	//		Returns the transition probability for transition from beginState
	//		to endState
	long  double transitionProbability(int beginState, int endState);

	// double logEmissionProbability(int state, char residue)
	//  Purpose: 
	//		Returns the log of the emission probability for the state and residue
	long  double logEmissionProbability(int state, char residue);

	// double logInitiationProbability(int state)
	//  Purpose: 
	//		Returns the log of the initiation probability for the state
	long  double logInitiationProbability(int state);

	// double logTransitionProbability(int beginState, int endState)
	//  Purpose: 
	//		Returns the log of the transition probability for transition from
	//		beginState to endState
	long double logTransitionProbability(int beginState, int endState);
	
	// setEmissionProbability(int state, char residue, double value)
	//  Purpose: 
	//		Sets the emission probability for the state and residue to value
	//	Postconditions:
	//		emissionProbabilites - value set for state/residue
	//		logEmissionProbabilites - value set for state/residue
	void setEmissionProbability(int state, char residue, long double value);

	// setInitiationProbability(int state, double value)
	//  Purpose: 
	//		Sets the initiation probability for the state to value
	//	Postconditions:
	//		initiationProbabilites - value set for state
	//		logInitiationProbabilites - value set for state
	void setInitiationProbability(int state, long double value);

	// setTransitionProbability(int beginState, int endState, double value)
	//  Purpose: 
	//		Sets the transition probability from the beginState to the endState
	//		to value
	//	Postconditions:
	//		transitionProbabilites - value set for beginState to endState
	//		logTransitionProbabilites - value set for beginState to endState
	void setTransitionProbability(int beginState, int endState, long double value);

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
	string probabilitiesResultsString();

	// string statesResultsString()
	//  Purpose:
	//		Returns a string representing the states
	//
	//		format:
	//			<result type="states">
	//				<<state1>>,<<state2>>,...
	//			</result>
	string statesResultsString();

	// string intitiationProbabiltiesResultsString()
	//  Purpose:
	//		Returns a string representing the initiation probablities
	//
	//		format:
	//			<result type="initiation_probabilites">
	//				<<state>>=<<initiation probability>>,
	//			</result>
	string intitiationProbabiltiesResultsString();

	// transitionProbablitiesResultsString(int state)
	//  Purpose:
	//		Returns a string representing the transition probablities for a state
	//
	//		format:
	//			<result type="transition_probabilites" state="<<state>>">
	//				<<to state>>=<<transition probability>>,
	//			</result>
	string transitionProbablitiesResultsString(int state);

	// emissionProbablitiesResultsString(int state)
	//  Purpose:
	//		Returns a string representing the emission probablities for a state
	//
	//		format:
	//			<result type="emission_probabilites" state="<<state>>">
	//				<<residue>>=<<emission probability>>,
	//			</result>
	string emissionProbablitiesResultsString(int state);

private:

	// Private Attributes
	// =============================================
	map<int, map<char, long double>> emissionProbabilities;
	map<int, map<char, long double>> logEmissionProbabilities;
	long double transitionProbabilities[2][2];
	long double logTransitionProbabilities[2][2];
	long double initiationProbabilities[2];
	long double logInitiationProbabilities[2];

};

#endif // HMMPROBABILITIES_H
