/*
 * driver.cpp
 *
 *	This is the driver file for creating an a hidden markov model from
 *  a fastafile.  Viterbi training is then exectued to try and locate
 *  GC rich portions of the sequence.
 *
 *	Typical use:
 *		hmm fastaFile numIterations
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "FastaFile.h"
#include "HiddenMarkovModel.h"
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

int main( int argc, char *argv[] ) {
/*
    // Check that file name and iterations were entered as arguments
    if (argc < 3) {
            cout << "Invalid # of arguments\n";
            cout << "usage: hmm fastaFile iterations\n";
            return -1;
    }

    cout << "Starting\n";

    // Get Parameters
    string fastaFileName = argv[1];
    int iterations = atoi(argv[2]);
*/
	// Set Fasta File names
	string fastaFileName = "c:/Users/kolart/Documents/Genome540/Assignment6/shigella_short.fna";
	int iterations = 1;

	// Create the fasta file object
	FastaFile* fastaFile = new FastaFile(fastaFileName);
	cout << "Fasta's done\n";

	cout << fastaFile->firstLineResultString();

	// Create the Hidden Markov Model
	HiddenMarkovModel hmm(fastaFile);
	hmm.viterbiTraining(iterations);

	cout << hmm.viterbiResultsString();
	cout << "Fred";
}
