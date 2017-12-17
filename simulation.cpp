#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <random>


using namespace std;


void staircase(string& read, default_random_engine& generator){
	normal_distribution<double> distribution((uint32_t)read.size()/2, 30);
	double size = distribution(generator);
	read = read.substr(0, size);
}







char randomNucleotide(){
	switch (rand() % 4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}


string addHomopolymer(char& nuc){
	uint32_t length(rand() % 10);
	string homopolymer(length, nuc);
	return homopolymer;
}



string randomSequence(const uint length){
	string result(length, 'A');
	for(uint i(0); i < length; ++i){
		result[i] = randomNucleotide();
	}	
	return result;
}

void insertion(double rate, string& result){
	uint dice(rand() % 100);
	if(dice < rate){
		char newNucleotide(randomNucleotide());
		result.push_back(newNucleotide);
		insertion(rate, result);
	}
}


string mutateSequence(const string& referenceSequence, uint maxMutRate=10, vector <double> ratioMutation={0.06,0.73,0.21}){
	string result, currentNuc, homopoly;
	result.reserve(5 * referenceSequence.size());
	
	for(uint i(0); i < referenceSequence.size(); ++i){
		uint mutRate(maxMutRate);
		double substitutionRate(mutRate * ratioMutation[0]);
		double insertionRate(mutRate * ratioMutation[1]);
		double deletionRate(mutRate * ratioMutation[2]);
		uint dice(rand() % 100);

		//// homopolymers ////
		if (currentNuc.size() < 5){
			if (i > 0){
				if (currentNuc.empty()){
					currentNuc.push_back(result.back());
				} else {
					if (result.back() == currentNuc.back()){
						currentNuc.push_back(result.back());
					} else {
						currentNuc.back() = result.back();
					}
				}
			}
		} else {
			uint32_t probaHomopoly(rand() % 100);
			if (probaHomopoly < 1){
				homopoly = addHomopolymer(currentNuc.back());
				result += homopoly;
			}
			currentNuc = {};
		}
		////
		if (dice <substitutionRate ){
			//SUBSTITUTION
			char newNucleotide(randomNucleotide());
			while(newNucleotide == referenceSequence[i]){
				newNucleotide = randomNucleotide();
			}
			result.push_back(newNucleotide);
			continue;
		} else if(dice < deletionRate+substitutionRate){
			//DELETION
			uint dice2(rand() % 100);
			while (dice2 < deletionRate+substitutionRate){ // deletions larger than 1
				++i;
				dice2 = rand() % 100;
			}
			continue;
		} else if (dice < deletionRate + substitutionRate + insertionRate){
			//INSERTION
			char newNucleotide(randomNucleotide());
			result.push_back(referenceSequence[i]);
			result.push_back(newNucleotide);
			insertion(deletionRate + substitutionRate + insertionRate, result); // larger than 1 insertions
			
			continue;
		} else {
		//NO ERROR
			result.push_back(referenceSequence[i]);
		}

	}
	return result;
}

vector<string> split(const string &s, char delim){
  stringstream ss(s);
  string item;
  vector<string> elems;
  while (getline(ss, item, delim)) {
    elems.push_back(move(item)); 
  }
  return elems;
}



vector<pair<string, uint32_t>> generateTranscriptReferences(){
	string sequence;
	vector<pair<string, uint32_t>> result;
	string transcript;
	pair<string, uint32_t> transcriptAndExpression;
	////////// for testing only, does not take the gtf information at the moment ///////////
	for (uint32_t i(0); i < 10; ++i){
		transcriptAndExpression = { randomSequence(2000), rand() % 10};
		result.push_back(transcriptAndExpression);
	}
	return result;
}





void generateReads(uint32_t coverage, const string& outFileName="simulatedReads.fa", const string& outRefFileName="simulatedPerfectSequences.fa"){
	default_random_engine generator;
	ofstream out(outFileName);
	ofstream outRef(outRefFileName);
	vector<pair<string, uint32_t>> referenceList(generateTranscriptReferences()); // transcripts and their expression levels (= number of copies that appear)
	string read, reference;

	uint32_t noRead(0);
	for(uint32_t i(0);i < referenceList.size(); ++i){
		for (uint32_t numberOfReads(0); numberOfReads < referenceList[i].second * coverage; ++numberOfReads){
			reference = referenceList[i].first;
			read = mutateSequence(reference); // here we must add homopolymers
			// here we need the staircase effect
			staircase(read, generator);
			out <<">" << noRead << "|referenceNumber:" << i << endl << read << endl;
			outRef <<">" << noRead << "|referenceNumber:" << i << endl << reference << endl;
			++noRead;
		}
		
	}
}



int main(int argc, char ** argv){
	if (argc > 1){
		uint32_t coverage(stoi(argv[1]));
		srand (time(NULL));
		auto startChrono = chrono::system_clock::now();
		generateReads(coverage);
		auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
		cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	} else {
		cout << "usage: ./theReadCreator coverage" << endl;
	}
	return 0;
}

