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



uint32_t SD_GAUSSIAN(30);
uint32_t SIZE_TO_START_HOMOPOL(5);
uint32_t SIZE_MAX_HOMOPOLY(10);
uint32_t SIZE_MIN_HOMOPOLY(3);

void staircase(string& read, default_random_engine& generator){
	normal_distribution<double> distribution((uint32_t)read.size()/2, SD_GAUSSIAN);
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
	uint32_t length(0);
	while (length < 1){
		length = rand() % SIZE_MAX_HOMOPOLY;
	}
	string homopolymer(length, nuc);
	return homopolymer;
}


void removeHomopolymer(uint64_t size, string& result){
	string tmpResult("");
	uint32_t length(0);
	while (length < 1){
		length = rand() % (size - SIZE_MIN_HOMOPOLY);
	}
	for (uint64_t i(0); i < (result.size() - length); ++i){
		tmpResult += result[i];
	}
	result = tmpResult;
}


string randomSequence(const uint32_t length){
	string result(length, 'A');
	for(uint32_t i(0); i < length; ++i){
		result[i] = randomNucleotide();
	}	
	return result;
}

void insertion(double rate, string& result){
	uint64_t dice(rand() % 100*100*100);
	if(dice < rate){
		char newNucleotide(randomNucleotide());
		result.push_back(newNucleotide);
		insertion(rate, result);
	}
}


string mutateSequence(const string& referenceSequence, 	unordered_map<string, double>& errorProfile){
//~ string mutateSequence(const string& referenceSequence, uint maxMutRate=10, vector <double> ratioMutation={0.06,0.73,0.21}){
	string result, currentNuc, homopoly;
	result.reserve(5 * referenceSequence.size());


	// get each error value
	// sort mism/ins/del
	double errors [] = {errorProfile["mismatches"], errorProfile["non-homopolymer_ins"], errorProfile["non-homopolymer_del"], errorProfile["homopolymer_ins"], errorProfile["non-homopolymer_del"]};
	double errorsHomo [] = {errorProfile["homopolymer_ins"], errorProfile["non-homopolymer_del"]};
	vector<double> errorOrder(errors, errors + 3);
	vector<double> errorOrderHomo(errorsHomo, errorsHomo + 2);
	sort(errorOrder.begin(), errorOrder.end());
	sort(errorOrderHomo.begin(), errorOrderHomo.end());
	map<string, double> errorProfileSorted;
	map<string, double> errorProfileHomoSorted;
	double prevValue(0);
	for (auto&& val : errorOrder){
		if (val == errorProfile["mismatches"]){
			errorProfileSorted.insert({"mismatches", val * 100*100*100 + prevValue});
		} else if (val == errorProfile["non-homopolymer_ins"]){
			errorProfileSorted.insert({"non-homopolymer_ins", val * 100*100*100 + prevValue});
		} else if (val == errorProfile["non-homopolymer_del"]){
			errorProfileSorted.insert({"non-homopolymer_del", val * 100*100*100 + prevValue});
		}
		prevValue += val;
	}
	prevValue = 0;
	for (auto&& val : errorOrderHomo){
		if (val == errorProfile["homopolymer_ins"]){
			errorProfileHomoSorted.insert({"homopolymer_ins", val * 100*100*100 + prevValue});
		} else if (val == errorProfile["homopolymer_del"]){
			errorProfileHomoSorted.insert({"homopolymer_del", val * 100*100*100 + prevValue});
		}
		prevValue += val;
	}

	uint32_t i(0);
	bool isError(false);

	while (i < referenceSequence.size()){
		
		isError = false;
		uint dice(rand() % 100*100*100);

		//// homopolymers ////
		if (currentNuc.size() <  SIZE_TO_START_HOMOPOL){
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
			
			uint32_t probaHomopoly(rand() % 100*100*100);
			for (auto it(errorProfileHomoSorted.begin()); it !=errorProfileHomoSorted.end(); ++it){
				if (probaHomopoly < it->second){
					if (it->first == "homopolymer_del"){
						removeHomopolymer(currentNuc.size(), result);
						isError = true;
						break;
					} else if (it->first == "homopolymer_ins"){
						homopoly = addHomopolymer(currentNuc.back());
						result += homopoly;
						isError = true;
						break;
					}
				}
			}
			currentNuc = {};
		}
		if (not isError){  // if not already added an error in homopolymer at this base
			for (auto it(errorProfileSorted.begin()); it !=errorProfileSorted.end(); ++it){
				if (dice < it->second){
					if (it->first == "mismatches"){
						//SUBSTITUTION
						char newNucleotide(randomNucleotide());
						while(newNucleotide == referenceSequence[i]){
							newNucleotide = randomNucleotide();
						}
						result.push_back(newNucleotide);
						isError = true;
						break;
					} else if (it->first == "non-homopolymer_del"){
						// DELETION
						uint dice2(rand() % 100*100*100);
						while (dice2 < it->second and i < referenceSequence.size()){ // deletions larger than 1
							++i;
							dice2 = rand() % 100*100*100;
						}
						isError = true;
						break;
					} else if (it->first == "non-homopolymer_ins"){
						// INSERTION
						char newNucleotide(randomNucleotide());
						result.push_back(referenceSequence[i]);
						result.push_back(newNucleotide);
						insertion(it->second, result); // larger than 1 insertions
						isError = true;
						break;
					}
				} 
			}
			if (not isError){
				//NO ERROR
				result.push_back(referenceSequence[i]);
			}
		}
		++i;
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


//// this function needs to be updated with the output filled with sequences from the ref transcriptome and expressions levels from FS ////
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





void generateReads(uint32_t coverage, unordered_map<string, double>& errorProfile, const string& outFileName="simulatedReads.fa", const string& outRefFileName="simulatedPerfectSequences.fa"){
	default_random_engine generator;
	ofstream out(outFileName);
	ofstream outRef(outRefFileName);
	vector<pair<string, uint32_t>> referenceList(generateTranscriptReferences()); // transcripts and their expression levels (= number of copies that appear)
	string read, reference;

	uint32_t noRead(0);
	for(uint32_t i(0);i < referenceList.size(); ++i){
		for (uint32_t numberOfReads(0); numberOfReads < referenceList[i].second * coverage; ++numberOfReads){
			reference = referenceList[i].first;
			read = mutateSequence(reference, errorProfile); // here we must add homopolymers
			staircase(read, generator);
			out <<">" << noRead << "|referenceNumber:" << i << endl << read << endl;
			outRef <<">" << noRead << "|referenceNumber:" << i << endl << reference << endl;
			++noRead;
		}
		
	}
}



void getErrorProfile(string& errorProfileFileName, unordered_map<string, double>& errorProfile){
	ifstream in(errorProfileFileName);
	string errorLine;
	vector<string> errorType;
	double rate;
	while (not in.eof()){
		getline(in, errorLine);
		errorType = split(errorLine, ' ');
		if (errorType.size() == 2){
			rate = stod(errorType[1]);
			errorProfile.insert({errorType[0], rate});
		}
	}
}

int main(int argc, char ** argv){
	uint32_t coverage(0);
	string errorProfileFileName("");
	srand (time(NULL));
	auto startChrono = chrono::system_clock::now();
	int32_t c(0);
	while ((c = getopt (argc, argv, "c:e:")) != -1){
		switch(c){
			case 'c':
				coverage=stoi(optarg);
				break;
			case 'e':
				errorProfileFileName=optarg;
				break;
		}
	}
	if(coverage == 0 or errorProfileFileName ==""){
		cout
		<<"Mandatory parameters"<<endl
		<<"-c coverage"<<endl
		<<"-e error profile file"<<endl
		<< "usage: ./theReadCreator -c coverage -e ERROR_PROFILE_FILE" << endl;
		return 0;
	}
	unordered_map<string, double> errorProfile;
	getErrorProfile(errorProfileFileName, errorProfile);
	generateReads(coverage, errorProfile);
	auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
	cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;



	return 0;
}
