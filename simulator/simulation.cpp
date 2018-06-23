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


//// truncate reads for the staircase effect ////
void staircase(string& read, default_random_engine& generator){
	normal_distribution<double> distribution((uint32_t)read.size()/2, SD_GAUSSIAN);
	double size = distribution(generator);
	read = read.substr(0, size);
}



//// get a random nucleotide ////
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


//// create an insertion in homopolymers ////
string addHomopolymer(char& nuc){
	uint32_t length(0);
	while (length < 1){
		length = rand() % SIZE_MAX_HOMOPOLY;
	}
	string homopolymer(length, nuc);
	return homopolymer;
}


//// create a deletion in homopolymers ////
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


//// generate a random seq of a given length ////
string randomSequence(const uint32_t length){
	string result(length, 'A');
	for(uint32_t i(0); i < length; ++i){
		result[i] = randomNucleotide();
	}	
	return result;
}


//// add insertion ////
void insertion(double rate, string& result){
	uint64_t dice(rand() % 100*100*100);
	if(dice < rate){
		char newNucleotide(randomNucleotide());
		result.push_back(newNucleotide);
		insertion(rate, result);
	}
}


//// add errors in reads usinf the error profile ////
string mutateSequence(string& referenceSequence, 	unordered_map<string, double>& errorProfile){
	
	uint INSERT(0),DELETION(0),SUB(0);
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
		prevValue += (val * 100*100*100);
	}
	prevValue = 0;
	for (auto&& val : errorOrderHomo){
		if (val == errorProfile["homopolymer_ins"]){
			errorProfileHomoSorted.insert({"homopolymer_ins", val * 100*100*100 + prevValue});
		} else if (val == errorProfile["homopolymer_del"]){
			errorProfileHomoSorted.insert({"homopolymer_del", val * 100*100*100 + prevValue});
		}
		prevValue += (val * 100*100*100);
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
			for (auto it(errorProfileHomoSorted.rbegin()); it !=errorProfileHomoSorted.rend(); ++it){
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
			for (auto it(errorProfileSorted.rbegin()); it !=errorProfileSorted.rend(); ++it){
				if (dice < it->second){
					if (it->first == "mismatches"){
						//SUBSTITUTION
						SUB++;
						char newNucleotide(randomNucleotide());
						while(newNucleotide == referenceSequence[i]){
							newNucleotide = randomNucleotide();
						}
						result.push_back(newNucleotide);
						isError = true;
						break;
					} else if (it->first == "non-homopolymer_del"){
						// DELETION
						DELETION++;
						uint dice2(rand() % 100*100*100);
						while (dice2 < it->second and i < referenceSequence.size()){ // deletions larger than 1
							++i;
							dice2 = rand() % 100*100*100;
						}
						isError = true;
						break;
					} else if (it->first == "non-homopolymer_ins"){
						// INSERTION
						INSERT++;
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


//// split a string ////
vector<string> split(const string &s, char delim){
  stringstream ss(s);
  string item;
  vector<string> elems;
  while (getline(ss, item, delim)) {
    elems.push_back(move(item)); 
  }
  return elems;
}


//// read multiline faste and get one read and its header ////
pair<string, string> getRead(ifstream* readFile){
	pair<string, string> reads;
	string header, read, inter;
	char c;
	getline(*readFile,header);
	getline(*readFile,read);
	point:
	c=readFile->peek();
	if(c=='>'){
		if(read.size()>0){
			bool fail(false);
			for(uint j(0);(j)<read.size();++j){
				if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
					fail=true;
					break;
				}
			}
			if(!fail){
				reads = {header,read};
			}
		}
		read="";
	}else{
		if(!readFile->eof()){
			getline(*readFile,inter);
			read+=inter;
			goto point;
		}else{
			if(read.size()>0){
				bool fail(false);
				for(uint j(0);(j)<read.size();++j){
					if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
						fail=true;
						break;
					}
				}
				if(!fail){
					reads = {header,read};
				}
			}
		}
	}
	return reads;
}





//// get reference transcripts from fasta file derived from gtf, and expression from FS .pro file, associate header/sequence/expression ////
void generateTranscriptReferences(string& expressionFileName, string& transcriptsFileName, unordered_map <string, pair<string, uint32_t>>& transcriptSeqAndExpression){
	ifstream expression(expressionFileName);
	string expressionLine, transcript;
	uint expr;
	vector<string> transcriptExpression;
	// read the .pro file and extract transcript + expression
	while (not expression.eof()){
		getline(expression, expressionLine);
		transcriptExpression = split(expressionLine, '\t');
		if (transcriptExpression.size() > 5){
			expr = stoi(transcriptExpression[5])/10;
			transcript = transcriptExpression[1];
			transcriptSeqAndExpression.insert({transcript, {"", expr}});
		}
	}
	//~ // read the fasta file and get sequences
	ifstream sequences(transcriptsFileName);
	pair<string, string> read;
	while (not sequences.eof()){
		read = getRead(&sequences);
		transcript = split(read.first.substr(1), ' ')[0];
		if (transcriptSeqAndExpression.count(transcript)){
			transcriptSeqAndExpression[transcript].first = read.second;
		}
	}
}




//// compute final reads ////
void generateReads(uint32_t coverage, unordered_map<string, double>& errorProfile, unordered_map <string, pair<string, uint32_t>>& transcriptSeqAndExpression, const string& outFileName="simulatedReads.fa", const string& outRefFileName="simulatedPerfectSequences.fa"){
	default_random_engine generator;
	ofstream out(outFileName);
	ofstream outRef(outRefFileName);
	string read, reference;
	uint32_t noRead(0);
	for (auto it(transcriptSeqAndExpression.begin()); it != transcriptSeqAndExpression.end(); ++it){
		for (uint32_t numberOfReads(0); numberOfReads < it->second.second * coverage; ++numberOfReads){
			reference = it->second.first;
			read = mutateSequence(reference, errorProfile); // here we must add homopolymers
			staircase(read, generator);
			out <<">" << noRead << "|reference:" << it->first << endl << read << endl;
			outRef <<">" << noRead << "|reference:" << it->first << endl << reference << endl;
			++noRead;
		}
		
	}
}


//// read error profile from file ////
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
	string errorProfileFileName(""), transcriptsFileName(""), expressionFileName("");
	srand (time(NULL));
	auto startChrono = chrono::system_clock::now();
	int32_t c(0);
	while ((c = getopt (argc, argv, "c:e:t:p:")) != -1){
		switch(c){
			case 'c':
				coverage=stoi(optarg);
				break;
			case 'e':
				errorProfileFileName=optarg;
				break;
			case 't':
				transcriptsFileName=optarg;
				break;
			case 'p':
				expressionFileName=optarg;
				break;
		}
	}
	if(coverage == 0 or errorProfileFileName == "" or transcriptsFileName == "" or expressionFileName == ""){
		cout
		<<"Mandatory parameters"<<endl
		<<"-c coverage"<<endl
		<<"-e error profile file"<<endl
		<<"-t transcripts fasta file"<<endl
		<<"-p transcripts expression file"<<endl
		<< "usage: ./theReadCreator -c coverage -e ERROR_PROFILE_FILE -t TRANSCRIPTS_FILE -p EXPRESSION_PRO_FILE" << endl;
		return 0;
	}
	unordered_map<string, double> errorProfile;
	cout << "Reading error profile..." << endl;
	getErrorProfile(errorProfileFileName, errorProfile);
	unordered_map <string, pair<string, uint32_t>> transcriptSeqAndExpression;
	cout << "Generating reference transcripts and expressions..." << endl;
	generateTranscriptReferences(expressionFileName, transcriptsFileName, transcriptSeqAndExpression);
	cout << "Computing reads." << endl;
	generateReads(coverage, errorProfile, transcriptSeqAndExpression);
	auto end = chrono::system_clock::now(); auto waitedFor = end - startChrono;
	cout << "Time  in ms : " << (chrono::duration_cast<chrono::milliseconds>(waitedFor).count()) << endl;
	return 0;
}

