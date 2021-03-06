#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace std;

namespace synth {


std::default_random_engine generator;
std::uniform_int_distribution<int> uint_dist(0, RAND_MAX);

double doubleRand(double mx) {
	//OLD RANDOM
//   return double(rand()) / (double(RAND_MAX) + 1.0) * mx;
   return double(uint_dist(generator)) / (double(RAND_MAX) + 1.0) * mx;
}

int intRand(int mx) {
	//OLD RANDOM
//   return (rand() * 1.0) / (RAND_MAX) * mx;
	return (uint_dist(generator) * 1.0) / (RAND_MAX) * mx;
}


const int MAXN = 10000;
const int MAXLEN = 2000;
const int MAX_SAMPLE = 1000;

int n;
vector<double> adv, advCumSum;
vector<int> parent, parentDistance;
vector<vector<int>> sequence;

//double adv[MAXN], advCumSum[MAXN];
//int seq[MAXN][MAXLEN];
//int parent[MAXN], parentDistance[MAXN];

int locusCount, step;
// , potentialLocusCount;
double incAdvProb, decAdvProb, keepAdvProb;
double advIncStep, advDecStep;
double missingValueRate, zeroToOneRate, oneToZeroRate;

double stepMutationRate;

int sampleWithAdvantage() {
	double prolR = doubleRand(advCumSum[n-1]);
	vector<double>::iterator prolIt = lower_bound(advCumSum.begin(), advCumSum.end(), prolR);
	int prolIdx = prolIt - advCumSum.begin();
	return prolIdx;
}

void simulate() {
	assert(step < MAXN);
	//everything is not on vectors, so no more need to check the sizes here
	//assert(locusCount < MAXLEN);
	//init:
	sequence.clear();
	sequence.push_back(vector<int>(locusCount, 0));
	// for (int i=0; i<locusCount; i++)
	// 	seq.push_back(0);
		// seq[0][i] = 0;
	adv.clear();
	adv.push_back(1);
	// adv[0] = 1;
	advCumSum.clear();
	advCumSum.push_back(adv[0]);
	// advCumSum[0] = adv[0];
	parent.clear();
	parent.push_back(-1);
	// parent[0] = -1;
	parentDistance.clear();
	parentDistance.push_back(0);
	n = 1;

	std::poisson_distribution<int> distribution(stepMutationRate);

	
	for (int s = 0; s < step; s++) {
		int prolIdx = sampleWithAdvantage();
		// cerr << "  s=" << s << " prolIdx=" << prolIdx << endl;

		assert(prolIdx < n);
		assert(0 <= prolIdx);

		// copy(seq[prolIdx], seq[prolIdx] + locusCount, seq[n]);
		sequence.push_back(sequence[prolIdx]);
		int cnt = distribution(generator);
		for (int j=0; j<cnt; j++) {
			int l = intRand(locusCount);
			// cerr << "  " << l << " " << locusCount << endl;
			sequence[n][l] = 1-sequence[n][l];
		}

		double advR = doubleRand(1);
		if (advR < incAdvProb) 
			adv.push_back(adv[prolIdx] * (1 + advIncStep));
			// adv[n] = adv[prolIdx] * (1 + advIncStep);
		else if (advR < decAdvProb)
			adv.push_back(adv[prolIdx] * (1 - advDecStep));
			// adv[n] = adv[prolIdx] * (1 - advDecStep);
		else
			adv.push_back(adv[prolIdx]);
			// adv[n] = adv[prolIdx];
		// advCumSum[n] = advCumSum[n-1] + adv[n];
		advCumSum.push_back(advCumSum[n-1] + adv[n]);

		// parent[n] = prolIdx;
		parent.push_back(prolIdx);
		// parentDistance[n] = 1;
		parentDistance.push_back(1);
		n++;
	}

	// cerr << "All Sequences: " << endl;
	// for (int i=0; i<n; i++) {
	// 	for (int j=0; j<locusCount; j++) {
	// 		cerr << seq[i][j];
	// 	}
	// 	cerr << endl;
	// }

	// cerr << "Simulation done" << endl;
}


struct Output {
    int sampleCount;
	// vector<int> cloneSamples[MAXN];
	// int sampleClone[MAX_SAMPLE];
	// int output[MAXLEN][MAX_SAMPLE];
	// int parentCompressed[MAXN],
	// 	parentCompressedDistance[MAXN];
	vector<vector<int>> cloneSamples;
	vector<int> sampleClone;
	vector<vector<int>> output;
	vector<int> parentCompressed,
		parentCompressedDistance;


    Output(int _sampleCount) : sampleCount(_sampleCount) {
		cloneSamples = vector<vector<int>>(n, vector<int>());
        // for (int i=0; i<n; i++) {
        //     cloneSamples[i].clear();
        // }

		parentCompressedDistance = vector<int>(n, -1);
	// for (int i=0; i<n; i++) {
	// 	o.parentCompressedDistance[i] = -1;
	// }
		parentCompressed = vector<int>(n,-1);

    }

	void compressParent(int v) {
		if (parentCompressedDistance[v] != -1)
			return;
		if (parent[v] == -1) {
			parentCompressed[v] = -1;
			parentCompressedDistance[v] = 0;
			return;
		}
		compressParent(parent[v]);
		if (parent[parent[v]] != -1 && cloneSamples[parent[v]].size() == 0) {
			parentCompressed[v] = parentCompressed[parent[v]];
			parentCompressedDistance[v] = parentCompressedDistance[parent[v]] + 1;
		} else {
			parentCompressed[v] = parent[v];
			parentCompressedDistance[v] = 1;
		}
	}

};


Output sample(int sampleCount) {
	assert(sampleCount < MAX_SAMPLE);

	Output o(sampleCount);

	for (int i=0; i<sampleCount; i++) {
		// double prolR = doubleRand(advCumSum[n-1]);
		// double* prolIt = lower_bound(advCumSum, advCumSum+n, prolR);
		// int prolIdx = prolIt - advCumSum;

		int prolIdx = sampleWithAdvantage();

		o.cloneSamples[prolIdx].push_back(i);
		o.sampleClone.push_back(prolIdx);
		// o.sampleClone[i] = prolIdx;
	}


	std::uniform_real_distribution<double> distribution(0.0,1.0);

	o.output = vector<vector<int>>(locusCount, vector<int>(sampleCount, -1));

	for (int i=0; i<sampleCount; i++) {
		for (int j=0; j<locusCount; j++) {
			int v = 3;
			if (distribution(generator) < missingValueRate) {
				v = 3;
			} else {
				int lv = sequence[o.sampleClone[i]][j];
				if (lv == 0) {
					if (distribution(generator) < zeroToOneRate)
						lv = 1;
				} else {
					if (distribution(generator) < oneToZeroRate)
						lv = 0;
				}
				v = lv;
			}
			o.output[j][i] = v;
		}
	}

	for (int i=0; i<n; i++) {
		o.compressParent(i);
	}
	return o;
}

void printSample(const Output& o, string cloneFileName, string seqFileName, string treeFileName, string trueSeqFileName) {

	ofstream fclone(cloneFileName);
	for (int i=0; i<n; i++) {
		fclone << i << " ";
		for (auto j: o.cloneSamples[i]) {
			fclone << j+1 << " ";
		}
		fclone << endl;
	}

	ofstream fseq(seqFileName);
	for (int j=0; j<locusCount; j++) {
		for (int i=0; i<o.sampleCount; i++) {
			fseq << o.output[j][i] << " ";
		}
		fseq << endl;
	}

	ofstream ftree(treeFileName);
	for (int i=0; i<n; i++) {
		if (parent[i] == -1 || o.cloneSamples[i].size() > 0)
			ftree << i << " ";
	}

	ftree << endl;
	for (int i=0; i<n; i++) {
		if (parent[i] != -1 && o.cloneSamples[i].size() > 0)
			ftree << i << "->" << o.parentCompressed[i] << " " << o.parentCompressedDistance[i] << endl;
		// cerr << i << " " << parent[i] << " " << adv[i] << " " << advCumSum[i] << endl;
	}

	if (trueSeqFileName != "") {
		ofstream fseq(trueSeqFileName);
		for (int j=0; j<locusCount; j++) {
			for (int i=0; i<o.sampleCount; i++) {
				fseq << sequence[o.sampleClone[i]][j] << " ";
			}
			fseq << endl;
		}
	}

}


}

