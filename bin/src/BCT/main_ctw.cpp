#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>  
#include <iomanip>    
#include <vector>
#include <map>

using namespace std;

// Global variables required by Functions.h
const int m = 4; // DNA alphabet size (ACGT)
const int D = 10; // max depth (matching CTW.R)
long double beta = 1.0 / m; // Default beta
const long double alpha = pow((1.0 - beta), (1.0 / (m - 1.0)));
const short k_max = 1; // Not used for CTW but required by Functions.h

// Include the original BCT functions
#include "Functions.h"

// Helper to convert DNA to 0-3
short base_to_int(char c) {
    switch(toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1; // Unknown/N
    }
}

// Helper to read FASTA
struct Chromosome {
    string name;
    vector<short> sequence;
};

vector<Chromosome> read_fasta(string filename) {
    vector<Chromosome> chromosomes;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    string line;
    Chromosome current_chr;
    current_chr.name = "";

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (current_chr.name != "") {
                chromosomes.push_back(current_chr);
            }
            current_chr.name = line.substr(1);
            // Remove description if present (keep only ID up to first space)
            size_t space_pos = current_chr.name.find(' ');
            if (space_pos != string::npos) {
                current_chr.name = current_chr.name.substr(0, space_pos);
            }
            current_chr.sequence.clear();
        } else {
            for (char c : line) {
                short val = base_to_int(c);
                if (val != -1) {
                    current_chr.sequence.push_back(val);
                }
            }
        }
    }
    if (current_chr.name != "") {
        chromosomes.push_back(current_chr);
    }
    return chromosomes;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_fasta> <output_csv> [window_size]" << endl;
        return 1;
    }

    string input_fasta = argv[1];
    string output_csv = argv[2];
    int window_size = 100000; // Default from CTW.R
    if (argc > 3) {
        window_size = atoi(argv[3]);
    }

    cout << "Reading FASTA: " << input_fasta << endl;
    vector<Chromosome> chromosomes = read_fasta(input_fasta);

    ofstream out(output_csv);
    out << "chromosome,bin_mid,bin_value" << endl;

    for (const auto& chr : chromosomes) {
        cout << "Processing " << chr.name << " (" << chr.sequence.size() << " bp)" << endl;
        
        int len = chr.sequence.size();
        int bins = max(1, (int)round((double)len / window_size));
        double bin_size_real = (double)len / bins;

        for (int i = 0; i < bins; i++) {
            int start = (int)round(i * bin_size_real);
            int end = (int)round((i + 1) * bin_size_real) - 1;
            if (end >= len) end = len - 1;
            
            int width = end - start + 1;
            if (width < D + 1) continue; // Too short for context tree

            // Extract window sequence
            vector<short> window_seq;
            window_seq.reserve(width);
            for (int k = start; k <= end; k++) {
                window_seq.push_back(chr.sequence[k]);
            }

            // Calculate CTW
            // Logic adapted from main_bct.cpp
            tree T;
            init_tree(T);
            
            // Update tree with sequence
            for (int k = D; k < window_seq.size(); k++) {
                short s = window_seq[k];
                vector<short> ct(D);
                for (int j = 0; j < D; j++) {
                    ct[j] = window_seq[k - j - 1];
                }
                update(T, s, ct);
            }

            long double log_prob = ctw(T);
            
            // Normalize (match R logic: -log2(P) / (N-D))
            // ctw(T) returns log2(P)
            // R: CTW.values = CTW.values * -log2e  <-- This converts log2 to -ln
            // But wait, let's output the raw log2 probability and let R/Python handle normalization
            // OR, replicate the exact metric here.
            // R metric: (-ln(P)) / (size - 10) * 100
            // ctw(T) is log2(P). ln(P) = log2(P) * ln(2).
            // So -ln(P) = -ctw(T) * ln(2).
            
            double log2e = log(2.0); // ln(2)
            double neg_ln_P = -log_prob * log2e;
            double metric = neg_ln_P / (width - D);
            
            // Clamp 0-1
            if (metric > 1.0) metric = 1.0;
            if (metric < 0.0) metric = 0.0;
            
            metric *= 100.0;

            double mid = start + 1 + (double)width / 2.0; // 1-based coordinates for R compatibility
            
            out << chr.name << "," << fixed << setprecision(2) << mid << "," << metric << endl;
            
            // Clean up tree memory? 
            // The original code doesn't have a clear destructor for the tree structure, 
            // and it uses 'new' extensively. 
            // For a long running process, this will leak memory.
            // We should implement a simple cleanup or just accept it for now (it's per window).
            // Actually, 'tree' is vector<vector<node*>>.
            // We need to delete nodes.
            for(auto& level : T) {
                for(auto& n : level) {
                    delete n;
                }
            }
        }
    }

    out.close();
    cout << "Done. Results saved to " << output_csv << endl;
    return 0;
}
