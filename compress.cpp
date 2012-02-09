#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

class DATCompression {
private:
    vector <int> identity_compress(vector <int> data);
    vector <int> identity_decompress(vector <int> data);
    
public:
    vector <int> compress(vector <int> data) {
        return identity_compress(data);
    }
    vector <int> decompress(vector <int> compressed) {
        return identity_decompress(compressed);
    }

    int init() {}
};


vector<int> DATCompression::identity_compress(vector<int> data) {
    vector<int> compressed;
    compressed.push_back(data[0]);
    compressed.push_back(data[1]);
    compressed.push_back(data[2]);
    
    vector<int>::iterator viter = data.begin();
    viter++;
    viter++;
    viter++;
    for (viter; viter != data.end(); ++viter) {
        int val = *viter;
        
        compressed.push_back(val / 256);
        compressed.push_back(val % 256);
    }
    
    return compressed;
}
vector<int> DATCompression::identity_decompress(vector<int> compressed) {
    vector<int> decompressed;
    decompressed.push_back(compressed[0]);
    decompressed.push_back(compressed[1]);
    decompressed.push_back(compressed[2]);
    
    vector<int>::iterator viter = compressed.begin();
    viter++;
    viter++;
    viter++;
    for (viter; viter != compressed.end(); ++viter) {
        int val = *viter;
        val = 256 * val + *(++viter);
        
        decompressed.push_back(val);
    }
    
    return decompressed;
}










vector<int> read_data(string filename) {
    ifstream data_file;
    string data_line;
    vector<int> data;
    
    data_file.open(filename.c_str());
    
    int x, y, L, valCnt, lineCnt;
    string pHval;
    L = 0;
    lineCnt = 0;
    while (getline(data_file, data_line)) {
        lineCnt++;
        stringstream line_stream(data_line);
        line_stream >> x >> y;
        valCnt = 0;
        while (getline(line_stream, pHval, '\t')) {
            data.push_back(atoi(pHval.c_str()));
            valCnt++;
        }
        if (L > 0 && L != valCnt) {
            printf("ERROR: wrong number of values on this line! Got %d values, expected %d values on line %d\n", valCnt, L, lineCnt);
            exit(1);
        }
        L = valCnt;
    }
    data_file.close();
    
    data.insert(data.begin(), L);
    data.insert(data.begin(), y);
    data.insert(data.begin(), x);    
    
    printf("Read in %ld pH values with (x, y, L) = (%d, %d, %d)\n", data.size() - 3, x, y, L);
    
    return data;
}

int main() {
    DATCompression *dat = new DATCompression();
    
    dat->init();
    
    // vector<int> data(110*110*75);
    vector<int> data = read_data("data/B28-39_100_100_acq_0007.tab");
    for (int i=0; i<data.size(); ++i) {
        if (data[i] < 0 || data[i] > 16383) {
            cout << "ERROR: initial data value is outside of [0, 16383] range!" << endl;
            return 1;
        }
    }
        
    vector<int> compressed = dat->compress(data);
    for (int i=0; i<compressed.size(); ++i) {
        if (compressed[i] < 0 || compressed[i] > 255) {
            cout << "ERROR: compressed value is outside of [0, 255] range!" << endl;
            return 1;
        }
    }
    cout << "Verified that all " << compressed.size() << " compressed ints are in within [0, 255] legal range" << endl;
    
    vector<int> decompressed = dat->decompress(compressed);
    
    if (data.size() != decompressed.size()) {
        cout << "ERROR: decompressed into the wrong size vector!" << endl;
        return 1;
    }
    for (int i=0; i<data.size(); ++i) {
        if (data[i] != decompressed[i]) {
            cout << "ERROR: decompression doesn't equal compression!" << endl;
            return 1;
        }
    }
    
    printf("Correctly recovered the original %ld elements of the data\n", data.size());
    printf("TopCoder Compression ratio is %f\n", 2 * data.size() / (float) compressed.size());
    
    return 0;
}


