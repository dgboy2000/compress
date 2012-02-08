#include <iostream>
#include <vector>

using namespace std;

class DATCompression {
public:
    vector <int> compress(vector <int> data) {
        return data;
    }
    vector <int> decompress(vector <int> data) {
        return data;
    }

    int init() {}
};


int main() {
    DATCompression *dat = new DATCompression();
    
    dat->init();
    
    vector<int> data(110*110*75);
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
    
    cout << "Correctly recovered the original " << data.size() << " elements of the data" << endl;
    cout << "TopCoder Compression ratio is " << 2 * data.size() / (float) compressed.size() << endl;
    
    return 0;
}


