#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#define PH_MAX 16383

using namespace std;

class DATCompression {
private:
    vector <int> identity_compress(vector <int> data);
    vector <int> identity_decompress(vector <int> data);
    
    vector <int> diff_compress(vector <int> data);
    vector <int> diff_decompress(vector <int> data);
    
    vector <int> compute_diffs(vector <int> data);
public:
    vector <int> compress(vector <int> data) { return diff_compress(data); }
    vector <int> decompress(vector <int> compressed) { return diff_decompress(compressed); }

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


vector<int> DATCompression::diff_compress(vector<int> data) {
    vector<int> compressed;
    int x, y, L;
    
    x = data[0];
    y = data[1];
    L = data[2];
    
    compressed.push_back(x);
    compressed.push_back(y);
    compressed.push_back(L);
    
    for (int x_ind=0; x_ind < x; ++x_ind) {
        for (int y_ind=0; y_ind < y; ++y_ind) {
            int val = data[3 + L * (x_ind*y + y_ind)];

            compressed.push_back(val / 256);
            compressed.push_back(val % 256);
            
            for (int L_ind=1; L_ind < L; ++L_ind) {
                int diff = data[3 + L * (x_ind*y + y_ind) + L_ind] - val;
                val = diff + val;
                if (abs(diff) < 128) {
                    compressed.push_back(diff + 127);
                    continue;
                }
                
                compressed.push_back(255);
                compressed.push_back(val / 256);
                compressed.push_back(val % 256);
            }
        }
    }
    
    return compressed;
}
vector<int> DATCompression::diff_decompress(vector<int> compressed) {
    vector<int> decompressed;
    int x, y, L;
    
    x = compressed[0];
    y = compressed[1];
    L = compressed[2];
    
    decompressed.push_back(x);
    decompressed.push_back(y);
    decompressed.push_back(L);
    
    vector<int>::iterator viter = compressed.begin();
    viter++;
    viter++;
    viter++;
    
    for (int x_ind=0; x_ind < x; ++x_ind) {
        for (int y_ind=0; y_ind < y; ++y_ind) {
            int val = *(viter++);
            val = 256 * val + *(viter++);

            decompressed.push_back(val);
            
            for (int L_ind=1; L_ind < L; ++L_ind) {
                int diff = *(viter++);
                if (diff < 255) {
                    val = val + diff - 127;
                    decompressed.push_back(val);
                    continue;
                }
                
                val = *(viter++);
                val = 256 * val + *(viter++);
                decompressed.push_back(val);
            }
        }
    }
    
    if (viter != compressed.end()) {
        printf("ERROR: finished decompressing before the end of the compressed string!\n");
        exit(1);
    }
    
    return decompressed;
}










/*********** Compression Utilities ***************/
// vector<int> DATCompression::compute_diffs(vector<int> data) {
//     vector<
// }








/*********** Testing Utilities ***************/

vector<int> read_data(string filename) {
    ifstream data_file;
    string data_line;
    vector<int> data;
    
    data_file.open(filename.c_str());
    
    int x, y, x_ind, y_ind, L, valCnt, lineCnt;
    string pHval;
    L = 0;
    lineCnt = 0;
    while (getline(data_file, data_line)) {
        lineCnt++;
        stringstream line_stream(data_line);
        
        string x_str, y_str;
        getline(line_stream, x_str, '\t');
        getline(line_stream, y_str, '\t');
        x_ind = atoi(x_str.c_str());
        y_ind = atoi(y_str.c_str());
        x = max(x, x_ind + 1);
        y = max(y, y_ind + 1);
        
        valCnt = 0;
        while (getline(line_stream, pHval, '\t')) {
            int val = atoi(pHval.c_str());
            data.push_back(val);
            valCnt++;
        }
        if (L > 0 && L != valCnt) {
            printf("ERROR: wrong number of values on this line! Got %d values, expected %d values on line %d\n", valCnt, L, lineCnt);
            exit(1);
        }
        L = valCnt;
    }
    data_file.close();
    
    if (data.size() != x * y  * L) {
        printf("ERROR: read in wrong amount of data! Expected %d*%d*%d = %d pH values, but read %ld\n", x, y, L, x*y*L, data.size());
        exit(1);
    }
    
    data.insert(data.begin(), L);
    data.insert(data.begin(), y);
    data.insert(data.begin(), x);    
    
    printf("Read in %ld pH values with (x, y, L) = (%d, %d, %d)\n", data.size() - 3, x, y, L);
    
    return data;
}

void test_compression_on_file(string filename) {
    DATCompression dat;
    dat.init();
    
    vector<int> data = read_data(filename);
    for (int i=0; i<data.size(); ++i) {
        if (data[i] < 0 || data[i] > PH_MAX) {
            printf("ERROR: initial data value is outside of [0, %d] range!\n", PH_MAX);
            exit(1);
        }
    }
        
    vector<int> compressed = dat.compress(data);
    for (int i=0; i<compressed.size(); ++i) {
        if (compressed[i] < 0 || compressed[i] > 255) {
            cout << "ERROR: compressed value is outside of [0, 255] range!" << endl;
            exit(1);
        }
    }
    cout << "Verified that all " << compressed.size() << " compressed ints are in within [0, 255] legal range" << endl;
    
    vector<int> decompressed = dat.decompress(compressed);
    
    if (data.size() != decompressed.size()) {
        printf("ERROR: decompressed %ld elements but expected %ld\n", decompressed.size(), data.size());
        exit(1);
    }
    for (int i=0; i<data.size(); ++i) {
        if (data[i] != decompressed[i]) {
            cout << "ERROR: decompression doesn't equal compression!" << endl;
            exit(1);
        }
    }
    
    printf("Correctly recovered the original %ld elements of the data\n", data.size());
    printf("TopCoder Compression ratio is %f\n", 2 * data.size() / (float) compressed.size());
}

int main() {
    DATCompression *dat = new DATCompression();
    
    dat->init();
    
    // vector<int> data(110*110*75);
    test_compression_on_file("data/B28-39_100_100_acq_0007.tab");
    test_compression_on_file("data/B28-39_100_100_acq_0400.tab");
    test_compression_on_file("data/B28-39_1600_1000_acq_0007.tab");
    test_compression_on_file("data/B28-39_1600_1000_acq_0400.tab");
    
    return 0;
}



