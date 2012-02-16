#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#define PH_MAX 16383

using namespace std;

class DATCompression {
private:
    int golomb_n; // par for golomb encoding, power of 2

    vector <int> identity_compress(vector <int> data);
    vector <int> identity_decompress(vector <int> data);

    // golomb encoding suitable for the case that the majority of data are small nonnegative numbers
    vector <int> golomb_compress(vector <int> data);
    vector <int> golomb_decompress(vector <int> data);

    vector <int> huffman_compress(vector <int> data);
    vector <int> huffman_decompress(vector <int> data);
    
    vector <int> diff_compress(vector <int> data);
    vector <int> diff_decompress(vector <int> data);

    vector <int> frame_diff_compress(vector <int> data);
    vector <int> frame_diff_decompress(vector <int> data);
    
    vector <int> compute_diffs(vector <int> data);
    
    vector <int> encode_positive(vector <int> nonpositive);
    vector <int> decode_positive(vector <int> positive);
    
    vector <int> burrows_wheeler_encode(vector <int> data);
    vector <int> burrows_wheeler_decode(vector <int> data);
public:
    vector <int> compress(vector <int> data) { return golomb_compress(data); }
    vector <int> decompress(vector <int> compressed) { return golomb_decompress(compressed); }

    vector <int> compute_diff_statistics(vector <int> data);

    int init() {
        golomb_n = 4;
    }
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

vector<int> DATCompression::huffman_compress(vector<int> data)
{
  vector<int> compressed;
  

  return compressed;
}

vector<int> DATCompression::huffman_decompress(vector<int> compressed)
{
  vector<int> decompressed;

  return decompressed;
}

vector<int> DATCompression::golomb_compress(vector<int> data)
{
  vector<int> compressed;
  int q, r, temp = 0, t = 0;

  for(int i = 0; i < data.size(); i++)
  {
    if(data[i] < 0)
    {
      cout << "Error: golomb_compress do not accept negative data" << endl;
      exit(0);
    }
    q = data[i] / golomb_n, r = data[i] % golomb_n;
    for(int j = 0; j < q; j++)
    {
      temp <<= 1;
      temp |= 0x1;
      t++;
      if(t >= 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }
    temp <<= 1;
    t++;
    if(t >= 8)
    {
      compressed.push_back(temp);
      temp = 0, t = 0;
    }
    for(int j = 0; j < int(log(golomb_n) / log(2.0) + 1e-8); j++)
    {
      temp <<= 1;
      if(r % 2)
        temp |= 0x1;
      t++;
      if(t >= 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
      r /= 2;
    }
  }
  
  if(t)
  {
    for(int i = 0; i < 8 - t; i++)
    {
      temp <<= 1;
      temp |= 0x1;
    }
    compressed.push_back(temp);
  }

  return compressed;
}

vector<int> DATCompression::golomb_decompress(vector<int> compressed)
{
  vector<int> decompressed;
  int temp, q = 0, r = 0, power = 1, nr = 0;
  char status = 'q';

  for(int i = 0; i < compressed.size(); i++)
  {
    temp = compressed[i];
    for(int j = 0; j < 8; j++)
    {
      if(status == 'q')
      {
        if((temp & 0x80) == 0x80)
          q++;
        else
          status = 'r', r = 0, power = 1, nr = 0;
      }
      else
      {
        if((temp & 0x80) == 0x80)
          r += power;
        nr++, power *= 2;
        if(nr == int(log(golomb_n) / log(2.0) + 1e-8))
        {
          status = 'q';
          decompressed.push_back(q * golomb_n + r);
          q = 0;
        }
      }
      temp <<= 1;
    }
  }

  return decompressed;
}

vector<int> DATCompression::frame_diff_decompress(vector<int> compressed)
{
  vector<int> decompressed;
  int x = compressed[0], y = compressed[1], l = compressed[2], q, r;

  decompressed.push_back(x);
  decompressed.push_back(y);
  decompressed.push_back(l);

  for(int x_ind = 0; x_ind < x; x_ind++)
    for(int y_ind = 0; y_ind < y; y_ind++)
    {
      q = compressed[3 + 2 * (x_ind * y + y_ind)];
      r = compressed[3 + 2 * (x_ind * y + y_ind) + 1];
      decompressed.push_back(q * 256 + r);
    }

  return decompressed;
}

vector<int> DATCompression::frame_diff_compress(vector<int> data)
{
  vector<int> compressed;
  double total;
  int x = data[0], y = data[1], l = data[2], d;

  compressed.push_back(x);
  compressed.push_back(y);
  compressed.push_back(l);

  // for the first frame, l_ind = 0
  for(int x_ind = 0; x_ind < x; x_ind++)
    for(int y_ind = 0; y_ind < y; y_ind++)
    {
      compressed.push_back(data[3 + l * (x_ind * y + y_ind)] / 256);
      compressed.push_back(data[3 + l * (x_ind * y + y_ind)] % 256);
    }

  // for the following frames
  for(int l_ind = 0; l_ind < l - 1; l_ind++)
  {
    // calculate the different between the l_ind frame and (l_ind + 1) frame
    total = 0.0;
    for(int x_ind = 0; x_ind < x; x_ind++)
      for(int y_ind = 0; y_ind < y; y_ind++)
        total += data[3 + l * (x_ind * y + y_ind) + l_ind] - data[3 + l * (x_ind * y + y_ind) + l_ind + 1];
    total /= (x * y);
    total = floor(total + 0.5 + 1e-8);

    if(total < 0)
      compressed.push_back(-2 * total);
    else if(total > 0)
      compressed.push_back(2 * total - 1);
    else if(!total)
      compressed.push_back(0);
    cout << total << " ";
    for(int x_ind = 0; x_ind < x; x_ind++)
      for(int y_ind = 0; y_ind < y; y_ind++)
      {
        d = data[3 + l * (x_ind * y + y_ind) + l_ind] - data[3 + l * (x_ind * y + y_ind) + l_ind + 1] - total;
        if(d < 0)
          d = -2 * d;
        else if(d > 0)
          d = 2 * d - 1;
        compressed.push_back(d);
        cout << d << " ";
      }
    cout << endl;
  }

  return compressed;
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
vector<int> DATCompression::compute_diffs(vector<int> data) {
    vector<int> diffs;
    
    int x = data[0];
    int y = data[1];
    int L = data[2];
    
    vector<int>::iterator viter = data.begin();
    viter++;
    viter++;
    viter++;
    
    for (int x_ind=0; x_ind < x; ++x_ind) {
        for (int y_ind=0; y_ind < y; ++y_ind) {
            int val = data[3 + L * (x_ind*y + y_ind)];
                        
            for (int L_ind=1; L_ind < L; ++L_ind) {
                int diff = data[3 + L * (x_ind*y + y_ind) + L_ind] - val;
                val = diff + val;
                diffs.push_back(diff);
            }
        }
    }
    
    return diffs;
}


vector<int> DATCompression::encode_positive(vector<int> nonpositive) {
    vector<int> positive;
    positive.reserve(nonpositive.size());
    
    for (vector<int>::iterator viter = nonpositive.begin(); viter != nonpositive.end(); ++viter) {
        int val = *viter;
        if (val == 0) positive.push_back(0);
        else if (val > 0) positive.push_back(2*val - 1);
        else positive.push_back(-2*val);
    }
    return positive;
}

vector<int> DATCompression::decode_positive(vector<int> positive) {
    vector<int> nonpositive;
    nonpositive.reserve(positive.size());
    
    for (vector<int>::iterator viter = positive.begin(); viter != positive.end(); ++viter) {
        int val = *viter;
        if (val % 2 == 1) nonpositive.push_back((val+1) / 2);
        else nonpositive.push_back(-val/2);
    }
    return nonpositive;
}

// Computes the BWT of the specified data. The returned vector contains as it's first element
// the 0-based index in the BWT matrix of the first row that contains the original string
vector<int> DATCompression::burrows_wheeler_encode(vector<int> nonpositive) {
    
}


vector<int> DATCompression::compute_diff_statistics(vector<int> data) {
    vector<int> diffs = compute_diffs(data);
    map<int, int> diff_counts;

    for (vector<int>::iterator viter = diffs.begin(); viter < diffs.end(); ++viter) {
        int val = *viter;
        if (diff_counts.count(val) == 0) diff_counts[val] = 0;
        ++diff_counts[val];
    }
    
    cout << "Counts of each difference in the data:" << endl;
    for (map<int, int>::iterator miter = diff_counts.begin(); miter != diff_counts.end(); ++miter) {
        printf("%d: %d\n", (*miter).first, (*miter).second);
    }
    
    return diffs;
}


class SuffixArray {
    inline bool leq(int a1, int a2, int b1, int b2) // lexicographic order
        {return(a1 < b1 || a1 == b1 && a2 <= b2); } // for pairs
    inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3)
        {return(a1 < b1 || a1 == b1 && leq(a2,a3, b2,b3)); } // and triples
    
    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
    void radixPass(int* a, int* b, int* r, int n, int K);
    
    // find the suffix array SA of s[0..n-1] in {1..K}ˆn
    // require s[n]=s[n+1]=s[n+2]=0, n>=2
    void doSuffixArrayComputation(int* s, int* SA, int n, int K);
    
    void prepareSuffixArray();
    
    int n;
    int *sa, *s;
public:
    SuffixArray(vector<int> data);
    SuffixArray(int *data, int n);
    ~SuffixArray() { free(sa); free(s); }
    
    int *getSuffixArray();
};

SuffixArray::SuffixArray(vector<int> data) {
    n = data.size();
    
    s = (int *) calloc(n + 3, sizeof(int));
    int *cur_data_elt = s;
    
    vector<int>::iterator data_end = data.end();
    for (vector<int>::iterator viter=data.begin(); viter!=data_end; ++viter) {
        *cur_data_elt = *viter;
        ++cur_data_elt;
    }
    
    prepareSuffixArray();
}

SuffixArray::SuffixArray(int *data, int _n) {
    n = _n;
    
    s = (int *) calloc(n + 3, sizeof(int));
    memcpy(s, data, n);
    
    prepareSuffixArray();
}

void SuffixArray::prepareSuffixArray() {
    // Add 3 trailing zeroes for skew algorithm
    int *cur_data_elt = s + (sizeof(int)) * n;
    *cur_data_elt = 0; ++cur_data_elt;
    *cur_data_elt = 0; ++cur_data_elt;
    *cur_data_elt = 0; ++cur_data_elt;
    
    sa = (int *) calloc(n, sizeof(int));
    doSuffixArrayComputation(s, sa, n, 16);
}

int * SuffixArray::getSuffixArray() {
    return sa;
}



// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
void SuffixArray::radixPass(int* a, int* b, int* r, int n, int K)
{// count occurrences
    int* c = new int[K + 1]; // counter array
    for (int i = 0; i <= K; i++) c[i] = 0; // reset counters
    for (int i = 0; i < n; i++) c[r[a[i]]]++; // count occurrences
    for (int i = 0, sum = 0; i <= K; i++) // exclusive prefix sums
    {int t = c[i]; c[i] = sum; sum += t; }
    for (int i = 0; i < n; i++) b[c[r[a[i]]]++] = a[i]; // sort
    delete [] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}ˆn
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void SuffixArray::doSuffixArrayComputation(int* s, int* SA, int n, int K) {
    int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
    int* s12 = new int[n02 + 3]; s12[n02]= s12[n02+1]= s12[n02+2]=0;
    int* SA12 = new int[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
    int* s0 = new int[n0];
    int* SA0 = new int[n0];
// generate positions of mod 1 and mod 2 suffixes
// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    for (int i=0, j=0; i < n+(n0-n1); i++) if (i%3 != 0) s12[j++] = i;
// lsb radix sort the mod 1 and mod 2 triples
    radixPass(s12 , SA12, s+2, n02, K);
    radixPass(SA12, s12 , s+1, n02, K);
    radixPass(s12 , SA12, s , n02, K);
// find lexicographic names of triples
    int name = 0, c0 = -1, c1 = -1, c2 = -1;
    for (int i = 0; i < n02; i++) {
        if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2)
            {name++; c0 = s[SA12[i]]; c1 = s[SA12[i]+1]; c2 = s[SA12[i]+2]; }
        if (SA12[i] % 3 == 1) { s12[SA12[i]/3] = name; } // left half
        else {s12[SA12[i]/3 + n0] = name; } // right half
    }
// recurse if names are not yet unique
    if (name < n02) {
        doSuffixArrayComputation(s12, SA12, n02, name);
// store unique names in s12 using the suffix array
        for (int i = 0; i < n02; i++) s12[SA12[i]] = i + 1;
    } else // generate the suffix array of s12 directly
        for (int i = 0; i < n02; i++) SA12[s12[i] - 1] = i;
// stably sort the mod 0 suffixes from SA12 by their first character
    for (int i=0, j=0; i < n02; i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
    radixPass(s0, SA0, s, n0, K);
// merge sorted SA0 suffixes and sorted SA12 suffixes
    for (int p=0, t=n0-n1, k=0; k < n; k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
        int i = GetI(); // pos of current offset 12 suffix
        int j = SA0[p]; // pos of current offset 0 suffix
        if (SA12[t] < n0 ? // different compares for mod 1 and mod 2 suffixes
            leq(s[i], s12[SA12[t] + n0], s[j], s12[j/3]) :
            leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
        {// suffix from SA12 is smaller
            SA[k] = i; t++;
            if (t == n02) // done --- only SA0 suffixes left
                for (k++; p < n0; p++, k++) SA[k] = SA0[p];
        } else {// suffix from SA0 is smaller
            SA[k] = j; p++;
            if (p == n0) // done --- only SA12 suffixes left
                for (k++; t < n02; t++, k++) SA[k] = GetI();
        }
    }
    delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}






/*********** Testing Utilities ***************/

void test_suffix_arrays() {
    int len = 9;
    int word[] = {1, 2, 3, 4, 5, 4, 3, 2, 1};
    int true_suffixes[] = {8, 0, 7, 1, 6, 2, 5, 3, 4};
    
    SuffixArray sa(word, len);
    int *suffixes = sa.getSuffixArray();
    
    for (int i=0; i<len; ++i) {
        if (suffixes[i] == true_suffixes[i]) printf("CORRECT: True suffix SA[%d] = %d\n", i, suffixes[i]);
        else printf("ERROR: True suffix SA[%d] is %d but found %d\n", i, true_suffixes[i], suffixes[i]);
    }
}

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

void investigate_file(string filename) {
    DATCompression dat;
    dat.init();
    
    printf("Investigating file %s\n", filename.c_str());
    
    vector<int> data = read_data(filename);
    dat.compute_diff_statistics(data);
}

int main() {
    DATCompression *dat = new DATCompression();
    
    dat->init();
    
    test_suffix_arrays();
    
    //investigate_file("data/B28-39_100_100_acq_0007.tab");
    
    // test_compression_on_file("data/test.tab");

    //test_compression_on_file("data/B28-39_100_100_acq_0007.tab");
    //test_compression_on_file("data/B28-39_100_100_acq_0400.tab");
    //test_compression_on_file("data/B28-39_1600_1000_acq_0007.tab");
    //test_compression_on_file("data/B28-39_1600_1000_acq_0400.tab");
    
    return 0;
}



