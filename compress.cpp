#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <utility>
#include <vector>

#define PH_MAX 16383
#define HUFFMAN_MAX 16
#define RUN_MIN 4
#define RUN_MAX 259

using namespace std;


void test_compression_on_file(string filename);


class DATCompression {
private:
    int golomb_n; // par for golomb encoding, power of 2

    vector <int> identity_compress(vector <int> data);
    vector <int> identity_decompress(vector <int> data);

    vector <int> compression(vector <int> data);
    vector <int> decompression(vector <int> data);

    // golomb encoding suitable for the case that the majority of data are small nonnegative numbers
    vector <int> golomb_compress(vector <int> data);
    vector <int> golomb_decompress(vector <int> data);

    // the input data of huffman encoding must be in range of [0..255]
    vector <int> huffman_compress(vector <int> data);
    vector <int> huffman_decompress(vector <int> data);
    
    vector <int> diff_compress(vector <int> data);
    vector <int> diff_decompress(vector <int> data);

    vector <int> frame_diff_compress(vector <int> data);
    vector <int> frame_diff_decompress(vector <int> data);
    
    vector <int> compute_diffs(vector <int> data);
    
    vector <int> encode_positive(vector <int> nonpositive);
    vector <int> decode_positive(vector <int> positive);
public:
    int *burrows_wheeler_encode(int *data, int n);
    int *burrows_wheeler_decode(int *data, int n);
    
    vector<int> burrows_wheeler_encode(vector<int> data);
    vector<int> burrows_wheeler_decode(vector<int> data);
    
    int *run_length_encode(int *data, int n);
    int *run_length_decode(int *data);
    
    vector <int> run_length_encode(vector <int> data);
    vector <int> run_length_decode(vector <int> data);

    vector <int> compress(vector <int> data) { return compression(data); }
    vector <int> decompress(vector <int> compressed) { return decompression(compressed); }

    vector <int> compute_diff_statistics(vector <int> data);

    int init() {
        golomb_n = 4;
    }
};


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
    SuffixArray(vector<int> *data);
    SuffixArray(int *data, int n);
    ~SuffixArray() { free(sa); free(s); }
    
    inline int operator[](int index) { return sa[index]; }
    inline long int size() { return n; }
    
    void print();
    
    int *getSuffixArray();
    int *getString() { return s; }
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

vector <int> DATCompression::compression(vector<int> data)
{
    vector <int> compressed = data;
    int i, n=30;
  
    // // Simple compression stack
    // compressed = diff_compress(data);
    // compressed = huffman_compress(compressed);

    // bzip-inspired compression stack (but simpler)
    compressed = diff_compress(compressed);
    compressed = run_length_encode(compressed);
    compressed = burrows_wheeler_encode(compressed);
    compressed = run_length_encode(compressed);
    compressed = huffman_compress(compressed);

    return compressed;
}

vector <int> DATCompression::decompression(vector<int> compressed)
{
  vector <int> decompressed = compressed;
  int i,n=30;
  
  // // Simple decompression stack
  // decompressed = huffman_decompress(decompressed);
  // decompressed = diff_decompress(decompressed);

  // bzip-inspired decompression stack (but simpler)
  decompressed = huffman_decompress(decompressed);
  decompressed = run_length_decode(decompressed);
  decompressed = burrows_wheeler_decode(decompressed);
  decompressed = run_length_decode(decompressed);
  decompressed = diff_decompress(decompressed);

  return decompressed;
}

vector<int> DATCompression::huffman_compress(vector<int> data)
{
  vector<int> compressed;
  int freq[HUFFMAN_MAX], in_set[HUFFMAN_MAX], masked[HUFFMAN_MAX];
  int max_freq, second_max_freq, max_ind, second_max_ind, s, temp = 0, code_size, t = 0, q;
  vector<int> code[HUFFMAN_MAX];
  
  memset(freq, 0, sizeof(freq));
  memset(masked, 0, sizeof(masked));
  for(int i = 0; i < data.size(); i++)
  {
    if(data[i] < 0 || data[i] >= HUFFMAN_MAX * HUFFMAN_MAX)
    {
      cout << "Error: huffman_compress do not accept data out of range [0..255]" << endl;
      exit(0);
    }
    freq[data[i] / HUFFMAN_MAX]++;
    freq[data[i] % HUFFMAN_MAX]++;
  }

  for(int i = 0; i < HUFFMAN_MAX; i++)
  {
    code[i].clear();
    in_set[i] = i;
  }

  // generate code dictionary
  for(int i = 0; i < HUFFMAN_MAX - 1; i++)
  {
    // find out the minimum and second minimum freq
    max_freq = -1, second_max_freq = -1;
    for(int j = 0; j < HUFFMAN_MAX; j++)
      if(!masked[j])
      {
        if(freq[j] < max_freq || max_freq == -1)
          second_max_freq = max_freq, second_max_ind = max_ind, max_freq = freq[j], max_ind = j;
        else if(freq[j] < second_max_freq || second_max_freq == -1)
          second_max_freq = freq[j], second_max_ind = j;
      }

    // merge them
    freq[max_ind] += freq[second_max_ind];
    masked[second_max_ind] = 1;

    for(int j = 0; j < HUFFMAN_MAX; j++)
    {
      if(in_set[j] == in_set[second_max_ind])
        code[j].push_back(1);
      else if(in_set[j] == in_set[max_ind])
        code[j].push_back(0);
    }
    
    s = in_set[second_max_ind];
    for(int j = 0; j < HUFFMAN_MAX; j++)
      if(in_set[j] == s)
        in_set[j] = in_set[max_ind];
    in_set[second_max_ind] = in_set[max_ind];

  }

  for(int i = 0; i < HUFFMAN_MAX; i++)
  {
    code_size = code[i].size();
    for(int j = 0; j < 4; j++)
    {
      temp <<= 1;
      if(code_size % 2)
        temp |= 1;
      code_size /= 2, t++;
      if(t == 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }

    for(int j = code[i].size() - 1; j >= 0; j--)
    {
      temp <<= 1;
      if(code[i][j])
        temp |= 1;
      t++;
      if(t == 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }
  }

  /*for(int i = 0; i < HUFFMAN_MAX; i++)
  {
    for(int j = code[i].size() - 1; j >= 0; j--)
      cout << code[i][j] << " ";
    cout << endl;
  }*/

  for(int i = 0; i < data.size(); i++)
  {
    for(int k = 0; k < 2; k++)
    {
      if(k == 0)
        q = data[i] / HUFFMAN_MAX;
      else
        q = data[i] % HUFFMAN_MAX;
      for(int j = code[q].size() - 1; j >= 0; j--)
      {
        temp <<= 1;
        if(code[q][j])
          temp |= 1;
        t++;
        if(t == 8)
        {
          compressed.push_back(temp);
          temp = 0, t = 0;
        }
      }
    }
  }

  if(t)
  {
    for(int i = 0; i < 8 - t; i++)
      temp <<= 1;
    compressed.push_back(temp);
  }

  if(!t)
    compressed.push_back(0);
  else
    compressed.push_back(8 - t);

  return compressed;
}

vector<int> DATCompression::huffman_decompress(vector<int> compressed)
{
  vector<int> decompressed;
  vector<int> code[HUFFMAN_MAX];
  int code_space[65536], position;
  int number, d = 0, power = 1, code_ind = 0, c, temp = 0, q, r, g = 0;
  int useless_bits = compressed[compressed.size() - 1];
  char status = 'd';

  memset(code_space, -1, sizeof(code_space));
  for(int i = 0; i < compressed.size() - 1; i++)
  {
    number = compressed[i];
    for(int j = 0; j < 8; j++)
    {
      if(status == 'd')
      {
        if((number & 0x80) == 0x80)
          temp += power;
        number <<= 1;
        power *= 2, d++;
        if(d == 4)
          status = 'c', c = temp;
      }
      else if(status == 'c')
      {
        if((number & 0x80) == 0x80)
          code[code_ind].push_back(1);
        else
          code[code_ind].push_back(0);
        number <<= 1;
        c--;
        if(!c)
          status = 'd', code_ind++, d = 0, power = 1, temp = 0;
        if(code_ind == HUFFMAN_MAX)
        {
          status = 't';
          // generate code space
          for(int k = 0; k < HUFFMAN_MAX; k++)
          {
            position = 0;
            for(int l = 0; l < code[k].size(); l++)
            {
              if(!code[k][l])
                position = position * 2 + 1;
              else
                position = position * 2 + 2;
            }
            code_space[position] = k;
          }
          position = 0;
        }
      }
      else if(status == 't')
      {
        if(i == compressed.size() - 2 && 8 - j == useless_bits)
          break;

        if((number & 0x80) == 0x80)
          position = position * 2 + 2;
        else
          position = position * 2 + 1;
        if(code_space[position] >= 0)
        {
          if(!g)
            q = code_space[position], g = 1;
          else
          {
            r = code_space[position], g = 0;
            decompressed.push_back(q * HUFFMAN_MAX + r);
          }
          position = 0;
        }
        number <<= 1;
      }
    }
  }

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

// Returns a vector where every run of length RUN_MIN or higher has been
// encoded as A{RUN_MIN}(run_count - RUN_MIN). The first entry is the length
// of the original uncompressed vector
int *DATCompression::run_length_encode(int *data, int n) {
    vector<int> encoded;
    int *rle;
    encoded.reserve(n); // Heuristic; reserve enough space for the data
    
    int i=0;
    while (i < n) {
        int run_int = data[i];
        int run_count = 1;
        ++i;
        while (i < n && data[i] == run_int && run_count < RUN_MAX) {
            ++run_count;
            ++i;
        }
        for (int j=0; j<run_count && j<RUN_MIN; ++j) encoded.push_back(run_int);
        if (run_count >= RUN_MIN) encoded.push_back(run_count - RUN_MIN);
    }
    
    rle = (int *) calloc(encoded.size()+1, sizeof(encoded));
    rle[0] = n;
    copy(encoded.begin(), encoded.end(), rle+1);
    
    return rle;
}
// Decode a vector in the format produced by run_length_encode
int *DATCompression::run_length_decode(int *encoded) {
    int n = encoded[0];
    int *decoded = (int *) calloc(n, sizeof(int));

    int cur_ind = 1;
    int total_count = 0;
    while (total_count < n) {
        int run_int = encoded[cur_ind++];
        int run_count = 1;
        while (encoded[cur_ind] == run_int && run_count < RUN_MIN) {
            ++cur_ind;
            ++run_count;
        }
        if (run_count == RUN_MIN) run_count += encoded[cur_ind++];
        for (int i=0; i<run_count; ++i) decoded[total_count++] = run_int;
    }
    
    return decoded;
}

// Returns a vector where every run of length RUN_MIN or higher has been
// encoded as A{RUN_MIN}(run_count - RUN_MIN).
vector<int> DATCompression::run_length_encode(vector<int> data) {
    int n = data.size();
    vector<int> encoded;
    encoded.reserve(n); // Heuristic; reserve enough space for the data
    
    int i=0;
    while (i < n) {
        int run_int = data[i];
        int run_count = 1;
        ++i;
        while (i < n && data[i] == run_int && run_count < RUN_MAX) {
            ++run_count;
            ++i;
        }
        for (int j=0; j<run_count && j<RUN_MIN; ++j) encoded.push_back(run_int);
        if (run_count >= RUN_MIN) encoded.push_back(run_count - RUN_MIN);
    }
    
    return encoded;
}
// Decode a vector in the format produced by run_length_encode
vector<int> DATCompression::run_length_decode(vector<int> encoded) {
    int n = encoded.size();
    vector<int> decoded;
    decoded.reserve(n); // Heuristic; reserve some space for the data

    vector<int>::iterator viter = encoded.begin(), data_end = encoded.end();
    while (viter != data_end) {
        int run_int = *(viter++);
        int run_count = 1;
        while (*viter == run_int && run_count < RUN_MIN && viter != data_end) {
            ++viter;
            ++run_count;
        }
        if (run_count == RUN_MIN) run_count += *(viter++);
        for (int i=0; i<run_count; ++i) decoded.push_back(run_int);
    }
    
    return decoded;
}

// Computes the BWT of the specified data. The returned vector contains as it's last 3 elements
// the 0-based index in the BWT matrix of the first row that contains the original string.
// This index can be very large, so it is represented as a 3 byte integer.
// The corresponding '$' character has been replaced with int min in the bwt sequence that follows.
// The return vector is therefore 4 longer than the input vector.
int *DATCompression::burrows_wheeler_encode(int *data, int n) {
    SuffixArray sa(data, n);
    
    int *bwt = (int *) calloc(n+4, sizeof(int));
    long int dollar_ind;
    bwt[0] = data[n-1];
    
    for (int i=0; i<n; ++i) {
        if (sa[i] == 0) {
            bwt[i+1] = INT_MIN;
            dollar_ind = i + 1;
        }
        else bwt[i+1] = data[sa[i] - 1];
    }
    
    bwt[n+1] = dollar_ind / (256*256);
    bwt[n+2] = (dollar_ind / 256) % 256;
    bwt[n+3] = dollar_ind % 256;
    
    return bwt;
}

// Computes the BWT of the specified data. The returned vector contains as it's last 3 elements
// the 0-based index in the BWT matrix of the first row that contains the original string.
// This index can be very large, so it is represented as a 3 byte integer.
// The corresponding '$' character has been replaced with int min in the bwt sequence that follows.
// The return vector is therefore 4 longer than the input vector.
vector<int> DATCompression::burrows_wheeler_encode(vector<int> data) {
    int n = data.size();
    SuffixArray sa(&data);
    
    vector<int> bwt;
    bwt.resize(n+4);
    vector<int>::iterator bwt_iter = bwt.begin();
    
    long int dollar_ind;
    *(bwt_iter++) = data.back();
    
    for (int i=0; i<n; ++i) {
        if (sa[i] == 0) {
            *(bwt_iter++) = 0;
            dollar_ind = i + 1;
        }
        else *(bwt_iter++) = data[sa[i] - 1];
    }
    
    *(bwt_iter++) = dollar_ind / (256*256);
    *(bwt_iter++) = (dollar_ind / 256) % 256;
    *(bwt_iter++) = dollar_ind % 256;
    
    return bwt;
}

bool int_comparator(int i,int j) { return (i<j); }
bool pair_comparator(pair<int,int> a, pair<int,int> b) { return (a.first < b.first); } // Compare pairs by first elt
// Takes a BWT vector of the format output by burrows_wheeler_encode and returns the original
// string, without the '$' character and without the I integers on the end.
// The return vector is therefore 4 shorter than the input vector.
int *DATCompression::burrows_wheeler_decode(int *data, int n) {
    long int dollar_pos = 256*256*data[n+1] + 256*data[n+2] + data[n+3];
    int bwt_size = n + 1;
    
    pair<int,int> *f = (pair<int,int> *) calloc(bwt_size, sizeof(pair<int,int>));
    int *t = (int *) calloc(bwt_size, sizeof(int));
    int *s = (int *) calloc(n, sizeof(int));
    
    for (int i=0; i<bwt_size; ++i) f[i] = pair<int,int>(data[i], i);
    stable_sort(f, f + bwt_size, pair_comparator);
    
    for (int i=0; i<bwt_size; ++i) t[f[i].second] = i;
    free(f);
    
    int last_t = dollar_pos;
    for (int i=1; i<=n; ++i) {
        last_t = t[last_t];
        s[n-i] = data[last_t];
    }
    free(t);
    
    return s;
}
// Takes a BWT vector of the format output by burrows_wheeler_encode and returns the original
// string, without the '$' character and without the I integers on the end.
// The return vector is therefore 4 shorter than the input vector.
vector<int> DATCompression::burrows_wheeler_decode(vector<int> data) {
    int n = data.size() - 4;
    long int dollar_pos = 256*256*data[n+1] + 256*data[n+2] + data[n+3];
    int bwt_size = n + 1;
    
    pair<int,int> *f = (pair<int,int> *) calloc(bwt_size, sizeof(pair<int,int>));
    int *t = (int *) calloc(bwt_size, sizeof(int));
    vector<int> s;
    s.resize(n);
    
    for (int i=0; i<bwt_size; ++i) f[i] = pair<int,int>(data[i], i);
    f[dollar_pos].first = INT_MIN;
    stable_sort(f, f + bwt_size, pair_comparator);
    
    for (int i=0; i<bwt_size; ++i) t[f[i].second] = i;
    free(f);
    
    int last_t = dollar_pos;
    for (vector<int>::reverse_iterator viter = s.rbegin(); viter != s.rend(); ++viter) {
        last_t = t[last_t];
        *viter = data[last_t];
    }
    free(t);
    
    return s;
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


SuffixArray::SuffixArray(vector<int> *data) {
    n = data->size();    
    s = (int *) calloc(n + 3, sizeof(int));

    copy(data->begin(), data->end(), s);
    
    prepareSuffixArray();
}

SuffixArray::SuffixArray(int *data, int _n) {
    n = _n;
    
    s = (int *) calloc(n + 3, sizeof(int));
    memcpy(s, data, n * (sizeof(int)));
    
    prepareSuffixArray();
}

void SuffixArray::prepareSuffixArray() {
    // Add 3 trailing zeroes for skew algorithm
    int i;
    int *cur_data_elt = s + n;
    *cur_data_elt = 0; ++cur_data_elt;
    *cur_data_elt = 0; ++cur_data_elt;
    *cur_data_elt = 0; ++cur_data_elt;
    
    sa = (int *) calloc(n, sizeof(int));
    for (i=0; i<n; ++i) ++s[i];
    doSuffixArrayComputation(s, sa, n, 256);
    for (i=0; i<n; ++i) --s[i];
}

int * SuffixArray::getSuffixArray() {
    return sa;
}


void SuffixArray::print() {
    for (int i=0; i<(n+3); ++i) {
        printf("s[%d] = %d\n", i, s[i]);
    }
    for (int i=0; i<n; ++i) {
        printf("sa[%d] = %d\n", i, sa[i]);
    }
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
        if (suffixes[i] == true_suffixes[i]) {
            //printf("CORRECT: True suffix SA[%d] = %d\n", i, suffixes[i]);
        } else {
            printf("FAIL: True suffix SA[%d] is %d but found %d\n", i, true_suffixes[i], suffixes[i]);
            return;
        }
    }
    printf("PASS: suffix array is correct\n");
}

void test_bwt() {
    int i, n = 6;
    int s[] = {0,1,3,0,2,0};
    int bwt_expected[] = {0,2,INT_MIN,3,0,0,1,0,0,2};
    int *bwt, *decoded;
    
    vector<int> vec_s, vec_bwt, vec_bwt_expected, vec_decoded;
    
    DATCompression dat;
    dat.init();
    
    bwt = dat.burrows_wheeler_encode(s, n);
    decoded = dat.burrows_wheeler_decode(bwt, n);
    
    for (i=0; i<n+4; ++i) {
        if (bwt[i] != bwt_expected[i]) {
            printf("FAIL: bwt[%d] = %d != bwt_expected[%d] = %d\n", i, bwt[i], i, bwt_expected[i]);
            goto test_btw_cleanup;
        }
    }
    printf("PASS: burrows_wheeler_encode is correct\n");
    
    for (i=0; i<n; ++i) {
        if (s[i] != decoded[i]) {
            printf("FAIL: s[%d] = %d != decoded[%d] = %d\n", i, s[i], i, decoded[i]);
            goto test_btw_cleanup;
        }
    }
    printf("PASS: burrows_wheeler_decode is correct\n");
    
    vec_s.reserve(n);
    vec_bwt_expected.reserve(n+4);
    copy(s, s+n, back_inserter(vec_s));
    copy(bwt_expected, bwt_expected+n+4, back_inserter(vec_bwt_expected));
    vec_bwt_expected[vec_bwt_expected.back()] = 0; // In vector bwt, $ is '0' instead of 'INT_MIN'
    
    vec_bwt = dat.burrows_wheeler_encode(vec_s);
    vec_decoded = dat.burrows_wheeler_decode(vec_bwt);
    
    for (i=0; i<n+4; ++i) {
        if (vec_bwt[i] != vec_bwt_expected[i]) {
            printf("FAIL: vec_bwt[%d] = %d != vec_bwt_expected[%d] = %d\n", i, bwt[i], i, bwt_expected[i]);
            goto test_btw_cleanup;
        }
    }
    printf("PASS: vector burrows_wheeler_encode is correct\n");
    
    for (i=0; i<n; ++i) {
        if (vec_s[i] != vec_decoded[i]) {
            printf("FAIL: vec_s[%d] = %d != vec_decoded[%d] = %d\n", i, vec_s[i], i, vec_decoded[i]);
            goto test_btw_cleanup;
        }
    }
    printf("PASS: vector burrows_wheeler_decode is correct\n");
    
    n=1000;
    vec_s.clear();
    vec_s.reserve(n);
    for (i=0; i<n; ++i) vec_s.push_back(rand() % 256);
    vec_bwt = dat.burrows_wheeler_encode(vec_s);
    vec_decoded = dat.burrows_wheeler_decode(vec_bwt);
    for (i=0; i<n; ++i) {
        if (vec_s[i] != vec_decoded[i]) {
            printf("FAIL: random vec_s[%d] = %d != vec_decoded[%d] = %d\n", i, vec_s[i], i, vec_decoded[i]);
            goto test_btw_cleanup;
        }
    }
    printf("PASS: random burrows_wheeler_decode is correct\n");        
    
test_btw_cleanup:
    free(bwt);
    free(decoded);
}

void test_rle() {
    int n = 14;
    int s[] = {0,0,0,0,0,1,1,1,1,0,1,1,0,0};
    int rle_expected[] = {n,0,0,0,0,1,1,1,1,1,0,0,1,1,0,0};

    int *encoded, *decoded;

    vector<int> vec_s, vec_encoded, vec_decoded;

    DATCompression dat;
    dat.init();

    encoded = dat.run_length_encode(s, n);
    decoded = dat.run_length_decode(encoded);

    for (int i=0; i<n; ++i) {
        if (s[i] != decoded[i]) {
            printf("FAIL: s[%d] = %d != decoded[%d] = %d\n", i, s[i], i, decoded[i]);
            goto test_rle_cleanup;
        }
    }
    printf("PASS: rle is correct\n");
    
    copy(s, s+n, back_inserter(vec_s));
    vec_encoded = dat.run_length_encode(vec_s);
    vec_decoded = dat.run_length_decode(vec_encoded);
    for (int i=0; i<n; ++i) {
        if (vec_s[i] != vec_decoded[i]) {
            printf("FAIL: vec_s[%d] = %d != vec_decoded[%d] = %d\n", i, vec_s[i], i, vec_decoded[i]);
            goto test_rle_cleanup;
        }
    }
    printf("PASS: vector rle is correct\n");
    
test_rle_cleanup:
    free(encoded);
    free(decoded);
}

vector<int> read_data(string filename) {
    ifstream data_file;
    string data_line;
    vector<int> data;
    
    data_file.open(filename.c_str());
    
    int x, y, x_ind, y_ind, L, valCnt, lineCnt;
    string pHval;
    x = 0;
    y = 0;
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
    
    //printf("Read in %ld pH values with (x, y, L) = (%d, %d, %d)\n", data.size() - 3, x, y, L);
    
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
    //cout << "Verified that all " << compressed.size() << " compressed ints are in within [0, 255] legal range" << endl;
    
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
    test_suffix_arrays();
    test_bwt();
    test_rle();
    
    // investigate_file("data/B28-39_100_100_acq_0007.tab");
    
    // test_compression_on_file("data/test.tab");

    test_compression_on_file("data/B28-39_100_100_acq_0007.tab");
    // test_compression_on_file("data/B28-39_100_100_acq_0400.tab");
    //test_compression_on_file("data/B28-39_1600_1000_acq_0007.tab");
    // test_compression_on_file("data/B28-39_1600_1000_acq_0400.tab");
    
    return 0;
}



