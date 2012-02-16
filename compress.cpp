#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <cstring>
#include <vector>

#define PH_MAX 16383
#define HUFFMAN_MAX 16

using namespace std;


class DATCompression {
private:
    int golomb_n; // par for golomb encoding, power of 2

    vector <int> identity_compress(vector <int> data);
    vector <int> identity_decompress(vector <int> data);

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
    vector <int> compress(vector <int> data) { return huffman_compress(data); }
    vector <int> decompress(vector <int> compressed) { return huffman_decompress(compressed); }

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
  int freq[HUFFMAN_MAX], in_set[HUFFMAN_MAX], masked[HUFFMAN_MAX];
  int max_freq, second_max_freq, max_ind, second_max_ind, s, temp = 0, code_size, t = 0, q, r;
  vector<int> code[HUFFMAN_MAX];
  
  memset(freq, 0, sizeof(freq));
  memset(masked, 0, sizeof(masked));
  for(int i = 0; i < data.size(); i++)
  {
    if(data[i] < 0 || data[i] >= HUFFMAN_MAX * HUFFMAN_MAX)
    {
      cout << "Error: huffman_compress do not accept data out of range [0..255]" << endl;
      //exit(0);
    }
    freq[data[i] / HUFFMAN_MAX]++;
    freq[data[i] % HUFFMAN_MAX]++;
  }

  for(int i = 0; i < HUFFMAN_MAX; i++)
  {
    code[i].clear();
    in_set[i] = i;
  }

  for(int i = 0; i < HUFFMAN_MAX - 1; i++)
  {
    // find out the minimum and second miniimum freq
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
      temp = temp * 2 + (code_size % 2);
      code_size /= 2, t++;
      if(t == 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }

    for(int j = code[i].size() - 1; j >= 0; j--)
    {
      temp = temp * 2 + code[i][j], t++;
      if(t == 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }
  }

  for(int i = 0; i < data.size(); i++)
  {
    q = data[i] / HUFFMAN_MAX;
    r = data[i] % HUFFMAN_MAX;
    for(int j = code[q].size() - 1; j >= 0; j--)
    {
      temp = temp * 2 + code[q][j], t++;
      if(t == 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }
    for(int j = code[r].size() - 1; j >= 0; j--)
    {
      temp = temp * 2 + code[r][j], t++;
      if(t == 8)
      {
        compressed.push_back(temp);
        temp = 0, t = 0;
      }
    }
  }

  if(t)
  {
    for(int i = 0; i < 8 - t; i++)
      temp = temp * 2;
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
  int code_pointer[HUFFMAN_MAX];
  int number, d = 0, power = 1, code_ind = 0, c, temp = 0, t = 0, cur_bit, q, r, g = 0;
  int useless_bits = compressed[compressed.size() - 1];
  char status = 'd';

  memset(code_pointer, 0, sizeof(code_pointer));
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
          status = 't';
      }
      else if(status == 't')
      {
        if(i == compressed.size() - 2 && 8 - i == useless_bits)
          break;

        if((number & 0x80) == 0x80)
          cur_bit = 1;
        else
          cur_bit = 0;
        number <<= 1;
        for(int k = 0; k < HUFFMAN_MAX; k++)
          if(!code_pointer[k])
            if(code[k][t] == cur_bit && t == code[k].size() - 1)
            {
              t = -1;
              memset(code_pointer, 0, sizeof(code_pointer));
              if(!g)
                q = k, g = 1;
              else
              {
                r = k, g = 0;
                decompressed.push_back(q * HUFFMAN_MAX + r);
              }
              break;
            }
            else if(code[k][t] != cur_bit)
              code_pointer[k] = 1;
        t++;
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
    
    //investigate_file("data/B28-39_100_100_acq_0007.tab");
    
    test_compression_on_file("data/test.tab");

    //test_compression_on_file("data/B28-39_100_100_acq_0007.tab");
    //test_compression_on_file("data/B28-39_100_100_acq_0400.tab");
    //test_compression_on_file("data/B28-39_1600_1000_acq_0007.tab");
    //test_compression_on_file("data/B28-39_1600_1000_acq_0400.tab");
    
    return 0;
}



