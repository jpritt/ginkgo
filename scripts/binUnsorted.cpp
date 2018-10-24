#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <stdlib.h>

using namespace std;

struct char201
{
  char val[201];

  char201() {
    memset(val, '\0', 201);
  }

  char201(const char201& oth) {
    memmove(val, oth.val, 201);
  }

  char201(const char * str) {
    memset(val, '\0', 201);
    if (str)
    {
      strncpy(val, str, 200);
    }
  }


  bool operator<(const char201& oth) const {
    return strcmp(val , oth.val) < 0; 
  }

};

int main(int argc, char *argv[]){
   
  FILE * bin_file;
  FILE * reads_file;
  //FILE * outfile;
  bin_file = fopen(argv[1], "r");
  reads_file = fopen(argv[3], "r");
  //outfile = fopen(argv[5], "w");
 
  if(argc < 7)
  {
    cout << "USAGE: binFile #bins readsFile nameReadsFile outfile num_samples [cells]" << endl;
    return 1;
  }

  FILE * cells_file;
  int num_samples = atoi(argv[6]);
  if (num_samples > 1)
  {
    if(argc < 8)
    {
      cout << "USAGE: binFile #bins readsFile nameReadsFile outfile num_samples [cells]" << endl;
      return 1;
    }
    cells_file = fopen(argv[7], "r");
  }

  
  if(bin_file == NULL || reads_file == NULL)
  {
    cout << "Unable to open input files: " << argv[1] << " or " << argv[3] << " " << argv[3] << endl;
    return 1;
  }

  string* cells;
  if (num_samples == 1 || cells_file == NULL)
  {
    cout << "No valid cells file, reading file as a single sample" << endl;
  }
  else
  {
    cout << "Splitting file into samples listed in " << argv[7] << endl;
    cells = new string[num_samples];
    char line[100];
    size_t buffer = 1000;
    int id = 0;
    ssize_t read;
    while (fgets(line, 100, cells_file) != NULL)
    {
      line[strlen(line) - 1] = '\0';
      cells[id] = line;
      id++;
    }
    fclose(cells_file);
  }

  pair<int, int> bound;
  map<char201, pair<int,int> > bin_map;
  char chrom[201];
  char prev_chr[201];
  char new_chr[201];
  char dump[201];
  int len = atoi(argv[2]) - 1;
  int *loc = new int[len];
  int *mappable = new int[len];
  float *good_count = new float[len];  

  //int *bins = new int[num_samples][len];
  int bins[num_samples][len];
  memset(bins, 0, num_samples*len*sizeof(int));
  int low_bound = 0;
  int high_bound = 0;
  bool skip = true;  // If this is the first line, don't compare chromosome to previous value

  for (int i = 0; i < num_samples; i++)
    for (int j = 0; j < len; j++)
      bins[i][j] = 0;

  //Creates a map where:
  //Key = Chromosome
  //Value = Pair of integers containing bin boundaries for that chromosome
  //fscanf(bin_file, "%201s%201s", &dump, &dump);  // Ignore header line

  //Count columns in bin file
  char * line;
  size_t buffer = 1000;
  ssize_t header;
  int cols = 1;
  header = getline(&line, &buffer, bin_file);
  for (int i=0; i<header; i++)
  {
    if (line[i] == '\t')
      cols++;
  }

  if (cols == 2) {
    fscanf(bin_file, "%201s%201s", &dump, &dump);
    while (fscanf(bin_file, "%201s%i", &new_chr, &loc[high_bound]) != EOF)
    {
      mappable[high_bound] = 1;
      good_count[high_bound] = 1;
      if (strcmp(prev_chr, new_chr) != 0)
      {
        if (skip == false)
        {
          bound.first = low_bound;
          bound.second = high_bound-1;
          bin_map[prev_chr]= bound;
          if (strncmp(prev_chr, "chr", 3) == 0) {
              bin_map[prev_chr + 3]= bound;
          }
          low_bound = high_bound;
        }
      skip=false;
      }
      strcpy(prev_chr, new_chr);
      high_bound++;
    }
  }
  else if (cols >= 4) {
    while (fscanf(bin_file, "%201s%i%i%f", &new_chr, &loc[high_bound], &mappable[high_bound], &good_count[high_bound]) != EOF)
    {
      if (strcmp(prev_chr, new_chr) != 0)
      {
        if (skip == false)
        {
          bound.first = low_bound;
          bound.second = high_bound-1;
          bin_map[prev_chr]= bound;
          if (strncmp(prev_chr, "chr", 3) == 0) {
              bin_map[prev_chr + 3]= bound;
          }
          low_bound = high_bound;
        }
      skip=false;
      }
      strcpy(prev_chr, new_chr);
      high_bound++;
    }
  }
  else {
    cout << "Error! Bin file must contain either 2 or 4 columns." << endl;
    return 1;
  }

  bound.first = low_bound;
  bound.second = high_bound-2;
  bin_map[new_chr] = bound;
  if (strncmp(new_chr, "chr", 3) == 0) {
    bin_map[new_chr + 3]= bound;
  }

  //char * line;
  //size_t buffer = 1000;
  ssize_t read;
  int tabs = 0;

  // Count the columns; only care about the first 2, and the last if we are splitting by cell
  read = getline(&line, &buffer, reads_file);
  for (int i=0; i<read; i++)
  {
    if (line[i] == '\t')
      tabs++;
  }
  if (num_samples > 1)
    tabs -= 2;
  else
    tabs --;
  rewind(reads_file);
  
  int pos, low, mid, high;
  char cell[201];
  int cnt = 0;

  //Uses binary search within chrom boundaries calculated above
  //to map reads into bins
  if (num_samples == 1)
  {
    while (fscanf(reads_file, "%201s%i", &chrom, &pos) != EOF)
    {
   
      // Skip extra columns
      for (int i=0; i<tabs; i++) {
        fscanf(reads_file, "%201s", &dump);
      }

      if (bin_map.count(chrom) == 0) {
        continue;
      }
   
      low = (bin_map[chrom]).first;
      high = (bin_map[chrom]).second;
      while (low <= high)
      {
        mid = (low + high) / 2;
        if (pos < loc[mid])
	  high = mid - 1;
        else if (pos > loc[mid])
          low = mid + 1;
        else
          low = high +1;
      }
      
      if (pos > loc[mid])
        mid++;
    
      bins[0][mid]++; 
    }
  }
  else
  {
    while (fscanf(reads_file, "%201s%i", &chrom, &pos) != EOF)
    {

      // Skip extra columns
      for (int i=0; i<tabs; i++) {
        fscanf(reads_file, "%201s", &dump);
      }

      // Read cell name
      fscanf(reads_file, "%201s", &cell);
      int sample_id = -1;

      for (int i = 0; i < num_samples; i++)
      {
        if (strcmp(cells[i].c_str(), cell) == 0)
        {
          sample_id = i;
          break;
        }
      }

      if (sample_id < 0)
        continue;

      if (bin_map.count(chrom) == 0) {
        continue;
      }
   
      low = (bin_map[chrom]).first;
      high = (bin_map[chrom]).second;
      while (low <= high)
      {
        mid = (low + high) / 2;
        if (pos < loc[mid])
	  high = mid - 1;
        else if (pos > loc[mid])
          low = mid + 1;
        else
          low = high +1;
      }
      
      if (pos > loc[mid])
        mid++;
    
      bins[sample_id][mid]++; 
    }
  }

  int wgt = 1;
  if (num_samples == 1) {
    FILE * outfile;
    outfile = fopen(argv[5], "w");

    fprintf(outfile, "%s\n", argv[4]);
    for(int i=0; i < len; i++) {
        //cout << bins[0][i] << ", " << mappable[i] << ", " << good_count[i] << endl;
        if (mappable[i] == 1) {
          fprintf(outfile, "%f\n", (wgt * bins[0][i] / float(good_count[i])));
        }
        else {
          //if (i == 0 || !mappable[i-1])
          //  fprintf(outfile, "%f\n", (wgt * bins[0][i+1] / float(good_count[i+1])));
          //else if (i == (len-1) || !mappable[i+1])
          //  fprintf(outfile, "%f\n", (wgt * bins[0][i-1] / float(good_count[i-1])));
          //else
          //  fprintf(outfile, "%f\n", (wgt * bins[0][i-1] / (2*float(good_count[i-1])) + wgt * bins[0][i+1] / (2*float(good_count[i+1]))));
          fprintf(outfile, "-1\n");
        }
        //fprintf(outfile, "%i\n", bins[0][i]);
    }
  }
  else {
    for (int s=0; s < num_samples; s++) {
      FILE * outfile;
      std::string str = argv[5]+cells[s];
      //std::string str = cells[s];
      outfile = fopen(str.c_str(), "w");

      //std::string name = argv[4]+'.'+cells[s]+'\n';
      //fprintf(outfile, name.c_str());
      fprintf(outfile, "%s.%s\n", argv[4], cells[s].c_str());
      for(int i=0; i < len; i++) {
        if (mappable[i]) {
          fprintf(outfile, "%f\n", (wgt * bins[s][i] / float(good_count[i])));
        }
        else {
          fprintf(outfile, "%f\n", -1);
        }
        //fprintf(outfile, "%i\n", bins[s][i]);
      }
      fclose(outfile);
    }
  }

  return 0;
}
