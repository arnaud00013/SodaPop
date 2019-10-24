#include "global.h"

VectStr PrimordialAASeq;
double fold_DG = 0;
double bind_DG = 0;

int Total_Cell_Count = 0;
int dummy = 0;
double frame_time = 0;
std::string outPath;
char buffer[200];

double matrix[gene_number][res_number][20];
double matrix_supp[gene_number][res_number][20];

/******* GENETIC CODE MAPPINGS *******/
// these const mappings are hard-coded and populated at compile-time

struct codon_to_num{
    static std::map<std::string,int> create_map()
        {
          std::map<std::string,int> m;
          m["GCA"] = 1;
          m["GCC"] = 1;
          m["GCG"] = 1;
          m["GCT"] = 1;
          m["TGC"] = 2;
          m["TGT"] = 2;
          m["GAC"] = 3;
          m["GAT"] = 3;
          m["GAA"] = 4;
          m["GAG"] = 4;
          m["TTC"] = 5;
          m["TTT"] = 5;
          m["GGA"] = 6;
          m["GGC"] = 6;
          m["GGG"] = 6;
          m["GGT"] = 6;
          m["CAC"] = 7;
          m["CAT"] = 7;
          m["ATA"] = 8;
          m["ATC"] = 8;
          m["ATT"] = 8;
          m["AAA"] = 9;
          m["AAG"] = 9;
          m["CTA"] = 10;
          m["CTC"] = 10;
          m["CTG"] = 10;
          m["CTT"] = 10;
          m["TTA"] = 10;
          m["TTG"] = 10;
          m["ATG"] = 11;
          m["AAC"] = 12;
          m["AAT"] = 12;
          m["CCA"] = 13;
          m["CCC"] = 13;
          m["CCG"] = 13;
          m["CCT"] = 13;
          m["CAA"] = 14;
          m["CAG"] = 14;
          m["AGA"] = 15;
          m["AGG"] = 15;
          m["CGA"] = 15;
          m["CGC"] = 15;
          m["CGG"] = 15;
          m["CGT"] = 15;
          m["AGC"] = 16;
          m["AGT"] = 16;
          m["TCA"] = 16;
          m["TCC"] = 16;
          m["TCG"] = 16;
          m["TCT"] = 16;
          m["ACA"] = 17;
          m["ACC"] = 17;
          m["ACG"] = 17;
          m["ACT"] = 17;
          m["GTA"] = 18;
          m["GTC"] = 18;
          m["GTG"] = 18;
          m["GTT"] = 18;
          m["TGG"] = 19;
          m["TAC"] = 20;
          m["TAT"] = 20;
          m["TAA"] = 21;
          m["TAG"] = 21;
          m["TGA"] = 21;
          return m;
        }
    static const std::map<std::string,int> cnum;
};

struct codon_to_prot{
    static std::map<std::string,char> create_map()
        {
          std::map<std::string,char> m;
          m["AAA"] = 'K';
          m["AAC"] = 'N';
          m["AAG"] = 'K';
          m["AAT"] = 'N';
          m["ACA"] = 'T';
          m["ACC"] = 'T';
          m["ACG"] = 'T';
          m["ACT"] = 'T';
          m["AGA"] = 'R';
          m["AGC"] = 'S';
          m["AGG"] = 'R';
          m["AGT"] = 'S';
          m["ATA"] = 'I';
          m["ATC"] = 'I';
          m["ATG"] = 'M';
          m["ATT"] = 'I';
          m["CAA"] = 'Q';
          m["CAC"] = 'H';
          m["CAG"] = 'Q';
          m["CAT"] = 'H';
          m["CCA"] = 'P';
          m["CCC"] = 'P';
          m["CCG"] = 'P';
          m["CCT"] = 'P';
          m["CGA"] = 'R';
          m["CGC"] = 'R';
          m["CGG"] = 'R';
          m["CGT"] = 'R';
          m["CTA"] = 'L';
          m["CTC"] = 'L';
          m["CTG"] = 'L';
          m["CTT"] = 'L';
          m["GAA"] = 'E';
          m["GAC"] = 'D';
          m["GAG"] = 'E';
          m["GAT"] = 'D';
          m["GCA"] = 'A';
          m["GCC"] = 'A';
          m["GCG"] = 'A';
          m["GCT"] = 'A';
          m["GGA"] = 'G';
          m["GGC"] = 'G';
          m["GGG"] = 'G';
          m["GGT"] = 'G';
          m["GTA"] = 'V';
          m["GTC"] = 'V';
          m["GTG"] = 'V';
          m["GTT"] = 'V';
          //STOP
          m["TAA"] = 'X';
          m["TAC"] = 'Y';
          //STOP
          m["TAG"] = 'X';
          m["TAT"] = 'Y';
          m["TCA"] = 'S';
          m["TCC"] = 'S';
          m["TCG"] = 'S';
          m["TCT"] = 'S';
          //STOP
          m["TGA"] = 'X';
          m["TGC"] = 'C';
          m["TGG"] = 'W';
          m["TGT"] = 'C';
          m["TTA"] = 'L';
          m["TTC"] = 'F';
          m["TTG"] = 'L';
          m["TTT"] = 'F';
          return m;
        }
    static const std::map<std::string,char> cprot;
};

struct prot_to_num{
    static std::map<char,int> create_map()
        {
          std::map<char,int> m;
          m['A'] = 1;
          m['C'] = 2;
          m['D'] = 3;
          m['E'] = 4;
          m['F'] = 5;
          m['G'] = 6;
          m['H'] = 7;
          m['I'] = 8;
          m['K'] = 9;
          m['L'] = 10;
          m['M'] = 11;
          m['N'] = 12;
          m['P'] = 13;
          m['Q'] = 14;
          m['R'] = 15;
          m['S'] = 16;
          m['T'] = 17;
          m['V'] = 18;
          m['W'] = 19;
          m['Y'] = 20;
          m['X'] = 21;//STOP
          return m;
        }
    static std::map<char,int> const pnum;
};

// populate genetic code mappings
std::map<std::string,int> const codon_to_num::cnum = codon_to_num::create_map();
std::map<char,int> const prot_to_num::pnum = prot_to_num::create_map();
std::map<std::string,char> const codon_to_prot::cprot = codon_to_prot::create_map();

// string hash function from Github user rioki
constexpr
unsigned int hashString(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (hashString(str, h+1)*33) ^ str[h];
}

Encoding_Type intToEncoding_Type(int type){
    switch (type){
        case 0:
            return Encoding_Type::by_default;
        case 1:
            return Encoding_Type::no_sequence;
        case 2:
            return Encoding_Type::full;
        case 3:
            return Encoding_Type::other;
        default:
            return Encoding_Type::by_default;
    }
}

Init_Pop intToPop_Type(int type){
    switch (type){
        case 0:
            return Init_Pop::from_snapFile;
        case 1:
            return Init_Pop::from_cellFile;
        default:
            return Init_Pop::from_snapFile;
    }
}

Input_Type stringToInput_Type(const std::string& type){
    switch (hashString(type.c_str())){
        case hashString("s"):
            return Input_Type::selection_coefficient;
          break;
        case hashString("stability"):
            return Input_Type::stability;
          break;
        default:
            return Input_Type::selection_coefficient;
          break;
    }
}

void openCommandLog(std::ofstream& cmdlog, std::string outDir, char *argv[], int argc){
    sprintf(buffer,"out/%s/command.log",outDir.c_str());
    cmdlog = std::ofstream(buffer, std::ios::out | std::ios::trunc);

    if (cmdlog.is_open()){
        // file was opened successfully
        std::cout << "-> Command log was opened successfully ..." << std::endl;
        std::string args;
        std::for_each( argv + 1, argv + argc , [&]( const char* c_str ){ args += std::string ( c_str ) + " "; } );
        cmdlog << "sodapop " << args << std::endl;
        cmdlog << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open command log.");
    }
}

void openMutationLog(std::ofstream& mutlog, std::string outDir){
    sprintf(buffer, "out/%s/MUTATION_LOG",outDir.c_str());
    mutlog = std::ofstream(buffer, std::ios::out | std::ios::trunc);

    if (mutlog.is_open()){
        // file was opened successfully
        std::cout << "-> Mutation log was opened successfully ..." << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open mutation log.");
    }
}

void openPevLog(std::ofstream& pevlog, std::string outDir){
	// Open pevlog
	sprintf(buffer, "out/%s/PANGENOMES_EVOLUTION_LOG.txt",outDir.c_str());
	pevlog.open(buffer);
	if ( !pevlog.is_open() ) {
		std::cerr << "Pangenomes evolution log file could not be opened";
		exit(1);
	}else{
		std::cout << "Opening pangenomes evolution log ..." << std::endl;
		pevlog <<"Generation_ctr"<<"\t"<<"cell_ID"<<"\t"<<"x"<<"\t"<<"r_x"<<"\t"<<"Beta_x"<<"\t"<<"Alpha_x"<<"\t"<<"fitness"<<std::endl;
	}
}

void openCellsgenecontentLog(std::ofstream& cells_gene_content_log, std::string outDir){
	// Open cells_gene_content_log
	sprintf(buffer, "out/%s/CELL_GENE_CONTENT_LOG.txt",outDir.c_str());
	cells_gene_content_log.open(buffer);
	if ( !cells_gene_content_log.is_open() ) {
		std::cerr << "Cell gene content log file could not be opened";
		exit(1);
	}else{
		std::cout << "Opening Cell gene content log ..." << std::endl;
		cells_gene_content_log <<"Generation_ctr"<<"\t"<<"cell_ID"<<"\t"<<"gene_ID"<<"\t"<<"nucl_sequence"<<"\t"<<"cai"<<std::endl;
	}
}

void openGainLog(std::ofstream& gainlog, std::string outDir){
	// Open gainlog
	sprintf(buffer, "out/%s/GENE_GAIN_EVENTS_LOG.txt",outDir.c_str());
	gainlog.open(buffer);
	if ( !gainlog.is_open() ) {
		std::cerr << "Gene gain events log file could not be opened";
		exit(1);
	}else{
		std::cout << "Opening gene gain events log ..." << std::endl;
		//d_cell is the donor cell and r_cell is the receiving cell
		gainlog <<"gain_event_ID"<<"\t"<<"Generation_ctr"<<"\t"<<"d_cell_ID"<<"\t"<<"d_cell_barcode"<<"\t"<<"r_cell_ID"<<"\t"<<"r_cell_barcode"<<"\t"<<"gene_ID"<<std::endl;
	}
}

void openLossLog(std::ofstream& losslog, std::string outDir){
	// Open losslog
	sprintf(buffer, "out/%s/GENE_LOSS_EVENTS_LOG.txt",outDir.c_str());
	losslog.open(buffer);
	if ( !losslog.is_open() ) {
		std::cerr << "Gene loss events log file could not be opened";
		exit(1);
	}else{
		std::cout << "Opening gene loss events log ..." << std::endl;
		//t_cell is the target cell
		losslog <<"loss_event_ID"<<"\t"<<"Generation_ctr"<<"\t"<<"t_cell_ID"<<"\t"<<"t_cell_barcode"<<"\t"<<"gene_ID"<<std::endl;
	}
}


void openStartingPop(std::string filePath, std::ifstream& fileStream){
    std::cout << "Opening starting population snapshot ..." << std::endl;
    fileStream = std::ifstream(filePath.c_str(),std::ios::in|std::ios::binary);
    if (fileStream.is_open()){
        // file was opened successfully
        std::cout << "-> File was opened successfully ..." << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open starting population snapshot.");
    }
}

void readSnapshotHeader(std::ifstream& snapshot)
{
    snapshot.read((char*)(&frame_time),sizeof(double));
    //read number of cells in file
    snapshot.read((char*)(&Total_Cell_Count),sizeof(int));
    //read file outputEncoding
    snapshot.read((char*)(&dummy),sizeof(int));
}

void createOutputDir(std::string dirName){
    sprintf(buffer,"out/%s/snapshots",dirName.c_str());
    outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... " << (makePath(outPath) ? "OK" : "failed") << std::endl;
}

/******* MAPPING FUNCTIONS *******/

// input:  3 nucleotide codon as string
// output: index number of codon
int GetIndexFromCodon(std::string in_codon)
{
    auto it = codon_to_num::cnum.find(in_codon);
    if (it == codon_to_num::cnum.end()){
        std::cerr << "Invalid codon: "<< in_codon << std::endl;
        std::cerr << "Nucleotide sequence must not contain STOP codons."<< std::endl;
        exit(2);
    }
    return it->second;
}

// input:  3 nucleotide codon sequence as string
// output: amino acid sequence as string
std::string GetProtFromNuc(std::string in_seq)
{
    const int ln = in_seq.length();
    if ((ln % 3) != 0){
        std::cerr << "Invalid length for nucleotide sequence: " << ln << std::endl;
        std::cerr << "Nucleotide sequence length must be divisible by 3." << std::endl;
        exit(2);
    }
    const int la=ln/3;
    std::string AA="";
    for (int i=0; i<la;++i){
        std::string temp=in_seq.substr(i*3,3);
        
        auto it = codon_to_prot::cprot.find(temp);
        if (it == codon_to_prot::cprot.end()){
          std::cerr << "Invalid codon: "<< temp << std::endl;
          std::cerr << "Nucleotide sequence must not contain STOP codons."<< std::endl;
          exit(2);
        }   
        AA.push_back(codon_to_prot::cprot.at(temp));
    }
    return AA;
}

// input:  amino acid letter as string
// output: index number of amino acid
int GetIndexFromAA(char aa)
{  
    //check for valid amino acid
    auto it = prot_to_num::pnum.find(aa);
    if (it == prot_to_num::pnum.end()){
        std::cerr << "Invalid amino acid: "<< aa << std::endl;
        exit(2);
    }
    return it->second;
}

// returns the next or second to next bp from a given nucleotide
const char AdjacentBP(char a, int j){
 
    if ( j > 2 ){
        std::cerr << "Invalid bp distance. Error in AdjacentBP(). "<< j << std::endl;
        exit(2);
    }

    std::map <char, int> A;
    A['A'] = 0;
    A['T'] = 1;
    A['G'] = 2;
    A['C'] = 3;

    std::map <int, char> B;
    B[0] = 'A';
    B[1] = 'T';
    B[2] = 'G';
    B[3] = 'C';

    int x = (A[a] + j + 1) % 4; 	//get (j+1) nucleotide from a
    return B[x];
}

//Checks if codon b when mutated resulted in a STOP codon a.
//If b is a STOP codon, it is replaced by 
//Input: 	a, new codon
//		    b, old codon
//		    i, mutation site in codon (< 3)
std::string n3_to_n3(std::string a, std::string b, int i){
  //double r = RandomNumber();
  double r = randomNumber();
  double l;

  //change statement to switch
  if ( a == "TAA"){
      switch (i)
      {
          case 0:
              l =  (2*r);

              if (b == "CAA" ){
                if  ( l<1 )
                  a = "GAA";
                else
                  a = "AAA";
              }else if (b == "GAA" ){
                if  ( l<1 )
                  a = "CAA";
                else
                  a = "AAA";
              }else if (b == "AAA" ){
                if  ( l<1 )
                  a = "GAA";
                else
                  a = "CAA";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 1:
              if 	(b == "TTA")
                a = "TCA";
              else if (b == "TCA")
                a = "TTA";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

            break;

          case 2:
              if 	(b == "TAT")
                a = "TAC";
              else if (b == "TAC")
                a = "TAT";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

            break;

           default:
              std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
      }
  }
  else if (a == "TAG")
  {
      switch (i)
      {
          case 0:
              l = (2*r);

              if (b == "AAG" ){
                if  ( l<1 )
                  a = "GAG";
                else
                  a = "CAG";
              }else if (b == "GAG" ){
                if  ( l<1 )
                  a = "AAG";
                else
                  a = "CAG";
              }else if (b == "CAG" ){
                if  ( l<1 )
                  a = "GAG";
                else
                  a = "AAG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 1:
              l = (2*r);

              if (b == "TTG" ){
                if  ( l<1 )	a = "TGG";
                else 		a = "TGC";
              }else if (b == "TGG" ){
                if  ( l<1 )	a = "TTG";
                else 		a = "TGC";
              }else if (b == "TCG" ){
                 if  ( l<1 )	a = "TTG";
                else 		a = "TGG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 2:
              if 	(b == "TAT") a = "TAC";
              else if (b == "TAC") a = "TAT";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

              break;

             default:
                std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
       }
   }
   else if (a == "TGA")
   {
      switch (i)
      {
          case 0:
              l =  (2*r);

              if (b == "AGA" ){
                if  ( l<1 )	a = "GGA";
                else 		a = "CGA";
              }else if (b == "GGA" ){
                if  ( l<1 )	a = "AGA";
                else 		a = "CGA";
              }else if (b == "CGA" ){
                 if  ( l<1 )	a = "AGA";
                else 		a = "GGA";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

         case 1:
              if 	(b == "TTA") a = "TCA";
              else if (b == "TCA") a = "TTA";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

              break;

          case 2:
              l =  (2*r);

              if (b == "TGT" ){
                if  ( l<1 )	a = "TGG";
                else 		a = "TGC";
              }else if (b == "TGG" ){
                if  ( l<1 )	a = "TGT";
                else 		a = "TGC";
              }else if (b == "TGC" ){
                 if  ( l<1 )	a = "TGT";
                else 		a = "TGG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

           default:
              std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
      }  
  }
  return a;
}

// generates a random, 15nt barcode
std::string getBarcode()
{
    char seq [16];
    for (int i = 0; i < 15; ++i){
        if (randomNumber()<0.5){
            if (randomNumber()<0.5){
                seq[i] = 'G';
            }
            else seq[i] = 'C';
        }
        else{
            if (randomNumber()<0.5){
                seq[i] = 'A';
            }
            else seq[i] = 'T';
        }
    }
    seq[15] = '\0';
    return std::string(seq);
}

// initializes the 3D matrix for DDG values
void InitMatrix()
{
    std::cout << "Initializing matrix ..." << std::endl;
    for (int i = 0; i != gene_number; ++i)
      for (int j = 0; j != res_number; ++j)
        for (int k = 0; k != 20; ++k){
            matrix[i][j][k] = 1; //i.e., exp(-0/kT) = 1
            matrix_supp[i][j][k] = 1; // maybe it would make more sense here to set it to inf
        }
}

// extracts values from the DDG file and stores them in the matrix
double ExtractDDGMatrix(std::string filepath, Matrix_Type m)
{
    std::ifstream temp(filepath);
    if (!temp.is_open()){
        std::cerr << "File could not be opened: "<< filepath << std::endl;
        exit(2);
    }
    std::cout << "Extracting DDG matrix ..." << std::endl;
    std::string line;
    int gene_num = 0;
    double sum = 0;
    int idx = 0;
    while(!temp.eof()){
        std::string word;
        getline(temp,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if ( word == "DDG"){
            iss >> word;
            //residue index
            const int i = atoi(word.c_str());
            for (int j = 0; iss>>word; ++j){
                //extract DDG values
                const double x = atof(word.c_str());
                sum +=x;
                ++idx;
                //if(m){
                    //matrix_supp[gene_num][i-1][j] = exp(-x/kT);
                //}
                //else{
                    matrix[gene_num][i-1][j] = exp(-x/kT);
                //}
            }
        }
        else if ( word == "rCat"){
            iss >> word;
            //residue index
            const int i = atoi(word.c_str());
            for (int j = 0; iss>>word; ++j){
                //extract DDG values
                const double x = atof(word.c_str());
                sum +=x;
                ++idx;
                matrix_supp[gene_num][i-1][j] = x;
            }
        }
        else if ( word == "Gene_NUM"){
            iss >> word;
            gene_num = atoi(word.c_str());
        }
    }
    temp.close();
    return sum/idx;
}

void ExtractDMSMatrix(std::string filepath)
{
    std::ifstream temp(filepath);
    if (!temp.is_open()){
        std::cerr << "File could not be open: "<< filepath << std::endl;
        exit(2);
    }
    std::cout << "Extracting DMS matrix ..." << std::endl;
    std::string line;
    int gene_num = 0;
    while (!temp.eof()){
        std::string word;
        getline(temp,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if ( word == "DMS"){
            iss >> word;
            //residue index
            const int i = atoi(word.c_str());
            for (int j = 0; iss>>word; ++j){
                //extract DMS values
                const double x = atof(word.c_str());
                matrix[gene_num][i-1][j] = x;
            }
        }
        else if ( word == "Gene_NUM"){
            iss >> word;
            gene_num = atoi(word.c_str());
        }
    }
    temp.close();
}

double Ran_Gaussian(double const mean, double const sigma)
{
    double x, y, r2;
    do{
        // choose x,y in uniform square [-1,+1]
        x = -1 + 2 * randomNumber();
        y = -1 + 2 * randomNumber();
        // check if it is in the unit circle
        r2 = x * x + y * y;
    }while (r2 > 1.0 || r2 == 0); 
    // Box-Muller transform
    return mean + sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// Loads primordial genes in a VectStr
int LoadPrimordialGenes(const std::string& genelistfile, const std::string& genesPath)
{  
    std::ifstream genelistIN (genelistfile.c_str());
    if (!genelistIN.is_open()){
        std::cerr << "File could not be open: "<< genelistfile <<std::endl;
        exit(2);
    }
    std::cout << "Loading primordial genes file ..." << std::endl;
    int flag_AASeq = 0; 
    int gc = 0;
    while(!genelistIN.eof()){
        std::string word, line;
        getline(genelistIN,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word; 
        if ( word=="G" ){
            iss >> word;
            word = genesPath + word;
            std::ifstream genefileIN (word.c_str(), std::ifstream::in | std::ifstream::out);
            if (!genefileIN.is_open()){
                std::cerr << "File could not be open: " << word << std::endl;
                exit(2);
            }
            assert ( gc > 0);

            int gn = -1;
            while(!genefileIN.eof()){
                  std::string w, l;
                  getline(genefileIN,l);
                  std::istringstream iss(l, std::istringstream::in);
                  iss >> w;
                  if (w=="Gene_NUM"){
                      iss >> w;
                      gn = atoi(w.c_str());
                  }
                  else if ( w=="N_Seq" ){
                      assert( gn >= 0);
                      iss >> w;
                      std::string aaseq=GetProtFromNuc(w);
                      //check stop codons in midsequence
                      size_t loc = aaseq.find('X', 0);
                      assert( loc == std::string::npos ); // no match
                      auto iter = PrimordialAASeq.begin();
                      PrimordialAASeq.insert(iter+gn, aaseq); 
                      flag_AASeq += 1;
                  }
            }
        }
        else if ( word=="Gene_Count" ){
            iss >> word;
            gc = atoi(word.c_str());
            PrimordialAASeq.reserve(gc);
        }
    }
    assert( flag_AASeq==gc );
    return gc;
}

// Reads a unit cell stored in binary format using Cell::dumpShort()
void qread_Cell(std::ifstream& IN, std::ofstream& OUT)
{
    char mybuffer[140];
    int na = 0;
    int ns = 0;
    double f = 0;
    std::string barcode;

    int l = 0;
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&na),sizeof(int));
    IN.read((char*)(&ns),sizeof(int));
    IN.read((char*)(&f),sizeof(double));
    
    sprintf(mybuffer,"\t%d\t%d\t%e", na, ns, f);
    OUT << mybuffer;
}

// Reads a unit cell stored in binary format using Cell::dumpSeq()
void seqread_Cell(std::ifstream& IN, std::ofstream& OUT)
{
    char mybuffer[140];
    int cell_id(0);
    int cell_index(0);
    int gene_size(0);
    double f(0);
    double m(0);
    std::string barcode;

    IN.read((char*)(&cell_index),sizeof(int));
    IN.read((char*)(&cell_id),sizeof(int));

    int l(0);
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&f),sizeof(double));
    IN.read((char*)(&m),sizeof(double));
    IN.read((char*)(&gene_size),sizeof(int));
    
    sprintf(mybuffer,"\t%d\t%e\t%e\t", cell_index, f, m);

    OUT << mybuffer << std::endl;

    for (int j=0; j<gene_size; ++j){
        std::string DNAsequence;   
        int Na(0);
        int Ns(0);

        IN.read((char*)(&Na),sizeof(int));
        IN.read((char*)(&Ns),sizeof(int));

        //read DNA sequence
        int nl(0);
        IN.read((char*)&nl, sizeof(int));
        std::vector<char> buff(nl);
        IN.read(&buff[0], nl);  
        DNAsequence.assign(buff.begin(), buff.end());

        sprintf(mybuffer,"%d\tG\t%d\t%d\t",j,Ns,Na);
        OUT << mybuffer << std::endl;
        OUT << DNAsequence << std::endl;
    }
}

void read_Parent(std::ifstream& IN, std::ofstream& OUT)
{
    uint32_t a;

    IN.read((char*)(&a),sizeof(a));
    
    OUT << a;
}

// Reads a unit cell stored in binary format using Cell::dump()
void read_Cell(std::ifstream& IN, std::ofstream& OUT, bool DNA)
{
    char mybuffer[140];
    int cell_id(0);
    int cell_index(0);
    int gene_size(0);
    double f(0);
    double m(0);
    std::string barcode;

    IN.read((char*)(&cell_index),sizeof(int));
    IN.read((char*)(&cell_id),sizeof(int));

    int l(0);
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);  
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&f),sizeof(double));
    IN.read((char*)(&m),sizeof(double));
    IN.read((char*)(&gene_size),sizeof(int));
    
    sprintf(mybuffer,"\t%d\t%.9f\t%e\t", cell_index, f, m);
    OUT << mybuffer << std::endl;

    //read gene info
    double e(0);
    double c(0);
    double dg(0);
    double eff(0);

    int gene_nid(0);
    int Ns(0);
    int Na(0);

    bool is_mob(false);

    for (int j=0; j<gene_size; ++j){
        std::string DNAsequence;

        IN.read((char*)(&gene_nid),sizeof(int));   
        IN.read((char*)(&e),sizeof(double));
        IN.read((char*)(&c),sizeof(double));
        IN.read((char*)(&eff),sizeof(double));
        IN.read((char*)(&dg),sizeof(double));
        IN.read((char*)(&f),sizeof(double));

        IN.read((char*)(&Na),sizeof(int));
        IN.read((char*)(&Ns),sizeof(int));
        IN.read((char*)(&is_mob),sizeof(bool));
        
        //read DNA sequence
        int nl(0);
        IN.read((char*)&nl, sizeof(int));
        std::vector<char> buff(nl);
        IN.read(&buff[0], nl);  
        DNAsequence.assign(buff.begin(), buff.end());

        sprintf(mybuffer,"%d\tG\t%e\t%.8f\t%.9f\t%d\t%d\t",j, c, dg, f, Na, Ns);
        OUT << mybuffer << std::endl;
        if (DNA) OUT << DNAsequence << std::endl;
        else OUT << GetProtFromNuc(DNAsequence) << std::endl;
    } 
}

// Counts the number of dissimilar characters in 2 strings
int StringDiff(const std::string& A, const std::string& B)
{
    unsigned int L = A.length();
    assert(L == B.length());
    int ctr = 0;  
    for (unsigned int i =0; i<L; ++i){
        if ( A.at(i) != B.at(i))
          ctr+=1;
    } 
    return ctr;
}

std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first){
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

// checks if given directory exists
bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & _S_IFDIR) != 0;
#else 
    struct stat info;
    if (stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
#endif
}

// creates directory from specified path
bool makePath(const std::string& path)
{
#if defined(_WIN32)
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
        return true;

    switch (errno){
    case ENOENT:
        // parent didn't exist, try to create it
        {
            std::size_t pos = path.find_last_of('/');
            if (pos == std::string::npos)
#if defined(_WIN32)
                pos = path.find_last_of('\\');
            if (pos == std::string::npos)
#endif
                return false;
            if (!makePath( path.substr(0, pos) ))
                return false;
        }
        // try to create again
#if defined(_WIN32)
        return 0 == _mkdir(path.c_str());
#else 
        return 0 == mkdir(path.c_str(), mode);
#endif
    case EEXIST:
        // done
        return isDirExist(path);

    default:
        return false;
    }
}

void printProgress (double progress)
{
    std::cout << '[';
    int pos = PBWidth * progress;
    for (int i = 0; i < PBWidth; ++i){
        if (i < pos)
          std::cout << '=';
        else if (i == pos)
          std::cout << '+';
        else
          std::cout << ' ';
    }

    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

std::map<std::string, double> get_codon_usage_map(const int the_cell_id, const std::string& the_cell_file_workspace){
    int nb_codons_input = 0; //initialization : should be equal to 61 at the end!
    bool b_f_does_exit = false;
    std::stringstream sstm;
    std::string line;
    sstm << the_cell_id << ".cell";
    sprintf(buffer,"%s/%s",the_cell_file_workspace.c_str(),sstm.str().c_str());
    std::cout << buffer <<std::endl;
    if (FILE *file = fopen(buffer, "r")) {
        b_f_does_exit = true;
        fclose(file);
    }
    //verify if file exists
    if (b_f_does_exit){
        // read codon usage data and build map 
        std::ifstream the_cell_ifs (buffer, std::ifstream::in);
        std::map<std::string, double> out_codon_map;        
        out_codon_map = {{"UUU", 0},{"UCU", 0},{"UAU", 0},{"UGU", 0},{"UUC", 0},{"UCC", 0},{"UAC", 0},{"UGC", 0},{"UUA", 0},{"UCA", 0},{"UUG", 0},{"UCG", 0},{"UGG", 0},{"CUU", 0},{"CCU", 0},{"CAU", 0},{"CGU", 0},{"CUC", 0},{"CCC", 0},{"CAC", 0},{"CGC", 0},{"CUA", 0},{"CCA", 0},{"CAA", 0},{"CGA", 0},{"CUG", 0},{"CCG", 0},{"CAG", 0},{"CGG", 0},{"AUU", 0},{"ACU", 0},{"AAU", 0},{"AGU", 0},{"AUC", 0},{"ACC", 0},{"AAC", 0},{"AGC", 0},{"AUA", 0},{"ACA", 0},{"AAA", 0},{"AGA", 0},{"AUG", 0},{"ACG", 0},{"AAG", 0},{"AGG", 0},{"GUU", 0},{"GCU", 0},{"GAU", 0},{"GGU", 0},{"GUC", 0},{"GCC", 0},{"GAC", 0},{"GGC", 0},{"GUA", 0},{"GCA", 0},{"GAA", 0},{"GGA", 0},{"GUG", 0},{"GCG", 0},{"GAG", 0},{"GGG", 0}};
        
        while (!the_cell_ifs.eof()) {
            getline(the_cell_ifs, line);
            std::string word;
            std::istringstream iss(line, std::istringstream:: in );
            iss >> word;
            if (word == "UUU"){
                iss >> word;
                out_codon_map["UUU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UCU"){
                iss >> word;
                out_codon_map["UCU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UAU"){
                iss >> word;
                out_codon_map["UAU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UGU"){
                iss >> word;
                out_codon_map["UGU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UUC"){
                iss >> word;
                out_codon_map["UUC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UCC"){
                iss >> word;
                out_codon_map["UCC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UAC"){
                iss >> word;
                out_codon_map["UAC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UGC"){
                iss >> word;
                out_codon_map["UGC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UUA"){
                iss >> word;
                out_codon_map["UUA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UCA"){
                iss >> word;
                out_codon_map["UCA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UUG"){
                iss >> word;
                out_codon_map["UUG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UCG"){
                iss >> word;
                out_codon_map["UCG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "UGG"){
                iss >> word;
                out_codon_map["UGG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CUU"){
                iss >> word;
                out_codon_map["CUU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CCU"){
                iss >> word;
                out_codon_map["CCU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CAU"){
                iss >> word;
                out_codon_map["CAU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CGU"){
                iss >> word;
                out_codon_map["CGU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CUC"){
                iss >> word;
                out_codon_map["CUC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CCC"){
                iss >> word;
                out_codon_map["CCC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CAC"){
                iss >> word;
                out_codon_map["CAC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CGC"){
                iss >> word;
                out_codon_map["CGC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CUA"){
                iss >> word;
                out_codon_map["CUA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CCA"){
                iss >> word;
                out_codon_map["CCA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CAA"){
                iss >> word;
                out_codon_map["CAA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CGA"){
                iss >> word;
                out_codon_map["CGA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CUG"){
                iss >> word;
                out_codon_map["CUG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CCG"){
                iss >> word;
                out_codon_map["CCG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CAG"){
                iss >> word;
                out_codon_map["CAG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "CGG"){
                iss >> word;
                out_codon_map["CGG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AUU"){
                iss >> word;
                out_codon_map["AUU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "ACU"){
                iss >> word;
                out_codon_map["ACU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AAU"){
                iss >> word;
                out_codon_map["AAU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AGU"){
                iss >> word;
                out_codon_map["AGU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AUC"){
                iss >> word;
                out_codon_map["AUC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "ACC"){
                iss >> word;
                out_codon_map["ACC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AAC"){
                iss >> word;
                out_codon_map["AAC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AGC"){
                iss >> word;
                out_codon_map["AGC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AUA"){
                iss >> word;
                out_codon_map["AUA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "ACA"){
                iss >> word;
                out_codon_map["ACA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AAA"){
                iss >> word;
                out_codon_map["AAA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AGA"){
                iss >> word;
                out_codon_map["AGA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AUG"){
                iss >> word;
                out_codon_map["AUG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "ACG"){
                iss >> word;
                out_codon_map["ACG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AAG"){
                iss >> word;
                out_codon_map["AAG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "AGG"){
                iss >> word;
                out_codon_map["AGG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GUU"){
                iss >> word;
                out_codon_map["GUU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GCU"){
                iss >> word;
                out_codon_map["GCU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GAU"){
                iss >> word;
                out_codon_map["GAU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GGU"){
                iss >> word;
                out_codon_map["GGU"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GUC"){
                iss >> word;
                out_codon_map["GUC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GCC"){
                iss >> word;
                out_codon_map["GCC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GAC"){
                iss >> word;
                out_codon_map["GAC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GGC"){
                iss >> word;
                out_codon_map["GGC"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GUA"){
                iss >> word;
                out_codon_map["GUA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GCA"){
                iss >> word;
                out_codon_map["GCA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GAA"){
                iss >> word;
                out_codon_map["GAA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GGA"){
                iss >> word;
                out_codon_map["GGA"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GUG"){
                iss >> word;
                out_codon_map["GUG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GCG"){
                iss >> word;
                out_codon_map["GCG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GAG"){
                iss >> word;
                out_codon_map["GAG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            } else if (word == "GGG"){
                iss >> word;
                out_codon_map["GGG"] = atof(word.c_str());
                nb_codons_input = nb_codons_input +1;
            }
        }
        the_cell_ifs.close();

        //make sure the user gave the 61 required codon usage data
        if (nb_codons_input!=61){
            std::cerr << "Error : You need to input codon usage data for all 61 possible sense codons. Please look at the file " << the_cell_id << ".cell!" << std::endl;
            exit(2);
        }

        return out_codon_map;
    }
    else{
        // error file does not exist -> throw exception
        std::string msg_err = "Error! The following .cell file does not exist : "+buffer;
        throw std::runtime_error(msg_err);
    }

}



void normalize_CUF(const int& p_species_id, std::map<int, std::map<std::string, double>>& p_map_of_codon_usag_map){
    //Copy the original codon usage frequencies map; 
    std::map<std::string, double> old_codon_usage_freq_map(p_map_of_codon_usag_map[p_species_id]);
    //Calculate normalized RSCUs of synonymous codons (Equation 2 Sharp & Li 1987)
    for (auto it = p_map_of_codon_usag_map[p_species_id].begin();it!=p_map_of_codon_usag_map[p_species_id].end();++it){
        std::string current_codon = it->first;
        std::vector<double> vec_cufs_syn_codons_of_current_codon;
        if (current_codon == "UUU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UUU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UCU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UCU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UAU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UAU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UAC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UAU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UGU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UGC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UGU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UUC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UUC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UCC"){
                    //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UCC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UAC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UAC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UAU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UAC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UGC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UGU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UGC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UUA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UUA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UCA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UCA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UUG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UUG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UCG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["UCG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "UGG"){
            //UGG has no synonymous codons
            it->second = 1;
        } else if (current_codon == "CUU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CUU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CCU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CCU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CAU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CAU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CGU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CGU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CUC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CUC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CCC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CCC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CAC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CAC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CGC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CGC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CUA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUG"]); 
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CUA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CCA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CCA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CAA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CAA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CGA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CGA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CUG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CUA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CUG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CCG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CCA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CCG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CAG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CAA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CAG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "CGG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["CGG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AUU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AUU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "ACU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["ACU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AAU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAC"]);
            it->second = (old_codon_usage_freq_map["AAU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AGU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AGU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AUC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AUC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "ACC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["ACC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AAC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AAC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AGC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["UCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AGC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AUA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AUC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AUA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "ACA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["ACA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AAA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AAA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AGA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AGA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AUG"){
           //AUG has no synonymous codons
           it->second = 1;
        } else if (current_codon == "ACG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["ACA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["ACG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AAG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AAA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AAG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "AGG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["AGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["CGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["AGG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GUU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GUU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GCU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GCU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GAU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAC"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GAU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GGU"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GGU"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GUC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GUC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GCC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GCC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GAC"){
             //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAU"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GAC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GGC"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GGC"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GUA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GUA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GCA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GCA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GAA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GAA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GGA"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGA"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGG"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GGA"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GUG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GUA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GUG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GCG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GCA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GCG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GAG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GAA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GAG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        } else if (current_codon == "GGG"){
            //add all the CUFs of the synonymous codons of the current codon to vec_cufs_syn_codons_of_current_codon
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGG"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGU"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGC"]);
            vec_cufs_syn_codons_of_current_codon.push_back(old_codon_usage_freq_map["GGA"]);
            //NORMALIZE CUF
            it->second = (old_codon_usage_freq_map["GGG"])/(*max_element(vec_cufs_syn_codons_of_current_codon.begin(),vec_cufs_syn_codons_of_current_codon.end()));
        }
    }
    
    std::cout << "#########Original map########" << std::endl;//to remove after test
    for(auto const &elem : old_codon_usage_freq_map){//to remove after test
        std::cout << elem.first << " " << elem.second << ";";//to remove after test
    }//to remove after test
    std::cout << std::endl;//to remove after test
    std::cout << "#########New map########" << std::endl;//to remove after test
    for(auto const &elem : p_map_of_codon_usag_map[p_species_id]){//to remove after test
        std::cout << elem.first << " " << elem.second << ";";//to remove after test
    }//to remove after test
    std::cout << std::endl;//to remove after test
}     

//returns a vector of species ID based on the .cell files
std::vector<int> get_unique_species(){
    std::vector<int> vec_out_species_id;
    std::map<int,int> ref_subpops_abundance_map = get_ref_map_subpops_abundance(); 
    for (auto const &it : ref_subpops_abundance_map.begin(); it != ref_subpops_abundance_map.end();++it){
        vec_out_species_id.push_back(it->first);
    }

    return vec_out_species_id;
}

std::map<int,int> get_ref_map_subpops_abundance(){
    bool b_f_does_exit = false;
    std::string line;
    buffer = "files/start/population.dat"
    if (FILE *file = fopen(buffer, "r")) {
        b_f_does_exit = true;
        fclose(file);
    }
    //verify if file exists
    if (b_f_does_exit){
        // read subpopulation abundance data and add it to map 
        std::ifstream the_pop_ifs (buffer, std::ifstream::in);
        std::map<int,int> out_subpops_abundance_map;        
        while (!the_pop_ifs.eof()) {
            getline(the_pop_ifs, line);
            std::string word;
            std::istringstream iss(line, std::istringstream:: in );
            iss >> word;
            if (word == "C"){
                iss >> word; //read abundance of current cell type / subpopulation
                int count_current_subpop = atoi(word.c_str());
                iss >> word; //read .cell file path
                //convert .cell file path to cell ID
                std::string path_cell_file = word;
                std::string::size_type i = path_cell_file.find(".cell");
                if (i != std::string::npos){//erase .cell from string
                    path_cell_file.erase(i, 5); //erase .cell from string
                }
                std::string::size_type i = path_cell_file.find("files/start/");
                if (i != std::string::npos){//erase files/start/ from string
                    path_cell_file.erase(i, 12); //erase files/start/ from string 
                }
                int id_current_subpop = atoi(path_cell_file.c_str());
                out_subpops_abundance_map[id_current_subpop] = count_current_subpop;
            }
    }else{
        // error file does not exist -> throw exception
        std::string msg_err = "Error! The following .dat file does not exist : "+buffer;
        throw std::runtime_error(msg_err);
    }
    return out_subpops_abundance_map;
}

