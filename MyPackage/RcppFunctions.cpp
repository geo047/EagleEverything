// Rcpp code for:
//   1.  the calculation of M %*% t(M)
//   2.  checking the genotypes of M are 0,1, & 2


// Author:   Andrew W. George
// Purpose: to calculate M %*% t(M) when M may not fit into memory
// Outline: 
//          1. read data from ascii file. Assuming genotypes are 0,1,2. 
//             An error is produced if other values are found. 
//          2. convert genotypes into their binary values.
//          3. pack binary values into unsigned long int (could be 32 bits or 64 bits 
//             depending upon the system.
//          4. write packed longs to a new file in binary format.
//          5. read blocks of binarys to form submatrices of M.
//          6. Perform M %*% t(M) as a block multiplication. 
//
// Inputs from R:
//     1. file name of ASCII file with genotypes (not handling marker names yet, or different formats like cvs)
//     2. max memory size in bytes



#include <RcppEigen.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <bitset>
#include <string>
using namespace std;
using namespace Rcpp;
using Eigen::MatrixXi;
using Eigen::MatrixXd;  
using Eigen::Lower;
using Eigen::Map;   // maps rather than copies


// #ifdef _OPENMP
//#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
// #endif


const size_t bits_in_double = std::numeric_limits<long double>::digits;
const size_t bits_in_ulong = std::numeric_limits<unsigned long int>::digits;
const size_t bits_in_int = std::numeric_limits<unsigned int>::digits;



// internal function to remove a row from a dynamic matrix
void removeRow(MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

// internal function to remove a column from a dynamic matrix
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}









// check genotypes in file are 0,1, & 2 only
// [[Rcpp::export]]
void  checkGenotypes(CharacterVector f_name,
                     int AA,
                     int AB,
                     int BB)
{
std::string 
     fname = Rcpp::as<std::string>(f_name),
     line;
 int 
    genoval=-1,
    linenum=0;

 ostringstream 
      os;



 // open file and check for its existence. 
 std::ifstream fileIN(fname.c_str());
 if(!fileIN.good()) {
      os << "\n\nERROR: Could not open  " << fname << "\n\n" << std::endl;
      Rcpp::stop(os.str() );
 }

 // Determine number of rows in file
 Rprintf("\n\n Checking genotype file for incorrect genotypes ... ");
 while(fileIN.good()){
      while(getline(fileIN, line)){
           Rprintf(".");
           linenum++;
           istringstream streamA(line);
           while(streamA >> genoval){
           //  if(genoval !=0 & genoval !=1 & genoval != 2){
             if(genoval !=AA & genoval !=AB & genoval != BB){
               os << "\n\nERROR: File " << fname << " contains genotypes other than " << AA << "," << 
                        AB << ", and " << BB << " For example genotype " << genoval << " has been found on line.  " << linenum << "\n\n";
               Rcpp::stop(os.str() );
             }
           }
      }
 }



}





//get number of rows and columns in genotype file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// [[Rcpp::export]]
std::vector<long>   getRowColumn(std::string fname)
{
  // Purpose:  to open the ascii file where the genotypes are kept
  //           an error will be produced if the file cannot be found.
  //           I am assuming no row or column names
 int 
   genoval;

 std::string
   line;

 ostringstream 
      os;


 std::vector<long> dimen(2)  ;  // dim[0] row number
                               // dim[1] col number 

 // open file and check for its existence. 
 std::ifstream fileIN(fname.c_str());
 if(!fileIN.good()) {
      os << "\n\n ERROR: Could not open  " << fname << "\n\n" << std::endl;
      Rcpp::stop(os.str() );
 }


 // Determine number of rows in file
 while(fileIN.good()){
      while(getline(fileIN, line)){
         dimen[0]++;
      }
 }


 // Determine number of columns in file
 fileIN.clear(); // returns to beginning of line
 fileIN.seekg(0, ios::beg);

 getline(fileIN, line);
 istringstream streamA(line);
 while(streamA >> genoval){
          dimen[1]++;
 }
 fileIN.close();
 if(fileIN.bad())
 {
    os << "\n\nERROR:  There was a problem with reading the genotype file - possibly strange ASCII characters.\n\n";
    Rcpp::stop(os.str() );
 } else {
     return dimen;
 }

}





// recode ascii as packed binary file
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void  CreatePackedBinary(std::string fname, std::string binfname, std::vector<long> dims,
                         int AA, 
                         int AB, 
                         int BB)
{
int 
   indx_packed = 0,
   indx_packed_long_vec = 0,
   colindx = 0, 
   n_of_long = dims[1]/(bits_in_ulong/2),
   n_extra=0,
   n_total=0;
short 
    rowvec[dims[1]]; // holds entire row worth of genotypes from ascii file

std::string
   line;


 ostringstream 
      os;


std::bitset <bits_in_ulong>
     packed(0);


// check that number of bits to long is even
if ( (bits_in_ulong % 2)!=0){
  os << "\n\nERROR: Number of bits to a ulong is not even.\n\n";
  Rcpp::stop(os.str() );
}

// check if number of columns in file will fill n longs completely
if(dims[1] % (bits_in_ulong/2) != 0){
   n_extra = 1 ;  // an extra long is required to the extra columns
}

// number of longs needed to store a complete row of the ascii file
// where genotypes are being packed into 2 bits. 
n_total = n_of_long + n_extra;

// Vector that is to be packed 
std::vector<unsigned long int> packed_long_vec (n_total);



// open ascii genotype  file
std::ifstream fileIN(fname.c_str());
if(!fileIN.good()) {
  os << "\n\nERROR: ASCii genotype file could not be opened with filename  " << fname << "\n\n" << std::endl;
  Rcpp::stop(os.str() );
}

// open binary file that is to hold packed genotype data
std::ofstream fileOUT(binfname.c_str(), ios::binary );




while(getline(fileIN, line))
{
  istringstream streamA(line);
  indx_packed = 0;
   indx_packed_long_vec = 0;
  packed.reset();
 // Here, BB is coded into 2 when bit packed, 
 //       AB is coded into 1, 
 //       AA is coded into 0. 
  for(int i=0; i< dims[1] ; i++){
     streamA >> rowvec[i];
     if(rowvec[i] == BB){
          packed[indx_packed*2+1] = 1;
          packed[indx_packed*2] = 0;
     } else if (rowvec[i] == AB) {
          packed[indx_packed*2+1] = 0;
          packed[indx_packed*2] = 1;
     } else if (rowvec[i] == AA) {
          packed[indx_packed*2+1] = 0;
          packed[indx_packed*2] = 0;
     } else {
          os << "\n\nERROR: Genotype file contains genotypes that are not 0,1, or 2. For example " << rowvec[i] << "\n\n";
          Rcpp::stop(os.str() );
     } 

     if(  ( ((indx_packed+1)  % ( bits_in_ulong/2))==0) | (dims[1]-1) == i ) { 
        indx_packed = 0;
        packed_long_vec[indx_packed_long_vec]  =  packed.to_ulong();
        indx_packed_long_vec++;
        packed.reset();  // set bits back to 0
     } else  {
       indx_packed++;
    } 



  }


    // want to begin with a fresh long when we read in a new line
    // writing binary values to disk.
    fileOUT.write((char *)(&packed_long_vec[0]), packed_long_vec.size() * sizeof(unsigned long int));
  }


// close files
fileIN.close();
fileOUT.close();


}


Eigen::MatrixXi  ReadBlock(std::string binfname, 
                           int start_row,
                           int numcols,
                           int numrows_in_block) 

{
 // reads in packed data from binary file of longs
 // to form M Eign interger matrix 

int 
  coli = 0, 
  rowi = 0;

long 
    igeno;

double 
     geno;

Eigen::MatrixXi
      M(numrows_in_block, numcols) ;

const size_t bits_in_ulong = std::numeric_limits<unsigned long int>::digits;
std::bitset <bits_in_ulong> 
    packed(0), 
    geno_bitset(0),
    mask(3);

// Open binary file
   std::ifstream fileBIN(binfname.c_str(), ios::in | ios::binary );



// Determine size (in bytes) of block
   int number_of_longs_in_row  =  (int) numcols/(bits_in_ulong/2);
   if (numcols % (bits_in_ulong/2) !=0) 
        number_of_longs_in_row++;
    int size_in_bytes_of_block  = numrows_in_block * number_of_longs_in_row  * bits_in_ulong/8;
   fileBIN.seekg(start_row*number_of_longs_in_row*bits_in_ulong/8, std::ios_base::beg);
// create float vector to store block of binary results. This is going 
// to be a chunk of data that is subrows x colnum in size To Do. 
   std::vector<unsigned long int> v(size_in_bytes_of_block/(bits_in_ulong/8) );

// Load the data
fileBIN.read((char*)&v[0] , size_in_bytes_of_block );  // reads a block of bytes of size. 


// Close the binary file
   fileBIN.close();

// Convert integers into bitsets
packed.reset(); //to initialize bitset
Rprintf(" Beginning conversion of binary data into an integer matrix... \n");
for(int i=0;i < v.size(); i++)
{
  std::bitset <bits_in_ulong> packed(v[i]);
  for(int j=0; j< (bits_in_ulong/2); j++)
  {
     geno_bitset.reset();
     geno_bitset =  ((packed & (mask << (j*2)))) >> (j*2);
     igeno = geno_bitset.to_ulong();
     //if(igeno>1)
     //     igeno = -1;
     // it's igeno - 1 so that 0,1,2 map onto -1, 0, 1
     M(rowi, coli) = (double) igeno - 1; // converted to double but okay, its safe because igeno always small 
     coli++;
     if(coli == numcols)
     {
        coli=0;
        rowi++;
        break; //exit the for loop
     }
  }
}


 return M;

}




Eigen::MatrixXi  createMmat(std::string binfname, std::vector<int> dims)
{
 // reads in packed data from binary file of longs
 // to form M arma matrix of doubles for further analysis. 

int 
  coli = 0, 
  rowi = 0;

long 
    igeno;

double 
     geno;

Eigen::MatrixXi
      M(dims[0], dims[1]) ;

const size_t bits_in_ulong = std::numeric_limits<unsigned long int>::digits;
std::bitset <bits_in_ulong> 
    packed(0), 
    geno_bitset(0),
    mask(3);


// Read in binary file
   std::ifstream fileBIN(binfname.c_str(), ios::in | ios::binary );

// Determine the file length
   fileBIN.seekg(0, std::ios_base::end);
   std::size_t size=fileBIN.tellg();   // size of binary file in bytes
   fileBIN.seekg(0, std::ios_base::beg);
// create float vector to store block of binary results. This is going 
// to be a chunk of data that is subrows x colnum in size To Do. 
   std::vector<unsigned long int> v(size/(bits_in_ulong/8) );

// Load the data
fileBIN.read((char*)&v[0] , size );  // reads a block of bytes of size. 


// Close the binary file
   fileBIN.close();

// Convert integers into bitsets
packed.reset(); //to initialize bitset
Rprintf( " Beginning conversion of binary data into an integer matrix." );
for(int i=0;i < v.size(); i++)
{
  std::bitset <bits_in_ulong> packed(v[i]);
  for(int j=0; j< (bits_in_ulong/2); j++)
  {
     geno_bitset.reset();
     geno_bitset =  ((packed & (mask << (j*2)))) >> (j*2);
     igeno = geno_bitset.to_ulong();
     // if(igeno>1)
      //     igeno = -1;
      // it's igeno -1 so that 0,1,2 maps onto -1, 0, 1
     M(rowi, coli) = (double) igeno - 1; // converted to double but okay, its safe because igeno always small 
     coli++;
     if(coli == dims[1])
     {
        coli=0;
        rowi++;
        break; //exit the for loop
     }
  }
}



 return M;


}















//-----------------------------------------
//               Main 
//-----------------------------------------





// [[Rcpp::export]]
void  createMt_rcpp(CharacterVector f_name, CharacterVector f_name_bin, 
                              int AA, 
                              int AB, 
                              int BB,
                              double  max_memory_in_Gbytes,  std::vector <long> dims )
{

std::string
   line;

const size_t bits_in_ulong = std::numeric_limits<unsigned long int>::digits;

ostringstream
      os;



int
   genoval;

std::string
     fname = Rcpp::as<std::string>(f_name),
     fnamebin = Rcpp::as<std::string>(f_name_bin);

short 
  rowvec[dims[1]];

// Calculate number of packed longs ints needed for a single column of ASCII data
long
  n_extra = 0,
  n_total = 0,
  n_of_long = dims[0]/(bits_in_ulong/2);
  if(  (dims[0] % (bits_in_ulong/2)) != 0)
     n_extra = 1;

  n_total = n_of_long + n_extra;


// Calculate number of columns that can be read in as a block with XGb of
// memory. 
   // Amount of memory (in bytes) needed to store a single column of data
   // in packed binary form. 

  double mem_bytes = n_total * (bits_in_ulong/8);



  // calculate number of columns that can be read into XGb
  int n_of_cols_to_be_read = (max_memory_in_Gbytes * 1000000/mem_bytes) * (bits_in_ulong/2);


 // open binary output file
std::ofstream fileOUTbin(fnamebin.c_str(), ios::binary );

//-----------------------------------------------------------------------
//  Two situations
//   1.  memory X is sufficient to read all data into memory and transpose
//   2.  memory X is insufficient to read all data into memory. 
//------------------------------------------------------------------------


if(n_of_cols_to_be_read > dims[1]){


// Situation 1
//-------------

  // want packed-block object that is a matrix of bitset values. 
  std::vector< std::vector < std::bitset <bits_in_ulong> > > 
        packed_block( dims[1] , std::vector<bitset <bits_in_ulong> > (n_total, 0) );


 // initialize the packed 2D array to all 0's.
  for(int i=0; i < dims[1]; i++)
    for(int j=0; j < n_total; j++)
        packed_block[i][j].reset();

int
     indx_packed_within = 0,
     indx_packed_across = 0;

 // open ASCII file and check for its existence. 
 std::ifstream fileIN(fname.c_str());
 if(!fileIN.good()) {
      os << "\n\nERROR: Could not open  " << fname << "\n\n" << std::endl;
      Rcpp::stop(os.str() );

 }

 while(getline(fileIN, line))
 {

   // read a line of data from ASCII file
   istringstream streamA(line);



   for(int i=0; i < dims[1]; i++){
     streamA >> rowvec[i];
 
   // Here, BB is coded as 2 when bit packed,
   //       AB is coded as 1, 
   //       AA is coded as 0. 
   if(rowvec[i] == BB){
      packed_block[i][indx_packed_across][indx_packed_within*2 + 1] = 1;
      packed_block[i][indx_packed_across][indx_packed_within*2 ] = 0;
   } else if (rowvec[i] == AB) {
      packed_block[i][indx_packed_across][indx_packed_within*2 + 1] = 0;
      packed_block[i][indx_packed_across][indx_packed_within*2 ] = 1;
  } else if (rowvec[i] == AA) {
      packed_block[i][indx_packed_across][indx_packed_within*2 + 1] = 0;
      packed_block[i][indx_packed_across][indx_packed_within*2 ] = 0;
  } else {
      os << "Genotype file contains genotypes that are not " << AA << "," << AB << ", or " << BB << " For example " << rowvec[i] << "\n\n";
      Rcpp::stop(os.str() );
  }
    // Rcpp::Rcout <<  i << " " << rowvec[i] << " " << indx_packed_across << " " << packed_block[i][indx_packed_across] << std::endl;
  }  // end for int i


  if(  ( ((indx_packed_within + 1)  % ( bits_in_ulong/2))==0)   ) {
        indx_packed_within = 0;
        indx_packed_across++;
  } else  {
       indx_packed_within++;
  }



 }  // end while

 // write packed binary file to disc
  for(int i=0; i < dims[1]; i++){
    for(int j=0; j < n_total; j++){
           fileOUTbin.write((char *)(&packed_block[i][j]), sizeof(unsigned long int));
   }}


} else {
//  Situation 2 
//  Block approach needed due to lack of memory

 Rcpp::Rcout << " A block transpose is being performed due to lack of memory.  "  << std::endl;
 Rcpp::Rcout << " Memory parameter workingmemGb is set to " << max_memory_in_Gbytes << "Gbytes" << std::endl;
 Rcpp::Rcout << " If possible, increase workingmemGb parameter. " << std::endl;

  // Calculate number of blocks needed
  int n_blocks = dims[1]/n_of_cols_to_be_read;
  if (dims[1] % n_of_cols_to_be_read != 0)
      n_blocks++;

     Rcpp::Rcout  << " Block Tranpose of ASCII genotype file beginning ... " << std::endl;

  // Block read and transpose - requires n_blocks passes through the 
  // ASCII input file which could be slow if file is large and memory low
   for(int b=0; b < n_blocks; b++){
  //for(int b=1; b < 2; b++){
     Rcpp::Rcout << " Processing block ... " << b << " of a total number of blocks of " << n_blocks << std::endl;


     // want packed-block object that is a matrix of bitset values. 
     std::vector< std::vector < std::bitset <bits_in_ulong> > >
           packed_block( n_of_cols_to_be_read , std::vector<bitset <bits_in_ulong> > (n_total, 0) );


    // initialize the packed 2D array to all 0's.
     for(int i=0; i < n_of_cols_to_be_read ; i++)
       for(int j=0; j < n_total; j++)
           packed_block[i][j].reset();

   int
        indx_packed_within = 0,
        indx_packed_across = 0;


    int
        start_val = b * n_of_cols_to_be_read,
        end_val   = (b+1) * n_of_cols_to_be_read;


     if (end_val > dims[1])
        end_val = dims[1];




    // open ASCII file and check for its existence. 
    std::ifstream fileIN(fname.c_str());
    if(!fileIN.good()) {
      os << "ERROR: Could not open  " << fname << std::endl;
      Rcpp::stop(os.str() );
     }
    int counter = 0;
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << std::endl;
    while(getline(fileIN, line))
    {
       
      // read a line of data from ASCII file
      istringstream streamA(line);




      for(int i=0; i < dims[1] ; i++){
        streamA >> rowvec[i];

       if(i>= start_val & i <  end_val){
        //  if(counter==0)   Rcpp::Rcout << rowvec[i] << " " ;
          int iindx = i % n_of_cols_to_be_read; // converts it back to an 
                                                // index between 0 and n_of_cols_to_be_read
         if(rowvec[i] == BB){
            packed_block[iindx][indx_packed_across][indx_packed_within*2 + 1] = 1;
            packed_block[iindx][indx_packed_across][indx_packed_within*2 ] = 0;
         } else if (rowvec[i] == AB) {
            packed_block[iindx][indx_packed_across][indx_packed_within*2 + 1] = 0;
            packed_block[iindx][indx_packed_across][indx_packed_within*2 ] = 1;
        } else if (rowvec[i] == AA) {
            packed_block[iindx][indx_packed_across][indx_packed_within*2 + 1] = 0;
            packed_block[iindx][indx_packed_across][indx_packed_within*2 ] = 0;
        } else {
            os  << "Genotype file contains genotypes that are not 0,1, or 2. For example " << rowvec[i] << "\n\n";
            Rcpp::stop(os.str() );
        }
       } // end if
     }  // end for int i


     if(  ( ((indx_packed_within + 1)  % ( bits_in_ulong/2))==0)   ) {
           indx_packed_within = 0;
           indx_packed_across++;
     } else  {
          indx_packed_within++;
     }


      counter++;

    }     // end while
   // close ASCII file because I have read the entire file
   fileIN.close();

 // write packed binary file to disc
  for(int i=0; i < n_of_cols_to_be_read; i++){
    for(int j=0; j < n_total; j++){
           fileOUTbin.write((char *)(&packed_block[i][j]), sizeof(unsigned long int));
   }}



  } // end for block




}




// close files
fileOUTbin.close();

//
//
// // Returning transpose of matrix to R as a check
// MatrixXi
//   genoMat;
//
//  genoMat = ReadBlock(fnamebin, 0, dims[0], dims[1]);
// // genoMat = ReadBlock(fnamebin, 0, dims[0], 17857);
// 
//   Rcpp::Rcout << genoMat(0, 0) << " " << genoMat(1, 0) << " " << genoMat(2, 0) << std::endl;
//   return genoMat;
}





//--------------------------------------------
// Calculation of transformed blup a values
//--------------------------------------------
// [[Rcpp::export]]
MatrixXd calculate_reduced_a_rcpp ( CharacterVector f_name_bin, double varG, 
                                           Map<MatrixXd> P,
                                           Map<MatrixXd>  y,
                                           double max_memory_in_Gbytes,  
                                           std::vector <long> dims,
                                           Rcpp::NumericVector  selected_loci)
{
  // function to calculate the BLUPs for the dimension reduced model. 
  // It is being performed in Rcpp because it makes use of Mt. 
  // Args
  // f_name_bin    path + file name of Mt.bin
  // varG          variance of polygenic component
  // P             calculate in R
  // y             response/trait  but read in as a row matrix
  // max_memory_in_Gbytes  working memory in Gbytes
  // dims          dimension (row, column), of M.

std::string
     fnamebin = Rcpp::as<std::string>(f_name_bin);

Eigen::MatrixXd
      ar(dims[1],1);  // column vector



const size_t bits_in_double = std::numeric_limits<double>::digits;


   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt))
double mem_bytes_needed =   ( dims[0]*dims[1] + dims[0]*dims[0] + dims[0] ) *  ( bits_in_double/(8.0 * 1000000000));


Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
Rprintf("Max memory (Gbytes) available is: %f \n", max_memory_in_Gbytes);

if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory

  Eigen::MatrixXi
                  Mt;


  if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(int ii=0; ii < selected_loci.size() ; ii++)
          Mt.row(selected_loci(ii)).setZero();
   }

   Mt = ReadBlock(fnamebin, 0, dims[0], dims[1]);
   ar  =    varG * Mt.cast<double>() *  P   * y ;

   Rcpp::Rcout << "testing ... " << std::endl;
   Rcpp::Rcout << ar.rows() << std::endl;
   Rcpp::Rcout << ar.cols() << std::endl;
   Rcpp::Rcout << "ar(0,0) = " << ar(0,0) << std::endl;
   Rcpp::Rcout << "ar(1,0) = " << ar(1,0) << std::endl;
   Rcpp::Rcout << "ar(2,0) = " << ar(2,0) << std::endl;

} else {

      // calculation being processed in block form
      Rprintf(" Increasing maxmemGb would improve performance... \n");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algrebra but it is comes to the following
      int num_rows_in_block = (max_memory_in_Gbytes /  ( bits_in_double/(8.0 * 1000000000)) - dims[0] * dims[0] - dims[0])/dims[0] ;


      // blockwise multiplication

      // find out the number of blocks needed
      int num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;


      Rprintf(" Maximum memory has been set to %f Gb\n", max_memory_in_Gbytes);
      Rprintf(" Block multiplication necessary. \n");
      Rprintf(" Number of blocks needing block multiplication is ... % d \n", num_blocks);

      for(int i=0; i < num_blocks; i++){
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[1])
            num_rows_in_block1 = dims[1] - start_row1;

          Eigen::MatrixXi
                  Mt;

         Mt = ReadBlock(fnamebin, start_row1, dims[0], num_rows_in_block1) ;

         Eigen::MatrixXd
             ar_tmp;

         if(!R_IsNA(selected_loci(0))){
         // setting columns (or row when Mt) to 0
            for(int ii=0; ii < selected_loci.size() ; ii++)
            {
            // since we are now dealing with Mt, and blocking on columns, 
            // because columns are rows in Mt, then we have to be careful
            // that we do not select loci outside the block bounds. Also 
            // the values have to be adjusted based on the block number
                if(selected_loci(ii) >= start_row1 & selected_loci(ii) < start_row1 + num_rows_in_block1 )
                {   // selected loci index is in block 
                int block_selected_loci = selected_loci(ii) - start_row1;
                Mt.row(block_selected_loci).setZero();
                }
             }   
         }

         ar_tmp  =  varG * Mt.cast<double>() *  P  * y ;

          


         // assign block vector results to final vector (ar) of results
         int  counter = 0;
         for(int j=start_row1; j < start_row1 + num_rows_in_block1; j++){
              ar(j,0) = ar_tmp(counter,0);
              counter++;
         }

        Rcpp::Rcout << "block done ... " << std::endl;
      } // end for int




}  // end if mem_bytes_needed

  return(ar);

} // end function 













// ------------------------------------------------------
//    Calculation of untransformed BLUP a values 
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List   calculate_a_and_vara_rcpp(  CharacterVector f_name_bin,  
                                    Rcpp::NumericVector  selected_loci,
                                    Map<MatrixXd> inv_MMt_sqrt,
                                    Map<MatrixXd> dim_reduced_vara,
                                    double  max_memory_in_Gbytes,  
                                    std::vector <long> dims,
                                    Eigen::VectorXd  a  )
{
// Purpose: to calculate the untransformed BLUP (a) values from the 
//          dimension reduced BLUP value estimates. 
//          It is neccessary to have a block multiplication form of this function. 
//          Also, since the matrix multiplications are reliant upon the BLAS library, only 
//          double precision matrix multiplication is possible. This means, the Mt matrix must 
//          be converted into a douple precision matrix which has a large memory cost.  
// Note:
//      1. dims is the row, column dimension of the Mt matrix



std::string
     fnamebin = Rcpp::as<std::string>(f_name_bin);


Eigen::MatrixXd
      ans(dims[0],1),
      var_ans(dims[0],1);


const size_t bits_in_double = std::numeric_limits<double>::digits;


   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt))
Rcpp::Rcout << ( dims[0]   ) << std::endl;
Rcpp::Rcout << ( dims[1]   ) << std::endl;
Rcpp::Rcout << ((dims[1] * bits_in_double/(8.0 * 1000000000))) << std::endl;
double mem_bytes_needed =   ( dims[0]   +  dims[1]   + dims[1] + dims[1] + dims[1] ) *  (dims[1] * bits_in_double/(8.0 * 1000000000));


Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
Rprintf("Max memory (Gbytes) available is: %f \n", max_memory_in_Gbytes);

if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory


   Eigen::MatrixXi
                  Mt;

   Mt = ReadBlock(fnamebin, 0, dims[1], dims[0]);

   if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(int ii=0; ii < selected_loci.size() ; ii++)
          Mt.row(selected_loci(ii)).setZero();
   }

   // calculate untransformed BLUP values
   Rprintf(" calculating raw a values ... \n");
   ans =    Mt.cast<double>() *  inv_MMt_sqrt  * a ;


   // calculate untransformed variances of BLUP values
   Rprintf(" calculating raw var(a) values ... \n");
   Eigen::MatrixXd
      var_ans_tmp;

   var_ans_tmp =  Mt.cast<double>() *  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;

   for(int i=0; i< dims[0]; i++){
         var_ans(i,0) =   var_ans_tmp.row(i)   * ((Mt.row(i)).transpose()).cast<double>() ;
    }


} else {
      // calculation being processed in block form
      Rprintf(" Increasing maxmemGb would improve performance... \n");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algrebra but it is comes to the following
      int num_rows_in_block =  max_memory_in_Gbytes / ( dims[1] * bits_in_double/(8.0 * 1000000000)) - dims[1] - 1 ;



      // blockwise multiplication

      // find out the number of blocks needed
      int num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;


      Rprintf(" Maximum memory has been set to %f Gb\n", max_memory_in_Gbytes);
      Rprintf(" Block multiplication necessary. \n");
      Rprintf(" Number of blocks needing block multiplication is ... % d \n", num_blocks);

      for(int i=0; i < num_blocks; i++){
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[0])
            num_rows_in_block1 = dims[0] - start_row1;

          Eigen::MatrixXi
                  Mt;

         Mt = ReadBlock(fnamebin, start_row1, dims[1], num_rows_in_block1) ;

         Eigen::MatrixXd
             vt,
             ans_tmp,
             var_ans_tmp(num_rows_in_block1,1);

         if(!R_IsNA(selected_loci(0))){
         // setting columns (or row when Mt) to 0
            for(int ii=0; ii < selected_loci.size() ; ii++)
            {
            // since we are now dealing with Mt, and blocking on columns, 
            // because columns are rows in Mt, then we have to be careful
            // that we do not select loci outside the block bounds. Also 
            // the values have to be adjusted based on the block number
                if(selected_loci(ii) >= start_row1 & selected_loci(ii) < start_row1 + num_rows_in_block1 )
                {   // selected loci index is in block 
                int block_selected_loci = selected_loci(ii) - start_row1;
                Mt.row(block_selected_loci).setZero();
                }
             }   
         }

         ans_tmp  =  Mt.cast<double>() *  inv_MMt_sqrt  * a ;

          

         vt =  Mt.cast<double>() *  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;
         for(int j=0; j < num_rows_in_block1; j++)
              var_ans_tmp(j,0)  =   vt.row(j)  * ((Mt.row(j)).transpose()).cast<double>() ;




         // assign block vector results to final vector (ans) of results
         int  counter = 0;
         for(int j=start_row1; j < start_row1 + num_rows_in_block1; j++){
              ans(j,0) = ans_tmp(counter,0);
              var_ans(j,0) = var_ans_tmp(counter,0);
              counter++;
         }

        Rcpp::Rcout << "block done ... " << std::endl;
      } // end for int



}

  return Rcpp::List::create(Rcpp::Named("a")=ans,
                            Rcpp::Named("vara") = var_ans);


}



// [[Rcpp::export]]
void createM_rcpp(CharacterVector f_name, CharacterVector f_name_bin, 
                  int AA,
                  int AB, 
                  int BB,
                  double  max_memory_in_Gbytes,  std::vector <long> dims) 
{
  // Rcpp function to create binary packed file of ASCII marker genotype file.


std::string 
   line; 


ofstream
   fileOUT;

int 
   genoval,
   rown,
   coln;

std::string 
     fname = Rcpp::as<std::string>(f_name),
     fnamebin = Rcpp::as<std::string>(f_name_bin);



//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double 
  memory_needed_in_Gb =  (dims[0] *  dims[1] *   bits_in_int)/( (double) 8000000000) ;


//--------------------------------------
// Summary of Genotype File
//--------------------------------------

Rcpp::Rcout <<  " +-------------------------------------------------------------------+ " << std::endl;
Rcpp::Rcout <<  " |                  Summary of Genotype File                         | " << std::endl;
Rcpp::Rcout <<  " |                                                                   | " << std::endl;
Rcpp::Rcout <<  " | Number of rows (individuals) in file                " ;
Rcpp::Rcout.width(13); 
Rcpp::Rcout << right << dims[0] << " | " << std::endl;
Rcpp::Rcout <<  " | Number of columns (markers) in file                 " ;
Rcpp::Rcout.width(13);
Rcpp::Rcout  << right << dims[1] << " | " << std::endl;
Rcpp::Rcout <<  " | Amount of memory (in Gbytes) needed to                            | " << std::endl;
Rcpp::Rcout <<  " |    load entire file                                 ";
Rcpp::Rcout.precision(2);
Rcpp::Rcout.width(10);
Rcpp::Rcout  << right << memory_needed_in_Gb << " Gb | " << std::endl;
Rcpp::Rcout <<  " | Maximum memory (in Gbytes) has been set to          ";
Rcpp::Rcout.width(10);
Rcpp::Rcout << right << max_memory_in_Gbytes << " Gb | " << std::endl; 
Rcpp::Rcout <<  " |                                                                   | " << std::endl;
Rcpp::Rcout <<  " +-------------------------------------------------------------------+ " << std::endl;


//-------------------------------------------
// convert ascii file into packed binary file
//-----------------------------------------
// Here, we do not need to worry about the amount of memory because 
// we are processing a line of the file at a time. This is not the case when 
// creating a binary packed Mt because we have to read in blocks before we can 
// transpose. 
Rprintf( " Converting ascii file to packed binary file.\n " );
CreatePackedBinary(fname, fnamebin, dims, AA, AB, BB);

}




// [[Rcpp::export]]
Eigen::VectorXi  extract_geno_rcpp(CharacterVector f_name_bin, double  max_memory_in_Gbytes, 
                                    int selected_locus, std::vector<long> dims)
{
  std::string
     fnamebin = Rcpp::as<std::string>(f_name_bin);

//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double
  memory_needed_in_Gb =  (dims[0] *  dims[1] *   bits_in_int)/( (double) 8000000000) ;


Eigen::VectorXi
   column_of_genos(dims[0]);


if(max_memory_in_Gbytes > memory_needed_in_Gb ){
   // reading entire data file into memory
   Eigen::MatrixXi genoMat =  ReadBlock(fnamebin,  0, dims[1], dims[0]);
   column_of_genos = genoMat.col(selected_locus);



}  else {
    long num_rows_in_block = (max_memory_in_Gbytes  * (double) 8000000000 )/(bits_in_int * dims[1]);

         int num_blocks = dims[0]/num_rows_in_block;
          if (dims[0] % num_rows_in_block)
                 num_blocks++;


          for(int i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1;
              Rcpp::Rcout << num_rows_in_block1 << " num rows in block 1 " << std::endl;

              Eigen::MatrixXi    genoMat_block1 ( ReadBlock(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;

              int blockj=0;
              for(int j=start_row1; j< start_row1+num_rows_in_block1 ; j++){
                  column_of_genos(j) = genoMat_block1.col(selected_locus)(blockj);
                  blockj++;
              } // end for j
              
          } // end for  i


} // end if max_memory

return(column_of_genos);

}










// [[Rcpp::export]]
Eigen::MatrixXi  calculateMMt_rcpp(CharacterVector f_name_bin, 
                                   double  max_memory_in_Gbytes, int num_cores,
                                   Rcpp::NumericVector  selected_loci , std::vector<long> dims)
{
// set multiple cores
Eigen::initParallel();
omp_set_num_threads(num_cores);
Eigen::setNbThreads(num_cores);



std::string 
   line; 


ofstream
   fileOUT;

int 
   genoval,
   rown,
   coln;

std::string 
     fnamebin = Rcpp::as<std::string>(f_name_bin);


MatrixXi 
    MMt(MatrixXi(dims[0], dims[0]).setZero()),
    genoMat;





//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double 
  memory_needed_in_Gb =  (dims[0] *  dims[1] *   bits_in_int)/( (double) 8000000000) ;





//-------------------------
// Perform MMt calculation
//-------------------------
Rprintf( " Performing MMt calculation ... \n" );
if(max_memory_in_Gbytes > memory_needed_in_Gb ){
   // reading entire data file into memory
   genoMat =  ReadBlock(fnamebin,  0, dims[1], dims[0]);

   if(!R_IsNA(selected_loci(0))){
     // setting columns to 0
     for(int ii=0; ii < selected_loci.size() ; ii++)
       genoMat.col(selected_loci(ii)).setZero();
   }
   MMt = genoMat * genoMat.transpose(); 

} else {
    // based on user defined memory. Doing MMt via blockwise multiplication
    long num_rows_in_block = (max_memory_in_Gbytes  * (double) 8000000000 )/(bits_in_int * dims[1]);
//    if(num_rows_in_block > dims[0]){
//         genoMat =  ReadBlock(fnamebin,  0, dims[1], dims[0]);

//        // Efficient calculation of MMt
//       if(!R_IsNA(selected_loci(0) )){
//       // setting columns to 0
//       for(int ii=0; ii < selected_loci.size() ; ii++)
//          genoMat.col(selected_loci(ii)).setZero();
//       }
//       MMt = genoMat * genoMat.transpose(); 



//    } else {
           // blockwise multiplication

          // find out the number of blocks needed
          int num_blocks = dims[0]/num_rows_in_block;
          if (dims[0] % num_rows_in_block)
                 num_blocks++;


          for(int i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1;
              Rcpp::Rcout << num_rows_in_block1 << " num rows in block 1 " << std::endl;

              Eigen::MatrixXi    genoMat_block1 ( ReadBlock(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;
              Eigen::MatrixXi    MMtsub(MatrixXi(num_rows_in_block1, num_rows_in_block1).setZero());

             if(!R_IsNA(selected_loci(0) )){
             // setting columns to 0
             for(int ii=0; ii < selected_loci.size() ; ii++)
                genoMat_block1.col(selected_loci(ii)).setZero();
             }
              Rcpp::Rcout << "  Block 1  "  << std::endl;
              Rcpp::Rcout << genoMat_block1.rows() << std::endl;
              Rcpp::Rcout << genoMat_block1.cols() << std::endl;

              MMtsub = genoMat_block1 * genoMat_block1.transpose(); 

              //          i            j            num rows               num   cols
              MMt.block(start_row1, start_row1, num_rows_in_block1, num_rows_in_block1) = MMtsub;


              for(int j=i+1;j<num_blocks; j++){
                   long start_row2 = j * num_rows_in_block;
                   long num_rows_in_block2 = num_rows_in_block;
                   if ((start_row2 + num_rows_in_block2) > dims[0])
                          num_rows_in_block2 = dims[0] - start_row2;
                    Eigen::MatrixXi    genoMat_block2 ( ReadBlock(fnamebin,  start_row2, dims[1], num_rows_in_block2)) ;

                   Eigen::MatrixXi    MMtsub(MatrixXi(num_rows_in_block1, num_rows_in_block2).setZero());

                  if(!R_IsNA(selected_loci(0) )){
                   // setting columns to 0
                   for(int jj=0; jj < selected_loci.size() ; jj++)
                      genoMat_block2.col(selected_loci(jj)).setZero();
                   }
              Rcpp::Rcout << " Block 2 " << std::endl;
              Rcpp::Rcout << genoMat_block2.rows() << std::endl;
              Rcpp::Rcout << genoMat_block2.cols() << std::endl;
                   MMtsub = genoMat_block1 * genoMat_block2.transpose(); 
                   //          i,        j,     num rows,              num cols
                   MMt.block(start_row1, start_row2, num_rows_in_block1, num_rows_in_block2) = MMtsub;
                   // and its symmetric block
                   MMt.block(start_row2, start_row1, num_rows_in_block2, num_rows_in_block1) = MMtsub.transpose();


            }  // end for int j





          } // end for int


 // }  // end inner if else

}  // end outer if else


// Now working out square root of MMt via SVD
// Rprintf( " Performing SVD...\n " );
//  Eigen::MatrixXf Ohm (dims[0], dims[0]);
//  Eigen::MatrixXf L (dims[0], dims[0]);

// #pragma omp parallel 
// { 
 // Eigen::JacobiSVD<Eigen::MatrixXf,Eigen::NoQRPreconditioner> svd(MMt.cast<float>(), Eigen::ComputeThinU | Eigen::ComputeThinV);
// Ohm  = svd.singularValues().asDiagonal();
// Rcpp::Rcout << svd.singularValues() << std::endl;
// }




// Rcpp::Rcout << " Writing subset of results " << std::endl;
// for(int i=0; i < 2; i++)
// {
//  for(int j=0; j < dims[0] ; j++){
//   Rcpp::Rcout <<  MMt(i,j) << " "  ;
//  }
// Rcpp::Rcout << std::endl;
//}

  return MMt;

//  return Rcpp::List::create(Rcpp::Named("nummrks")=dims[1],
//                           Rcpp::Named("MMt") = MMt);



}






