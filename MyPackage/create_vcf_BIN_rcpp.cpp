// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include  "createM_BIN_rcpp.h"


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif



 Rcpp::NumericVector dimOfFile(     std::string fname,
                                  Rcpp::Function message)
 {

  // Internal function
  // Purpose: check that vcf has correct first row and 
  //          calculate the dimension of M (not Mt) 
  // Return:  array[2] with col, row 


  std::string token, line;
  Rcpp::NumericVector dim_of_M(2) ;  // row x col of M matrix


  // open vcf marker text  file
  std::ifstream fileIN(fname.c_str());

  if(!fileIN.good()) {
      message("ERROR: Vcf file could not be opened with filename  " , fname , "\n" );
      return 0;
  }

  // check that vcf file is a vcf file
    getline(fileIN, line );  // read first line of file
    std::istringstream streamA(line);
    streamA >> token;  // first field of first line


  if (token.rfind("##fileformat=VC", 0) != 0) {
     message("ERROR: This does not appear to be a vcf file because the first row is not beginning with  ##fileformat=VCF... ");
     return 0;
  }

  // Determine number of individuals and snps 
    long row_count = 0;
    long col_count = 0;

  while(getline(fileIN, line)){
     std::istringstream streamA(line);
     streamA >> token;  // first field 
     if (token.rfind("##", 0) != 0){
        // line is not preamble
        if (token.rfind("#CHROM", 0) != 0){
            // data line
            row_count++;  // number of snps
            if (row_count == 1){
               while(streamA >> token){
                 col_count++;  // number of individuals
               } // end while
            } // if row_count
        } // end if
     } // end if
  }  // end while getline

   col_count = col_count - 9 + 1 ; // removing cols that do not contain genogtype data but adjusting for
                                 // counts begin from 1, not 0. 
   Rcpp::Rcout << " Testing " << std::endl;
   Rcpp::Rcout << "number rows " << row_count << " number of cols " << col_count << std::endl;
   dim_of_M[0] = col_count;  // number of individuals
   dim_of_M[1] = row_count;  // number of columns

  // closing file
  fileIN.close();


   // return dim of M (not Mt)
   return dim_of_M;

  }



 Rcpp::NumericVector   createBINfiles_withinmemory(std::string fname, std::string fnamebinMt,  std::string fnamebinM,
                                                  Rcpp::NumericVector dim_of_M, bool quiet, Rcpp::Function message)
 {
  // Purpose:  create the Mt and M binary files of the recoded vcf file. 
  //           This is done within memory.  
  // Return:   boolean vector of whether the snp was removed due to being monomorphic, of low MAF, or multiallelic. 

  std::string token, line;



 // open output bin file that is to hold  Mt vcf genotype data
 // std::ofstream fileOUT(fnamebinMt.c_str(), std::ios::out | std::ios::binary );
 std::ofstream fileOUT(fnamebinMt.c_str(), std::ios::out  );
  if(!fileOUT.good()) {
      message("ERROR: Temporary binary file could not be opened with filename  " , fname , "\n" );
      return 0;
  }




 Rcpp::Rcout << " File name = " << fnamebinMt.c_str()  << std::endl;
  if (!quiet ){
      message("");
      message(" Reading vcf File  ");
      message("");
  }

/*

    // genotypes will hold recoded genotypes (0,1,2) 
    //Create an array of pointers that points to more arrays
    char** genotypes = new char*[ (int) dim_of_M[1] ];
    for (int i = 0; i < (int) dim_of_M[1] ; ++i) {
        genotypes[i] = new char[ (int) dim_of_M[0] ];
    }
*/


    // genotypes will hold recoded genotypes (0,1,2) 
    //Create an array of pointers that points to more arrays
    char** genotypesTrans = new char*[ (int) dim_of_M[0] ];
    for (int i = 0; i < (int) dim_of_M[0] ; ++i) {
        genotypesTrans[i] = new char[ (int) dim_of_M[1] ];
    }




  // open vcf marker text  file
  std::ifstream fileIN(fname.c_str());
  if(!fileIN.good()) {
      message("ERROR: Vcf file could not be opened with filename  " , fname , "\n" );
      return 0;
  }

  // skip over preamble 
  bool go = true;
  while(go){
     getline(fileIN, line);
     std::istringstream streamA(line);
     streamA >> token;  // first field of first line  
     if (token.rfind("##", 0) != 0) 
         go = false;  // finish passing over preamble
  }



  // check that this vcf file is a genotype file by looking for the FORMAT column name in the 9th column of the header row
  //getline(fileIN, line);
  std::istringstream streamB(line);
  for(int i=0; i<9; i++){
     streamB >> token; 
   }


  if (token.rfind("FORMAT", 0) != 0){
       message("ERROR: This does not appear to be a vcf file containing genotype calls because the 9th header column did not contain ");
       message("       the FORMAT field (which is required for genotype data). ");
       return 0;
   }



// Reading in the genotype calls and writing to Mt.bin line by line 
long ind_count = 0;  // column count
long snp_count = 0;  // row count

 // initializing input line 
std::vector<char> rowinfile( (int) dim_of_M[0] );
Rcpp::NumericVector  shouldsnpbremoved( (int) dim_of_M[1]  ) ;




 int counter_col,  counter_row = 0,
    al=0, 
    ii=0;
 bool geno=false,  noignore;
 char previous;
 while( getline(fileIN, line))
 {
   previous = line[0];  // first character in line
 
   counter_col = 0;  // genotype counter_col
    shouldsnpbremoved[counter_row] = 0;  // set to 0 initially
   for(int i=1; i < line.size(); i++)
   {

     if (previous == '\t' && line[i] == '0') {
        al = 0;
     }  else if (previous == '\t' && line[i] == '1') {
        al = 1;
     }  else if (previous == '\t' && line[i] == '2') {
        shouldsnpbremoved[counter_row] = 1;
    }  else if (previous == '\t' && line[i] == '.') {
        al = -9;
        geno=true;
    }  else if ( al != -9 && !shouldsnpbremoved[counter_row]  && (previous == '/' || previous == '|') && line[i] == '0' ){
      al = al + 0;
        geno=true;
    } else if ( al != -9 && !shouldsnpbremoved[counter_row]  && (previous == '/' || previous == '|') && line[i] == '1' ){
      al = al + 1;
        geno=true;
    } else if  ( al != -9 && !shouldsnpbremoved[counter_row]  && (previous == '/' || previous == '|') && line[i] == '2' ){
        shouldsnpbremoved[counter_row] = 1;
   }  else if (  al != -9 && !shouldsnpbremoved[counter_row]  && (previous == '/' || previous == '|') && line[i] == '.' ){
      al = -9;
        geno=true;
   } 
   if ( shouldsnpbremoved[counter_row] ){
      geno = false;
  }   



  if (geno ){

     // check that counter_row is not greater than 
     if (counter_col ==  dim_of_M[0]  ){
      message("\n");
      message("Error:  a problem has occurred when reading in the data from the vcf file. ");
      message("        The error has occurred for snp number " , counter_row+1, ".");
      message("        This snp may contain more than ", dim_of_M[0] , " columns of data or the snp may contain tabs at the end of the line. \n ");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
     }


   geno = false;
   if (al == 0)
      rowinfile[counter_col] = '0';
   if (al == 1 || al == -9)
      rowinfile[counter_col] = '1';
   if (al == 2)
      rowinfile[counter_col] = '2';


   // genotypes[counter_row][counter_col] = rowinfile[counter_col];
   genotypesTrans[counter_col][counter_row] = rowinfile[counter_col];
   counter_col++; 
  }

  previous = line[i];

  }  // end for



 if (counter_col != dim_of_M[0]  ){
      message("\n");
      message("Error:  The vcf file contains an unequal number of snp (or columns) per individual (or row).  ");
      message("        The error has occurred for snp number " , counter_row+1 , ". It contains " , counter_col , " columns of data.");
      message("        It should contain " , dim_of_M[0] , " columns of data. \n");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
}
 
   if ( !shouldsnpbremoved[counter_row] ){
      // writing vector to binary file
       fileOUT.write( (char *) &rowinfile[0], rowinfile.size() * sizeof(char));
   }
   counter_row++;

} // end while getline

fileOUT.close();
fileIN.close();

// WRite out transpose of Mt (M) to binary file

// open output bin file that is to hold  Mt vcf genotype data
 // std::ofstream fileOUT_M(fnamebinM.c_str(), std::ios::out | std::ios::binary );
 std::ofstream fileOUT_M(fnamebinM.c_str(), std::ios::out  );
 if(!fileOUT_M.good()) {
      message("ERROR: Temporary binary file could not be opened with filename  " , fname , "\n" );
      return 0;
  }
 Rcpp::Rcout << " File name of M file = " << fnamebinM.c_str()  << std::endl;

 Rcpp::Rcout << "about to do transpose .... " << std::endl;

/*

 for(int i=0; i < dim_of_M[0]; i++){
   for(int j=0; j < dim_of_M[1]; j++){
   // fileOUT_M.write( (char *) &genotypesTrans[i][j],  sizeof(char));
   Rcpp::Rcout << "i = " << i << " j = " << j << genotypesTrans[i][j] << std::endl;
   fileOUT_M.write( (char *) &genotypesTrans[i][j],  sizeof(char));
   }  

 }  // end for int

*/


      for(int i=0; i< dim_of_M[0]; i++){
       fileOUT_M.write( (char *) &genotypesTrans[i][0] , dim_of_M[1] * sizeof(char));
      }

fileOUT_M.close();





 //Free each sub-array
  for (int i = 0; i < (int) dim_of_M[0] ; ++i) {
        delete[] genotypesTrans[i];   
    }
 //Free the array of pointers
    delete[] genotypesTrans;

 

  return shouldsnpbremoved;



}




// [[Rcpp::export]]
Rcpp::NumericVector   create_vcf_BIN_rcpp(Rcpp::CharacterVector f_name, Rcpp::CharacterVector f_name_bin_M,
                           Rcpp::CharacterVector f_name_bin_Mt, double  max_memory_in_Gbytes,  bool quiet, Rcpp::Function message) 
{
  // Rcpp function to create binary file from vcf formatted input files



  double
     max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;

  std::string token, line;

  std::string
       fname = Rcpp::as<std::string>(f_name),
       fnamebinM = Rcpp::as<std::string>(f_name_bin_M),
       fnamebinMt = Rcpp::as<std::string>(f_name_bin_Mt);

  Rcpp::NumericVector dim_of_M(2) ;  // row x col of M matrix





  // Check if file exists and find its dimension
  dim_of_M = dimOfFile(fname, message);
  if(dim_of_M.length() == 1){
    // problem with file
    return 0;
  }




  // check if data can be held in memory as char matrix 
  double mem_bytes = (1.2 * dim_of_M[1] * dim_of_M[0] * (sizeof(char) * CHAR_BIT) )/8.0;  

 Rcpp::NumericVector shouldsnpbremoved  = createBINfiles_withinmemory(fname, fnamebinMt,  fnamebinM,  dim_of_M, quiet, message);







return 0;





/*

while(getline(fileIN, line))

    std::istringstream streamC(line);
    for(int i=0; i<9; i++){
       // ignore first 9 columns of data
       streamC >> token;
    }
   ind_count = 0;
   while (streamC >> token){
       std::string delimiter = ":";
       std::string allele  = token.substr(0, token.find(delimiter)); //
       if (allele == "0|0" || allele=="0/0"){
          rowinfile[ind_count] = '0';
       } else if (allele == "0/1" || allele == "1/0" || allele == "0|1" || allele == "1|0"){
          rowinfile[ind_count] = '1';
       } else if (allele == "1/1" || allele == "1|1"){
          rowinfile[ind_count] = '2';
      }  else if (allele == "." || "./0" || "0/." || ".|0" || "0|." || "./1" || "1/." || "1|." || ".|1" ){
          rowinfile[ind_count ] = '1';
      } else {
        // not a SNP marker and shoud be ignored and removed from the map
        shouldsnpbremoved[snp_count] = true;
      } // end if 
     ind_count++;
                      
    }  // end while streamA

 // ----> the correct one if (ind_count != dim_of_M[0]  ){
 if (ind_count != dim_of_M[0]  ){
      message("\n");
      message("Error:  The vcf file contains an unequal number of columns per row.  ");
      message("        The error has occurred for the " , snp_count+1 , " snp which contains " , ind_count , " but ");
      message("        it should contain " , dim_of_M[0] , " columns of data. ");
      message(" ReadVCF has terminated with errors");
      message("\n");
      message("\n");
      return 0 ;
  }

   // writing vector to binary file
   fileOUT.write( (char *) &rowinfile[0], rowinfile.size() * sizeof(char));

  snp_count++;


} // end while getline


return 0;

*/










  return dim_of_M;


}




