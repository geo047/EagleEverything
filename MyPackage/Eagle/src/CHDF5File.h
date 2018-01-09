/**
 *
 * Author: Josh Bowden
 * CSIRO
 * Filename: CHDFFile.h
 * Description: HDF5 file manipulation class definition. The specified type <T> will be the output data type of matrix values
 *  - string, table and atomic datasets are also allowed,
 *  - dataset attributes can be specified as a string or of type <T>.
 *
 *
 * Example usage:
 *
  CHDF5File<float>* tempHDF5File = (CHDF5File<float>*) new CHDF5File<float>(filenameIn.Filename , H5F_ACC_RDWR ) ;  // H5F_ACC_RDONLY
  hid_t datatype = H5T_NATIVE_FLOAT ;  //  H5T_NATIVE_DOUBLE ;

  tempHDF5File->CreateStringListDataset(filenameIn.Dataset + "/log_dataPCA", vector<strings> vs ) ;
  tempHDF5File->CreateDataset(filenameIn.Dataset + "/scoresPCA/x",ReturnMem(),numRows,numRows*numCols, datatype) ;
  tempHDF5File->CreateDataset(filenameIn.Dataset + "/scoresPCA/y",ReturnMem(),numRows,numRows*numCols, datatype) ;
  tempHDF5File->ExtendExistingDataset(filenameIn.Dataset + "/orig_and_final_VariancePCA/y",ReturnMem(),numRows,numRows*numCols, datatype) ;

  delete tempHDF5File ;
 *
 *
 *
 *
 *
 */
#ifndef CHDF5File_
#define CHDF5File_

#include <string>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <stack>
#include <exception>
#include <cstdlib>

using namespace std;

#define DIRECTORYSTRINGLENGTH 1024

#ifdef __unix
    //   The Windows version uses _findfirst, _findnext and _findclose methods
    //    for accessing filesystem directory entries, whereas the UNIX version
    //    uses opendir, readdir, closedir and stat, S_ISDIR and S_ISREG methods
    extern "C"
    {
    #include <sys/types.h> 
    #include <sys/stat.h>
    #include <dirent.h>
    #include <unistd.h>
    #include <errno.h>
    #include <time.h>
    #include <fcntl.h>
   // #include <stdlib.h>  // for system() call in FileMove()
    }
    #define SYSTEMDIRDELIM "/"
    #define SYSTEMDIRDELIMCHAR '/'

    enum t_IsDIRORFILE { isDIR = 16, isFILE = 32 } ;

#elif defined (_WIN32) || (_WIN64)
    #include <Windows.h>
    extern "C"
    {
    #include <time.h>
    #include <io.h>  // used for:  _findfirst, _findnext, _findclose
    }

    // enum t_IsDIRORFILE { isDIR = 16, isFILE = 32 } ;
    // _A_ARCH = 32 ; _A_NORMAL = 0
    enum t_IsDIRORFILE { isDIR = _A_SUBDIR, isFILE = (_A_NORMAL | _A_ARCH) } ;

    #define SYSTEMDIRDELIM "\\"
    #define SYSTEMDIRDELIMCHAR '\\'
#endif


using namespace std;
extern "C"
{
   #include "hdf5.h"
}

// Remember to define _HDF5USEDLL_


typedef enum ROWSORCOLS {rows, cols} ROWSORCOLS ;



struct FilenameAndDataset
{
    string Filename ;       // the filename
    string Dataset  ;       // the data set name - in a HDF file only
    string SampleRange ;    // this is the requested sample range
    size_t FileSampleRange ;// this is tested for on a file by file basis (HDF files only at present)
    size_t RangeStart ;
    size_t RangeEnd   ;
} ;
  



//#include "CBasicFunctions.h"

#if !defined ALIGNMENT_START_BYTES_
#define ALIGNMENT_START_BYTES_ (4096)
//#define ALIGNMENT_ROWS_BYTES_ (64)
#endif


template <class T>
class CHDF5File
{
public: 

    CHDF5File<T> ( string fileName, unsigned int access ) ;  // constructor - opens a hdf5 file ; access = H5F_ACC_RDONLY or H5F_ACC_RDWR
    ~CHDF5File( void ) ;                                     // destructor - closes the hdf5 file handle (hdf_file)
    void FreeDatasetData( void ) ;
    
    // Path creation
    bool    CreateFullHDFGroupPath(  string datasetPathAndNameIn, bool IsPathOnlyIn ) ;  // if IsPathOnlyIn == true then only a path is given (no dataset name)
    
    // Dataset (non-compound table type) functions 
    size_t ReturnDatasetDimensionSize( int dimensionIN ) ;
    hid_t ReturnTemplateVarType( void ) ;  // returns the HDF hid_t type of the instantiated class type <T>
    
    // This function creates a dataset of type <T> in the HDF file, and should transform input data in dataBuf from hid_tDataTypeIn on the fly to type T
    bool CreateDataset( string datasetNameIn, void * dataBuf, size_t numRows, size_t  numValsIn, hid_t hid_tDataTypeIn  ) ; // N.B. dataTypeIn should probably be equal to hid_t_MemAttr. datasetNameIn is full path and name of dataset to create; numRows = number of rows or if = 0 then 1 dimensional array; dataBuf is a pointer to the memory containing the data to place in the HDF5 file, numValsIn is the number of elements in the buffer; Could also use sizeof(dataBuf) / H5TGet_size( hid_tDataTypeIn ) to get the number of elements
    bool CreateDatasetNoCompression( string datasetNameIn, void * dataBuf, size_t numRows, size_t  numValsIn, hid_t hid_tDataTypeIn, hid_t hid_tDataWantedOnDisk  ) ;
    size_t SetDatasetDimensionSizesOnly( string datasetNameIn) ; // returns dataset size in bytes and sets arrays iDatasetRank and pDatasetDims[]
    hsize_t GetDatasetDimension(int dimension_wanted) ; 
    hsize_t GetDatasetDimension(ROWSORCOLS dimension_wanted) ;  // can now specify rows | cols
  
    // Data will be of same type as the data present on disk (e.g. is stored as ints then return data will be ints)
    // returns the pointer to the dataset in pDatasetData and returns the size of the data in bytes
    size_t GetFullDataset( string datasetNameIn ) ;
    // GetFullDatasetReturnInPtr() returns the pointer to the dataset in memStreamIn pointer and returns the size of the data in bytes. Disk datatype is converted to hid_t_DataTypeWantedIn type data
    size_t GetFullDatasetReturnInPtr( string datasetNameIn, T *& memStreamIn,  hid_t  hid_t_DataTypeWantedIn  ) ;
    size_t GetSectionOfDatasetReturnInPtr( string datasetNameIn,   T *&  memStreamIn, size_t rowStartIn, size_t numRowsIn,  hid_t  hid_t_DataTypeWantedIn ) ;
    size_t GetSectionOfDatasetReturnInPreallocatedPtr( string datasetNameIn,   T *  memStreamIn, size_t rowStartIn, size_t numRowsIn,  hid_t  hid_t_DataTypeWantedIn ) ;
    bool WriteRowsToExistingDataset(  string datasetNameIn, void * dataBuf, size_t numRows, size_t  numValsIn, size_t startRow, hid_t memTypeIn ) ; // will fail to write anything if any of the number of rows to be updated are not present. startRow is zero based
    bool SetMemDatatype( string datasetNameIn ) ; // sets the this->hid_t_MemAttr value to the datatype of the input requested dataset, also sets this->pDatasetDims[] and this->pDatasetMaxDims[] and if dataset is a table/compound datatype it will set this->compound_type_id


    // filled using GetFullDataset() and freed using FreeDatasetData() 
    void *  pDatasetData  ;    // raw dataset data. filled using GetFullDataset() and GetFullDatasetReturnInPtr() and freed using FreeDatasetData() 
    hsize_t * pDatasetDims  ;    // dimension size of array dimensions. filled using GetFullDataset(), SetMemDatatype(), CreateDataset() and ExtendExistingDataset() and  and freed using FreeDatasetData() 
    hsize_t * pDatasetMaxDims ;   // N.B. if a value in this array == H5S_UNLIMITED (-1), the max dimension size is unlimited
    hid_t  hid_t_MemAttr ;      // filled using GetFullDataset() and freed using FreeDatasetData(). This is a HDF5 data type that can be used for storing the data from memory back to a HDF5 file using CreateDataset()
    int   iDatasetRank  ;      // number of array dimensions

    string sHDFfilename ;
    hid_t hdf_file_id ;       // this is the file that was opened by the constructor
   // CBasicFunctions bb ;
    void   MyFreeAligned( void * _ptr) ;
    void * MyMallocAligned(  long long int  size_in, long long int alignment_in) ;
    
    // this function returns 0 when substring is not present
    // and returns the "1 based" position of the string
    // i.e same as delphi pos(substr,string) function
    // not the "0 based" C variety
    size_t  pos( char*   subStr, std::string  STR)
    {
        size_t retVal = STR.find(subStr) ;
        retVal += 1 ; // (as delphi version returns 0 if subStr does not exist and not -1 like C++ .find() function)
        return retVal ; 
    }

    // this function returns 0 when substring is not present
    // and returns the "1 based" position of the string
    // i.e same as delphi pos(substr,string) function
    // not the "0 based" C variety
    size_t  pos( std::string  subStr, std::string  STR) // returns the "1 based" position of the string (also returns 1 if subStr is an empty string)
    {
        size_t retVal = STR.find(subStr) ;
        retVal += 1 ; // (as delphi version returns 0 if subStr does not exist and not -1 like C++ .find() function)
        return retVal ; 
    }

    // finds last 
    size_t  poslast( std::string  subStr, std::string  STR) // returns the "1 based" position of the string (also returns 1 if subStr is an empty string)
    {
        size_t  posLast = STR.find_last_of( subStr ) ;
        posLast += 1 ; // (as delphi version returns 0 if subStr does not exist and not -1 like C++ .find() function)
        return posLast ;
    }

    string copy ( std::string STR, long long start,size_t number)
    {
        std::string resStr = "" ;
        if (number > 0)
        {
            if (start-1 > STR.length())  // TODO: check if this does not ignore last character
                resStr = "" ;
            else
                resStr = STR.substr(start-1,number) ;
        }
        return resStr ; 
    }
    
    int length( std::string   tStr )
    {
        int lengthStr = 0 ;
        lengthStr = tStr.length() ;
        return (lengthStr) ;
    }
} ;

// These functions are required to match any other dynamic memory creation allocations
template <class T>
void CHDF5File<T>::MyFreeAligned( void * _ptr)
{ 
    if (_ptr != NULL)
    {
    // These macros may require <features.h> _POSIX_C_SOURCE is a POSIX spec, not a C++ spec.
    #if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600)
            free(_ptr) ;
    #elif (__cplusplus > 199711L)
        #if defined(__INTEL_COMPILER)
            _mm_free (_ptr);  // beware that gnu and intel compilers use different header for this function
        #else
            free(_ptr) ;  // <cstdlib>
        #endif
    #elif defined(__INTEL_COMPILER)
            _mm_free (_ptr);  // beware that gnu and intel compilers use different header for this function
    #else
         free(_ptr) ;
    #endif
    }
}
// alignment_in is the alignment wanted in bytes
// size_in is the size of the the memory required in bytes.
template <class T>
void * CHDF5File<T>::MyMallocAligned(  long long int  size_in, long long int alignment_in)
{
    void * _resptr = NULL ;
// These macros may require <features.h> _POSIX_C_SOURCE is a POSIX spec, not a C++ spec.
#if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600)
    int ret = posix_memalign((void **) &_resptr, alignment_in, size_in) ;   // <cstdlib>
    if (ret != 0)
        std::cerr << "CHDF5File<T>::MyMallocAligned() failed to allocate " << size_in << " bytes base memory" << std::endl ;
#elif (__cplusplus > 199711L)
    _resptr =  aligned_alloc (alignment_in, size_in) ;  // <cstdlib>
#elif defined(__INTEL_COMPILER)
    _resptr = _mm_malloc( size_in, alignment_in ) ;  // beware that gnu and intel compilers use different header for this function
#else
    _resptr =  malloc(size_in) ;
#endif
    if (_resptr == NULL)
    {
        std::cerr << "CHDF5File<T>::MyMallocAligned Error: failed to allocate " << size_in << " bytes memory" << std::endl ;
        exit(-1) ;
    }
}
 
 

// input:
// fileName     "c:\path\hdffilename.h5"
// attributName "/wood_diff/DiffStepSize = 0.2"
// access = H5F_ACC_RDONLY or H5F_ACC_RDWR
template <class T>
CHDF5File<T>::CHDF5File( string fileName, unsigned int access)
{

//   Exception::dontPrint();
         herr_t ret_value = H5Eset_auto2( H5E_DEFAULT, NULL, NULL );
         if( ret_value < 0 )
              cerr << "EH5Eset_auto2() H5Eset_auto failed"  << endl ;

         if ( H5Fis_hdf5(fileName.c_str()) > 0)
         {
             // open the filename ;
            // if(CHDF5File<T>::ALLOWACCESSHDF == )
             hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
             plist_id = H5P_DEFAULT ;

             hdf_file_id = H5Fopen( fileName.c_str() , access, plist_id );
             H5Pclose(plist_id);
         }
         else
         {
            // otherwise create the file
            hdf_file_id = H5Fcreate( fileName.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT) ;  // Fails if file already exists.
         }

         if (hdf_file_id < 0)
        {
            cerr << "CHDF5File() constructor error: failed to create/open file: " << fileName << endl ;
        }

        
        sHDFfilename = fileName ;
       
        pDatasetData = NULL ;  // dataset data
        iDatasetRank = 0 ;
        pDatasetDims = NULL;
        pDatasetMaxDims = NULL ;
        this->hid_t_MemAttr = 0 ;

}




//  hsize_t * pDatasetDims  ;    // dimension size of array dimensions. filled using GetFullDataset(), SetMemDatatype(), CreateDataset() and ExtendExistingDataset() and  and freed using FreeDatasetData() 
// int   iDatasetRank  ;      // number of array dimensions
template <class T>
hsize_t CHDF5File<T>::GetDatasetDimension(int dimension_wanted) 
{
    if (dimension_wanted < this->iDatasetRank)
    {
      return pDatasetDims[dimension_wanted] ;
    } 
    else
    {
      return 0 ; // indicates the dimension does not exist (size is zero)
    }
  
}
template <class T>
hsize_t CHDF5File<T>::GetDatasetDimension( ROWSORCOLS dimension_wanted) 
{
   // return GetDatasetDimension( dimension_wanted) ;  
   if (this->iDatasetRank > 1 )
    {
     if (dimension_wanted == cols)
        return pDatasetDims[1] ;
     else if (dimension_wanted == rows)
        return pDatasetDims[0] ;
    } 
    else
    {
      return 0 ; // indicates the dimension does not exist (size is zero)
    }
}




// destructor
template <class T>
CHDF5File<T>::~CHDF5File( void )
{
    herr_t err = 0;
    if ( hdf_file_id > 0)
    {
        herr_t err =  H5Fclose( this->hdf_file_id ) ;
        if (err < 0 )
            cerr << "HDF5 file close error1" << endl;
        this->hdf_file_id = 0 ;
    }

    this->FreeDatasetData() ;
}




template <class T>
void CHDF5File<T>::FreeDatasetData( void )
{
    if (this->pDatasetData != NULL)
        delete [] ((T*)this->pDatasetData) ;
    pDatasetData = NULL ;
    if (this->pDatasetDims != NULL)
        delete [] (this->pDatasetDims) ;
    pDatasetDims = NULL ;
    if (this->pDatasetMaxDims != NULL)
        delete [] (this->pDatasetMaxDims) ;
    pDatasetMaxDims = NULL ;

    if (this->hid_t_MemAttr > 0)
        H5Tclose(this->hid_t_MemAttr) ;
    this->hid_t_MemAttr = 0 ;
}



template <class T>
bool CHDF5File<T>::CreateFullHDFGroupPath(  string datasetPathAndNameIn, bool IsPathOnlyIn )
{
    string t_sRemainingDatasetPath = "";
    string t_sCurrentGroupFullPath = "";
    int t_pos = 0 ;
    hid_t t_group ;


    t_sRemainingDatasetPath = datasetPathAndNameIn ;
    if (IsPathOnlyIn == false)  // remove the dataset name
    {
        // removal of the last item (which is the dataset name) but leave the last "/"
        t_pos = t_sRemainingDatasetPath.find_last_of( "/" ) ;
        if ((t_pos != -1) && (t_pos != t_sRemainingDatasetPath.length()-1 ))
        {
            t_sRemainingDatasetPath = t_sRemainingDatasetPath.substr(0,t_pos+1) ;
        }
    }
    else // path only is input so make sure it ends with '/'
    {
        if ( (t_sRemainingDatasetPath.find_last_of( "/" )+1) != t_sRemainingDatasetPath.length() )  // add last '/'
        {
            t_sRemainingDatasetPath = t_sRemainingDatasetPath + "/" ;
        }
    }


    t_pos = pos((string)"/", t_sRemainingDatasetPath) ;
    if (t_pos == 1)// remove "/" if it is in the first position
    {
        t_sRemainingDatasetPath  = copy(t_sRemainingDatasetPath, 2 , length(t_sRemainingDatasetPath) - t_pos) ;
        t_pos = pos((string)"/", t_sRemainingDatasetPath) ;
        t_sCurrentGroupFullPath  = "/" + copy(t_sRemainingDatasetPath,1,t_pos-1) ;
        t_sRemainingDatasetPath  = copy(t_sRemainingDatasetPath,t_pos+1 , length(t_sRemainingDatasetPath) - t_pos) ;

    }
    else
    {
        t_sCurrentGroupFullPath  = "/" + copy(t_sRemainingDatasetPath,1,t_pos) ;
        t_sRemainingDatasetPath  = copy(t_sRemainingDatasetPath,t_pos +  1, length(t_sRemainingDatasetPath) - t_pos) ;
    }


    while ( t_sRemainingDatasetPath.length() > 0 )
    {
        t_group = H5Gopen(this->hdf_file_id, t_sCurrentGroupFullPath.c_str(), H5P_DEFAULT ) ;
        if (t_group < 0)
        {
            t_group = H5Gcreate2(this->hdf_file_id, t_sCurrentGroupFullPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) ;

        }

        H5Gclose(t_group) ;

        t_pos = pos((string)"/", t_sRemainingDatasetPath) ;
        t_sCurrentGroupFullPath  +=  "/" + copy(t_sRemainingDatasetPath,1,t_pos) ;
        t_sRemainingDatasetPath  = copy(t_sRemainingDatasetPath,t_pos +  1, length(t_sRemainingDatasetPath) - t_pos) ;
    }

    // create the last group if it does not exist
    t_group = H5Gopen(this->hdf_file_id, t_sCurrentGroupFullPath.c_str(), H5P_DEFAULT ) ;
    if (t_group < 0)
    {
        t_group = H5Gcreate2(this->hdf_file_id, t_sCurrentGroupFullPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) ;
        if (t_group < 0)
            cerr << "Error: CHDF5File<T>::CreateFullHDFGroupPath(): H5Gcreate2: could not create group " <<  t_sCurrentGroupFullPath << endl;
    }
    H5Gclose(t_group) ;


    return true ;
}



// sets the this->hid_t_MemAttr value to the datatype of the input requested dataset, also sets this->pDatasetDims[]
// and this->pDatasetMaxDims[] and if dataset is a table/compound datatype it will set this->compound_type_id.
template <class T>
bool CHDF5File<T>::SetMemDatatype( string datasetNameIn )

{
        hid_t  t_dataset = -1 ;
        hid_t  t_file_attr  ;
        hid_t  t_dataspace  ;

        t_dataset  =  H5Dopen(this->hdf_file_id, datasetNameIn.c_str(), H5P_DEFAULT) ;

        if ( (t_dataset > 0) )  // dataset exists
        {
            // Discover datatype in the file
            t_file_attr     =  H5Dget_type( t_dataset ) ;

            // Find corresponding memory datatype
            if (this->hid_t_MemAttr != 0 ) H5Tclose( this->hid_t_MemAttr ) ;
            this->hid_t_MemAttr   = 0 ;
            this->hid_t_MemAttr   =  H5Tget_native_type( t_file_attr, H5T_DIR_DEFAULT );

            if (this->IsATableDataset(datasetNameIn))
            {
                if (this->compound_type_id != 0 ) H5Tclose( this->compound_type_id ) ;
                    this->compound_type_id = H5Tcopy(this->hid_t_MemAttr) ;
            }

            delete [] this->pDatasetDims  ;
            delete [] this->pDatasetMaxDims  ;
            t_dataspace    =  H5Dget_space( t_dataset  ) ;
            this->iDatasetRank  =  H5Sget_simple_extent_ndims( t_dataspace ) ;
            this->pDatasetDims  =  new hsize_t[this->iDatasetRank] ;  // pDatasetDims
            this->pDatasetMaxDims =  new hsize_t[this->iDatasetRank] ;
            H5Sget_simple_extent_dims( t_dataspace, this->pDatasetDims , this->pDatasetMaxDims ) ;
            // returns the number of dimensions (i.e. the rank) which we already have


            H5Dclose(t_dataset) ;
            H5Tclose(t_file_attr) ;
            H5Sclose(t_dataspace) ;

            return true ;
        }
        else
            return false ;

}



template <class T>
hid_t CHDF5File<T>::ReturnTemplateVarType( void )
{
    hid_t ret_hid ;

    if (typeid(T) == typeid(float))
        ret_hid = H5T_NATIVE_FLOAT ;
    else if (typeid(T) == typeid(double))
        ret_hid = H5T_NATIVE_DOUBLE ;
    else if (typeid(T) == typeid(long double))
            ret_hid = H5T_NATIVE_LDOUBLE ;
    else if (typeid(T) == typeid(char))
            ret_hid = H5T_NATIVE_CHAR ;
    else if (typeid(T) == typeid(unsigned char))
            ret_hid = H5T_NATIVE_UCHAR ;
    else if (typeid(T) == typeid(short))
            ret_hid = H5T_NATIVE_SHORT ;
    else if (typeid(T) == typeid(unsigned short))
            ret_hid = H5T_NATIVE_USHORT ;
    else if (typeid(T) == typeid(int))
            ret_hid = H5T_NATIVE_INT;
    else if (typeid(T) == typeid(unsigned int))
            ret_hid = H5T_NATIVE_UINT ;
    else if (typeid(T) == typeid(long))
            ret_hid = H5T_NATIVE_LONG ;
    else if (typeid(T) == typeid(unsigned long))
            ret_hid = H5T_NATIVE_ULONG ;
    else if (typeid(T) == typeid(long long))
            ret_hid = H5T_NATIVE_LLONG ;
    else if (typeid(T) == typeid(unsigned long long))
            ret_hid = H5T_NATIVE_ULLONG ;

    return ret_hid ;
}



// Creates a dataset of type <T> from input data of type hid_tDataTypeIn
// if numRows == 0 then create a 1D dataset of length numValsIn
// if numRows > 0 then creates a 2D dataset (an array) with an 'unlimited' number of rows and numcols == numValsIn / numRows
template <class T>
bool CHDF5File<T>::CreateDataset( string datasetNameIn, void * dataBuf, size_t numRows, size_t  numValsIn, hid_t hid_tDataTypeIn )
// datasetNameIn is full path and name of dataset to create;
// dataBuf is a pointer to the memory containing the data to place in the HDF5 file;
//  numRows: if == 0 then 1 dimensional array, otherwise 2D array and with numRows == the number of rows of the dataset
// numValsIn is the total number of elements in the buffer;
// Could also use sizeof(dataBuf) / H5TGet_size( hid_tDataTypeIn ) to get the number of elements
// dataTypeIn should probably be equal to hid_t_MemAttr
{

    herr_t err = 0   ;
    herr_t status   ;
    hid_t  t_dataspace_id ;
// hsize_t * dims   ;
// hsize_t * maxdims  ;
    hsize_t * chunk_size ;
    hid_t t_dataset_id ;
// hid_t t_datatype  ;
    hid_t   t_plist = 0  ;


    if (numRows == 0)
    {
        delete [] this->pDatasetDims;
        delete [] this->pDatasetMaxDims ;

        this->pDatasetDims  = new hsize_t[1] ;
        this->pDatasetMaxDims = new hsize_t[1] ;
        chunk_size    = new hsize_t[1] ;

        this->pDatasetDims[0]  = numValsIn ;  // number of rows
        this->pDatasetMaxDims[0] = H5S_UNLIMITED ; // unlimited size of array

        t_dataspace_id = H5Screate_simple(1, this->pDatasetDims,  this->pDatasetMaxDims ) ;

        t_plist  = H5Pcreate(H5P_DATASET_CREATE);

        // only chunk the dataset if it is greater than a certain amount
    // if ((this->pDatasetDims[0] >= 256))
        {
            if  (this->pDatasetDims[0] >=128 )
                chunk_size[0] = 128 ;
            else
                chunk_size[0] = this->pDatasetDims[0] ;

            herr_t status = H5Pset_chunk (t_plist, 1, chunk_size);
            status = H5Pset_deflate( t_plist, 6);  // aggressive deflating extends file writing time so leave on very low
        }

    }
    else
    {
        delete [] this->pDatasetDims;
        delete [] this->pDatasetMaxDims ;

        this->pDatasetDims  = new hsize_t[2] ;
        this->pDatasetMaxDims = new hsize_t[2] ;
        chunk_size    = new hsize_t[2] ;

        this->pDatasetDims[0]  = numRows ;  // number of rows
        this->pDatasetMaxDims[0] = H5S_UNLIMITED ; // unlimited size of array
        this->pDatasetDims[1]  = numValsIn /numRows ;  // number of columns
        this->pDatasetMaxDims[1] = numValsIn /numRows ; // max size of array

        t_dataspace_id = H5Screate_simple(2, this->pDatasetDims,  this->pDatasetMaxDims ) ;


        t_plist  = H5Pcreate(H5P_DATASET_CREATE);


        //if ( (this->pDatasetDims[0] != 1) && (this->pDatasetDims[1] != 1)  )
        {

            if  (this->pDatasetDims[0] >=128 )
                chunk_size[0] = 128 ;
            else
                chunk_size[0] = this->pDatasetDims[0] ;

            if  (this->pDatasetDims[1] >=128 )
                chunk_size[1] = 128 ;
            else
                chunk_size[1] = this->pDatasetDims[1] ;

            status = H5Pset_chunk (t_plist, 2, chunk_size);
            status = H5Pset_deflate( t_plist, 6);  // aggressive deflating extends file writing time so leave on very low
        }
    }

    // create full group path
    CreateFullHDFGroupPath( datasetNameIn, false ) ;
    H5Tclose(this->hid_t_MemAttr) ;
    this->hid_t_MemAttr = H5Tcopy( hid_tDataTypeIn );  // = H5T_NATIVE_INT etc


    // SZIP stuff
/*
    unsigned szip_options_mask;
    unsigned szip_pixels_per_block;
    szip_options_mask=H5_SZIP_NN_OPTION_MASK;
    szip_pixels_per_block=128;

    if (H5Pset_szip (t_plist, szip_options_mask, szip_pixels_per_block) > 0)
        cout << "szip success" << endl;
    else
        cout << "szip fail" << endl;
*/
    hid_t  hid_t_DS_typeAttr ; // this is the type of the dataset to create,
    hid_t_DS_typeAttr = this->ReturnTemplateVarType() ;  // the dataset will be the type of the T template function

    // t_dataset_id   = H5Dcreate2(this->hdf_file_id, datasetNameIn.c_str(), this->hid_t_MemAttr, t_dataspace_id, H5P_DEFAULT, t_plist, H5P_DEFAULT) ;
     t_dataset_id   = H5Dcreate2(this->hdf_file_id, datasetNameIn.c_str(), hid_t_DS_typeAttr, t_dataspace_id, H5P_DEFAULT, t_plist, H5P_DEFAULT) ;
    if (t_dataset_id < 0)
    {
        if (t_dataspace_id > 0) err = H5Sclose(t_dataspace_id) ;
    // if (t_datatype > 0)  err = H5Tclose(t_datatype) ;
        if (t_plist > 0)   err = H5Pclose(t_plist) ;
        cerr << "Error: CHDF5File<T>::CreateDataset() could not create dataset (" << datasetNameIn << ")" << ". Dataset may already exist"  << endl;
        delete [] chunk_size ;
        return false ;
    }

    // this should convert the input data (in dataBuff) to the dataset type specified in H5Dcreate2()
    err =  H5Dwrite( t_dataset_id, this->hid_t_MemAttr , H5S_ALL, H5S_ALL, H5P_DEFAULT, dataBuf ) ;
    if (err < 0)
    {
        if (t_dataspace_id > 0) err = H5Sclose(t_dataspace_id) ;
        if (t_dataset_id > 0) err = H5Dclose(t_dataset_id) ;
//  if (t_datatype > 0)  err = H5Dclose(t_datatype) ;
        if (t_plist > 0)   err = H5Pclose(t_plist) ;

        cerr << "Error: CHDF5File<T>::CreateDataset() could not write dataset (" << datasetNameIn << ")" << endl;
        delete [] chunk_size ;
        return false ;
    }

    err =  H5Fflush(t_dataset_id, H5F_SCOPE_GLOBAL) ;
    if (err < 0)
    {
        cerr << " CHDF5File<T>::CreateDataset() failed to flush dataset: " << datasetNameIn << endl;
    }

    delete [] chunk_size ;

    if (t_dataspace_id > 0) err += H5Sclose(t_dataspace_id) ;
    if (t_dataset_id > 0)   err += H5Dclose(t_dataset_id) ;
// if (t_datatype > 0)  err += H5Tclose(t_datatype) ;
    if (t_plist > 0)  err += H5Pclose(t_plist) ;


    if (err == 0)
        return true ;
    else
    {
        cerr << " CreateDataset() failed when closing HDF structures" << endl;
        return false ;
    }
}


/*!  datasetNameIn is full path and name of dataset to create;
\param dataBuf is a pointer to the memory containing the data to place in the HDF5 file;  If setas NULL then dataset is created but no data written.
\param  numRows: if == 0 then 1 dimensional array, otherwise 2D array and with numRows == the number of rows of the dataset
\param numValsIn is the total number of elements in the buffer;
          Could also use sizeof(dataBuf) / H5TGet_size( hid_tDataTypeIn ) to get the number of elements
\param dataTypeIn should probably be equal to hid_t_MemAttr
*/
template <class T>
bool CHDF5File<T>::CreateDatasetNoCompression( string datasetNameIn, void * dataBuf, size_t numRows, size_t  numValsIn, hid_t hid_tDataTypeIn, hid_t hid_tDataWantedOnDisk  )
{

    herr_t err = 0   ;
    herr_t status   ;
    hid_t  t_dataspace_id ;
    hid_t t_dataset_id ;
    hid_t   t_plist = 0  ;

    if (numRows == 0)
    {
        delete [] this->pDatasetDims;
        delete [] this->pDatasetMaxDims ;

        this->pDatasetDims  = new hsize_t[1] ;
        this->pDatasetMaxDims = new hsize_t[1] ;

        this->pDatasetDims[0]  = numValsIn ;  // number of rows
        this->pDatasetMaxDims[0] = H5S_UNLIMITED ; // unlimited size of array

        t_dataspace_id = H5Screate_simple(1, this->pDatasetDims,  this->pDatasetMaxDims ) ;

        t_plist  = H5Pcreate(H5P_DATASET_CREATE);
    }
    else
    {
        delete [] this->pDatasetDims;
        delete [] this->pDatasetMaxDims ;
        this->pDatasetDims  = new hsize_t[2] ;
        this->pDatasetMaxDims = new hsize_t[2] ;


        this->pDatasetDims[0]  = numRows ;  // number of rows
        this->pDatasetMaxDims[0] = numRows ; // H5S_UNLIMITED ; // unlimited size of array
        this->pDatasetDims[1]  = numValsIn / numRows ;  // number of columns
        this->pDatasetMaxDims[1] = numValsIn / numRows ; // max number of columns  of array

        t_dataspace_id = H5Screate_simple(2, this->pDatasetDims,  this->pDatasetMaxDims ) ;
        if (t_dataspace_id < 0)
        {
            cerr << "Error: CHDF5File<T>::CreateDataset() H5Screate_simple(); returned error: " << t_dataspace_id<<  " while creating dataset: " << datasetNameIn << ")"   << endl;
            return false ;
        }

        t_plist  = H5Pcreate(H5P_DATASET_CREATE);
        if (t_plist < 0)
        {
            H5Sclose(t_dataspace_id) ;
            cerr << "Error: CHDF5File<T>::CreateDataset()  H5Pcreate(H5P_DATASET_CREATE); returned in error (" << datasetNameIn << ")"   << endl;
            return false ;
        }
    }

    // create full group path
    CreateFullHDFGroupPath( datasetNameIn, false ) ;
    H5Tclose(this->hid_t_MemAttr) ;  // close any previously used hid_t
    this->hid_t_MemAttr = H5Tcopy( hid_tDataTypeIn );  // = H5T_NATIVE_INT etc

    hid_t  hid_t_DS_typeAttr ; // this is the type of the dataset to create,
    if (hid_tDataWantedOnDisk == 0)
        hid_t_DS_typeAttr = this->ReturnTemplateVarType() ;  // the dataset will be the type of the T template function
    else
        hid_t_DS_typeAttr = hid_tDataWantedOnDisk ;


     t_dataset_id   = H5Dcreate2(this->hdf_file_id, datasetNameIn.c_str(), hid_t_DS_typeAttr, t_dataspace_id, H5P_DEFAULT, t_plist, H5P_DEFAULT) ;
    if (t_dataset_id < 0)
    {
        cerr << "Error: CHDF5File<T>::CreateDataset() could not create dataset: " << datasetNameIn  << endl;
        cerr << "Dataset dimension: rows : " << this->pDatasetDims[0] << " columns: "<< this->pDatasetDims[1] << " Datatype size: " << sizeof(T) << " bytes" <<  endl;
        if (t_dataspace_id > 0) err = H5Sclose(t_dataspace_id) ;
        if (t_plist > 0)  err = H5Pclose(t_plist) ;
        return false ;
    }

    if (dataBuf != NULL)  // only write data to data set if buffer input is not NULL
    {
        err =  H5Dwrite( t_dataset_id, this->hid_t_MemAttr , H5S_ALL, H5S_ALL, H5P_DEFAULT, dataBuf ) ;
        if (err < 0)
        {
            if (t_dataspace_id > 0) err = H5Sclose(t_dataspace_id) ;
            if (t_dataset_id > 0) err = H5Dclose(t_dataset_id) ;
            if (t_plist > 0)  err = H5Pclose(t_plist) ;

            cerr << "Error: CHDF5File<T>::CreateDataset() could not write dataset (" << datasetNameIn << ")" << endl;
            return false ;
        }
    }

    err =  H5Fflush(t_dataset_id, H5F_SCOPE_GLOBAL) ;
    if (err < 0)
    {
        cerr << " CHDF5File<T>::CreateDataset() failed to flush dataset: " << datasetNameIn << endl;
    }


    if (t_dataspace_id > 0) err += H5Sclose(t_dataspace_id) ;
    if (t_dataset_id > 0)   err += H5Dclose(t_dataset_id) ;
    if (t_plist > 0)  err += H5Pclose(t_plist) ;


    if (err == 0)
        return true ;
    else
    {
        cerr << " CreateDataset() failed when closing HDF structures" << endl;
        return false ;
    }
}





// This only works for 2D datasets
template <class T>
bool CHDF5File<T>::WriteRowsToExistingDataset(  string datasetNameIn, void * dataBuf, size_t numRows, size_t  numValsIn, size_t startRow, hid_t memTypeIn )
// numSamples == the number of rows of the 2D dataset or if == 0 then 1D
{
    herr_t ret ;
    hid_t t_dataset ;
    hid_t t_file_dataspace ;
    hid_t t_mem_space ;
    int rank, status_n   ;

     hsize_t t_size[3]   ;
     hsize_t t_start[2]  ;
     hsize_t t_count[2]  ;
     hsize_t mem_dims[2] ;  // needed for H5Screate_simple()

    // 1/ open the Dataset
    t_dataset = H5Dopen2( this->hdf_file_id, datasetNameIn.c_str() , H5P_DEFAULT ) ;

/* if ( t_dataset > 0)
    {
        herr_t err =  H5Fflush( t_dataset, H5F_SCOPE_LOCAL ) ;  // H5F_SCOPE_LOCAL H5F_SCOPE_GLOBAL
        if (err < 0 )
            cerr << "HDF5 file H5Fflush error on dataset: " << datasetNameIn <<  endl;
    }*/

    // 2/ this retrieves the file dataspace
    t_file_dataspace = H5Dget_space(t_dataset);
    rank      = H5Sget_simple_extent_ndims (t_file_dataspace);
    if (rank != 2){
        printf("HDF error: WriteRowsToExistingDataset() Only works for 2D datasets\n ") ;
        printf("HDF error: Datset name: %s\n", datasetNameIn.c_str()) ;
        printf("HDF error: File: %s \n", this->sHDFfilename.c_str()  ) ;
        if (t_dataset > 0 )
            H5Dclose( t_dataset )  ;
        if (t_file_dataspace > 0 )
            H5Sclose( t_file_dataspace )  ;
        return false ;
    }

    // 3/ get the current size of the file dataspace (t_size[2])
    status_n  = H5Sget_simple_extent_dims (t_file_dataspace, t_size, NULL);
    if (t_size[0] < (startRow + numRows))  // There are not enough rows from the entry point row to fit the number of rows we want to write
    {
        printf("HDF error: There are not enough rows in the dataset. size: %d, startRow: %d, numRows: %d\n", t_size[0], startRow, numRows ) ;
        printf("HDF error: Datset name: %s\n", datasetNameIn.c_str()) ;
        printf("HDF error: File: %s \n", this->sHDFfilename.c_str()  ) ;
        if (t_dataset > 0 )
            H5Dclose( t_dataset )  ;
        if (t_file_dataspace > 0 )
            H5Sclose( t_file_dataspace )  ;
        return false ;
    }

    // 4/ set up hyperslab
    t_start[0] = startRow  ;
    t_start[1] = 0 ;    // this may have to be +1

    t_count[0]= numRows ;  // number of rows to add
    t_count[1]= numValsIn /  numRows  ; // use this if the data is shorter for some unknown reason // t_size[1] ; // number of columns per row (same as original matrix)

    // 5/ this selects a portion of the file dataspace
    ret = H5Sselect_hyperslab(t_file_dataspace, H5S_SELECT_SET, t_start, NULL, t_count, NULL);

    // 6/ this creates a memory dataspace
    mem_dims[0]  = numRows ;
    mem_dims[1]  = numValsIn /  numRows ;
    t_mem_space  = H5Screate_simple(2, mem_dims, NULL);


    // 7/ the type of the data to write
    hid_t new_memtype = H5Tcopy(memTypeIn) ;

    // 8/ Write selection from the vector buffer to the dataset in the file
    ret += H5Dwrite(t_dataset, new_memtype , t_mem_space, t_file_dataspace, H5P_DEFAULT, dataBuf) ;

    if ( t_dataset > 0)
    {
        herr_t err =  H5Fflush( t_dataset, H5F_SCOPE_GLOBAL ) ;  // H5F_SCOPE_LOCAL H5F_SCOPE_GLOBAL
        if (err < 0 )
            cerr << "HDF5 file H5Fflush error on dataset: " << datasetNameIn <<  endl;
    }

    if (t_dataset > 0 )
        ret +=  H5Dclose( t_dataset )  ;
    if (t_file_dataspace > 0 )
        ret +=  H5Sclose( t_file_dataspace )  ;
    if (t_mem_space > 0 )
        ret +=  H5Sclose( t_mem_space )  ;

    if (new_memtype > 0 )
        ret +=  H5Tclose( new_memtype )  ;


    if (ret == 0)
        return true ;
    else
        return false ;
    return true ;
}




template <class T>
size_t CHDF5File<T>::ReturnDatasetDimensionSize( int dimensionIN )
{
    size_t dimSize = 1 ;
    if (this->iDatasetRank > 0)
    {
        dimSize = pDatasetDims[dimensionIN] ;  // this is the size in terms of number of full rows
    }
    else
        return 0 ;

    return dimSize ;

}



// returns dataset size in bytes and sets arrays iDatasetRank and pDatasetDims[]
template <class T>
size_t CHDF5File<T>::SetDatasetDimensionSizesOnly( string datasetNameIn)
{

        string sDatasetPathAndName ;
        string sDatasetWantedVal   ;

        // if: attributNameAndWantedValueIn = "/wood_diff/DiffStepSize = >= 0.2"
        sDatasetPathAndName      = datasetNameIn ;
        
        hid_t  t_dataset = -1 ;
        hid_t  t_dataspace = -1 ;
        herr_t  ret;
        hid_t  t_file_attr  ;
        H5T_class_t dataset_class ;
        size_t  t_iSizeOfDiskElementBytes  = 0 ;
        size_t  t_iSizeOfMemoryElementInBytes = 0;
        size_t  t_iNumberofElementsInDataset  = 1 ;
        size_t  t_iReturnDataSizeInBytes = 0  ;

        this->FreeDatasetData() ;  // frees pDatasetData, this->pDatasetMaxDims and this->pDatasetDims

        t_dataset    =  H5Dopen(this->hdf_file_id, sDatasetPathAndName.c_str(), H5P_DEFAULT) ;
        t_dataspace  =  H5Dget_space( t_dataset  ) ;

        this->iDatasetRank  =  H5Sget_simple_extent_ndims( t_dataspace ) ;
        this->pDatasetDims  =  new hsize_t[this->iDatasetRank] ;  // pDatasetDims
        this->pDatasetMaxDims =  new hsize_t[this->iDatasetRank] ;
        H5Sget_simple_extent_dims( t_dataspace, this->pDatasetDims , this->pDatasetMaxDims ) ; // returns the numebr of dimensions (i.e. the rank) which we already have


        if ( (t_dataset > 0) && (t_dataspace > 0) )  // dataset exists and is a single value
        {
            for (int t1 = 0 ; t1 < this->iDatasetRank; t1++)
            {
                t_iNumberofElementsInDataset *= this->pDatasetDims[t1] ;  // calculate how many elements in the array
            }


            /* Discover datatype in the file */
            t_file_attr     =  H5Dget_type( t_dataset ) ;
            t_iSizeOfDiskElementBytes =  H5Tget_size( t_file_attr ) ;
        /// data_size_inbytes   *= datatype_size_inbytes ;  // = number of elements * size of element
            /* Find corresponding memory datatype */
            this->hid_t_MemAttr   =  H5Tget_native_type( t_file_attr, H5T_DIR_DEFAULT );
            dataset_class    =  H5Tget_class( this->hid_t_MemAttr ) ;

            /// determine size in memory
            t_iSizeOfMemoryElementInBytes =  H5Tget_size( this->hid_t_MemAttr ) ;
            t_iReturnDataSizeInBytes  =  t_iNumberofElementsInDataset * t_iSizeOfMemoryElementInBytes ;  // = number of elements * size of element

            ret =  H5Tclose(t_file_attr) ;
            //ret =  H5Tclose(mem_attr_t) ;
        }

        if (t_dataset > 0 )  ret =  H5Dclose(t_dataset) ;
        if (t_dataspace > 0 ) ret =  H5Sclose(t_dataspace) ;

//  delete pDims ;
//  delete pMaxDims ;

        return t_iReturnDataSizeInBytes ;

}


// returns dataset size in bytes
// Data will be of same type as the data present on disk (e.g. is stored as ints then return data will be ints)
// If you want to convert data then use: GetFullDatasetReturnInPtr()
template <class T>
size_t CHDF5File<T>::GetFullDataset( string datasetNameIn)
{

        string sDatasetPathAndName ;
        string sDatasetWantedVal   ;

        // if: attributNameAndWantedValueIn = "/wood_diff/DiffStepSize = >= 0.2"
        sDatasetPathAndName      =  datasetNameIn; 

        hid_t  t_dataset = -1 ;
        hid_t  t_dataspace = -1 ;
        herr_t  ret;
        hid_t  t_file_attr  ;
        H5T_class_t dataset_class ;
        size_t  t_iSizeOfDiskElementBytes ;
        size_t  t_iSizeOfMemoryElementInBytes ;
        size_t  t_iNumberofElementsInDataset  ;
        size_t  t_iReturnDataSizeInBytes = 0  ;

        t_iSizeOfDiskElementBytes  = 0 ;
        t_iSizeOfMemoryElementInBytes = 0 ;
        t_iNumberofElementsInDataset = 1 ;

        this->FreeDatasetData() ;  // frees pDatasetData, this->pDatasetMaxDims and this->pDatasetDims

        t_dataset  =  H5Dopen(this->hdf_file_id, sDatasetPathAndName.c_str(), H5P_DEFAULT) ;
        t_dataspace  =  H5Dget_space( t_dataset  ) ;

        this->iDatasetRank  =  H5Sget_simple_extent_ndims( t_dataspace ) ;
        this->pDatasetDims  =  new hsize_t[this->iDatasetRank] ;  // pDatasetDims
        this->pDatasetMaxDims =  new hsize_t[this->iDatasetRank] ;
        H5Sget_simple_extent_dims( t_dataspace, this->pDatasetDims , this->pDatasetMaxDims ) ; // returns the number of dimensions (i.e. the rank) which we already have


        if ( (t_dataset > 0) && (t_dataspace > 0) )  // dataset exists and is a single value
        {
            for (int t1 = 0 ; t1 < this->iDatasetRank; t1++)
            {
                t_iNumberofElementsInDataset *= this->pDatasetDims[t1] ;  // calculate how many elements in the array
            }


            /* Discover datatype in the file */
            t_file_attr     =  H5Dget_type( t_dataset ) ;
            t_iSizeOfDiskElementBytes =  H5Tget_size( t_file_attr ) ;
        /// data_size_inbytes   *= datatype_size_inbytes ;  // = number of elements * size of element
            /* Find corresponding memory datatype */
            this->hid_t_MemAttr   =  H5Tget_native_type( t_file_attr, H5T_DIR_DEFAULT );
        // this->hid_t_MemAttr   =  this->ReturnTemplateVarType() ;  // returns the datatype of the template class == <T> and data will be converted on the fly to from disk to this datatype
            dataset_class    =  H5Tget_class( this->hid_t_MemAttr ) ;

            /// determine size in memory
            t_iSizeOfMemoryElementInBytes =  H5Tget_size( this->hid_t_MemAttr ) ;
            t_iReturnDataSizeInBytes  =  t_iNumberofElementsInDataset * t_iSizeOfMemoryElementInBytes ;  // = number of elements * size of element


            // Save attribute data into a char array for the moment.
            // Cast it to correct data type when we have more information, or keep it as string data if that is the case
            //returnDatasetPointerInOut = (char  *) malloc( data_size_inbytes ) ;
            //this->FreeDatasetData() ; have to do this at start of method as it deletes this->hid_t_MemAttr
            if ( this->pDatasetData == NULL )
            {
                this->pDatasetData = (char  *) new char[ t_iReturnDataSizeInBytes ] ;
                ret = H5Dread (t_dataset ,this->hid_t_MemAttr , t_dataspace, H5S_ALL, H5P_DEFAULT, (void *)this->pDatasetData);
            }
            else
            {
                t_iReturnDataSizeInBytes = 0 ;  // retruns 0 if problem with allocation
            }


            ret =  H5Tclose(t_file_attr) ;
            //ret =  H5Tclose(mem_attr_t) ;

        }

        if (t_dataset > 0 )  ret =  H5Dclose(t_dataset) ;
        if (t_dataspace > 0 ) ret =  H5Sclose(t_dataspace) ;

//  delete pDims ;
//  delete pMaxDims ;

        return t_iReturnDataSizeInBytes ;

}


// User can specify what return type the data should be converted to using hid_t_DataTypeWantedIn
// Or if hid_t_DataTypeWantedIn == 0 then output the same datatype as is present on disk.
// returns dataset size in bytes
template <class T>
size_t CHDF5File<T>::GetFullDatasetReturnInPtr( string datasetNameIn, T *& memStreamIn,  hid_t  hid_t_DataTypeWantedIn )
// if  hid_t_DataTypeWantedIn == 0 then use the systems memory version of the file data type
// else convert the file data type to the desired memory data type on the fly by HDF library
{
        string sDatasetPathAndName ;
        string sDatasetWantedVal   ;

        // if: attributNameAndWantedValueIn = "/wood_diff/DiffStepSize = >= 0.2"
        sDatasetPathAndName      = datasetNameIn ; 

        hid_t  t_dataset = -1 ;
        hid_t  t_dataspace = -1 ;
        herr_t  ret;
        hid_t  t_file_attr  ;
        H5T_class_t dataset_class ;
        size_t  t_iSizeOfDiskElementBytes ;
        size_t  t_iSizeOfMemoryElementInBytes ;
        size_t  t_iNumberofElementsInDataset  ;
        size_t  t_iReturnDataSizeInBytes = 0  ;

        t_iSizeOfDiskElementBytes  = 0 ;
        t_iSizeOfMemoryElementInBytes = 0 ;
        t_iNumberofElementsInDataset = 1 ;

        this->FreeDatasetData() ;  // frees pDatasetData, this->pDatasetMaxDims and this->pDatasetDims

        t_dataset  =  H5Dopen(this->hdf_file_id, sDatasetPathAndName.c_str(), H5P_DEFAULT) ;
        if (t_dataset <= 0) 
        {
            std::cerr << "CHDF5File<T>::GetFullDatasetReturnInPtr() Error: Dataset does not exist: " << sDatasetPathAndName << " in file: " << this->sHDFfilename <<  std::endl ;
            memStreamIn = NULL ;
            return 0 ;
        }
        t_dataspace  =  H5Dget_space( t_dataset  ) ;

        this->iDatasetRank  =  H5Sget_simple_extent_ndims( t_dataspace ) ;
        this->pDatasetDims  =  new hsize_t[this->iDatasetRank] ;  // pDatasetDims
        this->pDatasetMaxDims =  new hsize_t[this->iDatasetRank] ;
        H5Sget_simple_extent_dims( t_dataspace, this->pDatasetDims , this->pDatasetMaxDims ) ; // returns the number of dimensions (i.e. the rank) which we already have


        if ( (t_dataset > 0) && (this->iDatasetRank > 0) )  // dataset exists
        {
            for (int t1 = 0 ; t1 < this->iDatasetRank; t1++)
            {
                t_iNumberofElementsInDataset *= this->pDatasetDims[t1] ;  // calculate how many elements in the array
            }


            /* Discover datatype in the file */
            t_file_attr     =  H5Dget_type( t_dataset ) ;
            t_iSizeOfDiskElementBytes =  H5Tget_size( t_file_attr ) ;
        /// data_size_inbytes   *= datatype_size_inbytes ;  // = number of elements * size of element
            /* Find corresponding memory datatype */
            if (hid_t_DataTypeWantedIn == 0)
                this->hid_t_MemAttr   =  H5Tget_native_type( t_file_attr, H5T_DIR_DEFAULT );
            else
                this->hid_t_MemAttr   =  hid_t_DataTypeWantedIn ;  // convert the file datatype to the desired memory datatype on the fly by HDF library

            dataset_class    =  H5Tget_class( this->hid_t_MemAttr ) ;

            /// determine size in memory
            t_iSizeOfMemoryElementInBytes =  H5Tget_size( this->hid_t_MemAttr ) ;
            t_iReturnDataSizeInBytes  =  t_iNumberofElementsInDataset * t_iSizeOfMemoryElementInBytes ;  // = number of elements * size of element


            // Save attribute data into a char array for the moment.
            // Cast it to correct data type when we have more information, or keep it as string data if that is the case
            //returnDatasetPointerInOut = (char  *) malloc( data_size_inbytes ) ;
            //this->FreeDatasetData() ; have to do this at start of method as it deletes this->hid_t_MemAttr

            // memStreamIn.SetSize(t_iReturnDataSizeInBytes) ; // allocate the memory required in bytes
            memStreamIn = (T*) this->MyMallocAligned(t_iNumberofElementsInDataset  *  sizeof(T), ALIGNMENT_START_BYTES_ ) ;
            // memStreamIn = (T*) _mm_malloc(t_iNumberofElementsInDataset  *  sizeof(T), ALIGNMENT_START_BYTES_ ) ;
            ret = H5Dread (t_dataset ,this->hid_t_MemAttr, t_dataspace, H5S_ALL, H5P_DEFAULT, (void *) memStreamIn );

            ret =  H5Tclose(t_file_attr) ;
            //ret =  H5Tclose(mem_attr_t) ;

        }

        if (t_dataset > 0 )  ret =  H5Dclose(t_dataset) ;
        if (t_dataspace > 0 ) ret =  H5Sclose(t_dataspace) ;

//  delete pDims ;
//  delete pMaxDims ;

        return t_iReturnDataSizeInBytes ;

}


// returns dataset size in bytes
// rowStartIn is the 1 based number of the row
template <class T>
size_t CHDF5File<T>::GetSectionOfDatasetReturnInPtr( string datasetNameIn, T *& memStreamIn, size_t rowStartIn, size_t numRowsIn,   hid_t  hid_t_DataTypeWantedIn )
// if  hid_t_DataTypeWantedIn == 0 then use the systems memory version of the file data type
// else convert the file data type to the desired memory data type on the fly by HDF library
{
        string sDatasetPathAndName ;
        string sDatasetWantedVal   ;

        // if: attributNameAndWantedValueIn = "/wood_diff/DiffStepSize = >= 0.2"
        sDatasetPathAndName      = datasetNameIn ;  //  

        hid_t  t_dataset = -1 ;
        hid_t  t_dataspace = -1 ;
        herr_t  ret;
        hid_t  t_file_attr  ;
//  H5T_class_t dataset_class ;
        size_t  t_iSizeOfDiskElementBytes ;
        size_t  t_iSizeOfMemoryElementInBytes ;
        size_t  t_iNumberofElementsInDataset  ;
        size_t  t_iReturnDataSizeInBytes = 0  ;

        t_iSizeOfDiskElementBytes  = 0 ;
        t_iSizeOfMemoryElementInBytes = 0 ;
        t_iNumberofElementsInDataset = 1 ;

        this->FreeDatasetData() ;  // frees pDatasetData, this->pDatasetMaxDims and this->pDatasetDims

        t_dataset  =  H5Dopen(this->hdf_file_id, sDatasetPathAndName.c_str(), H5P_DEFAULT) ;
        t_dataspace  =  H5Dget_space( t_dataset  ) ;

        this->iDatasetRank  =  H5Sget_simple_extent_ndims( t_dataspace ) ;
        this->pDatasetDims  =  new hsize_t[this->iDatasetRank] ;  // pDatasetDims
        this->pDatasetMaxDims =  new hsize_t[this->iDatasetRank] ;
        H5Sget_simple_extent_dims( t_dataspace, this->pDatasetDims , this->pDatasetMaxDims ) ; // returns the number of dimensions (i.e. the rank) which we already have



        if ( (t_dataset > 0) && (this->iDatasetRank == 2) )  // data set exists and is a standard 2D data set
        {

            /* Discover datatype in the file */
            t_file_attr     =  H5Dget_type( t_dataset ) ;
            t_iSizeOfDiskElementBytes =  H5Tget_size( t_file_attr ) ;

            /* Find corresponding memory datatype */
            if (hid_t_DataTypeWantedIn == 0)
                this->hid_t_MemAttr   =  H5Tget_native_type( t_file_attr, H5T_DIR_DEFAULT );
            else
                this->hid_t_MemAttr   =  hid_t_DataTypeWantedIn ;  // convert the file datatype to the desired memory datatype on the fly by HDF library

    //  dataset_class     =  H5Tget_class( this->hid_t_MemAttr ) ;

            // create the memory data space - only the size of the final data needed
            hid_t memSpace_id = H5Screate(H5S_SIMPLE) ;
            hsize_t currentsize[2] ;
            currentsize[0]  = numRowsIn ;
            currentsize[1]  = this->pDatasetDims[1] ;
            ret   = H5Sset_extent_simple(memSpace_id, 2,currentsize,currentsize ) ;
            for (int t1 = 0 ; t1 < this->iDatasetRank; t1++)
            {
                    t_iNumberofElementsInDataset *= currentsize[t1] ;  // calculate how many elements in the array
            }

            // create the file data space - offsets only the size of the final data needed
            // controls where to start reading data from and how much
            hid_t fileSpace_id = H5Dget_space(t_dataset) ;
            //hsize_t currentsize[2] ;
            hsize_t offset[2] ;
            offset[0] = rowStartIn - 1;
            offset[1] = 0 ;

            H5Sselect_hyperslab(fileSpace_id, H5S_SELECT_SET, offset, NULL, currentsize, NULL);
        // ret   = H5Sset_extent_simple(fileSpace_id, 2,currentsize,currentsize ) ;
        // ret   = H5Soffset_simple(fileSpace_id, offset) ;

            this->pDatasetDims[0] = currentsize[0] ;
            this->pDatasetDims[1] = currentsize[1] ;


            hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
            plist_id = H5P_DEFAULT ;


            /// determine size in memory
            t_iSizeOfMemoryElementInBytes =  H5Tget_size( this->hid_t_MemAttr ) ;
            t_iReturnDataSizeInBytes  =  t_iNumberofElementsInDataset * t_iSizeOfMemoryElementInBytes ;  // = number of elements * size of element

            memStreamIn = (T*) this->MyMallocAligned(t_iNumberofElementsInDataset  *  sizeof(T), ALIGNMENT_START_BYTES_ ) ;
            // memStreamIn = (T*) _mm_malloc(t_iNumberofElementsInDataset  *  sizeof(T), ALIGNMENT_START_BYTES_ )
            // memStreamIn.SetSize(t_iReturnDataSizeInBytes) ; // allocate the memory required in bytes

            ret = H5Dread (t_dataset ,this->hid_t_MemAttr , memSpace_id, fileSpace_id, plist_id, (void *) memStreamIn );

            ret =  H5Tclose(t_file_attr) ;
            if (memSpace_id > 0)   ret =   H5Sclose(memSpace_id) ;
            if (fileSpace_id > 0 )  ret = H5Sclose(fileSpace_id) ;
            if (plist_id > 0 )     H5Pclose(plist_id);
        }


        if (t_dataset > 0 )  ret =  H5Dclose(t_dataset) ;
        if (t_dataspace > 0 ) ret =  H5Sclose(t_dataspace) ;


        return t_iReturnDataSizeInBytes ;
}


// returns dataset size in bytes
// rowStartIn is the 1 based number of the row
// memStreamIn is a preallocated buffer area that will be filled with the data from the HDF file that has been converted to the desired data type (as per hid_t_DataTypeWantedIn)
template <class T>
size_t CHDF5File<T>::GetSectionOfDatasetReturnInPreallocatedPtr( string datasetNameIn, T * memStreamIn, size_t rowStartIn, size_t numRowsIn,   hid_t  hid_t_DataTypeWantedIn )
// if  hid_t_DataTypeWantedIn == 0 then use the systems memory version of the file data type
// else convert the file data type to the desired memory data type on the fly by HDF library
{
        string sDatasetPathAndName ;
        string sDatasetWantedVal   ;

        // if: attributNameAndWantedValueIn = "/wood_diff/DiffStepSize = >= 0.2"
        sDatasetPathAndName      = datasetNameIn ;  //  

        hid_t  t_dataset = -1 ;
        hid_t  t_dataspace = -1 ;
        herr_t  ret;
        hid_t  t_file_attr  ;
//  H5T_class_t dataset_class ;
        size_t  t_iSizeOfDiskElementBytes = 0  ;
        size_t  t_iSizeOfMemoryElementInBytes = 0  ;
        size_t  t_iNumberofElementsInDataset = 1 ;
        size_t  t_iReturnDataSizeInBytes = 0  ;

        this->FreeDatasetData() ;  // frees pDatasetData, this->pDatasetMaxDims and this->pDatasetDims

        t_dataset             =  H5Dopen(this->hdf_file_id, sDatasetPathAndName.c_str(), H5P_DEFAULT) ;
        t_dataspace           =  H5Dget_space( t_dataset  ) ;

        this->iDatasetRank    =  H5Sget_simple_extent_ndims( t_dataspace ) ;
        this->pDatasetDims    =  new hsize_t[this->iDatasetRank] ;  // pDatasetDims
        this->pDatasetMaxDims =  new hsize_t[this->iDatasetRank] ;
        H5Sget_simple_extent_dims( t_dataspace, this->pDatasetDims , this->pDatasetMaxDims ) ; // returns the number of dimensions (i.e. the rank) which we already have

        
        if ( (t_dataset > 0) && (this->iDatasetRank == 2) )  // data set exists and is a standard 2D data set
        {
            /* Discover datatype in the file */
            t_file_attr     =  H5Dget_type( t_dataset ) ;
            t_iSizeOfDiskElementBytes =  H5Tget_size( t_file_attr ) ;

            /* Find corresponding memory datatype */
            if (hid_t_DataTypeWantedIn == 0)
                this->hid_t_MemAttr   =  H5Tget_native_type( t_file_attr, H5T_DIR_DEFAULT );
            else
                this->hid_t_MemAttr   =  hid_t_DataTypeWantedIn ;  // convert the file datatype to the desired memory datatype on the fly by HDF library


            // create the memory data space - only the size of the final data needed
            hid_t memSpace_id = H5Screate(H5S_SIMPLE) ;
            hsize_t currentsize[2] ;
            currentsize[0]  = numRowsIn ;
            currentsize[1]  = this->pDatasetDims[1] ;
            ret   = H5Sset_extent_simple(memSpace_id, 2,currentsize,currentsize ) ;
            for (int t1 = 0 ; t1 < this->iDatasetRank; t1++)
            {
                    t_iNumberofElementsInDataset *= currentsize[t1] ;  // calculate how many elements in the array
            }

            // create the file data space - offsets only the size of the final data needed
            // controls where to start reading data from and how much
            hid_t fileSpace_id = H5Dget_space(t_dataset) ;
            //hsize_t currentsize[2] ;
            hsize_t offset[2] ;
            offset[0] = rowStartIn - 1;
            offset[1] = 0 ;

        //    H5Sselect_hyperslab(fileSpace_id, H5S_SELECT_SET, offset, NULL, currentsize, NULL);
            ret   = H5Sset_extent_simple(fileSpace_id, 2,currentsize,currentsize ) ;
            if (ret < 0) Rcpp::Rcout << "Error: HDF5file::GetSectionOfDatasetReturnInPreallocatedPtr(): H5Sset_extent_simple" << std::endl ;
            ret   = H5Soffset_simple(fileSpace_id, (const hssize_t*) offset) ;
            if (ret < 0) Rcpp::Rcout << "Error: HDF5file::GetSectionOfDatasetReturnInPreallocatedPtr(): H5Sset_extent_simple" << std::endl ;
            
            this->pDatasetDims[0] = currentsize[0] ;
            this->pDatasetDims[1] = currentsize[1] ;


            hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
            plist_id = H5P_DEFAULT ;


            /// determine size in memory
            t_iSizeOfMemoryElementInBytes =  H5Tget_size( this->hid_t_MemAttr ) ;
            t_iReturnDataSizeInBytes  =  t_iNumberofElementsInDataset * t_iSizeOfMemoryElementInBytes ;  // = number of elements * size of element

            // memStreamIn = (T*) this->MyMallocAligned(t_iNumberofElementsInDataset  *  sizeof(T), ALIGNMENT_START_BYTES_ ) ;
            // memStreamIn = (T*) _mm_malloc(t_iNumberofElementsInDataset  *  sizeof(T), ALIGNMENT_START_BYTES_ )
            // memStreamIn.SetSize(t_iReturnDataSizeInBytes) ; // allocate the memory required in bytes

            ret = H5Dread (t_dataset ,this->hid_t_MemAttr , memSpace_id, fileSpace_id, plist_id, (void *) memStreamIn );
            int nrowsp =  5;
            int ncolsp = 12;
            int counter = 0 ;
            Rcpp::Rcout << std::endl ;
            Rcpp::Rcout << " First " << nrowsp << " lines and " << ncolsp << " columns of the HDF data pointer" << std::endl ;

            while(counter < nrowsp)
            {   
               for(int i=0; i < ncolsp ; i++){
                   Rcpp::Rcout << memStreamIn[(counter*this->pDatasetDims[1])+i]+1 << " " ;
                }
                 Rcpp::Rcout << std::endl ;
                 Rcpp::Rcout << std::flush ;
                counter++;
              }  // end  while(getline(fileIN, line ))
          
            ret =  H5Tclose(t_file_attr) ;
            if (memSpace_id > 0)   ret =   H5Sclose(memSpace_id) ;
            if (fileSpace_id > 0 )  ret = H5Sclose(fileSpace_id) ;
            if (plist_id > 0 )     H5Pclose(plist_id);
        }


        if (t_dataset > 0 )  ret =  H5Dclose(t_dataset) ;
        if (t_dataspace > 0 ) ret =  H5Sclose(t_dataspace) ;


        return t_iReturnDataSizeInBytes ;
}

#endif
