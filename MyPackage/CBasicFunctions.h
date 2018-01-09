
/*
 *   CBasicFunctions.h
 *
 *      Author: Josh Bowden
 *      CSIRO
 *      class CBasicFunctions:  A class that manipulates input strings
 *
 */
 
#ifndef CBasicFunctions_
#define CBasicFunctions_


#include <string>


class CBasicFunctions {
public:
    // constructors  
    CBasicFunctions( ) ;
    // destructor 
    
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

#endif
