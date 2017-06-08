#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

//BEWARE, this program is not idiot-resistant!

template <typename T>
class gmatrix
{
public:
    gmatrix(std::vector< std::vector< T > > coef);
    gmatrix(const char* filename);
    gmatrix( size_t i );//construct i \times i diagonal gmatrix with 1 on diagonal;

    void exchangeRows( size_t i , size_t j );
    void multiplyRow( size_t row_no , T numberToMultiply );
    void addToOneRowAnotherRowMultipliedByAConstant( size_t toWhichRowAdd , size_t whichRowToAdd , T whatConstant );

    void exchangeColumns( size_t i , size_t j );
    void multiplyColumn( size_t row_no , T numberToMultiply );
    void addToOneColumnAnotherColumnMultipliedByAConstant( size_t toWhichRowAdd , size_t whichRowToAdd , T whatConstant );

    gmatrix transpose();

    size_t compute_rank();
	T Mat( size_t i , size_t j )
	{
		return this->mat.at(i).at(j);
	}
    size_t number_of_zero_rows()
    {
        size_t number_zero_rows = 0;
        for ( size_t i = 0 ; i != this->mat.size() ; ++i )
        {
            bool is_row_zero = true;
            for ( size_t j = 0 ; j != this->mat[i].size() ; ++j )
            {
                if ( this->mat[i][j] != 0 )
                {
                    is_row_zero = false;
                    break;
                }
            }
            if ( is_row_zero )++number_zero_rows;
        }
        return number_zero_rows;
    }

    int number_nonzero_rows()
    {
        return this->mat.size() - this->number_of_zero_rows();
    }

    friend ostream& operator << ( ostream& out , const gmatrix& m )
    {
        for ( size_t i = 0 ; i != m.mat.size() ; ++i )
        {
            for ( size_t j = 0 ; j != m.mat[i].size() ; ++j )
            {
                out << m.mat[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }

    double determinant()
    {
        double det = 1;
        for ( size_t i = 0 ; i != this->mat.size() ; ++i )
        {
            det *= this->mat[i][i];
        }
        return det;
    }

    size_t numberOfColumns(){return mat.size();}
    size_t numberOfRows(){return mat[0].size();}


    gmatrix gaussElimination();

    gmatrix gauss_elimination_over_reals();
    
    void write_to_file( const char* filename );
private:
    std::vector< std::vector< T > > mat;
};

template <typename T>
gmatrix<T>::gmatrix(std::vector< std::vector< T > > coef)
{
    std::vector< std::vector< T > > m( coef.size() );
    for ( size_t i = 0 ; i != coef.size() ; ++i )
    {
        std::vector< T > mi( coef[i].size() );
        for ( size_t j = 0 ; j != coef[i].size() ; ++j )
        {
            mi[j] = coef[i][j];
        }
        m[i] = mi;
    }
    this->mat = m;
}

template <typename T>
gmatrix<T>::gmatrix(const char* filename)
{
    ifstream in;
    in.open(filename);
    //we assume that the first two numbers in the file are the sizes of the gmatrix:
    int sizeX, sizeY;
    in >> sizeX >> sizeY;
    for ( size_t i = 0 ; i != (size_t)sizeX ; ++i )
    {
        std::vector< T > v;
        for ( size_t j = 0 ; j != (size_t)sizeY ; ++j )
        {
            T aa;
            in >> aa;
            v.push_back( aa );
        }
        this->mat.push_back(v);
    }
}

template <typename T>
gmatrix<T>::gmatrix( size_t sizee )
{
    for ( size_t i = 0 ; i != sizee ; ++i )
    {
        std::vector< T > v;
        for ( size_t j = 0 ; j != sizee ; ++j )
        {
            if ( i != j )
            {
                v.push_back(0);
            }
            else
            {
                v.push_back(1);
            }
        }
        this->mat.push_back(v);
    }
}

template <typename T>
void gmatrix<T>::exchangeRows( size_t i , size_t j )
{
    std::vector< T > rowI = this->mat[i];
    this->mat[i] = this->mat[j];
    this->mat[j] = rowI;
}

template <typename T>
void gmatrix<T>::multiplyRow( size_t rowNo , T numberToMultiply )
{
    for ( size_t i = 0 ; i != this->mat[rowNo].size() ; ++i )
    {
        this->mat[rowNo][i] *= numberToMultiply;
    }
}

template <typename T>
void gmatrix<T>::addToOneRowAnotherRowMultipliedByAConstant( size_t toWhichRowAdd , size_t whichRowToAdd , T whatConstant )
{
    for ( size_t i = 0 ; i != this->mat[toWhichRowAdd].size() ; ++i )
    {
        this->mat[toWhichRowAdd][i] += this->mat[whichRowToAdd][i]*whatConstant;
    }
}


template <typename T>
void gmatrix<T>::exchangeColumns( size_t i , size_t j )
{
    //TODO
}

template <typename T>
void gmatrix<T>::multiplyColumn( size_t row_no , T numberToMultiply )
{
    //TODO
}

template <typename T>
void gmatrix<T>::addToOneColumnAnotherColumnMultipliedByAConstant( size_t toWhichRowAdd , size_t whichRowToAdd , T whatConstant )
{
    //TODO
}


template <typename T>
gmatrix<T> gmatrix<T>::transpose()
{
    std::vector< std::vector<T> > newMatrix;
    for ( size_t i = 0 ; i != this->mat.size(); ++i )
    {
        std::vector<T> v;
        for ( size_t j = 0 ; j != this->mat.size() ; ++j )
        {
            v.push_back( this->mat[j][i] );
        }
        newMatrix.push_back(v);
    }
    return gmatrix( newMatrix );
}





//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//***********************************************************************************************************************************//
//And now here comes the Gaussian elimination

//we assume that M is a square gmatrix:
template <typename T>
gmatrix<T> gmatrix<T>::gaussElimination()
{
    bool dbg = false;
    gmatrix changeOfBasis( this->numberOfColumns() );
	unsigned current_rank = 0;
	
    for ( size_t colNo = 0 ; colNo != this->numberOfColumns() ; ++colNo )
    {
        while (true)
        {

            bool are_all_elements_below_diagonal_zero = true;
            for ( size_t i = current_rank ; i != this->numberOfRows() ; ++i )
            {
                if ( this->mat[i][colNo] != 0 ){are_all_elements_below_diagonal_zero = false;}
            }
            if ( are_all_elements_below_diagonal_zero )
            {
                break;
            }


            //find minimal element in the collumn colNo
            size_t minimal = current_rank;
            while ( (minimal != this->numberOfRows()) && (this->mat[minimal][colNo] == 0) )++minimal;
            for ( size_t i = minimal+1 ; i < this->numberOfRows() ; ++i )
            {
                if ( this->mat[i][colNo] == 0 )continue;//we are not interested in minimum being zero.
                if ( abs(this->mat[i][colNo]) < abs(this->mat[minimal][colNo]) )minimal = i;
            }
            
            this->exchangeRows(current_rank,minimal);
            changeOfBasis.exchangeRows(current_rank,minimal);
                         
            for ( size_t rowNo = current_rank+1 ; rowNo != this->numberOfRows() ; ++rowNo )
			{			
				if ( this->mat[rowNo][colNo] == 0 )continue;
				T number = -this->mat[rowNo][colNo]/this->mat[current_rank][colNo];   				
				this->addToOneRowAnotherRowMultipliedByAConstant(rowNo,current_rank,number);
				changeOfBasis.addToOneRowAnotherRowMultipliedByAConstant(rowNo,current_rank,number);
			}       
			current_rank++;   
        }
    }      
    return changeOfBasis;
}

//gmatrix[row][column]
template <typename T>
gmatrix<T> gmatrix<T>::gauss_elimination_over_reals()
{    
	bool dbg = false;
    gmatrix changeOfBasis( this->numberOfColumns() );
	unsigned current_rank = 0;
    for ( size_t colNo = 0 ; colNo != this->numberOfColumns() ; ++colNo )
    {
		if ( dbg )
		{
			std::cerr << "colNo : " << colNo << ", current_rank : " << current_rank << endl;
			cerr << "this->mat[current_rank][colNo] : " << this->mat[current_rank][colNo] << endl;
			getchar();
		}
        //if there is a zero element the position :colNo,current_rank?
        if ( this->mat[current_rank][colNo] == 0 )
        {
            //check if there is something nonzero down from here:
            for ( size_t rowNo = current_rank+1 ; rowNo != this->numberOfRows() ; ++rowNo )
            {
                if ( this->mat[rowNo][colNo] != 0 )
                {
                    this->exchangeRows( current_rank , rowNo );
                    changeOfBasis.exchangeRows( current_rank , rowNo );                
                    break;
                }
            }
        }
        if ( this->mat[current_rank][colNo] == 0 )
        {			
			continue; //in this case, this column is zero on the diagonal and below from it, so there is nothing left to be done.  
        }
              
		
		if  (dbg)cerr << "this->mat[current_rank][colNo];    : " << this->mat[current_rank][colNo]    << endl;
		
        //if we are there, that means that there may be some nonzero elements below the diagonal:
        for ( size_t rowNo = current_rank+1 ; rowNo != this->numberOfRows() ; ++rowNo )
        {			
			if ( this->mat[rowNo][colNo] == 0 )continue;
            T number = -this->mat[rowNo][colNo]/this->mat[current_rank][colNo];   
            if ( dbg )
            {
				cerr << "this->mat[rowNo][colNo] : " << this->mat[rowNo][colNo] << endl;
				cerr << "this->mat[current_rank][colNo] : " << this->mat[current_rank][colNo] << endl;
				cerr << "number : " << number << endl;
				getchar();      
			}
            this->addToOneRowAnotherRowMultipliedByAConstant(rowNo,current_rank,number);
            changeOfBasis.addToOneRowAnotherRowMultipliedByAConstant(rowNo,current_rank,number);
        }       
        current_rank++;
    }
    //cerr << "Rank after Gaussian elimination : " << current_rank << endl;
    
    return changeOfBasis;
}

template <typename T>
void gmatrix<T>::write_to_file(const char* filename )
{
    ofstream out;
    out.open( filename );
    out << *this;
    out.close();
}



template <typename T>
size_t gmatrix<T>::compute_rank()
{
    size_t rankk = 0;
    size_t col_no = 0;
    while ( col_no != this->numberOfColumns() )
    {
		while ( ( col_no != this->numberOfColumns() ) && ( this->mat[rankk][col_no] == 0 ) )++col_no;
		if ( col_no != this->numberOfColumns() )++rankk;
	}
	return rankk;
}
