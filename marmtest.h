
#ifndef marmtest_H
#define marmtest_H
#include <fstream>
void marmread(double *marmarray, const int ROWS, const int COLS, fstream & fin);
void marmread(double *marmarray, const int ROWS, const int COLS, fstream & fin)
{

    
   
    //fstream marmdat;
    int i=0;
    //marmdat.open(filename, ios::in);
    
    if (fin.is_open()) { // If file has correctly opened...
        // Output debug message
        cout << "File correctly opened" << endl;
        
        // Dynamically store data into array
        while (fin.good()) { // ... and while there are no errors,
            fin >> marmarray[i]; // fill the row with col elements
            i++;
        }
    }
    else cout << "Unable to open file" << endl;
    fin.close();
    /*for (i=0; i<ROWS; i++){
        for (int j=0; j<COLS; j++){
            cout << marmarray[j+i*COLS] << " " ;
        }
        cout << endl;
    }*/
    
   


    

    
    
    return;
    
    
}

#endif //marmtest_H
