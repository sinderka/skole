/////////////////////////////////////////////
///////////input/output handeling////////////
/////////////////////////////////////////////
/////////////////////////////////////////////

#include <iostream> // cout and cin



////////printing output to consoll ////////


void print(auto str) {
//INPUT: s , a string to print to the consolli

    std::cout << str << "\n";
}

void printArray(double* vec, long n, long m){
//INPUT: vec , an array of length n*m
// 		 n , number of rows to print vec in
//		 m , number of cols to print vec in
//OUTPUT: void
    std::string str;
    for (long ii = 0; ii < n; ii++) {

		for (long jj = 0; jj < m; jj++){

            //std::cout << vec[ii + jj*n] << "\t" ;
            printf("%.6f\t", vec[ii + jj*n] );
            std::cout << str;

		}

		std::cout << "\n";
    }
    std::cout << std::endl;
}


// Takes any string from user
std::string getStringFromUser() {
//INPUT: 
//OUTPUT: a string inputed by the user from consoll

    std::string input;
    std::getline(std::cin, input);

    return input;
}



void errorInInput(std::string type) {
//INPUT: type , a string with a datatype to be output to consoll in an error message
//OUTPUT: void


    std::cin.clear();
    std::cin.ignore(32767, '\n');

    std::cout << "Error in input, " << type << " was expected";
}

// Only real numbers
double getDoubleFromUser() {
//INPUT:
//OUTPUT: input , a double from consoll


    double input;

    std::cin >> input;

    while(true) {
        
        if (std::cin.fail()) {

            errorInInput("double");

        } else {

            std::cin.ignore(32767, '\n');

            return input;
        }
    }
}

// Only long
long getLongFromUser() {
//INPUT:
//OUTPUT: input , a long from consoll

    long input;

    std::cin >> input;

    while(true) {
        
        if (std::cin.fail()) {

            errorInInput("long");

        } else {

            std::cin.ignore(32767, '\n');

            return input;
        }
    }
}

// Only integers
int getIntFromUser() {
//INPUT:
//OUTPUT: input , an integer from consoll

    int input;

    std::cin >> input;

    while(true) {
        
        if (std::cin.fail()) {

            errorInInput("int");

        } else {

            std::cin.ignore(32767, '\n');

            return input;
        }
    }
}
