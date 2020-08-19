
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
//#include "f_bin.h"
#include <armadillo>

using namespace std;



struct parameter_types
{
    string delimiter = "::";
    double cash_payoff;
    double spot_price;
    string strikes;
    string option_type;
    string cash_or_asset_payoff;
    double interest_rate;
    double dividend_rate;
    double volatility;
    double maturity;
    double time_step;
    double confidence_level;
    string output;
    string image_processor;
    string image_viewer;
    string text_viewer;
} ;

parameter_types read_input(string argument, parameter_types par)
{

    string line;
    ifstream infile(argument);

    string delimiter = par.delimiter;

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.cash_payoff = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.spot_price = atof(line.c_str());

    getline(infile,line);
    par.strikes = line.substr(line.find(delimiter)+delimiter.length()) ;

    getline(infile,line);
    par.option_type = line.substr(line.find(delimiter)+delimiter.length()) ;

    getline(infile,line);
    par.cash_or_asset_payoff = line.substr(line.find(delimiter)+delimiter.length()) ;

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.interest_rate = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.dividend_rate = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.volatility = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.maturity = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.time_step = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(delimiter)+delimiter.length()) ;
    par.confidence_level = atof(line.c_str());

    getline(infile,line);
    par.output = line.substr(line.find(delimiter)+delimiter.length()) ;

    getline(infile,line);
    par.image_processor = line.substr(line.find(delimiter)+delimiter.length()) ;

    getline(infile,line);
    par.image_viewer = line.substr(line.find(delimiter)+delimiter.length()) ;

    getline(infile,line);
    par.text_viewer = line.substr(line.find(delimiter)+delimiter.length()) ;

    return par;
}

int main(int argc, char **argv)
{
    cout << "Hello structured!" << endl;

    parameter_types par;

    par = read_input(string(argv[1]), par);

    cout<<endl<<par.delimiter<<" "<<par.maturity<<" "<<par.time_step<<" "<<par.text_viewer;

    return 0;
}
