#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <cstdio>
#include <ctime>

#include <armadillo>

using namespace std;

struct parameter_types
{
    string arg_file_name;
    string arg_delim = "::";
    string delimiter = "\t";
    double cash_payoff;
    double spot_price;
    string strikes;
    arma::vec K1K2;
    string strike_type;
    string option_type;
    int eta;
    string cash_or_asset_payoff;
    double interest_rate;
    double dividend_rate;
    double volatility;
    double maturity;
    double time_step;
    int numb_tsteps;
    double u;
    double pu;
    double pd;
    double pm;
    double confidence_level;
    string output;
    string image_processor;
    string image_viewer;
    string text_viewer;
    string report;
} ;

struct variable_types
{
    arma::dmat pji;
    arma::dmat sji;
    arma::dmat vji;
    arma::dmat ap;
    arma::dmat ee_pfe;
};

parameter_types read_input(string argument, parameter_types par);

parameter_types setting_parameters(parameter_types par);

variable_types setting_variables( variable_types var, parameter_types par);

variable_types trinomial_tree(variable_types var, parameter_types par);

variable_types option_node_values(variable_types var, parameter_types par);

double BS_binary(int eta, double S, double K, double r, double q, double sigma, double tau, parameter_types par);

double cdf(double x);

double pfe_on_tree(arma::vec vj, arma::vec pj, int m, double c);

variable_types calc_ee_pfe(variable_types var, parameter_types par);

void save_ee_pfe(variable_types var, parameter_types par);

void plot_chart(parameter_types par, int cols, string xlabel, string ylabel, string title, string linesorpoints);

vector<string> split_string(string a_string, string delimiter);

arma::vec split_string_to_doubles(string s, string delimiter);

void show_result(string editor, string file);

string read_file(string file);
