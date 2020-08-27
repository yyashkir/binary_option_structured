//	Copyright © 2020 Yashkir Consulting

#include "f_binary.h"

using namespace std;

parameter_types read_input(string argument, parameter_types par)
{
    string line;
    string arg_delimiter = par.arg_delim;

    ifstream infile(argument);

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.cash_payoff = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.spot_price = atof(line.c_str());

    getline(infile,line);
    par.strikes = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    par.option_type = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    par.cash_or_asset_payoff = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.interest_rate = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.dividend_rate = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.volatility = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.maturity = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.time_step = atof(line.c_str());

    getline(infile,line);
    line = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;
    par.confidence_level = atof(line.c_str());

    getline(infile,line);
    par.output = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    par.image_processor = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    par.image_viewer = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    par.text_viewer = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    getline(infile,line);
    par.report = line.substr(line.find(arg_delimiter)+arg_delimiter.length()) ;

    return par;
}

parameter_types setting_parameters(parameter_types par)
{
    double r = par.interest_rate;
    double q = par.dividend_rate;

    if(par.option_type == "call")
        par.eta = 1;
    else
        par.eta = -1;

    par.numb_tsteps = (int)( par.maturity / par.time_step) + 1;
    par.time_step =  par.maturity / (par.numb_tsteps  - 1);

    par.u = exp(par.volatility * sqrt(2 * par.time_step));
    par.pu = pow( ( sqrt(par.u) * exp((r - q) * par.time_step / 2) - 1 )  / ( par.u - 1. ),2);
    par.pd = pow( ( sqrt(par.u) * exp((r - q) * par.time_step / 2) - par.u )  / ( par.u - 1. ),2);
    par.pm = 1. -  par.pu -  par.pd;

    par.K1K2.set_size(2);
    par.K1K2.fill(0);

    auto y = split_string_to_doubles(par.strikes, ";");
    if(y.size() == 1)
    {
        par.K1K2(0) = y(0);
        par.strike_type = "single";
    }
    else
    {
        par.K1K2(0) = y(0);
        par.K1K2(1) = y(1);
        par.strike_type = "double";
    }

    return par;
}

variable_types setting_variables(variable_types var, parameter_types par)
{
    int n = par.numb_tsteps;

    var.pji.set_size(n, 2*n-1);
    var.pji.fill(0);

    var.sji.set_size(n, 2*n-1);
    var.sji.fill(0);

    var.vji.set_size(n, 2*n-1);
    var.vji.fill(0);

    var.ap.set_size(n, 2*n-1);
    var.ap.fill(0);

    var.ee_pfe.set_size(n,4);
    var.ee_pfe.fill(0);

    return var;
}

variable_types trinomial_tree(variable_types var, parameter_types par)
{
    int i, j;

    int n = par.numb_tsteps;
    double S = par.spot_price;
    double u = par.u;
    double pu = par.pu;
    double pd = par.pd;
    double pm = par.pm;

    var.pji(0,0) = 1;
    var.pji(1,1) = pm;
    var.sji(0, 0) = S;
    for (j = 1; j < n; j++)
    {
        for (i = 0; i <= 2*j; i++)
        {
            var.sji(j, i) = S * pow(u, i - j);
            if(i == 0)
                var.pji(j,i) = var.pji(j-1,i) * pd;

            if(i == 1 && j > 1)
                var.pji(j,i) = var.pji(j-1,i) * pd + var.pji(j-1,i-1) * pm;


            if(i >= 2 && i <= 2*j-2 && j >1)
                var.pji(j,i) = var.pji(j-1,i-2) * pu + var.pji(j-1,i-1) * pm  + var.pji(j-1,i) * pd;


            if(i == 2*j-1 && j > 1 )
                var.pji(j,i) = var.pji(j-1,i-2) * pu + var.pji(j-1,i-1) * pm;

            if(i == 2*j )
                var.pji(j,i) = var.pji(j-1,i-2) * pu ;
        }
    }
    return var;
}

double cdf(double x)
{
	// cumulative distribution function
    double  A1 = 0.31938153,
            A2 = -0.356563782,
            A3 = 1.781477937,
            A4 = -1.821255978,
            A5 = 1.330274429,
            RSQRT2PI = 0.39894228040143267793994605993438;
	double  K = 1.0 / (1.0 + 0.2316419 * fabs(x));
	double  cnd = RSQRT2PI * exp(-0.5 * x * x) *
		(K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
	if (x > 0)
		cnd = 1.0 - cnd;
	return cnd;
}

double BS_binary(int eta, double S, double K, double r, double q,
                    double sigma, double tau, parameter_types par)
{
    double d1, d2, v = 0;

    if(tau == 0 )
    {
        if(eta * (S - K) > 0 )
        {
            if(par.cash_or_asset_payoff == "cash")
                return par.cash_payoff;
            if(par.cash_or_asset_payoff == "asset")
                return S;
        }
        else return 0;
    }
    d1 = (log(S/K) + (r - q + sigma * sigma/2.) * tau) / (sigma * sqrt(tau));
    d2 = d1 - sigma * sqrt(tau);
    if(par.cash_or_asset_payoff  == "cash")
        v = par.cash_payoff * exp(-r * tau) * cdf(eta * d2)  ;
    if(par.cash_or_asset_payoff == "asset")
        v = S * exp(-r * tau) * cdf(eta * d1)  ;
    return v;
}

variable_types option_node_values(variable_types var, parameter_types par)
{
    int j, i;
    double tau;
    int n = par.numb_tsteps;
    double dt = par.time_step;
    int eta = par.eta;
    double K1 = par.K1K2(0);
    double K2 = par.K1K2(1);
    double r = par.interest_rate;
    double q = par.dividend_rate;
    double vol = par.volatility;
//cout<<endl<<par.K1K2;//<<" "<<K1;//<<" "<<K2;exit(0);
    for (j = 0; j < n; j++)
    {
        tau = (n - j - 1) * dt;
        for (i = 0; i <= 2*j; i++)
        {
            if( par.strike_type == "single")
                var.vji(j,i) =         BS_binary(eta, var.sji(j,i), K1, r, q, vol, tau, par);
            if(par.strike_type == "double")
                var.vji(j,i) =   eta * BS_binary(eta, var.sji(j,i), K1, r, q, vol, tau, par)
                               - eta * BS_binary(eta, var.sji(j,i), K2, r, q, vol, tau, par);
        }
    }
    return var;
}

double pfe_on_tree(arma::vec vj, arma::vec pj, int m, double c)
{
    int i;
    double pfe = 0;
    arma::uvec ind(m);
    arma::vec p(m);

    //indexes of sorted vj (which is not yet sorted)
    ind = sort_index(vj);

    //sorted vj
    vj = sort(vj);

    //probabilities corresponding to sorted vj:
    for(i = 0; i< m; i++)
        p(i) = pj(ind(i));

    //cumulative probability:
    for(i = 1; i< m; i++)
        p(i) = p(i) + p(i-1);

    //looking for index i corresponding to conf level between p(i) and p(i+1)
    for(i = 0; i< m-1; i++)
        if(c > p(i) && c < p(i+1))
            break;
    // special case of low conf level:
    if(c < p(0))
        i = 0;

    if(i == 0)
        pfe = vj(i) ;
    else
        pfe = vj(i) + (vj(i+1) - vj(i)) * (c - p(i)) / (p(i+1) - p(i));

    return pfe;
}

variable_types calc_ee_pfe(variable_types var, parameter_types par)
{
    int j, i, m;
    double ee;      //expected exposure
    int n = par.numb_tsteps;
    double c = par.confidence_level;

    arma::vec vj;
    arma::vec pj;

    for(j=1; j < n; j++)
    {
        m = 2*j+1;      //number of elements of vji(j,i) for current j
        vj.set_size(m);
        pj.set_size(m);
        ee = 0;
        //node values and probabilities for time point j extracting to vector vj() and pj()
        for(i = 0; i < m; i++)
        {
            vj(i) = var.vji(j,i);
            pj(i) = var.pji(j,i);

            ee = ee + vj(i) * pj(i);        //expected value accumulated
        }

        var.ee_pfe(j,0) = j * par.time_step;
        var.ee_pfe(j,2) = ee;

        var.ee_pfe(j,1) = pfe_on_tree(vj, pj, m, 1 - c);
        var.ee_pfe(j,3) = pfe_on_tree(vj, pj, m, c);
    }
    var.ee_pfe(0,0) = 0;
    var.ee_pfe(0,1) = var.ee_pfe(0,2) = var.ee_pfe(0,3) = var.vji(0,0);

    return var;
}

void plot_chart(parameter_types par, int cols, string xlabel, string ylabel, string title,
                string linesorpoints)
{
    int k;
    string sep = par.delimiter;
	stringstream plot_string;
	string show_string;
	string data_file = par.output;

	string plot;
	string pngfile = data_file + ".png";

    plot_string
        << "set datafile separator '" <<sep<<"';"
        << endl <<"set terminal png;"
        << endl << "set output " << "'" << pngfile << "';"
        << endl << "set grid;"
        << endl << "set key autotitle columnhead;"
        << endl << "set xlabel '" << xlabel << "';"
        << endl << "set ylabel '" << ylabel << "';"
        << endl << "set title '" << title << "';"
        << endl << "plot ";

    for (k = 2; k <= cols; k++)
        plot_string << "'" << data_file << "' using " << 0 << ":" << k << "  with " << linesorpoints << ", ";

    plot = plot_string.str();
    ofstream plot_command_file;
    plot_command_file.open("plotcommand");
    plot_command_file << plot;
    plot_command_file.close();

    string command = par.image_processor + " plotcommand";
    const char *make_graph = command.c_str();



    int ret = system(make_graph);
    if (ret == 0)
    {
        show_string = par.image_viewer + " " + pngfile;
        const char *show_graph = show_string.c_str();
        int retp = system(show_graph);
        if (retp != 0)
            cout << endl << par.image_viewer + " failure" << endl;
    }
    else
    {
        cout << endl << par.image_processor +  " failure\n";
    }

}

void save_ee_pfe(variable_types var, parameter_types par)
{
    int j;
    int n = par.numb_tsteps;
    string sep = par.delimiter;

    ofstream rep(par.report);
    string report_text = "                Command parameters:\n";
    report_text += read_file(par.arg_file_name);
    rep << report_text;

    ofstream myfile;
    myfile.open (par.output);
    myfile << "time"+sep+"pfe-dn"+sep+"ee"+sep+"pfe-up";
    rep <<    "time"+sep+"pfe-dn"+sep+"ee"+sep+"pfe-up";
    for(j=0; j<n; j++)
    {
        myfile <<endl<<var.ee_pfe(j,0)<<sep<< var.ee_pfe(j,1)<<sep<< var.ee_pfe(j,2)<<sep<< var.ee_pfe(j,3);
        rep <<   endl<<var.ee_pfe(j,0)<<sep<< var.ee_pfe(j,1)<<sep<< var.ee_pfe(j,2)<<sep<< var.ee_pfe(j,3);
    }
    myfile.close();

    plot_chart(par, 4, "Time(y)", "Option price", "Binary Option ("+par.strike_type+" strike) pfe envelop ", "linespoints pt 7 lw 1");

    rep.close();
}

string read_file(string file)
{
    string line;
    string full_text = "";
    ifstream rf(file);
    do
    {
        getline(rf, line);
        full_text += line + "\n";
    }while(!rf.eof());

    return full_text;
}

void show_result(string editor, string file)
{
    string read_command = editor + " " + file;
    int retp = system(read_command.c_str());
	if (retp != 0)
        cout << endl << read_command  + " failure" << endl;
}

arma::vec split_string_to_doubles(string s, string delimiter)
{
	int i, n;
	char *end;
	vector<string> a = split_string(s, delimiter);
	n = int(a.size());
	arma::vec doubles(n);
	for (i = 0; i < n; i++)
		doubles(i) = std::strtod(a[i].c_str(), &end);
	return doubles;
}

vector<string> split_string(string a_string, string delimiter)
{
	vector<string> a;
	size_t pos = 0;
	std::string token;
	while ((pos = a_string.find(delimiter)) != std::string::npos) {
		token = a_string.substr(0, pos);
		if (token.size() != 0)
			a.push_back(token);
		a_string.erase(0, a_string.find(delimiter) + delimiter.length());
	}
	if (a_string.size() != 0)
		a.push_back(a_string);
	return a;
}
