
#include "f_binary.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello structured!" << endl;

    parameter_types par;
    variable_types var;

    par.arg_file_name = argv[1];

    par = read_input(par.arg_file_name, par);

    par = setting_parameters(par);

    var = setting_variables(var, par);

    var = trinomial_tree(var, par);

    var = option_node_values(var, par);

    var = calc_ee_pfe(var, par);

    save_ee_pfe(var, par);

    show_result(par.text_viewer, par.report);

    return 0;
}
