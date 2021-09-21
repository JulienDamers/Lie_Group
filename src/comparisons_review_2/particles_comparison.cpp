//
// Created by julien-damers on 23/08/2021.
//


#include "codac.h"
#include "codac-rob.h"
#include "codac-capd.h"
#include "tools.h"
#include "chrono"
#include "ctc_cn.h"



using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;

void example_2(bool lohner_done, bool capd_done)
{
    cout << "##########################" << endl;
    cout << "########Example 2#########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();
    codac::CtcEval ctc_eval;

    Interval domain_2_a(0, 35);
    Interval domain_2_x(0, 35); // Define domain of work on which we want to integrate
    double timestep_2 = 0.1;
    IntervalVector x0_2({{0.5, 0.5},
                         {0.5, 0.5}}); // Define initial condition
    Function f_2("x", "y", "(1;sin(x))"); // Evolution function to integrate

    TubeVector x_capd_2 = TubeVector(domain_2_x, timestep_2, 2);
    try
    {
        start = chrono::steady_clock::now();
        x_capd_2 = CAPD_integrateODE(domain_2_x, f_2, x0_2, timestep_2);
        stop = chrono::steady_clock::now();
        cout << "CAPD integration ex 2 processed in : "
             << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
        capd_done = true;
    }
    catch ( exception &e )
    {
        capd_done = false;
        cout << "\n\nException caught!\n" << e.what() << endl;
    }

    VIBesFigMap fig_map("Example 2 discrete");
    fig_map.set_properties(50,50,800,306);
    fig_map.add_tube(&x_capd_2,"reference");
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"g[g]");
    fig_map.axis_limits(x);
    fig_map.show();
    endDrawing();
}