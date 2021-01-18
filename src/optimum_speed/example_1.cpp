//
// Created by julien-damers on 11/12/2019.
//
#include <iostream>
#include <cstdlib>
#include <random>
#include "tools.h"
#include "tubex.h"
#include "tubex-capd.h"

#include "chrono"

using namespace std;
using namespace ibex;
using namespace tubex;
using namespace vibes;
using namespace pyibex;


int main (int argc, char* argv[])
{
    Tube::enable_syntheses();

    IntervalVector X0({{0,1},{2,3}}); // The uncertain initial condition

    IntervalVector x({{0,4},{-0.2,4}}); // The space to explore for the set inversion

    // Generate the separator for the forward reach set

    ibex::Function phi("x1","x2","t","(x1-t;x2*exp(t))");
    SepFwdBwd fullSep(phi,X0);
    IntervalVector proj(1);
    proj[0] = Interval(0,5);
    SepProj sepProj(fullSep,proj,0.01);
    double epsilon = 0.01;

    beginDrawing();
    // Visuals initialization
    VIBesFigMap fig_map("Example 1");
    fig_map.set_properties(50,50,800,800);

    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"g[g]");
    fig_map.axis_limits(x);
    endDrawing();
    return(0);
    

}
