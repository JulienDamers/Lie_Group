//
// Created by julien-damers on 11/12/2019.
//
#include <iostream>
#include <cstdlib>
#include <random>
#include "tools.h"
#include "codac.h"
#include "codac-capd.h"

#include "chrono"

using namespace std;
using namespace ibex;
using namespace codac;
using namespace vibes;
using namespace pyibex;


void example_1_continous()
{
    IntervalVector X0({{0,1},{2,3}}); // The uncertain initial condition

    IntervalVector x({{0,4},{-0.2,4}}); // The space to explore for the set inversion
    double epsilon = 0.1;

    // Generate the separator for the forward reach set

    ibex::Function phi("x1","x2","t","(x1-t;x2*exp(t))");
    SepFwdBwd fullSep(phi,X0);
    IntervalVector proj(1);
    proj[0] = Interval(0,5);
    SepProj sepProj(fullSep,proj,epsilon);

    beginDrawing();
    // Visuals initialization
    VIBesFigMap fig_map("Example 1 continuous");
    fig_map.set_properties(50,50,800,800);

    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig_map.axis_limits(x);
    endDrawing();
    return;
}





void example_1_discrete()
{
    IntervalVector X0({{0,1},{2,3}}); // The uncertain initial condition

    IntervalVector x({{0,4},{-0.2,4}}); // The space to explore for the set inversion

    double epsilon = 0.1;

    // Generate the separator for the forward reach set

    ibex::Function phi("x1","x2","t","(x1-t;x2*exp(t))");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);

    vector<float> projection_times{0.,1.,2.,3.};

    vector<Sep*> seps;
    vector<IntervalVector*> projections;
    for (size_t i=0;i<projection_times.size();i++)
    {
        IntervalVector *proj = new IntervalVector(1); // Defining interval to project
        (*proj)[0] = Interval(projection_times[i],projection_times[i]);
        projections.push_back(proj);
        SepProj *sepProj = new SepProj(*fullSep,*proj,epsilon);
        seps.push_back(sepProj);

    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);


    beginDrawing();
    // Visuals initialization
    VIBesFigMap fig_map("Example 1 discrete");
    fig_map.set_properties(50,50,800,800);

    auto start = chrono::steady_clock::now();
    sivia(x,usep,epsilon); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"g[g]");
    fig_map.axis_limits(x);
    endDrawing();

    int n = seps.size();
    for (int i =0; i<n; i++)
    {
        delete(projections[i]);
        delete(seps[i]);
    }

    return;
}


int main (int argc, char* argv[])
{
    Tube::enable_syntheses();

    example_1_continous();
    example_1_discrete();


    

}
