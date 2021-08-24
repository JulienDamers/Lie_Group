//
// Created by julien-damers on 12/12/2019.
//
#include <iostream>
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


void example_2_continuous()
{
    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd fullSep(phi,X0);
    double epsilon = 0.1;


    beginDrawing();

    VIBesFigMap fig_map("Example 2 continuous");
    fig_map.set_properties(50,50,800,306);



    IntervalVector proj(1);
    proj[0] = Interval(0,8);
    SepProj sepProj(fullSep,proj,epsilon);
    cout << "X: " << x << endl;
    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"g[g]");
    fig_map.axis_limits(x);
    endDrawing();

    return;
}


void example_2_discrete()
{
    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    double epsilon = 0.1;

    vector<float> projection_times{8.};

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

    VIBesFigMap fig_map("Example 2 discrete");
    fig_map.set_properties(50,50,800,306);
    auto start = chrono::steady_clock::now();
    sivia(x,usep,epsilon);
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

    //example_2_continuous();
    example_2_discrete();

    return(0);

}
