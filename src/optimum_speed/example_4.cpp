//
// Created by julien-damers on 06/01/2020.
//

#include <iostream>
#include <cstdlib>
#include "tools.h"
#include "codac.h"
#include "codac-capd.h"
#include "chrono"



using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;

void example_4_continuous()
{
    Interval domain(0,15);
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01;

    IntervalVector X0(3);
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);

    lie_group_ex4_separator fullSep(&a,&X0);

    IntervalVector x({{-4,4},{-4,4}});
    beginDrawing();
    VIBesFigMap fig("Example 4");
    fig.set_properties(50,50,800,800);


    IntervalVector proj({{-6,6},{0,15}});
    SepProj sepProj(fullSep,proj,0.001);

    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(x0[0].lb(),x0[0].ub(),x0[1].lb(),x0[1].ub(),"g[g]");
    ColorMap myColorMap(InterpolMode::RGB);
    vector<codac::rgb> colors;
    fig.add_tube(&a,"reference",0,1);
    fig.set_tube_color(&a,"k[k]");
    fig.set_tube_max_disp_slices(20000);
    fig.show();
    fig.axis_limits(x);
    endDrawing();
    return;
}


void example_4_discrete()
{
    Interval domain(0,15);
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01;

    IntervalVector X0(3);
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);

    lie_group_ex4_separator* fullSep;
    fullSep = new lie_group_ex4_separator(&a,&X0);

    IntervalVector x({{-4,4},{-4,4}});
    beginDrawing();
    VIBesFigMap fig("Example 4 continuous");
    fig.set_properties(50,50,800,800);

    vector<IntervalVector> projection_times{
        IntervalVector({{-6,6},{1,1}}),
        IntervalVector({{-6,6},{2,2}}),
        IntervalVector({{-6,6},{3,3}}),
        IntervalVector({{-6,6},{4,4}}),
        IntervalVector({{-6,6},{14,14}}),
        IntervalVector({{-6,6},{15,15}})
    };

    vector<Sep*> seps;
    vector<IntervalVector*> projections;
    for (size_t i=0;i<projection_times.size();i++)
    {
        IntervalVector *proj = new IntervalVector(2); // Defining interval to project
        (*proj) = projection_times[i];
        projections.push_back(proj);
        SepProj *sepProj = new SepProj(*fullSep,*proj,epsilon);
        seps.push_back(sepProj);

    }

    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);

    auto start = chrono::steady_clock::now();
    sivia(x,usep,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(x0[0].lb(),x0[0].ub(),x0[1].lb(),x0[1].ub(),"g[g]");
    ColorMap myColorMap(InterpolMode::RGB);
    vector<codac::rgb> colors;
    fig.add_tube(&a,"reference",0,1);
    fig.set_tube_color(&a,"k[k]");
    fig.set_tube_max_disp_slices(20000);
    fig.show();
    fig.axis_limits(x);
    endDrawing();
    return;

    int n = seps.size();
    for (int i =0; i<n; i++)
    {
        delete(projections[i]);
        delete(seps[i]);
    }
}

int main(int argc, char* argv[])
{

    Tube::enable_syntheses();

    example_4_continuous();
    example_4_discrete();



}