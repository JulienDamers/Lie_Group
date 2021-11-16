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
    // Generate reference
    Interval domain(0,15);
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    CtcLohner ctc_lohner(f);
    TubeVector a(domain, timestep, f.image_dim());
    a.set(x0,0.);
    ctc_lohner.contract(a);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01; // precision of SIVIA

    IntervalVector X0(3); // Large initial condition
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);
    IntervalVector x({{-4,4},{-4,4}}); // Space to be explored with SIVIA



    lie_group_ex4_separator fullSep(&a,&X0);
    IntervalVector proj({{-6,6},{-15,0}});
    SepProj sepProj(fullSep,proj,timestep);

    beginDrawing();
    VIBesFigMap fig("Example 4 continuous");
    fig.set_properties(50,50,800,800);
    fig.axis_limits(x);
    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig.add_tube(&a,"reference",0,1);
    fig.set_tube_color(&a,"k[k]");
    fig.show();
    fig.axis_limits(x);
    endDrawing();
    return;
}


void example_4_discrete()
{
    // Generate reference
    Interval domain(0,15);
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    CtcLohner ctc_lohner(f);
    TubeVector a(domain, timestep, f.image_dim());
    a.set(x0,0.);
    ctc_lohner.contract(a);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01;

    IntervalVector X0(3);
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);
    IntervalVector x({{-4,4},{-4,4}});


    lie_group_ex4_separator* fullSep;
    fullSep = new lie_group_ex4_separator(&a,&X0);

    vector<IntervalVector*> projections{
            new IntervalVector({{-6,6},{-1,-1}}),
            new IntervalVector({{-6,6},{-2,-2}}),
            new IntervalVector({{-6,6},{-3,-3}}),
            new IntervalVector({{-6,6},{-4,-4}}),
            new IntervalVector({{-6,6},{-14,-14}}),
            new IntervalVector({{-6,6},{-15,-15}})
    };
    vector<Sep*> seps;
    for (size_t i=0;i<projections.size();i++)
    {
        SepProj *sepProj = new SepProj(*fullSep,*(projections[i]),epsilon);
        seps.push_back(sepProj);
    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);


    beginDrawing();
    VIBesFigMap fig("Example 4 discrete");
    fig.set_properties(50,50,800,800);
    fig.axis_limits(x);

    auto start = chrono::steady_clock::now();
    sivia(x,usep,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig.add_tube(&a,"reference",0,1);
    fig.set_tube_color(&a,"k[k]");
    fig.show();
    fig.axis_limits(x);
    endDrawing();

    int n = seps.size();
    for (int i =0; i<n; i++)
    {
        delete(projections[i]);
        delete(seps[i]);
    }

    return;
}

int main(int argc, char* argv[])
{
    Tube::enable_syntheses();
    example_4_continuous();
    example_4_discrete();

}