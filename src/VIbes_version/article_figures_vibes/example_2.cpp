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
    // Generate reference
    Interval domain(0, 8);
    double timestep = 0.01;
    IntervalVector x0({{0., 0.},
                       {0., 0.}});
    Function f("x", "y", "(1;sin(x))");
    CtcLohner ctc_lohner(f);
    TubeVector a(domain, timestep,f.image_dim());
    a.set(x0,0.);
    ctc_lohner.contract(a);
    cout << "Reference generated " << a << endl;

    double epsilon = timestep;

    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1+t;x2+cos(x1)-cos(x1+t) )");
    SepFwdBwd fullSep(phi,X0);
    IntervalVector proj(1);
    proj[0] = Interval(-8,0);
    SepProj sepProj(fullSep,proj,epsilon);

    beginDrawing();

    VIBesFigMap fig_map("Example 2 continuous");
    fig_map.set_properties(50,50,800,368);
    fig_map.axis_limits(x);

    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon);
    auto stop = chrono::steady_clock::now();

    cout << "elapsed time test-case 2 continuous: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig_map.add_tube(&a,"reference",0,1);
    fig_map.set_tube_color(&a,"k[k]");
    fig_map.show();
    fig_map.axis_limits(x);
    endDrawing();

    return;
}


void example_2_discrete()
{
    Interval domain(0, 8);
    double timestep = 0.01;
    double epsilon = timestep;
    Function f("x", "y", "(1;sin(x))");
    CtcLohner ctc_lohner(f);

    // Trajectory of the corners of the initial condition (for drawing purposes)
    TubeVector a1(domain,timestep,f.image_dim());
    TubeVector a2(a1);
    TubeVector a3(a1);
    TubeVector a4(a1);
    a1.set(IntervalVector({{0,0},{0,0}}),0.);
    a2.set(IntervalVector({{1,1},{0,0}}),0.);
    a3.set(IntervalVector({{1,1},{1,1}}),0.);
    a4.set(IntervalVector({{0,0},{1,1}}),0.);
    ctc_lohner.contract(a1);
    ctc_lohner.contract(a2);
    ctc_lohner.contract(a3);
    ctc_lohner.contract(a4);


    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1+t;x2+cos(x1)-cos(x1+t) )");

    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    vector<Interval*> projections{
            new Interval(0),
            new Interval(-2),
            new Interval(-4),
            new Interval(-6),
            new Interval(-8)};

    vector<Sep*> seps;
    for (size_t i=0;i<projections.size();i++)
    {
        SepProj *sepProj = new SepProj(*fullSep,*(projections[i]),epsilon);
        seps.push_back(sepProj);

    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);


    beginDrawing();

    VIBesFigMap fig_map("Example 2 continuous");
    fig_map.set_properties(50,50,800,368);
    fig_map.axis_limits(x);
    auto start = chrono::steady_clock::now();
    sivia(x,usep,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time test-case 2 discrete: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig_map.add_tube(&a1,"reference",0,1);
    fig_map.set_tube_color(&a1,"k[k]");
    fig_map.add_tube(&a2,"a2",0,1);
    fig_map.set_tube_color(&a2,"grey[grey]");
    fig_map.add_tube(&a3,"a3",0,1);
    fig_map.set_tube_color(&a3,"grey[grey]");
    fig_map.add_tube(&a4,"a4",0,1);
    fig_map.set_tube_color(&a4,"grey[grey]");
    fig_map.show();
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

    example_2_continuous();
    example_2_discrete();

    return(0);

}
