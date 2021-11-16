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
    Interval domain(0,5); // Integration time of the reference, must be > proj
    double timestep = 0.01; // Time step for the creation of the reference
    IntervalVector x0({{0,0},{1,1}}); // Initial condition for reference
    Function f("x1","x2","(1; -x2)"); // Evolution function to integrate
    CtcLohner ctc_lohner(f);
    TubeVector a(domain,timestep,f.image_dim());
    a.set(x0,0.);
    ctc_lohner.contract(a);
    cout << "Reference generated " << a << endl;

    IntervalVector X0({{0,1},{2,3}}); // The uncertain initial condition

    IntervalVector x({{-0.1,6.5},{-0.2,3.5}}); // The space to explore for the set inversion

    double epsilon = timestep; // define accuracy of paving

    // Generate the separator for the forward reach set
    ibex::Function phi("x1","x2","t","(x1+t;x2/exp(t))"); // define transformation function
    SepFwdBwd fullSep(phi,X0); // Create the separator on Phi with X0 as constraint
    IntervalVector proj(1);
    proj[0] = Interval(-5,0); // Define the interval of time on which we want to integrate
    SepProj sepProj(fullSep,proj,epsilon); // Create the separator object

    beginDrawing();
    // Visuals initialization
    VIBesFigMap fig_map("Example 1 continuous");
    fig_map.set_properties(50,50,800,464);
    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time test-case 1 continuous: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig_map.add_tube(&a,"reference",0,1);
    fig_map.set_tube_color(&a,"k[k]");
    fig_map.show();
    fig_map.axis_limits(x);
    endDrawing();
    return;
}





void example_1_discrete()
{
    // Generate reference
    Interval domain(0,5);
    double timestep = 0.01;
    IntervalVector x0({{0,0},{1,1}});
    Function f("x1","x2","(1; -x2)");
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;


    IntervalVector X0({{0,1},{2,3}});
    IntervalVector x({{-0.1,6.5},{-0.2,3.5}});
    double epsilon = timestep;

    // Generate the separator for the forward reach set
    ibex::Function phi("x1","x2","t","(x1+t;x2/exp(t))");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);

    // Define the different times on which we want to integrate
    vector<Interval*> projections{
            new Interval(0.),
            new Interval(-1.),
            new Interval(-2.),
            new Interval(-3.),
            new Interval(-4.),
            new Interval(-5.)};


    vector<Sep*> seps;
    for (size_t i=0;i<projections.size();i++) // Generate the separator for each individual time
    {
        SepProj *sepProj = new SepProj(*fullSep,*(projections[i]),epsilon);
        seps.push_back(sepProj);
    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep); // Create the union of all separators


    beginDrawing();
    // Visuals initialization
    VIBesFigMap fig_map("Example 1 discrete");
    fig_map.set_properties(50,50,800,464);

    auto start = chrono::steady_clock::now();
    sivia(x,usep,epsilon); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time test-case 1 discrete: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A3[#00FF00A3]");
    fig_map.add_tube(&a,"reference",0,1);
    fig_map.set_tube_color(&a,"k[k]");
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

    example_1_continous();
    example_1_discrete();


    

}
