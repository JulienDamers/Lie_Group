//
// Created by julien-damers on 01/07/2021.
//

#include <iostream>
#include <cstdlib>
#include "tools.h"
#include "codac.h"
#include "codac-capd.h"
#include "chrono"
#include "ipegenerator.h"


using namespace std;
using namespace codac;
using namespace ibex;
using namespace pyibex;

void example_4_continuous_article()
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



    IntervalVector proj({{-6,6},{0,15}});
    SepProj sepProj(fullSep,proj,0.001);

    IntervalVector frame({{-4,4},{-4,4}});
    ipegenerator::Figure fig(frame,150,150);

    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    fig.set_opacity(30);
    fig.draw_box(X0.subvector(0,1),"green","green");
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,0,1,"black","black",ipegenerator::STROKE_AND_FILL);
    fig.set_graduation_parameters(-4,0.5,-4,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_4_continuous.ipe");
    fig.save_pdf("example_4_continuous.pdf");
    return;
}

void example_4_discrete_article()
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


    IntervalVector x({{-4,4},{-4,4}});

    IntervalVector frame({{-4,4},{-4,4}});
    ipegenerator::Figure fig(frame,150,150);

    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    fig.set_opacity(30);
    fig.draw_box(X0.subvector(0,1),"green","green");
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,0,1,"black","black",ipegenerator::STROKE_AND_FILL);
    fig.set_graduation_parameters(-4,0.5,-4,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_4_discrete.ipe");
    fig.save_pdf("example_4_discrete.pdf");

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
    example_4_continuous_article();
    example_4_discrete_article();


}