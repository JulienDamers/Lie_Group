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


void example_3_continuous_article()
{
    Interval domain(0, 4); // Define full integration time
    double timestep = 0.01; // The timestep used for our integration
    IntervalVector x0(2); // Initial condition for our reference
    x0[0] = Interval(0.5, 0.5);
    x0[1] = Interval(0, 0);
    Function f("x", "y", "(x^3+x*y^2-x+y; y^3+x^2*y-x-y)"); // Evolution function for our reference
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain, f, x0, timestep); // Generating reference trajectory
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01; // precision of the sivia

    IntervalVector X0({{0.4,  0.6},
                       {-0.1, 0.1}}); // Large initial box
    IntervalVector x({{-1.2, 1.2},
                      {-1.2, 1.2}}); // Space to be explored with SIVIA algorithm

    // Create our separator object
    lie_group_ex3_separator fullSep(&a, &X0);
    IntervalVector proj(1);
    proj[0] = Interval(0, 4);
    SepProj sepProj(fullSep, proj, 0.01);



    IntervalVector frame({{-1.2,1.2},{-1.2,1.2}});
    ipegenerator::Figure fig(frame,150,150);

    auto start = chrono::steady_clock::now();
    sivia_article(x, sepProj, epsilon,fig); // Perform the SIVIA algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;

    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_circle(0,0,1);
    fig.set_opacity(30);
    fig.draw_box(X0, "green","green");
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,"a",0,1,"black","black",ipegenerator::STROKE_AND_FILL);
    fig.set_graduation_parameters(-1.5,0.5,-1.5,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_3_continuous.ipe");
    fig.save_pdf("example_3_continuous.pdf");

    return;
}

void example_3_discrete_article()
{
    Interval domain(0, 4);
    double timestep = 0.01;
    IntervalVector x0(2);
    x0[0] = Interval(0.5, 0.5);
    x0[1] = Interval(0, 0);
    // Generate the reference  as we do not have a analytical expression for it
    Function f("x", "y", "(x^3+x*y^2-x+y; y^3+x^2*y-x-y)");
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain, f, x0, timestep);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01; // precision of the sivia

    IntervalVector X0({{0.4,  0.6},
                       {-0.1, 0.1}}); // Starting box constraint
    IntervalVector x({{-1.2, 1.2},
                      {-1.2, 1.2}}); // Space to be explored
    lie_group_ex3_separator* fullSep;
    fullSep = new lie_group_ex3_separator(&a, &X0);
    vector<float> projection_times{0.,1.,2.,3.,4.};

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

    IntervalVector frame({{-1.2,1.2},{-1.2,1.2}});
    ipegenerator::Figure fig(frame,150,150);


    auto start = chrono::steady_clock::now();
    sivia_article(x, usep, epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;

    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_circle(0,0,1);
    fig.set_opacity(30);
    fig.draw_box(X0, "green","green");
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,"a", 0,1,"black","black",ipegenerator::STROKE_AND_FILL);
    fig.set_graduation_parameters(-1.5,0.5,-1.5,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_3_discrete.ipe");
    fig.save_pdf("example_3_discrete.pdf");

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

    example_3_continuous_article();
    example_3_discrete_article();



}