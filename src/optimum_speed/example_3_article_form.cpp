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





int main(int argc, char* argv[])
{
    Tube::enable_syntheses();


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
    lie_group_ex3_separator fullSep(&a, &X0);
    IntervalVector proj(1);
    proj[0] = Interval(0, 4);
    SepProj sepProj(fullSep, proj, 0.01);

    IntervalVector frame({{-1.2,1.2},{-1.2,1.2}});
    ipegenerator::Figure fig(frame,150,150);


    auto start = chrono::steady_clock::now();
    sivia_article(x, sepProj, epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;

    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_circle(0,0,1);
    fig.set_opacity(30);
    fig.draw_box(X0, "green","green");
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,0,1,"black","black",ipegenerator::STROKE_AND_FILL);
    fig.set_graduation_parameters(-1.5,0.5,-1.5,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_3.ipe");
    fig.save_pdf("example_3.pdf");
}