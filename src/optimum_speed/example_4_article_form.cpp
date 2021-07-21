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

    Interval domain(0,15);
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.1;

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
    fig.draw_box(x0,"green","green");
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,0,1,"black","black",ipegenerator::STROKE_AND_FILL);
    fig.set_graduation_parameters(-4,0.5,-4,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_4.ipe");
    fig.save_pdf("example_4.pdf");
}