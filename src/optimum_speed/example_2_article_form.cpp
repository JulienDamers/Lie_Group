//
// Created by julien-damers on 30/06/2021.
//
#include <iostream>
#include <random>
#include "tools.h"
#include "codac.h"
#include "codac-capd.h"
#include "ipegenerator.h"

#include "chrono"

using namespace std;
using namespace ibex;
using namespace codac;
using namespace pyibex;





int main (int argc, char* argv[])
{
    Tube::enable_syntheses();

    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd fullSep(phi,X0);
    double epsilon = 0.1;

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,400,153);
    fig.set_graduation_parameters(-1,0.5,-1,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    IntervalVector proj(1);
    proj[0] = Interval(0,8);
    SepProj sepProj(fullSep,proj,epsilon);
    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    fig.set_opacity(30);
    fig.draw_box(X0,"green","green");
    fig.set_opacity(100);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("\\Huge\\mathbb{X}_0",0.4,0.5,true);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_2.ipe");
    fig.save_pdf("example_2.pdf");


    return(0);


}
