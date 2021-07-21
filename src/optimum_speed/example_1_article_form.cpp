//
// Created by julien-damers on 29/06/2021.
//

#include <iostream>
#include <cstdlib>
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

    IntervalVector X0({{0,1},{2,3}}); // The uncertain initial condition

    IntervalVector x({{0,4},{-0.2,4}}); // The space to explore for the set inversion

    // Generate the separator for the forward reach set

    ibex::Function phi("x1","x2","t","(x1-t;x2*exp(t))");
    SepFwdBwd fullSep(phi,X0);
    IntervalVector proj(1);
    proj[0] = Interval(0,5);
    SepProj sepProj(fullSep,proj,0.01);
    double epsilon = 0.1;

    IntervalVector frame({{-0.1,4},{-0.3,4}});
    ipegenerator::Figure fig(frame,300,300);
    fig.set_graduation_parameters(0,0.5,-0.2,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);

    // Visuals initialization

    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    fig.set_opacity(30);
    fig.draw_box(X0,"green","green");
    fig.set_opacity(100);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("\\mathbb{X}_0",0.4,2.5,true);
    fig.draw_axis("x1","x2");
    fig.save_ipe("example_1.ipe");
    fig.save_pdf("example_1.pdf");


    return(0);


}



// Visuals initialization

//ipegenerator::Figure fig(x,800,800);
//fig.set_number_digits_axis_x(0);
//fig.set_number_digits_axis_y(1);
//fig.set_graduation_parameters(0, 0.1, -0.2, 0.1);
//fig.draw_axis("x_1", "x_2");

//auto start = chrono::steady_clock::now();
//sivia_article(x,sepProj,epsilon,fig); // Perform the set inversion algorithm
//auto stop = chrono::steady_clock::now();
//cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

//fig.draw_box(X0,"green","green");
//fig.save_ipe("test.ipe");
//fig.save_pdf("test.pdf");