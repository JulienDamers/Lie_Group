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


void example_2_continuous_article()
{
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
    fig.save_ipe("example_2_continuous.ipe");
    fig.save_pdf("example_2_continuous.pdf");

    return;
}

void example_2_discrete_article()
{
    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    double epsilon = 0.1;

    vector<float> projection_times{0.,2.,4.,6.,8.};

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

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,400,153);
    fig.set_graduation_parameters(-1,0.5,-1,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig);
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
    fig.save_ipe("example_2_discrete.ipe");
    fig.save_pdf("example_2_discrete.pdf");

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

    example_2_continuous_article();
    example_2_discrete_article();

    return(0);


}
