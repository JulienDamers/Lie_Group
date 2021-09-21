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

void example_2_proof_reviewer_13(int i)
{
    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    double epsilon = 0.001;

    Interval domain(0, 8);
    double timestep = 0.001;
    // Generate the reference  as we do not have a analytical expression for it
    Function f("x", "y", "(1;sin(x))");
    // CAPD integration version
    TubeVector a1 = CAPD_integrateODE(domain, f, IntervalVector({{0.0},{0.0}}), timestep);
    TubeVector a2 = CAPD_integrateODE(domain, f, IntervalVector({{1,1},{1,1}}), timestep);
    TubeVector a3 = CAPD_integrateODE(domain, f, IntervalVector({{0.0},{1.0}}), timestep);
    TubeVector a4 = CAPD_integrateODE(domain, f, IntervalVector({{1.0},{0.0}}), timestep);

    float k = (float) i;
    IntervalVector proj({{k*8/1000}});
    SepProj sepProj = SepProj(*fullSep,proj,epsilon);

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,400,153);
    fig.set_graduation_parameters(-1,0.5,-1,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    fig.set_opacity(30);
    fig.draw_box(X0,"green","green");
    fig.set_opacity(100);
    fig.set_color_stroke("blue");
    fig.set_color_fill("blue");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a1,0,1);
    fig.set_color_stroke("green");
    fig.set_color_fill("green");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a2,0,1);
    fig.set_color_stroke("black");
    fig.set_color_fill("black");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a3,0,1);
    fig.set_color_stroke("red");
    fig.set_color_fill("red");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a4,0,1);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("\\Huge\\mathbb{X}_0",0.4,0.5,true);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x1","x2");

    string file_name = "img_vid/example_2_reply_rev13_" +  to_string(i);
    fig.save_ipe(file_name+".ipe");
    fig.save_pdf(file_name+".pdf");

    return;
}

void comparison_sivia_capd()
{
    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    double epsilon = 0.1;

    vector<float> projection_times{8.};

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
    fig.save_ipe("comparison_sivia_capd.ipe");
    fig.save_pdf("comparison_sivia_capd.pdf");

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

    //example_2_continuous_article();
    //example_2_discrete_article();
    example_2_proof_reviewer_13(1000);
    //comparison_sivia_capd();

    //for (int i=0; i< 1000;i++)
    //{
    //    example_2_proof_reviewer_13(i);
    //}

    return(0);


}
