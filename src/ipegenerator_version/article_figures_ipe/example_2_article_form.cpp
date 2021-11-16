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
    // Generate reference
    Interval domain(0, 8);
    double timestep = 0.01;
    IntervalVector x0({{0., 0.},
                         {0., 0.}});
    Function f("x", "y", "(1;sin(x))");
    TubeVector a = CAPD_integrateODE(domain, f, x0, timestep);

    double epsilon = timestep;

    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1+t;x2+cos(x1)-cos(x1+t) )");
    SepFwdBwd fullSep(phi,X0);
    IntervalVector proj(1);
    proj[0] = Interval(-8,0);
    SepProj sepProj(fullSep,proj,epsilon);

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,150,69);
    fig.set_graduation_parameters(-1,1,-1,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);



    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig);
    auto stop = chrono::steady_clock::now();

    cout << "elapsed time for test-case 2 continuous version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    fig.set_opacity(30);
    fig.draw_box(X0,"green","green");
    fig.set_opacity(100);
    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig.draw_tubeVector(&a,"reference", 0, 1,&colorMap_reference,nullptr,ipegenerator::STROKE_AND_FILL,true);

    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\small$\\mathbb{X}_0$}",0.45,0.45,false);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x_1","x_2");
    fig.save_ipe("example_2_continuous.ipe");
    fig.save_pdf("example_2_continuous.pdf");

    return;
}

void example_2_discrete_article()
{
    Interval domain(0, 8);
    double timestep = 0.01;
    double epsilon = timestep;
    Function f("x", "y", "(1;sin(x))");

    // Trajectory of the corners of the initial condition (for drawing purposes)
    TubeVector a1 = CAPD_integrateODE(domain, f, IntervalVector({{0.0},{0.0}}), timestep);
    TubeVector a2 = CAPD_integrateODE(domain, f, IntervalVector({{1,1},{1,1}}), timestep);
    TubeVector a3 = CAPD_integrateODE(domain, f, IntervalVector({{0.0},{1.0}}), timestep);
    TubeVector a4 = CAPD_integrateODE(domain, f, IntervalVector({{1.0},{0.0}}), timestep);



    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1+t;x2+cos(x1)-cos(x1+t) )");

    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    vector<float> projection_times{0.,-2.,-4.,-6.,-8.};

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
    ipegenerator::Figure fig(frame,150,69);
    fig.set_graduation_parameters(-1,1,-1,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);

    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig);
    auto stop = chrono::steady_clock::now();

    cout << "elapsed time or test-case 2 discrete version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    fig.set_opacity(30);
    fig.draw_box(X0,"green","green");
    fig.set_opacity(100);
    fig.set_color_stroke("black");
    fig.set_color_fill("black");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a1,"a1",0,1);
    fig.set_color_stroke("grey");
    fig.set_color_fill("grey");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a2,"a2",0,1);
    fig.set_color_stroke("grey");
    fig.set_color_fill("grey");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a3,"a3",0,1);
    fig.set_color_stroke("grey");
    fig.set_color_fill("grey");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a4,"a4",0,1);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\small$\\mathbb{X}_0$}",0.45,0.45,false);
    fig.draw_text("{\\small$\\mathbb{X}_2$}",2.45,2.2,false);
    fig.draw_text("{\\small$\\mathbb{X}_4$}",4.45,1.6,false);
    fig.draw_text("{\\small$\\mathbb{X}_6$}",6.45,0.40,false);
    fig.draw_text("{\\small$\\mathbb{X}_8$}",8.45,1.9,false);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x_1","x_2");
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



void example_2_discrete_video(int i)
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
    fig.draw_tubeVector(&a1,"a1",0,1);
    fig.set_color_stroke("green");
    fig.set_color_fill("green");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a2,"a2",0,1);
    fig.set_color_stroke("black");
    fig.set_color_fill("black");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a3,"a3",0,1);
    fig.set_color_stroke("red");
    fig.set_color_fill("red");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a4,"a4",0,1);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\Huge$\\mathbb{X}_0$}",0.45,0.45,false);
    fig.draw_text("{\\Huge$\\mathbb{X}_2$}",2.45,2.3,false);
    fig.draw_text("{\\Huge$\\mathbb{X}_4$}",4.45,1.7,false);
    fig.draw_text("{\\Huge$\\mathbb{X}_6$}",6.45,0.45,false);
    fig.draw_text("{\\Huge$\\mathbb{X}_8$}",8.45,2.2,false);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x_1","x_2");

    string file_name = "img_vid/example_2_reply_rev13_" +  to_string(i);
    fig.save_ipe(file_name+".ipe");
    fig.save_pdf(file_name+".pdf");

    return;
}



int main (int argc, char* argv[])
{
    Tube::enable_syntheses();

    example_2_continuous_article();
    example_2_discrete_article();


    // to create 1000 images to make a video
    //for (int i=0; i< 1000;i++)
    //{
    //    example_2_discrete_video(i);
    //}

    return(0);


}
