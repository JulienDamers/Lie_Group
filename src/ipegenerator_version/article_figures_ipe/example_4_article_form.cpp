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
    // Generate reference
    Interval domain(0,15);
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01; // precision of SIVIA

    IntervalVector X0(3); // Large initial condition
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);
    IntervalVector x({{-4,4},{-4,4}}); // Space to be explored with SIVIA


    lie_group_ex4_separator fullSep(&a,&X0);
    IntervalVector proj({{-6,6},{-15,0}});
    SepProj sepProj(fullSep,proj,timestep);

    IntervalVector frame({{-4,4},{-4,4}});
    ipegenerator::Figure fig(frame,150,150);

    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig); // Performing SIVIA algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time for test-case 4 continuous version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    fig.set_opacity(30);
    fig.draw_box(X0.subvector(0,1),"green","green");
    fig.set_opacity(100);
    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig.draw_tubeVector(&a,"reference", 0, 1,&colorMap_reference,nullptr,ipegenerator::STROKE_AND_FILL,true);
    fig.add_layer("Text");
    fig.set_current_layer("Text");
    fig.set_color_stroke("black");
    fig.set_color_fill("black");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_text("{$\\normalsize\\mathbb{X}_0$}",0,-0.07,false);
    fig.set_graduation_parameters(-4,1,-4,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x_1","x_2");
    fig.save_ipe("example_4_continuous.ipe");
    fig.save_pdf("example_4_continuous.pdf");
    return;
}

void example_4_discrete_article()
{
    // Generate reference
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
    IntervalVector x({{-4,4},{-4,4}});


    lie_group_ex4_separator* fullSep;
    fullSep = new lie_group_ex4_separator(&a,&X0);

    vector<IntervalVector*> projections{
            new IntervalVector({{-6,6},{-1,-1}}),
            new IntervalVector({{-6,6},{-2,-2}}),
            new IntervalVector({{-6,6},{-3,-3}}),
            new IntervalVector({{-6,6},{-4,-4}}),
            new IntervalVector({{-6,6},{-14,-14}}),
            new IntervalVector({{-6,6},{-15,-15}})
    };
    vector<Sep*> seps;
    for (size_t i=0;i<projections.size();i++)
    {
        SepProj *sepProj = new SepProj(*fullSep,*(projections[i]),epsilon);
        seps.push_back(sepProj);
    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);



    IntervalVector frame({{-4,4},{-4,4}});
    ipegenerator::Figure fig(frame,150,150);

    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    fig.set_opacity(30);
    fig.draw_box(X0.subvector(0,1),"green","green");
    fig.set_opacity(100);
    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig.draw_tubeVector(&a,"reference", 0, 1,&colorMap_reference,nullptr,ipegenerator::STROKE_AND_FILL,true);
    fig.add_layer("Text");
    fig.set_current_layer("Text");
    fig.set_color_stroke("black");
    fig.set_color_fill("black");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_text("{$\\normalsize\\mathbb{X}_0$}",0,-0.07,false);
    fig.draw_text("{$\\normalsize\\mathbb{X}_1$}",1.0,-0.2,false);
    fig.draw_text("{$\\normalsize\\mathbb{X}_2$}",1.95,-0.2,false);
    fig.draw_text("{$\\normalsize\\mathbb{X}_3$}",2.6,0.5,false);
    fig.draw_text("{$\\normalsize\\mathbb{X}_4$}",1.0,2.6,false);
    fig.draw_text("{$\\normalsize\\mathbb{X}_{14}$}",-1.2,-0.5,false);
    fig.draw_text("{$\\normalsize\\mathbb{X}_{15}$}",0.1,-0.5,false);
    fig.set_graduation_parameters(-4,1,-4,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);
    fig.draw_axis("x_1","x_2");
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