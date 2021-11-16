//
// Created by julien-damers on 29/06/2021.
//

#include <iostream>
#include <random>
#include "tools.h"
#include "codac-capd.h"
#include "ipegenerator.h"

#include "chrono"

using namespace std;
using namespace ibex;
using namespace codac;
using namespace pyibex;


/*
 * Continuous version of the first test-case
 */

void example_1_continuous_article()
{
    Interval domain(0,5); // Integration time of the reference, must be > proj
    double timestep = 0.01; // Time step for the creation of the reference
    IntervalVector x0({{0,0},{1,1}}); // Initial condition for reference
    Function f("x1","x2","(1; -x2)"); // Evolution function to integrate
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep); // The reference
    cout << "Reference generated " << a << endl;

    IntervalVector X0({{0,1},{2,3}}); // The uncertain initial condition

    IntervalVector x({{-0.1,6.5},{-0.2,3.5}}); // The space to explore for the set inversion

    double epsilon = timestep; // define accuracy of paving

    // Generate the separator for the forward reach set
    ibex::Function phi("x1","x2","t","(x1+t;x2/exp(t))"); // define transformation function
    SepFwdBwd fullSep(phi,X0); // Create the separator on Phi with X0 as constraint
    IntervalVector proj(1);
    proj[0] = Interval(-5,0); // Define the interval of time on which we want to integrate
    SepProj sepProj(fullSep,proj,epsilon); // Create the separator object


    // Visuals initialization
    IntervalVector frame({{-0.1,6.5},{-0.3,3.5}});
    ipegenerator::Figure fig(frame,150,87);
    fig.set_graduation_parameters(0,1,0,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);



    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time for test-case 1 continuous version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;


    fig.set_opacity(30);
    fig.draw_box(X0,"green","green"); // draw the initial condition
    fig.set_opacity(100);
    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig.draw_tubeVector(&a,"reference", 0, 1,&colorMap_reference,
                        nullptr,ipegenerator::STROKE_AND_FILL,true); // draws reference
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\Large$\\mathbb{X}_0$}",0.4,2.5,false);
    fig.draw_axis("x_1","x_2");
    fig.save_ipe("example_1_continuous.ipe");
    fig.save_pdf("example_1_continuous.pdf");

    return;
}

/*
 * Discrete version of the first test case
 */

void example_1_discrete_article()
{
    // Generate reference
    Interval domain(0,5);
    double timestep = 0.01;
    IntervalVector x0({{0,0},{1,1}});
    Function f("x1","x2","(1; -x2)");
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;


    IntervalVector X0({{0,1},{2,3}});
    IntervalVector x({{-0.1,6.5},{-0.2,3.5}});
    double epsilon = timestep;

    // Generate the separator for the forward reach set
    ibex::Function phi("x1","x2","t","(x1+t;x2/exp(t))");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);

    // Define the different times on which we want to integrate
    vector<float> projection_times{0.,-1.,-2.,-3.,-4.,-5.};


    vector<Sep*> seps;
    vector<IntervalVector*> projections;
    for (size_t i=0;i<projection_times.size();i++) // Generate the separator for each individual time
    {
        IntervalVector *proj = new IntervalVector(1); // Defining time interval for projection
        (*proj)[0] = Interval(projection_times[i],projection_times[i]);
        projections.push_back(proj);
        SepProj *sepProj = new SepProj(*fullSep,*proj,epsilon);
        seps.push_back(sepProj);

    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep); // Create the union of all separators


    // Visuals initialization
    IntervalVector frame({{-0.1,6.5},{-0.3,3.5}});
    ipegenerator::Figure fig(frame,150,87);
    fig.set_graduation_parameters(0,1,0,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);



    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();

    cout << "elapsed time for test-case 1 discrete version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
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
    fig.draw_text("{\\Large$\\mathbb{X}_0$}",0.4,2.5,false);
    fig.draw_text("{\\Large$\\mathbb{X}_1$}",1.45,0.9,false);
    fig.draw_text("{\\Large$\\mathbb{X}_2$}",2.45,0.3,false);
    fig.draw_text("{\\Large$\\mathbb{X}_3$}",3.45,0.1,false);
    fig.draw_text("{\\Large$\\mathbb{X}_4$}",4.45,0.05,false);
    fig.draw_text("{\\Large$\\mathbb{X}_5$}",5.45,0.0,false);
    fig.draw_axis("x_1","x_2");
    fig.save_ipe("example_1_discrete.ipe");
    fig.save_pdf("example_1_discrete.pdf");

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

    example_1_continuous_article();
    example_1_discrete_article();


}
