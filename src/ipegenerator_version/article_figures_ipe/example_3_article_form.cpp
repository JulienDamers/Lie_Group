//
// Created by julien-damers on 01/07/2021.
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


void example_3_continuous_article()
{
    Interval domain(0,6); // Integration time of the reference, must be > proj
    double timestep = 0.005; // Time step for the creation of the reference
    IntervalVector x0({{0.5,0.5},{0,0}}); // Initial condition for reference
    Function f("x","y","(-x^3-x*y^2+x-y; -y^3-x^2*y+x+y)"); // Evolution function to integrate
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep); // The reference
    cout << "Reference generated " << a << endl;

    IntervalVector x({{-2,2},{-2,2}}); // The space to explore for the set inversion

    double epsilon = 0.01; // define accuracy of paving

    Function phi("x10", "x20","t", "r1", "r2",
                   "((x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))-x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)))-r1;\
             (x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))+x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(  ((1-(1/(3*exp(-2*t)+1)))/0.75)+( ((1/(3*exp(-2*t)+1))*(4+(1/0.75)))-(4/3) ) * (x10^2+x20^2)  ))-r2)");


    Function f_dom("x10", "x20","t", "r1", "r2",
                     "((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)");

    Function f_circle ("i1","i2","( (i1-1.5)^2 + (i2-1.5)^2)");

    CtcFwdBwd ctc_phi(phi);
    CtcFwdBwd ctc_dom(f_dom, Interval::NEG_REALS);
    CtcUnion ctc_phi_int(ctc_phi,ctc_dom);

    CtcFwdBwd ctc_initial(f_circle, Interval(0, 0.04));
    CtcNotIn ctc_not_initial(f_circle, Interval(0, 0.04));


    ContractorNetwork cn_out;
    IntervalVectorVar box_out(3);
    IntervalVector& x_init_out = cn_out.create_interm_var(IntervalVector(2));
    cn_out.add(ctc_phi, {box_out, x_init_out});
    cn_out.add(ctc_initial, {x_init_out});
    CtcCn ctc_cn_out(&cn_out, &box_out);


    ContractorNetwork cn_in;
    IntervalVectorVar box_in(3);
    IntervalVector& x_init_in = cn_in.create_interm_var(IntervalVector(2));
    cn_in.add(ctc_phi_int,{box_in,x_init_in});
    cn_in.add(ctc_not_initial,{x_init_in});
    CtcCn ctc_cn_in(&cn_in,&box_in);



    SepCtcPair sep(ctc_cn_in,ctc_cn_out);
    SepProj sep_proj(sep,Interval(-6,0), epsilon);


    // Visuals initialization
    IntervalVector frame({{-2,2},{-2,2}});
    ipegenerator::Figure fig(frame,150,150);
    fig.set_graduation_parameters(-2,1,-2,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    auto start = chrono::steady_clock::now();
    sivia_article(x,sep_proj,epsilon,fig); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time for test-case 3continuous version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    fig.set_color_fill("green");
    fig.set_color_stroke("green");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.set_opacity(30);
    fig.draw_circle(1.5,1.5,0.2); // draw the initial condition
    fig.set_opacity(100);
    fig.draw_tubeVector(&a,"reference", 0, 1,"black","black",ipegenerator::STROKE_AND_FILL,true); // draws reference
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\Large$\\mathbb{X}_0$}",1.5,1.5,false);
    fig.draw_axis("x_1","x_2");
    fig.save_ipe("example_3_continuous.ipe");
    fig.save_pdf("example_3_continuous.pdf");

    return;


}

void example_3_discrete_article()
{
    Interval domain(0,6); // Integration time of the reference, must be > proj
    double timestep = 0.001; // Time step for the creation of the reference
    IntervalVector x0({{0.5,0.5},{0,0}}); // Initial condition for reference
    Function f("x","y","(-x^3-x*y^2+x-y; -y^3-x^2*y+x+y)"); // Evolution function to integrate
    CtcLohner ctc_lohner(f);
    TubeVector a(domain,timestep, f.image_dim()); // The reference
    a.set(x0,0.);
    ctc_lohner.contract(a);
    cout << "Reference generated " << a << endl;

    IntervalVector x({{-2,2},{-2,2}}); // The space to explore for the set inversion

    double epsilon = 0.01; // define accuracy of paving

    Function phi("x10", "x20","t", "r1", "r2",
                 "((x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))-x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)))-r1;\
             (x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))+x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(  ((1-(1/(3*exp(-2*t)+1)))/0.75)+( ((1/(3*exp(-2*t)+1))*(4+(1/0.75)))-(4/3) ) * (x10^2+x20^2)  ))-r2)");


    Function f_dom("x10", "x20","t", "r1", "r2",
                   "((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)");

    Function f_circle ("i1","i2","( (i1-1.5)^2 + (i2-1.5)^2)");

    CtcFwdBwd ctc_phi(phi);
    CtcFwdBwd ctc_dom(f_dom, Interval::NEG_REALS);
    CtcUnion ctc_phi_int(ctc_phi,ctc_dom);

    CtcFwdBwd ctc_initial(f_circle, Interval(0, 0.04));
    CtcNotIn ctc_not_initial(f_circle, Interval(0, 0.04));


    ContractorNetwork cn_out;
    IntervalVectorVar box_out(3);
    IntervalVector& x_init_out = cn_out.create_interm_var(IntervalVector(2));
    cn_out.add(ctc_phi, {box_out, x_init_out});
    cn_out.add(ctc_initial, {x_init_out});
    CtcCn ctc_cn_out(&cn_out, &box_out);


    ContractorNetwork cn_in;
    IntervalVectorVar box_in(3);
    IntervalVector& x_init_in = cn_in.create_interm_var(IntervalVector(2));
    cn_in.add(ctc_phi_int,{box_in,x_init_in});
    cn_in.add(ctc_not_initial,{x_init_in});
    CtcCn ctc_cn_in(&cn_in,&box_in);

    SepCtcPair sep(ctc_cn_in,ctc_cn_out);

    vector<Interval*> projections{
            new Interval(0),
            new Interval(-0.1),
            new Interval(-0.25),
            new Interval(-0.5),
            new Interval(-0.75),
            new Interval(-1.),
            new Interval(-2.),
            new Interval(-3.),
            new Interval(-4.),
            new Interval(-5.),
            new Interval(-6.)};
    vector<Sep*> seps;
    for (size_t i=0;i<projections.size();i++)
    {
        SepProj *sepProj = new SepProj(sep,*(projections[i]),epsilon);
        seps.push_back(sepProj);
    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);


    // Visuals initialization
    // Visuals initialization
    IntervalVector frame({{-2,2},{-2,2}});
    ipegenerator::Figure fig(frame,150,150);
    fig.set_graduation_parameters(-2,1,-2,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig); // Perform the set inversion algorithm
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time for test-case 3 discrete version: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    fig.set_color_fill("green");
    fig.set_color_stroke("green");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.set_opacity(30);
    fig.draw_circle(1.5,1.5,0.2); // draw the initial condition
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.set_opacity(100);
    fig.draw_circle(0,0,1);
    fig.draw_tubeVector(&a,"reference", 0, 1,"black","black",ipegenerator::STROKE_AND_FILL,true); // draws reference
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\Large$\\mathbb{X}_0$}",1.5,1.5,false);
    fig.draw_text("{\\Large$\\mathbb{X}_{0.1}$}",1,1.45,false);
    fig.draw_text("{\\Large$\\mathbb{X}_{0.25}$}",0.7,1.3,false);
    fig.draw_text("{\\Large$\\mathbb{X}_{0.5}$}",0.3,1.25,false);
    fig.draw_text("{\\Large$\\mathbb{X}_{0.75}$}",0.05,1.15,false);
    fig.draw_text("{\\Large$\\mathbb{X}_1$}",-0.35,1.05,false);
    fig.draw_text("{\\Large$\\mathbb{X}_2$}",-1.,0.6,false);
    fig.draw_text("{\\Large$\\mathbb{X}_3$}",-0.8,-0.9,false);
    fig.draw_text("{\\Large$\\mathbb{X}_4$}",0.,-1.2,false);
    fig.draw_text("{\\Large$\\mathbb{X}_5$}",0.9,-0.7,false);
    fig.draw_text("{\\Large$\\mathbb{X}_6$}",0.9,0.65,false);

    fig.draw_axis("x_1","x_2");
    fig.save_ipe("example_3_discrete.ipe");
    fig.save_pdf("example_3_discrete.pdf");
    return;


}


int main (int argc, char* argv[])
{
    Tube::enable_syntheses();

    example_3_continuous_article();
    example_3_discrete_article();


}