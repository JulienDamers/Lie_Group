//
// Created by julien-damers on 22/07/21.
//

#include "codac.h"
#include "codac-rob.h"
#include "codac-capd.h"
#include "tools.h"
#include "chrono"
#include "ctc_cn.h"



using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;

void example_2()
{
    cout << "##########################" << endl;
    cout << "########Example 2#########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();

    Interval domain_2(0, 15); // Define domain of work on which we want to integrate
    double timestep_2 = 0.1;
    IntervalVector x0_2({{0., 1.},
                         {0., 1.}}); // Define initial condition
    Function f_2("x", "y", "(1;sin(x))"); // Evolution function to integrate


    // Integration using CAPD
    TubeVector x_capd_2 = TubeVector(domain_2, timestep_2, 2);
    start = chrono::steady_clock::now();
    x_capd_2 = CAPD_integrateODE(domain_2, f_2, x0_2, timestep_2);
    stop = chrono::steady_clock::now();
    cout << "CAPD integration ex 2 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;



    // Integration using Lie
    TubeVector a_lie_2 = CAPD_integrateODE(domain_2, f_2, IntervalVector({{0.,0.},{0.,0}}), timestep_2);
    TubeVector x_lie_2(domain_2, timestep_2,2);
    CtcFunction ctc_sub_2(Function("x[3]","b","(x[1]-x[0]-b)"));
    Function phi_2("x[3]", "w[2]", "z[2]", "(w[0]; x[2] - z[1] + w[1])");
    CtcFunction ctc_phi_2(phi_2, x0_2);
    CtcFunction ctc_domain(Function("x","(x)"),a_lie_2.tdomain());
    CtcEval ctc_eval;
    ctc_eval.preserve_slicing();
    ctc_eval.set_fast_mode();


    // Complete version


    ContractorNetwork cn_out_2;
    IntervalVectorVar box_out_2(3);
    IntervalVector& z_out_2 = cn_out_2.create_interm_var(IntervalVector(2)); // z = a(x_1)
    IntervalVector& w_out_2 = cn_out_2.create_interm_var(IntervalVector(2)); // w = a(t+x_1)
    Interval& beta_out_2 = cn_out_2.create_interm_var(Interval()); // beta = t+x1
    cn_out_2.add(ctc_sub_2, {box_out_2, beta_out_2});
    cn_out_2.add(ctc_domain,{beta_out_2});
    cn_out_2.add(ctc_eval,{beta_out_2,w_out_2,a_lie_2});
    cn_out_2.add(ctc_eval,{box_out_2[0],z_out_2,a_lie_2});
    cn_out_2.add(ctc_phi_2, {box_out_2, w_out_2, z_out_2});
    cn_out_2.set_fixedpoint_ratio(0.);
    ctc_cn ctc_cn_out_2(&cn_out_2, &box_out_2);
    CtcStatic ctc_static_2(ctc_cn_out_2, true);

    start = chrono::steady_clock::now();
    ctc_static_2.contract(x_lie_2);
    stop = chrono::steady_clock::now();
    cout << "Lie integration ex 2 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    /*
     *
     *
     *                                Graphics
     *
     */

    IntervalVector frame_2({{-1, 20},
                            {-1, 4}});
    ipegenerator::Figure fig_2(frame_2, 400, 153);
    fig_2.set_graduation_parameters(-1, 5, -1, 0.5);
    fig_2.set_number_digits_axis_x(1);
    fig_2.set_number_digits_axis_y(1);

    // Visuals initialization
    fig_2.set_color_stroke("colorBlindOutStroke");
    fig_2.set_color_fill("colorBlindOutFill");
    fig_2.set_opacity(100);
    fig_2.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_2.draw_tubeVector(&x_capd_2,"capd", 0, 1);

    fig_2.set_color_stroke("colorBlindMaybeStroke");
    fig_2.set_color_fill("colorBlindMaybeFill");
    fig_2.set_opacity(75);
    fig_2.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_2.draw_tubeVector(&x_lie_2,"lie", 0, 1);

    fig_2.set_opacity(30);
    fig_2.draw_box(x0_2, "green", "green");


    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig_2.set_opacity(100);
    fig_2.draw_tubeVector(&a_lie_2,"reference",0,1,&colorMap_reference);
    fig_2.set_color_type(ipegenerator::STROKE_ONLY);
    fig_2.add_layer("text");
    fig_2.set_current_layer("text");
    fig_2.set_color_stroke("black");
    fig_2.set_color_type(ipegenerator::STROKE_ONLY);
    fig_2.draw_text("\\mathbb{X}_0", 0.4, 0.5, true);
    fig_2.draw_axis("x1", "x2");
    fig_2.save_ipe("comparison_ex_2.ipe");
    fig_2.save_pdf("comparison_ex_2.pdf");


}

int main(int argc, char* argv[])
{

    Tube::enable_syntheses();
    example_2();


}