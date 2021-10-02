//
// Created by julien-damers on 27/09/21.
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

void example_4(bool lohner_done, bool capd_done)
{
    cout << "##########################" << endl;
    cout << "########Example 4#########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();


    Interval domain_4(0, 15);
    double timestep_4 = 0.01;
    IntervalVector x0_4({{-0.1, 0.1},
                         {-0.1, 0.1},
                         {-0.4, 0.4},
                         {0,    0}});
    Function f_4("x", "y", "z", "t", "(cos(z); sin(z); sin(0.4*t); 1)");

    // Integration using CAPD

    TubeVector x_capd_4 = TubeVector(domain_4, timestep_4, 4);

    start = chrono::steady_clock::now();
    x_capd_4 = CAPD_integrateODE(domain_4, f_4, x0_4, timestep_4);
    stop = chrono::steady_clock::now();
    cout << "CAPD integration ex 4 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;



    // Integration using Lie Group

    IntervalVector x0_lie_4({{0.0, 0.0},
                             {0.,  0.},
                             {0.,  0.},
                             {0.,  0.}});
    TubeVector a_lie_4 = CAPD_integrateODE(domain_4, f_4, x0_lie_4, timestep_4);
    TubeVector reference(domain_4,timestep_4,3); // Generate a 3-dimension reference
    reference[0] = a_lie_4[0];
    reference[1] = a_lie_4[1];
    reference[2] = a_lie_4[2];

    TubeVector x_lie_4(Interval(0,15), timestep_4,3);


    IntervalVector X0_simplified_4({{-0.1, 0.1},
                                    {-0.1, 0.1},
                                    {-0.4, 0.4}});


    // Define all necessary contractors
    ibex::CtcFwdBwd ctc_sub(*new ibex::Function("x[4]", "y[3]", "z", "x[3]-y[2]-z"));
    ibex::Function phi("x[4]","w[3]","b","( x[1] + cos(b) * (-w[0]) - sin(b) * (-w[1]) - [-0.1,0.1]; x[2] + sin(b) * (-w[0]) + cos(b) * (-w[1]) - [-0.1,0.1]; b - [-0.4,0.4] )");
    CtcFunction ctc_phi(phi);
    codac::CtcEval ctc_eval;
    ctc_eval.preserve_slicing();
    ctc_eval.set_fast_mode();


    // simplified version
    ContractorNetwork cn_simplified_out_4;
    IntervalVectorVar box_simp_out_4(4);
    IntervalVector& w_simp_out_4 = cn_simplified_out_4.create_interm_var(IntervalVector(3));
    Interval& beta_simp_out_4 = cn_simplified_out_4.create_interm_var( Interval());
    cn_simplified_out_4.add(ctc_eval, {box_simp_out_4[0], w_simp_out_4, reference});
    cn_simplified_out_4.add(ctc_sub, {box_simp_out_4, w_simp_out_4, beta_simp_out_4});
    cn_simplified_out_4.add(ctc_phi, {box_simp_out_4, w_simp_out_4, beta_simp_out_4});
    ctc_cn ctc_cn_out_4(&cn_simplified_out_4, &box_simp_out_4);
    CtcStatic ctc_static(ctc_cn_out_4,true);

    start = chrono::steady_clock::now();
    ctc_static.contract(x_lie_4);
    stop = chrono::steady_clock::now();

    cout << x_lie_4(0) << endl;
    cout << "Lie integration ex 4 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    IntervalVector frame_4({{-4, 4},
                            {-4, 4}});
    ipegenerator::Figure fig_4(frame_4, 150, 150);
    fig_4.set_graduation_parameters(-4, 0.5, -4, 0.5);
    fig_4.set_number_digits_axis_x(1);
    fig_4.set_number_digits_axis_y(1);


    fig_4.set_color_stroke("colorBlindOutStroke");
    fig_4.set_color_fill("colorBlindOutFill");
    fig_4.set_opacity(30);
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_tubeVector(&x_capd_4,"x_capd", 0, 1);

    fig_4.set_opacity(30);
    fig_4.set_color_stroke("colorBlindMaybeStroke");
    fig_4.set_color_fill("colorBlindMaybeFill");
    fig_4.set_opacity(30);
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_tubeVector(&x_lie_4,"x_lie", 0, 1,true);

    fig_4.set_opacity(30);
    fig_4.draw_box(x0_4.subvector(0, 1), "green", "green");

    fig_4.set_opacity(100);
    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig_4.draw_tubeVector(&a_lie_4,"reference", 0, 1,&colorMap_reference,nullptr,ipegenerator::STROKE_AND_FILL,true);

    fig_4.add_layer("text");
    fig_4.set_current_layer("text");
    fig_4.set_color_stroke("black");
    fig_4.set_color_type(ipegenerator::STROKE_ONLY);
    //fig_4.draw_text("\\mathbb{X}_0", 0.4, 2.5, true);
    fig_4.draw_axis("x1", "x2");
    fig_4.save_ipe("comparison_ex_4.ipe");
    fig_4.save_pdf("comparison_ex_4.pdf");
}

int main(int argc, char* argv[])
{
    Tube::enable_syntheses();
    example_4(true,true);
}