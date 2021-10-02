//
// Created by julien-damers on 22/07/2021.
//

#include "codac.h"
#include "codac-rob.h"
#include "codac-capd.h"
#include "tools.h"
#include "chrono"


using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;


void example_1(bool lohner_done, bool capd_done)
{


    cout << "##########################" << endl;
    cout << "########Example 1#########" << endl;
    cout << "##########################" << endl << endl;


    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();



    Interval domain_1(0, 4);
    double timestep_1 = 0.01;
    IntervalVector x0_1({{0, 1},
                         {2,  3}});
    Function f_1("x", "y", "(1;-y)");

    // Integration using Lohner contractor
    TubeVector x_lohner_1(domain_1, timestep_1, 2);
    CtcLohner ctc_lohner_1(f_1);
    x_lohner_1.set(x0_1, 0.);
    start = chrono::steady_clock::now();
    ctc_lohner_1.contract(x_lohner_1);
    stop = chrono::steady_clock::now();
    cout << "Lohner integration ex 1 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms"
         << endl;


    // Integration using CAPD
    TubeVector x_capd_1(domain_1, timestep_1, 2);
    start = chrono::steady_clock::now();
    x_capd_1 = CAPD_integrateODE(domain_1, f_1, x0_1, timestep_1);
    stop = chrono::steady_clock::now();

    cout << "CAPD integration ex 1 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    //Intergration using Lie Groups

    IntervalVector x0_lie_1({{0, 0},
                             {1, 1}});
    TubeVector a_lie_1 = CAPD_integrateODE(domain_1, f_1, x0_lie_1, timestep_1);
    TubeVector x_lie_1(domain_1, timestep_1, 2);


    Function phi_1("x[3]", "a[2]", "(x[1] - a[0]; x[2]/a[1])");
    CtcFunction ctc_phi_1(phi_1, x0_1);
    codac::CtcEval ctc_eval;
    ctc_eval.preserve_slicing();
    ctc_eval.set_fast_mode();

    ContractorNetwork cn_out_1;
    IntervalVectorVar box_out_1(3);
    IntervalVector& z_out_1 = cn_out_1.create_interm_var(IntervalVector(2)); // a(t)
    cn_out_1.add(ctc_eval, {box_out_1[0], z_out_1, a_lie_1}); // t = x3
    cn_out_1.add(ctc_phi_1, {box_out_1, z_out_1});
    cn_out_1.set_fixedpoint_ratio(0);
    ctc_cn ctc_cn_out_1(&cn_out_1, &box_out_1);
    CtcStatic ctc_static_1(ctc_cn_out_1,true);

    start = chrono::steady_clock::now();
    ctc_static_1.contract(x_lie_1);
    stop = chrono::steady_clock::now();
    cout << "Lie integration ex 1 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    IntervalVector frame_1({{-1,   5},
                            {-0.5, 3}});
    ipegenerator::Figure fig_1(frame_1, 720, 300);
    fig_1.set_graduation_parameters(-1, 0.5, -0.5, 0.5);
    fig_1.set_number_digits_axis_x(1);
    fig_1.set_number_digits_axis_y(1);

    // Visuals initialization

    fig_1.set_color_stroke("colorBlindInStroke");
    fig_1.set_color_fill("colorBlindInFill");
    fig_1.set_opacity(100);
    fig_1.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_1.draw_tubeVector(&x_lohner_1, "lohner",0, 1);




    fig_1.set_color_stroke("colorBlindOutStroke");
    fig_1.set_color_fill("colorBlindOutFill");
    fig_1.set_opacity(50);
    fig_1.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_1.draw_tubeVector(&x_capd_1, "capd",0, 1);


    fig_1.set_color_stroke("colorBlindMaybeStroke");
    fig_1.set_color_fill("colorBlindMaybeFill");
    fig_1.set_opacity(10);
    fig_1.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_1.draw_tubeVector(&x_lie_1,"lie", 0, 1);

    fig_1.set_opacity(30);
    fig_1.draw_box(x0_1, "green", "green");
    fig_1.set_opacity(100);

    codac::ColorMap colorMap_reference(codac::InterpolMode::RGB);
    codac::rgb black= codac::make_rgb((float)0.,(float)0.,(float)0.);
    codac::rgb colorBlind= codac::make_rgb((float)0.,(float)0.619,(float)0.451);
    colorMap_reference.add_color_point(black,0);
    colorMap_reference.add_color_point(colorBlind,1);
    fig_1.draw_tubeVector(&a_lie_1,"reference", 0, 1,&colorMap_reference);


    fig_1.add_layer("text");
    fig_1.set_current_layer("text");
    fig_1.set_color_stroke("black");
    fig_1.set_color_type(ipegenerator::STROKE_ONLY);
    fig_1.draw_text("{\\Huge$\\mathbb{X}_0$}",0.4,2.5,false);
    fig_1.draw_axis("x1", "x2");
    fig_1.save_ipe("comparison_ex_1.ipe");
    fig_1.save_pdf("comparison_ex_1.pdf");
}


int main(int argc, char* argv[])
{
    Tube::enable_syntheses();
    example_1(true,true);

}