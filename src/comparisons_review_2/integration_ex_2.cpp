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

void example_2(bool lohner_done, bool capd_done)
{
    cout << "##########################" << endl;
    cout << "########Example 2#########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();
    codac::CtcEval ctc_eval;

    Interval domain_2_a(0, 35);
    Interval domain_2_x(0, 35); // Define domain of work on which we want to integrate
    double timestep_2 = 0.1;
    IntervalVector x0_2({{0., 1.},
                         {0., 1.}}); // Define initial condition
    Function f_2("x", "y", "(1;sin(x))"); // Evolution function to integrate


    // Integration using Lohner

    TubeVector x_lohner_2 = TubeVector(domain_2_a, timestep_2, 2);
    CtcLohner ctc_lohner_2(f_2);
    x_lohner_2.set(x0_2, 0.);
    try
    {
        start = chrono::steady_clock::now();
        ctc_lohner_2.contract(x_lohner_2);
        stop = chrono::steady_clock::now();
        cout << "Lohner integration ex 2 processed in : "
             << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms"
             << endl;
        lohner_done = true;

    }
    catch ( exception &e )
    {
        lohner_done = false;
        cout << "\n\nException caught!\n" << e.what() << endl;
    }

    // Integration using CAPD
    TubeVector x_capd_2 = TubeVector(domain_2_x, timestep_2, 2);
    try
    {
        start = chrono::steady_clock::now();
        x_capd_2 = CAPD_integrateODE(domain_2_x, f_2, x0_2, timestep_2);
        stop = chrono::steady_clock::now();
        cout << "CAPD integration ex 2 processed in : "
             << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
        capd_done = true;
    }
    catch ( exception &e )
    {
        capd_done = false;
        cout << "\n\nException caught!\n" << e.what() << endl;
    }


    // Integration using Lie
    TubeVector a_lie_2 = CAPD_integrateODE(domain_2_a, f_2, IntervalVector({{0.,0.},{0.,0}}), timestep_2);
    TubeVector ref_mid_2 = CAPD_integrateODE(domain_2_x, f_2, x0_2.mid(), timestep_2);

    TubeVector x_lie_2(domain_2_x, timestep_2,2);
    codac::CtcFunction ctc_sub_2(Function("a", "b", "c", "a-b-c"));
    codac::CtcFunction ctc_a_2(Function("t", "a[2]", "(t-a[0];1-cos(t)-a[1])"));
    ibex::CtcFixPoint ctc_fix_2(ctc_a_2);
    ibex::Function phi_2("x[3]", "w[2]", "z[2]", "(w[0]; x[1] - z[1] + w[1])");
    ibex::CtcFwdBwd ctc_phi_2(phi_2, x0_2);
    codac::CtcFunction ctc_domain(Function("x","(x)"),a_lie_2.tdomain());



    // Complete version


    ContractorNetwork cn_out_2;
    IntervalVectorVar box_out_2(3);
    IntervalVector& z_out_2 = cn_out_2.create_interm_var(IntervalVector(2)); // z = a(x_1)
    IntervalVector& w_out_2 = cn_out_2.create_interm_var(IntervalVector(2)); // w = a(t+x_1)
    Interval& beta_out_2 = cn_out_2.create_interm_var(Interval()); // beta = t+x1
    cn_out_2.add(ctc_sub_2, {box_out_2[0], box_out_2[2], beta_out_2});
    cn_out_2.add(ctc_domain,{beta_out_2});
    cn_out_2.add(ctc_eval,{beta_out_2,w_out_2,a_lie_2});
    cn_out_2.add(ctc_eval,{box_out_2[0],z_out_2,a_lie_2});
    //cn_out_2.add(ctc_fix_2,{beta_out_2,w_out_2});
    //cn_out_2.add(ctc_fix_2,{box_out_2[0],z_out_2});
    cn_out_2.add(ctc_phi_2, {box_out_2, w_out_2, z_out_2});
    cn_out_2.set_fixedpoint_ratio(0.);
    ctc_cn ctc_cn_out_2(&cn_out_2, &box_out_2);

    start = chrono::steady_clock::now();
    ctc_cn_out_2.contract(x_lie_2);
    stop = chrono::steady_clock::now();
    cout << "Lie integration ex 2 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    /*
     *
     *
     *                                Graphics
     *
     */

    IntervalVector frame_2({{-1, 30},
                            {-1, 4}});
    ipegenerator::Figure fig_2(frame_2, 400, 153);
    fig_2.set_graduation_parameters(-1, 5, -1, 0.5);
    fig_2.set_number_digits_axis_x(1);
    fig_2.set_number_digits_axis_y(1);

    // Visuals initialization
    if ( lohner_done )
    {
        fig_2.set_color_stroke("red");
        fig_2.set_color_fill("red");
        fig_2.set_opacity(30);
        fig_2.set_color_type(ipegenerator::STROKE_AND_FILL);
        fig_2.draw_tubeVector(&x_lohner_2, 0, 1);
        lohner_done = false;
    }
    if ( capd_done )
    {
        fig_2.set_color_stroke("blue");
        fig_2.set_color_fill("blue");
        fig_2.set_opacity(30);
        fig_2.set_color_type(ipegenerator::STROKE_AND_FILL);
        fig_2.draw_tubeVector(&x_capd_2, 0, 1);
        capd_done = false;
    }

    fig_2.set_color_stroke("yellow");
    fig_2.set_color_fill("yellow");
    fig_2.set_opacity(30);
    fig_2.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_2.draw_tubeVector(&x_lie_2, 0, 1);
    fig_2.set_opacity(30);
    fig_2.draw_box(x0_2, "green", "green");
    fig_2.set_color_stroke("black");
    fig_2.set_color_fill("black");
    fig_2.set_opacity(100);
    fig_2.draw_tubeVector(&a_lie_2,0,1);
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
    example_2(true,true);


}