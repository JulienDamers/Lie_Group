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


    // Integration using Lohner
    TubeVector x_lohner_4 = TubeVector(domain_4, timestep_4, 4);
    CtcLohner ctc_lohner_4(f_4);
    x_lohner_4.set(x0_4, 0.);
    try
    {
        start = chrono::steady_clock::now();
        ctc_lohner_4.contract(x_lohner_4);
        stop = chrono::steady_clock::now();
        cout << "Lohner integration ex 4 processed in : "
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

    TubeVector x_capd_4 = TubeVector(domain_4, timestep_4, 4);
    try
    {
        start = chrono::steady_clock::now();
        x_capd_4 = CAPD_integrateODE(domain_4, f_4, x0_4, timestep_4);
        stop = chrono::steady_clock::now();
        cout << "CAPD integration ex 3 processed in : "
             << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
        capd_done = true;
    }
    catch ( exception &e )
    {
        capd_done = false;
        cout << "\n\nException caught!\n" << e.what() << endl;

    }


    // Integration using Lie Group

    IntervalVector x0_lie_4({{0.0, 0.0},
                             {0.,  0.},
                             {0.,  0.},
                             {0.,  0.}});
    TubeVector a_lie_4 = CAPD_integrateODE(domain_4, f_4, x0_lie_4, timestep_4);
    TubeVector x_lie_4(Interval(0,15), timestep_4,3);


    IntervalVector X0_simplified_4({{-0.1, 0.1},
                                    {-0.1, 0.1},
                                    {-0.4, 0.4}});


    // Define all necessary contractors
    ibex::CtcFwdBwd ctc_sub(*new ibex::Function("x", "y", "z", "x-y-z"));


    ibex::Variable x_v(4, "x_simp"); // x_simp=(x1,x2,x3,x4)
    ibex::Variable w_v(4, "w"); // a(x4)
    ibex::Variable beta_v(1, "beta"); // x3-a3(x4)

    ibex::Function phi_simplified_4(x_v, w_v, beta_v,
                                    Return(x_v[0] + cos(beta_v[0]) * (-w_v[0]) - sin(beta_v[0]) * (-w_v[1]),
                                           x_v[1] + sin(beta_v[0]) * (-w_v[0]) + cos(beta_v[0]) * (-w_v[1]),
                                           beta_v[0]));
    ibex::CtcFwdBwd ctc_phi_simp_4(phi_simplified_4, X0_simplified_4);
    codac::CtcEval ctc_eval;


    // simplified version
    ContractorNetwork cn_simplified_out_4;
    IntervalVectorVar box_simp_out_4(4);
    IntervalVector& w_simp_out_4 = cn_simplified_out_4.create_interm_var(IntervalVector(4));
    Interval& beta_simp_out_4 = cn_simplified_out_4.create_interm_var( Interval());
    cn_simplified_out_4.add(ctc_eval, {box_simp_out_4[3], w_simp_out_4, a_lie_4});
    cn_simplified_out_4.add(ctc_sub, {box_simp_out_4[2], w_simp_out_4[2], beta_simp_out_4});
    cn_simplified_out_4.add(ctc_phi_simp_4, {box_simp_out_4, w_simp_out_4, beta_simp_out_4});
    ctc_cn ctc_cn_out_4(&cn_simplified_out_4, &box_simp_out_4);

    start = chrono::steady_clock::now();
    ctc_cn_out_4.contract(x_lie_4);
    stop = chrono::steady_clock::now();

    cout << x_lie_4(0) << endl;
    cout << "Lie integration ex 3 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    IntervalVector frame_4({{-4, 4},
                            {-4, 4}});
    ipegenerator::Figure fig_4(frame_4, 150, 150);
    fig_4.set_graduation_parameters(-4, 0.5, -4, 0.5);
    fig_4.set_number_digits_axis_x(1);
    fig_4.set_number_digits_axis_y(1);

    // Visuals initialization
    if ( lohner_done )
    {
        fig_4.set_color_stroke("red");
        fig_4.set_color_fill("red");
        fig_4.set_opacity(30);
        fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
        //fig_4.draw_tubeVector(&x_lohner_4, 0, 1);
        lohner_done = false;
    }
    if ( capd_done )
    {
        fig_4.set_color_stroke("blue");
        fig_4.set_color_fill("blue");
        fig_4.set_opacity(30);
        fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
        //fig_4.draw_tubeVector(&x_capd_4, 0, 1);
        capd_done = false;
    }
    fig_4.set_opacity(30);
    fig_4.set_color_stroke("yellow");
    fig_4.set_color_fill("yellow");
    fig_4.set_opacity(30);
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_tubeVector(&x_lie_4, 0, 1);
    fig_4.set_opacity(30);
    fig_4.draw_box(x0_4.subvector(0, 1), "green", "green");
    fig_4.set_opacity(100);
    fig_4.draw_circle(1.2,1.2,0.1);
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