//
// Created by julien-damers on 1/18/21.
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

void example_3(bool lohner_done, bool capd_done)
{


    // Exemple 3
    cout << "##########################" << endl;
    cout << "########Example 3#########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();
    codac::CtcEval ctc_eval;

    Interval domain_3(0, 8);
    double timestep_3 = 0.01;
    IntervalVector x0_3({{0.4,  0.6},
                         {-0.1, 0.1}});
    Function f_3("x", "y", "(-x^3-x*y^2+x-y; -y^3-x^2*y+x+y)");


    // Integration using Lohner
    TubeVector x_lohner_3(domain_3, timestep_3, 2);
    CtcLohner ctc_lohner_3(f_3);
    x_lohner_3.set(x0_3, 0.);
    try
    {
        start = chrono::steady_clock::now();
        ctc_lohner_3.contract(x_lohner_3);
        stop = chrono::steady_clock::now();
        cout << "Lohner integration ex 3 processed in : "
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

    TubeVector x_capd_3(domain_3, timestep_3, 2);
    try
    {
        start = chrono::steady_clock::now();
        x_capd_3 = CAPD_integrateODE(domain_3, f_3, x0_3, timestep_3);
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


    // Integration using Lie

    IntervalVector x0_lie_3({{0.5, 0.5},
                             {0.,  0.}});
    Function f_3_lie("x", "y", "(x^3+x*y^2-x+y; y^3+x^2*y-x-y)");
    TubeVector a_lie_3 = CAPD_integrateODE(domain_3, f_3_lie, x0_lie_3, timestep_3);
    TubeVector x_lie_3(domain_3, timestep_3, 2);

    ibex::Function f_dom_3("x1", "x2", "x3", "z1", "z2", "(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2))");
    ibex::Function phi_3("x1", "x2", "x3", "z1", "z2",
                         "((sqrt(3)*(x1*z1-x2*z2))/sqrt(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2)); \
                                 (sqrt(3)*(x2*z1+x1*z2))/sqrt(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2)))");
    ibex::CtcFwdBwd ctc_phi_3(phi_3, x0_3);
    ibex::CtcFwdBwd ctc_dom_verif_3(f_dom_3, Interval::POS_REALS);
    ibex::Function f_ref_dom_3("t","(t)");
    ibex::CtcFwdBwd ctc_ref_domain_3(f_ref_dom_3,Interval(a_lie_3.tdomain()));


    // Complete version

    ContractorNetwork cn_out_3;
    vector<IntervalVector *> intermediary_iv_out_3;
    IntervalVector box_out_3 = IntervalVector(3);
    IntervalVector z_out_3 = IntervalVector(2); // z = a(x_1)
    intermediary_iv_out_3.push_back(&z_out_3);
    cn_out_3.add(ctc_ref_domain_3,{box_out_3[2]});
    cn_out_3.add(ctc_eval, {box_out_3[2], z_out_3, a_lie_3});
    cn_out_3.add(ctc_dom_verif_3, {box_out_3, z_out_3});
    cn_out_3.add(ctc_phi_3, {box_out_3, z_out_3});
    ctc_cn ctc_cn_out_3(&cn_out_3, &box_out_3, &intermediary_iv_out_3);

    start = chrono::steady_clock::now();
    ctc_cn_out_3.contract(x_lie_3);
    stop = chrono::steady_clock::now();
    cout << "Lie integration ex 3 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    IntervalVector frame_3({{-1.5, 1.5},
                            {-1.5, 1.5}});
    ipegenerator::Figure fig_3(frame_3, 150, 150);
    fig_3.set_graduation_parameters(-1.5, 0.5, -1.5, 0.5);
    fig_3.set_number_digits_axis_x(1);
    fig_3.set_number_digits_axis_y(1);

    // Visuals initialization
    if ( lohner_done )
    {
        fig_3.set_color_stroke("red");
        fig_3.set_color_fill("red");
        fig_3.set_opacity(30);
        fig_3.set_color_type(ipegenerator::STROKE_AND_FILL);
        //fig_3.draw_tubeVector(&x_lohner_3, 0, 1);
        lohner_done = false;
    }
    if ( capd_done )
    {
        fig_3.set_color_stroke("blue");
        fig_3.set_color_fill("blue");
        fig_3.set_opacity(30);
        fig_3.set_color_type(ipegenerator::STROKE_AND_FILL);
        fig_3.draw_tubeVector(&x_capd_3, 0, 1);
        capd_done = false;
    }
    fig_3.set_color_stroke("yellow");
    fig_3.set_color_fill("yellow");
    fig_3.set_opacity(30);
    fig_3.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_3.draw_tubeVector(&x_lie_3, 0, 1);
    fig_3.set_opacity(30);
    fig_3.draw_box(x0_3, "green", "green");
    fig_3.set_color_stroke("black");
    fig_3.set_color_fill("black");
    fig_3.set_opacity(100);
    fig_3.draw_tubeVector(&a_lie_3, 0, 1);
    fig_3.draw_text("\\mathbb{X}_0", 0.4, 2.5, true);
    fig_3.draw_axis("x1", "x2");
    fig_3.save_ipe("comparison_ex_3.ipe");
    fig_3.save_pdf("comparison_ex_3.pdf");
}


int main(int argc, char* argv[])
{

    Tube::enable_syntheses();
    example_3(true,true);

}