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
    codac::CtcEval ctc_eval;


    Interval domain_1(0, 4);
    double timestep_1 = 0.01;
    IntervalVector x0_1({{0, 1},
                         {2,  3}});
    Function f_1("x", "y", "(1;-y)");

    // Integration using Lohner contractor
    TubeVector x_lohner_1(domain_1, timestep_1, 2);
    CtcLohner ctc_lohner_1(f_1);
    x_lohner_1.set(x0_1, 0.);
    try
    {
        start = chrono::steady_clock::now();
        ctc_lohner_1.contract(x_lohner_1);
        stop = chrono::steady_clock::now();
        cout << "Lohner integration ex 1 processed in : "
             << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms"
             << endl;
        lohner_done = true;
    }
    catch ( exception &e )
    {
        lohner_done = false;
        cout << "\n\nException caught!\n" << e.what() << endl;
    }

    TubeVector x_capd_1(domain_1, timestep_1, 2);
    // Integration using CAPD
    try
    {
        start = chrono::steady_clock::now();
        x_capd_1 = CAPD_integrateODE(domain_1, f_1, x0_1, timestep_1);
        stop = chrono::steady_clock::now();
        capd_done = true;
    }
    catch ( exception &e )
    {
        capd_done = false;
        cout << "\n\nException caught!\n" << e.what() << endl;

    }
    cout << "CAPD integration ex 1 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    //Intergration using Lie Groups

    IntervalVector x0_lie_1({{0, 0},
                             {1, 1}});
    TubeVector a_lie_1 = CAPD_integrateODE(domain_1, f_1, x0_lie_1, timestep_1);
    TubeVector x_lie_1(domain_1, timestep_1, 2);


    ibex::Function phi_1("x1", "x2", "t", "a1", "a2", "(x1 - a1; x2/a2)");
    ibex::CtcFwdBwd ctc_phi_1(phi_1, x0_1);


    ContractorNetwork cn_out_1;
    vector<IntervalVector *> intermediary_iv_out_1;
    IntervalVector box_out_1 = IntervalVector(3);
    IntervalVector z_out_1 = IntervalVector(2); // a(t)
    intermediary_iv_out_1.push_back(&z_out_1);
    cn_out_1.add(ctc_eval, {box_out_1[2], z_out_1, a_lie_1}); // t = x3
    cn_out_1.add(ctc_phi_1, {box_out_1, z_out_1});
    cn_out_1.set_fixedpoint_ratio(0);
    ctc_cn ctc_cn_out_1(&cn_out_1, &box_out_1, &intermediary_iv_out_1);

    start = chrono::steady_clock::now();
    ctc_cn_out_1.contract(x_lie_1);
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
    if ( lohner_done )
    {
        fig_1.set_color_stroke("red");
        fig_1.set_color_fill("red");
        fig_1.set_opacity(30);
        fig_1.set_color_type(ipegenerator::STROKE_AND_FILL);
        fig_1.draw_tubeVector(&x_lohner_1, 0, 1);
        lohner_done = false;
    }


    if ( capd_done )
    {
        fig_1.set_color_stroke("blue");
        fig_1.set_color_fill("blue");
        fig_1.set_opacity(30);
        fig_1.set_color_type(ipegenerator::STROKE_AND_FILL);
        fig_1.draw_tubeVector(&x_capd_1, 0, 1);
        capd_done = false;
    }

    fig_1.set_color_stroke("yellow");
    fig_1.set_color_fill("yellow");
    fig_1.set_opacity(30);
    fig_1.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_1.draw_tubeVector(&x_lie_1, 0, 1);



    fig_1.set_opacity(30);
    fig_1.draw_box(x0_1, "green", "green");
    fig_1.add_layer("text");
    fig_1.set_current_layer("text");
    fig_1.set_color_stroke("black");
    fig_1.set_color_type(ipegenerator::STROKE_ONLY);
    fig_1.draw_text("\\mathbb{X}_0", 0.4, 2.5, true);
    fig_1.draw_axis("x1", "x2");
    fig_1.save_ipe("comparison_ex_1.ipe");
    fig_1.save_pdf("comparison_ex_1.pdf");
}

void example_2(bool lohner_done, bool capd_done)
{
    cout << "##########################" << endl;
    cout << "########Example 2#########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();
    codac::CtcEval ctc_eval;

    Interval domain_2(0, 30); // Define domain of work on which we wnat to integrate
    double timestep_2 = 0.1;
    IntervalVector x0_2({{0., 0.2},
                               {0., 1.}}); // Define initial condition
    Function f_2("x", "y", "(1;sin(x))"); // Evolution function to integrate


    // Integration using Lohner

    TubeVector x_lohner_2 = TubeVector(domain_2, timestep_2, 2);
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
    TubeVector x_capd_2 = TubeVector(domain_2, timestep_2, 2);
    try
    {
        start = chrono::steady_clock::now();
        x_capd_2 = CAPD_integrateODE(domain_2, f_2, x0_2, timestep_2);
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
    TubeVector a_lie_2 = CAPD_integrateODE(domain_2, f_2, IntervalVector({{0.,0.},{0.,0}}), timestep_2);
    TubeVector ref_mid_2 = CAPD_integrateODE(domain_2, f_2, x0_2.mid(), timestep_2);

    TubeVector x_lie_2(domain_2, timestep_2,2);
    codac::CtcFunction ctc_sub_2(Function("a", "b", "c", "a-b-c"));
    codac::CtcFunction ctc_a_2(Function("t", "a[2]", "(t-a[0];1-cos(t)-a[1])"));
    ibex::CtcFixPoint ctc_fix_2(ctc_a_2);
    ibex::Function phi_2("x[3]", "w[2]", "z[2]", "(w[0]; x[1] - z[1] + w[1])");
    ibex::CtcFwdBwd ctc_phi_2(phi_2, x0_2);



    // Complete version


    ContractorNetwork cn_out_2;
    vector<IntervalVector *> intermediary_iv_out_2(2);
    vector<Interval *> intermediary_i_out_2(1);
    IntervalVector& box_out_2 = cn_out_2.create_dom(IntervalVector(3));
    IntervalVector& z_out_2 = cn_out_2.create_dom(IntervalVector(2)); // z = a(x_1)
    IntervalVector& w_out_2 = cn_out_2.create_dom(IntervalVector(2)); // w = a(t+x_1)
    Interval beta_out_2 = cn_out_2.create_dom(Interval()); // beta = t+x1
    intermediary_iv_out_2[0] = &z_out_2;
    intermediary_iv_out_2[1] = &w_out_2;
    intermediary_i_out_2[0] = &beta_out_2;
    cn_out_2.add(ctc_sub_2, {box_out_2[0], box_out_2[2], beta_out_2});
    cn_out_2.add(ctc_fix_2,{beta_out_2,w_out_2});
    cn_out_2.add(ctc_fix_2,{box_out_2[0],z_out_2});
    cn_out_2.add(ctc_phi_2, {box_out_2, w_out_2, z_out_2});
    ctc_cn ctc_cn_out_2(&cn_out_2, &box_out_2, &intermediary_iv_out_2, &intermediary_i_out_2);


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
    fig_2.draw_tubeVector(&ref_mid_2,0,1);
    fig_2.set_color_type(ipegenerator::STROKE_ONLY);
    fig_2.add_layer("text");
    fig_2.set_current_layer("text");
    fig_2.set_color_stroke("black");
    fig_2.set_color_type(ipegenerator::STROKE_ONLY);
    fig_2.draw_text("\\mathbb{X}_0", 0.4, 0.5, true);
    fig_2.draw_axis("x1", "x2");
    fig_2.save_ipe("comparison_ex_2.ipe");
    fig_2.save_pdf("comparison_ex_2.pdf");

    for(size_t i=0; i<intermediary_i_out_2.size();i++)
    {
        cout <<intermediary_i_out_2[i] << endl;
        cout <<*(intermediary_i_out_2[i]) << endl;
    }
}

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
    TubeVector x_lie_4(domain_4, timestep_4, 4);


    IntervalVector X0_simplified_4({{-0.1, 0.1},
                                    {-0.1, 0.1},
                                    {-0.4, 0.4}});


    // Define all necessary contractors
    ibex::CtcFwdBwd ctc_sub(*new ibex::Function("x", "y", "z", "x-y-z"));


    Variable x_v(4, "x_simp"); // x_simp=(x1,x2,x3,x4)
    Variable w_v(4, "w"); // a(x4)
    Variable beta_v(1, "beta"); // x3-a3(x4)

    ibex::Function phi_simplified_4(x_v, w_v, beta_v,
                                    Return(x_v[0] + cos(beta_v[0]) * (-w_v[0]) - sin(beta_v[0]) * (-w_v[1]),
                                           x_v[1] + sin(beta_v[0]) * (-w_v[0]) + cos(beta_v[0]) * (-w_v[1]),
                                           beta_v[0]));
    ibex::CtcFwdBwd ctc_phi_simp_4(phi_simplified_4, X0_simplified_4);
    codac::CtcEval ctc_eval;


    // simplified version
    ContractorNetwork cn_simplified_out_4;
    vector<IntervalVector *> intermediary_iv_out_4;
    vector<Interval *> intermediary_i_out_4;
    IntervalVector box_simp_out_4 = IntervalVector(4);
    IntervalVector w_simp_out_4 = IntervalVector(4);
    Interval beta_simp_out_4 = Interval();
    intermediary_iv_out_4.push_back(&w_simp_out_4);
    intermediary_i_out_4.push_back(&beta_simp_out_4);
    cn_simplified_out_4.add(ctc_eval, {box_simp_out_4[3], w_simp_out_4, a_lie_4});
    cn_simplified_out_4.add(ctc_sub, {box_simp_out_4[2], w_simp_out_4[2], beta_simp_out_4});
    cn_simplified_out_4.add(ctc_phi_simp_4, {box_simp_out_4, w_simp_out_4, beta_simp_out_4});
    ctc_cn ctc_cn_out_4(&cn_simplified_out_4, &box_simp_out_4, &intermediary_iv_out_4, &intermediary_i_out_4);

    start = chrono::steady_clock::now();
    ctc_cn_out_4.contract(x_lie_4);
    stop = chrono::steady_clock::now();
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
        fig_4.draw_tubeVector(&x_lohner_4, 0, 1);
        lohner_done = false;
    }
    if ( capd_done )
    {
        fig_4.set_color_stroke("blue");
        fig_4.set_color_fill("blue");
        fig_4.set_opacity(30);
        fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
        fig_4.draw_tubeVector(&x_capd_4, 0, 1);
        capd_done = false;
    }
    fig_4.set_color_stroke("yellow");
    fig_4.set_color_fill("yellow");
    fig_4.set_opacity(30);
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    //fig_4.draw_tubeVector(&x_lie_4, 0, 1);
    fig_4.set_opacity(30);
    fig_4.draw_box(x0_4.subvector(0, 1), "green", "green");
    fig_4.set_opacity(100);
    fig_4.add_layer("text");
    fig_4.set_current_layer("text");
    fig_4.set_color_stroke("black");
    fig_4.set_color_type(ipegenerator::STROKE_ONLY);
    fig_4.draw_text("\\mathbb{X}_0", 0.4, 2.5, true);
    fig_4.draw_axis("x1", "x2");
    fig_4.save_ipe("comparison_ex_4.ipe");
    fig_4.save_pdf("comparison_ex_4.pdf");
}

int main(int argc, char* argv[])
{

    bool run_ex_1 = true;
    bool run_ex_2 = true;
    bool run_ex_3 = false;
    bool run_ex_4 = false;

    Tube::enable_syntheses();




    if (run_ex_1)
    {
        // Example 1
        example_1(true,true);
    }
    if (run_ex_2)
    {
        example_2(true, true);
    }

    if (run_ex_3)
    {
        example_3(true, true);
    }

    if (run_ex_4)
    {
        example_4(true, true);
    }


    




}