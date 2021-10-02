//
// Created by julien-damers on 27/09/2021.
//

#include "codac.h"
#include "codac-rob.h"
#include "codac-capd.h"
#include "tools.h"
#include "chrono"
#include "ctc_cn.h"
#include "ipegenerator.h"



using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;

void reachability()
{
    cout << "##########################" << endl;
    cout << "########Example 4#########" << endl;
    cout << "##########################" << endl << endl;

    Interval domain_4(0,15);
    double timestep_4 = 0.01;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    Function f_4("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    // CAPD integration version
    TubeVector a_lie_4 = CAPD_integrateODE(domain_4,f_4,x0,timestep_4);
    cout << "Reference generated " << a_lie_4 << endl;
    TubeVector reference(domain_4,timestep_4,3); // Generate a 3-dimension reference
    reference[0] = a_lie_4[0];
    reference[1] = a_lie_4[1];
    reference[2] = a_lie_4[2];

    double epsilon = 0.01;

    IntervalVector X0(3);
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);

    lie_group_ex4_separator fullSep(&a_lie_4,&X0);

    IntervalVector x({{-4,4},{-4,4}});
    IntervalVector proj({{-6,6},{0,15}});
    SepProj sepProj(fullSep,proj,timestep_4);

    IntervalVector frame({{-4,4},{-4,4}});
    ipegenerator::Figure fig_4(frame,150,150);
    fig_4.set_graduation_parameters(-4, 0.5, -4, 0.5);
    fig_4.set_number_digits_axis_x(1);
    fig_4.set_number_digits_axis_y(1);

    auto start = chrono::steady_clock::now();
    sivia_article(x,sepProj,epsilon,fig_4);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;





    ibex::CtcFwdBwd ctc_sub(*new ibex::Function("x[4]", "y[3]", "z", "x[3]-y[2]-z"));
    ibex::Function phi("x[4]","w[3]","b","( x[1] + cos(b) * (-w[0]) - sin(b) * (-w[1]) - [-0.1,0.1]; x[2] + sin(b) * (-w[0]) + cos(b) * (-w[1]) - [-0.1,0.1]; b - [-0.4,0.4] )");
    CtcFunction ctc_phi(phi);
    codac::CtcEval ctc_eval;
    ctc_eval.preserve_slicing(true);
    ctc_eval.set_fast_mode();



    ContractorNetwork cn_simplified_out_4;
    IntervalVectorVar box_simp_out_4(4);
    IntervalVector& w_simp_out_4 = cn_simplified_out_4.create_interm_var(IntervalVector(3));
    Interval& beta_simp_out_4 = cn_simplified_out_4.create_interm_var( Interval());
    cn_simplified_out_4.add(ctc_eval, {box_simp_out_4[0], w_simp_out_4, reference});
    cn_simplified_out_4.add(ctc_sub, {box_simp_out_4, w_simp_out_4, beta_simp_out_4});
    cn_simplified_out_4.add(ctc_phi, {box_simp_out_4, w_simp_out_4, beta_simp_out_4});
    ctc_cn ctc_cn_out_4(&cn_simplified_out_4, &box_simp_out_4);
    CtcStatic ctc_static(ctc_cn_out_4,true);

    // Initialising output tube
    TubeVector x_total_4(Interval(0,15), timestep_4,3);
    TubeVector x_unreachable(Interval(0,15), timestep_4,3);
    TubeVector x_lie_4(Interval(0,15), timestep_4,3);

    // Generating derivative of output tube
    TFunction tf_4("(t;t; sin(0.4*t))");
    TubeVector v_lie_4(Interval(0,15),timestep_4, tf_4); // derivative of tune to contract
    v_lie_4[0] = cos(x_total_4[2]);
    v_lie_4[1] = sin(x_total_4[2]);
    TubeVector v_lie_4_unreachable(Interval(0,15),timestep_4, tf_4); // derivative of tune to contract




    Interval time_domain = x_total_4.tdomain();
    IntervalVector to_reach({Interval::ALL_REALS,Interval::ALL_REALS,Interval::ALL_REALS});
    Function f_reachable("x[2]","(sqrt((x[0]-0.6)^2+(x[1]-0.45)^2)-0.1)");
    //Function f_unreachable("x[2]","(sqrt((x[0]-0.6)^2+(x[1]-0.6)^2)-0.1)"); uncomment for an unreachable area
    CtcFunction ctc_circle(f_reachable);


    // Main contractor network reachability
    ContractorNetwork cn_total;
    cn_total.add(ctc_static,{x_total_4}); // lie constraint
    cn_total.add(ctc_circle,{to_reach[0],to_reach[1]}); // circle constraint
    cn_total.add(ctc_eval,{time_domain,to_reach,x_total_4,v_lie_4}); // box to reach



    start = chrono::steady_clock::now();
    cn_total.contract(true);
    stop = chrono::steady_clock::now();
    cout << "Reachability ex 4 processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    ctc_static.contract(x_lie_4);


    codac::ColorMap colorMap_lie(codac::InterpolMode::RGB);
    codac::rgb yellow= codac::make_rgb((float)1.,(float)1.,(float)0.);
    codac::rgb magenta= codac::make_rgb((float)1.,(float)0.,(float)1.);
    colorMap_lie.add_color_point(yellow,0);
    colorMap_lie.add_color_point(magenta,1);
    fig_4.set_opacity(5);
    fig_4.draw_tubeVector(&x_lie_4,"x_lie", 0, 1, &colorMap_lie);

    codac::ColorMap colorMap(codac::InterpolMode::RGB);
    codac::rgb red= codac::make_rgb((float)1.,(float)0.,(float)0.);
    codac::rgb green= codac::make_rgb((float)0.,(float)1.,(float)0.);
    colorMap.add_color_point(red,0);
    colorMap.add_color_point(green,1);
    fig_4.set_opacity(5);
    fig_4.draw_tubeVector(&x_total_4,"x_total", 0, 1, &colorMap);

    fig_4.add_layer("reached");
    fig_4.set_current_layer("reached");
    fig_4.set_opacity(100);
    fig_4.set_color_stroke("green");
    fig_4.set_color_fill("green");
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_circle(0.6,0.45,0.1);

    fig_4.draw_tubeVector(&a_lie_4,"a_lie",0,1,"black","black",ipegenerator::STROKE_AND_FILL);

    // To check with unreachable area
    fig_4.add_layer("unreachable");
    fig_4.set_current_layer("unreachable");
    fig_4.set_color_stroke("red");
    fig_4.set_color_fill("red");
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_circle(0.6,0.6,0.1);


    fig_4.set_current_layer("data");
    fig_4.set_opacity(30);
    fig_4.draw_box(X0.subvector(0, 1), "green", "green");
    fig_4.set_opacity(100);
    fig_4.add_layer("text");
    fig_4.set_current_layer("text");
    fig_4.set_color_stroke("black");
    fig_4.set_color_type(ipegenerator::STROKE_ONLY);
    fig_4.draw_text("{$\\mathbb{X}_0$}",0,-0.07,false);
    fig_4.draw_axis("x1", "x2");
    fig_4.save_ipe("reachability_static.ipe");
    fig_4.save_pdf("reachability_static.pdf");


}

int main(int argc, char* argv[])
{
    Tube::enable_syntheses();
    reachability();
}