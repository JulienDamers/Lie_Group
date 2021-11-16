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

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();

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

    double epsilon = timestep_4;

    IntervalVector X0(3);
    X0[0] = Interval(-0.1,0.1);
    X0[1] = Interval(-0.1,0.1);
    X0[2] = Interval(-0.4,0.4);

    lie_group_ex4_separator fullSep(&a_lie_4,&X0);

    IntervalVector x({{-4,4},{-4,4}});
    IntervalVector proj({{-6,6},{0,15}});
    SepProj sepProj(fullSep,proj,timestep_4);


    /*
     * Generating contractors related to Lie symmetries constraints
     */

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
    CtcStatic ctc_lie_static(ctc_cn_out_4,true);


    /*
     * Applying only Lie Constraint
     */
    TubeVector x_lie_4(Interval(0,15), timestep_4,3);
    ctc_lie_static.contract(x_lie_4);



    /*
     * Proving an area will be reached no matter our uncertainties
     */
    TubeVector x_reachable_4(Interval(0,15), timestep_4,3);

    Function f_reachable("x[4]","( (x[0]-0.1)^2 + (x[1]-1)^2 )");
    CtcNotIn ctc_reachable(f_reachable,Interval(0,0.75));
    CtcStatic ctc_reachable_static(ctc_reachable,true);

    // Main contractor network reachability in
    ContractorNetwork cn_reachable;
    cn_reachable.add(ctc_lie_static,{x_reachable_4}); // lie constraint
    cn_reachable.add(ctc_reachable_static,{x_reachable_4}); // circle constraint


    start = chrono::steady_clock::now();
    cn_reachable.contract(true); //Should return an empty tube (no trajectory not going through area)
    stop = chrono::steady_clock::now();
    cout << "Reachability for reachable area processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    /*
     * Proving an area will not be reached no matter our uncertainties
     */

    TubeVector x_unreachable_4(Interval(0,15), timestep_4,3);
    // Generating derivative of output tube x_unreachable_4
    TFunction tf_4("(t;t; sin(0.4*t))");
    TubeVector v_lie_4(Interval(0,15),timestep_4, tf_4); // derivative of tube to contract
    v_lie_4[0] = cos(x_unreachable_4[2]);
    v_lie_4[1] = sin(x_unreachable_4[2]);

    Function f_unreachable("x[3]","( (x[0]-1.2)^2 + (x[1]-1.3)^2)");
    CtcFunction ctc_unreachable(f_unreachable, Interval(0,0.01));

    IntervalVector to_reach(3);  // intermediate variable;
    Interval time_to_reach(x_unreachable_4.tdomain());

    ContractorNetwork cn_unreachable;
    cn_unreachable.add(ctc_lie_static,{x_unreachable_4}); // lie constraint
    cn_unreachable.add(ctc_unreachable,{to_reach}); // circle constraint
    cn_unreachable.add(ctc_eval,{time_to_reach, to_reach,x_unreachable_4,v_lie_4});

    start = chrono::steady_clock::now();
    cn_unreachable.contract(true); //Should return an empty tube (no trajectory going through area)
    stop = chrono::steady_clock::now();
    cout << "Reachability for an unreachable area processed in : "
         << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;


    // Loading Flow* integration
    cout << "Flow* integration reachability processed in : "
         <<  6000 << " ms"
         << endl;
    TubeVector x_flow_reachability("/home/julien-damers/Shared_Data/Perso/Thesis/CPP_WORKSPACE/SCHOLAR/Lie_Group/src/comparisons_review_2/flow_star_reachability.tube");



    /*
     *  Showing result of circle constraints
     */

    Tube area_reacheable_tube(domain_4,timestep_4);
    area_reacheable_tube = sqr(x_lie_4[0] - 0.1) + sqr(x_lie_4[1]-1.1) - 0.75;
    cout << "area_reachable_tube: "<< area_reacheable_tube << endl;

    Tube area_unreacheable_tube(domain_4,timestep_4);
    area_unreacheable_tube = sqr(x_lie_4[0] - 0.65) + sqr(x_lie_4[1]-0.63) - 0.01;
    cout << "area_unreachable_tube: "<< area_unreacheable_tube << endl;




    /*
     * Graphics
     */

    IntervalVector frame({{-4,4},{-4,4}});
    ipegenerator::Figure fig_4(frame,150,150);
    fig_4.set_graduation_parameters(-4, 0.5, -4, 0.5);
    fig_4.set_number_digits_axis_x(1);
    fig_4.set_number_digits_axis_y(1);

    fig_4.set_color_stroke("black");
    fig_4.set_color_fill("colorBlind1");
    fig_4.set_opacity(50);
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_tubeVector(&x_flow_reachability,"flow", 0, 1);

    // Drawing tube with Lie constraint
    codac::ColorMap colorMap_lie(codac::InterpolMode::RGB);
    codac::rgb yellow= codac::make_rgb((float)1.,(float)1.,(float)0.);
    codac::rgb magenta= codac::make_rgb((float)1.,(float)0.,(float)1.);
    colorMap_lie.add_color_point(yellow,0);
    colorMap_lie.add_color_point(magenta,1);
    fig_4.set_opacity(5);
    fig_4.draw_tubeVector(&x_lie_4,"x_lie", 0, 1, &colorMap_lie);


    // Drawing area we want to prove reachable
    fig_4.add_layer("reached");
    fig_4.set_current_layer("reached");
    fig_4.set_opacity(100);
    fig_4.set_color_stroke("green");
    fig_4.set_color_fill("green");
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_circle(0.1,1,sqrt(0.75));

    // Drawing forbidden area
    fig_4.add_layer("unreachable");
    fig_4.set_current_layer("unreachable");
    fig_4.set_color_stroke("red");
    fig_4.set_color_fill("red");
    fig_4.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig_4.draw_circle(1.2,1.3,0.1);


    // Draw  reference trajectory
    fig_4.draw_tubeVector(&a_lie_4,"a_lie",0,1,"black","black",ipegenerator::STROKE_AND_FILL);


    start = chrono::steady_clock::now();
    //sivia_article(x,sepProj,epsilon,fig_4);
    stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

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




    IntervalVector frame_reach({{0,15},{-5,55}});
    ipegenerator::Figure fig_reach(frame_reach,540,150);
    fig_reach.set_graduation_parameters(0, 1, -5, 5);
    fig_reach.set_number_digits_axis_x(1);
    fig_reach.set_number_digits_axis_y(1);
    fig_reach.draw_tube(&area_reacheable_tube, "area_reachable","black","black");
    IntervalVector threshold({{0,15},{0,0}});
    fig_reach.draw_box(threshold,"red","red");
    fig_reach.draw_axis("x1", "x2");
    fig_reach.save_ipe("area_reachable_tube.ipe");
    fig_reach.save_pdf("area_reachable_tube.pdf");

    IntervalVector frame_unreach({{0,15},{-5,25}});
    ipegenerator::Figure fig_unreach(frame_unreach,540,150);
    fig_unreach.set_graduation_parameters(0, 1, -5, 5);
    fig_unreach.set_number_digits_axis_x(1);
    fig_unreach.set_number_digits_axis_y(1);
    fig_unreach.draw_tube(&area_unreacheable_tube, "area_unreachable","black","black");
    fig_unreach.draw_box(threshold,"red","red");
    fig_unreach.draw_axis("x1", "x2");
    fig_unreach.save_ipe("area_unreachable_tube.ipe");
    fig_unreach.save_pdf("area_unreachable_tube.pdf");




}

int main(int argc, char* argv[])
{
    Tube::enable_syntheses();
    reachability();
}