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

int main(int argc, char* argv[])
{

    Tube::enable_syntheses();

    Interval domain(0, 4);
    codac::Function f("x","y","(x^3+x*y^2-x+y; y^3+x^2*y-x-y)");
    double timestep = 0.001;
    IntervalVector x0({{0.5, 0.5},
                       {0, 0}});


    // Create the reference curve using Lohner contractor
    TubeVector a(domain, timestep,2);
    CtcLohner ctc_a(f);
    a.set(x0,0.);
    ctc_a.contract(a);

    // Create the reference curve using CAPD
    //TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << a << endl;

    double epsilon = 0.1;

    IntervalVector X0({{0.4, 0.6},
                       {-0.1, 0.1}});
    ibex::Function f_dom("x1","x2","x3","z1","z2","(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2))");
    ibex::Function phi("x1","x2","x3","z1","z2",
                       "((sqrt(3)*(x1*z1-x2*z2))/sqrt(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2)); \
                                 (sqrt(3)*(x2*z1+x1*z2))/sqrt(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2)))");
    ibex::CtcFwdBwd ctc_phi(phi, X0);
    ibex::CtcNotIn ctc_phi_not(phi, X0);
    ibex::CtcFwdBwd ctc_dom_verif(f_dom, Interval::POS_REALS);
    ibex::CtcNotIn ctc_dom_verif_not(f_dom, Interval::POS_REALS);
    CtcUnion ctc_full(ctc_phi_not,ctc_dom_verif_not);
    codac::CtcEval ctc_eval;
    ctc_eval.preserve_slicing();
    ctc_eval.set_fast_mode();


    // Complete version

    ContractorNetwork cn_out;
    IntervalVectorVar box_out(3);
    IntervalVector& z_out = cn_out.create_interm_var(IntervalVector(2)); // z = a(x_1)
    cn_out.add(ctc_eval, {box_out[2], z_out, a});
    cn_out.add(ctc_dom_verif,{box_out, z_out});
    cn_out.add(ctc_phi, {box_out, z_out});
    ctc_cn ctcCn_out(&cn_out, &box_out);


    ContractorNetwork cn_in;
    IntervalVectorVar box_in(3);
    IntervalVector& z_in = cn_in.create_interm_var(IntervalVector(2));
    cn_in.add(ctc_eval, {box_in[2], z_in, a});
    cn_in.add(ctc_full,{box_in,z_in});
    ctc_cn ctcCn_in(&cn_in, &box_in);



    SepCtcPair sep(ctcCn_in, ctcCn_out);
    SepProj sep_proj(sep, Interval(0, 3), epsilon);



    IntervalVector map({{-1.2,   1.2},
                        {-1.2, 1.2}});


    beginDrawing();
    VIBesFigMap fig("ex_3_cn_method");
    fig.set_properties(50, 50, 800, 800);
    IntervalVector background_box({{-1.2,1.2},{-1.2,1.2}});
    fig.axis_limits(background_box);
    stack<IntervalVector> s;

    auto start = chrono::steady_clock::now();
    sivia(map,sep_proj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
    fig.add_tube(&a, "reference", 0, 1);
    drawBox(X0, "g[g]");
    fig.show();
    fig.axis_limits(background_box);
    endDrawing();

}