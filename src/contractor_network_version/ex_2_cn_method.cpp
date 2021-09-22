//
// Created by julien-damers on 1/14/21.
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

    double epsilon = 0.1;

    codac::CtcFunction ctc_sub(Function("a", "b", "c", "a-b-c"));
    codac::CtcFunction ctc_a(Function("t", "a[2]", "(t-a[0];1-cos(t)-a[1])"));
    ibex::CtcFixPoint ctc_fix(ctc_a);
    IntervalVector X0({{0, 1},
                       {0, 1}});
    ibex::Function phi("x[3]", "w[2]", "z[2]", "(w[0]; x[1] - z[1] + w[1])");
    ibex::CtcFwdBwd ctc_phi(phi, X0);
    ibex::CtcNotIn ctc_phi_not(phi, X0);
    codac::CtcEval ctc_eval;


    // Complete version

    ContractorNetwork cn_out;
    IntervalVectorVar box_out(3);
    IntervalVector& z_out = cn_out.create_interm_var(IntervalVector(2)); // z = a(x_1)
    IntervalVector& w_out = cn_out.create_interm_var(IntervalVector(2)); // w = a(t+x_1)
    Interval& beta_out = cn_out.create_interm_var(Interval()); // beta = t+x1
    cn_out.add(ctc_sub, {box_out[0], box_out[2], beta_out});
    cn_out.add(ctc_fix,{beta_out, w_out});   // cn_out.add(ctc_eval, {beta_out, w_out, a});
    cn_out.add(ctc_fix,{box_out[0],z_out});  // cn_out.add(ctc_eval, {box_out[0], z_out, a});
    cn_out.add(ctc_phi, {box_out, w_out, z_out});
    ctc_cn ctcCn_out(&cn_out, &box_out);

    ContractorNetwork cn_in;

    IntervalVectorVar box_in(3);
    IntervalVector& z_in = cn_in.create_interm_var(IntervalVector(2));
    IntervalVector& w_in = cn_in.create_interm_var(IntervalVector(2));
    Interval& beta_in = cn_in.create_interm_var(Interval());
    cn_in.add(ctc_sub, {box_in[0], box_in[2], beta_in});
    cn_in.add(ctc_fix,{beta_in,w_in});   // cn_in.add(ctc_eval, {beta_in, w_in, a});
    cn_in.add(ctc_fix,{box_in[0],z_in});  // cn_in.add(ctc_eval, {box_in[0], z_in, a}); // t = x3
    cn_in.add(ctc_phi_not, {box_in, w_in, z_in});
    ctc_cn ctcCn_in(&cn_in, &box_in);
    SepCtcPair sep(ctcCn_in, ctcCn_out);
    SepProj sep_proj(sep, Interval(0, 8), epsilon);

    IntervalVector map({{-1,   10},
                        {-1,3.2}});


    beginDrawing();
    VIBesFigMap fig("ex_2_cn_method");
    fig.set_properties(50, 50, 800, 306);
    IntervalVector background_box({{-1,10},{-1,3.2}});
    fig.axis_limits(background_box);
    stack<IntervalVector> s;

    cout << "X: " << map << endl;
    auto start = chrono::steady_clock::now();
    sivia(map,sep_proj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
    drawBox(X0, "g[]");
    fig.show();
    fig.axis_limits(background_box);
    endDrawing();

}