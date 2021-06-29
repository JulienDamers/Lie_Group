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

    codac::CtcFunction ctc_plus(Function("a", "b", "c", "a+b-c"));
    codac::CtcFunction ctc_a(Function("t", "a[2]", "(t-a[0];1-cos(t)-a[1])"));
    IntervalVector X0({{0, 1},
                       {0, 1}});
    ibex::Function phi("x[3]", "w[2]", "z[2]", "(w[0]; x[1] - z[1] + w[1])");
    ibex::CtcFwdBwd ctc_phi(phi, X0);
    ibex::CtcNotIn ctc_phi_not(phi, X0);
    codac::CtcEval ctc_eval;


    // Complete version

    ContractorNetwork cn_out;
    vector<IntervalVector *> intermediary_iv_out;
    vector<Interval *> intermediary_i_out;
    IntervalVector box_out = IntervalVector(3);
    IntervalVector z_out = IntervalVector(2); // z = a(x_1)
    IntervalVector w_out = IntervalVector(2); // w = a(t+x_1)
    Interval beta_out = Interval(); // beta = t+x1
    intermediary_iv_out.push_back(&z_out);
    intermediary_iv_out.push_back(&w_out);
    intermediary_i_out.push_back(&beta_out);
    cn_out.add(ctc_plus, {box_out[0], box_out[2], beta_out});
    cn_out.add(ctc_a,{beta_out,w_out});   // cn_out.add(ctc_eval, {beta_out, w_out, a});
    cn_out.add(ctc_a,{box_out[0],z_out});  // cn_out.add(ctc_eval, {box_out[0], z_out, a});
    cn_out.add(ctc_phi, {box_out, w_out, z_out});
    ctc_cn ctcCn_out(&cn_out, &box_out, &intermediary_iv_out, &intermediary_i_out);

    ContractorNetwork cn_in;
    vector<IntervalVector *> intermediary_iv_in;
    vector<Interval *> intermediary_i_in;
    IntervalVector box_in = IntervalVector(3);
    IntervalVector z_in = IntervalVector(2);
    IntervalVector w_in = IntervalVector(2);
    Interval beta_in = Interval();
    intermediary_iv_in.push_back(&z_in);
    intermediary_iv_in.push_back(&w_in);
    intermediary_i_in.push_back(&beta_in);
    cn_in.add(ctc_plus, {box_in[0], box_in[2], beta_in});
    cn_in.add(ctc_a,{beta_in,w_in});   // cn_in.add(ctc_eval, {beta_in, w_in, a});
    cn_in.add(ctc_a,{box_in[0],z_in});  // cn_in.add(ctc_eval, {box_in[0], z_in, a}); // t = x3
    cn_in.add(ctc_phi_not, {box_in, w_in, z_in});
    ctc_cn ctcCn_in(&cn_in, &box_in, &intermediary_iv_in, &intermediary_i_in);

    SepCtcPair sep(ctcCn_in, ctcCn_out);
    SepProj sep_proj(sep, Interval(-8, 0), epsilon);

    IntervalVector map({{-1,   10},
                        {-10, 10}});


    beginDrawing();
    VIBesFigMap fig("ex_2_cn_method");
    fig.set_properties(50, 50, 800, 306);
    IntervalVector background_box({{-1,10},{-1,3.2}});
    fig.axis_limits(background_box);
    stack<IntervalVector> s;

    auto start = chrono::steady_clock::now();
    sivia(map,sep_proj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
    drawBox(X0, "g[]");
    fig.show();
    fig.axis_limits(background_box);
    endDrawing();

}