//
// Created by julien-damers on 12/18/20.
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

    Interval domain(0,15);
    codac::TFunction f("x","y","z","x4","(cos(z); sin(z); sin(0.4*x4); 1)");
    double timestep = 0.001;
    IntervalVector x0(4);
    x0.init(Interval(0,0));
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;

    double epsilon=0.1;

    IntervalVector X0_simplified({{-0.1,0.1},{-0.1,0.1},{-0.4,0.4}});


    // Define all necessary contractors
    ibex::CtcFwdBwd ctc_add(*new ibex::Function("x","y","z","x+y-z"));
    ibex::CtcFwdBwd ctc_sub(*new ibex::Function("x","y","z","x-y-z"));


    Variable x_v(4,"x_simp"); // x_simp=(x1,x2,x3,x4)
    Variable w_v(4,"w"); // a(x4)
    Variable beta_v(1,"beta"); // x3-a3(x4)

    ibex::Function phi_simplified(x_v,w_v,beta_v,Return(x_v[0]+cos(beta_v[0])*(-w_v[0])-sin(beta_v[0])*(-w_v[1]),
                                                        x_v[1]+sin(beta_v[0])*(-w_v[0])+cos(beta_v[0])*(-w_v[1]),
                                                        beta_v[0]));
    ibex::CtcFwdBwd ctc_phi_simp(phi_simplified,X0_simplified);
    ibex::CtcNotIn ctc_phi_simp_not(phi_simplified,X0_simplified);

    codac::CtcEval ctc_eval;

    // simplified version
    ContractorNetwork cn_simplified_out;
    vector<IntervalVector*> intermediary_iv_out;
    vector<Interval*> intermediary_i_out;
    IntervalVector box_simp_out = IntervalVector(4);
    IntervalVector w_simp_out = IntervalVector(4);
    Interval beta_simp_out = Interval();
    intermediary_iv_out.push_back(&w_simp_out);
    intermediary_i_out.push_back(&beta_simp_out);
    cn_simplified_out.add(ctc_eval,{box_simp_out[3],w_simp_out,a});
    cn_simplified_out.add(ctc_sub,{box_simp_out[2],w_simp_out[2],beta_simp_out});
    cn_simplified_out.add(ctc_phi_simp,{box_simp_out,w_simp_out,beta_simp_out});
    ctc_cn ctc_cn_out(&cn_simplified_out, &box_simp_out, &intermediary_iv_out, &intermediary_i_out);


    ContractorNetwork cn_simplified_in;
    vector<IntervalVector*> intermediary_iv_in;
    vector<Interval*> intermediary_i_in;
    IntervalVector box_simp_in = IntervalVector(4);
    IntervalVector w_simp_in = IntervalVector(4);
    Interval beta_simp_in = Interval();
    intermediary_iv_in.push_back(&w_simp_in);
    intermediary_i_in.push_back(&beta_simp_in);
    cn_simplified_in.add(ctc_eval,{box_simp_in[3],w_simp_in,a});
    cn_simplified_in.add(ctc_sub,{box_simp_in[2],w_simp_in[2],beta_simp_in});
    cn_simplified_in.add(ctc_phi_simp_not,{box_simp_in,w_simp_in,beta_simp_in});
    ctc_cn ctc_cn_in(&cn_simplified_in, &box_simp_in, &intermediary_iv_in, &intermediary_i_in);


    SepCtcPair sep(ctc_cn_in,ctc_cn_out);
    IntervalVector proj({{-6,6},{1,1}});
    SepProj sep_proj(sep,proj, epsilon);

    IntervalVector map({{-4,4},{-4,4}});

    beginDrawing();
    VIBesFigMap fig("ex_4_cn_method");
    fig.set_properties(50, 50, 800, 800);
    fig.axis_limits(map.subvector(0,1));
    stack<IntervalVector> s;



    auto start = chrono::steady_clock::now();
    sivia(map,sep_proj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
    fig.add_tube(&a,"reference",0,1);
    drawBox(X0_simplified,"g[]");
    fig.show();
    fig.axis_limits(map.subvector(0,1));
    endDrawing();











}