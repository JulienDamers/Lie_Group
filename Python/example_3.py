from collections import deque

from pyibex import *
from codac import *
import time


class CtcCn(Ctc):

    def __init__(self, cn, box_to_contract_cn):
        Ctc.__init__(self, len(box_to_contract_cn))

        self.cn = cn
        self.box_to_contract_cn = box_to_contract_cn

    def contract(self, box):
        self.cn.reset_interm_vars()
        self.cn.contract({self.box_to_contract_cn: box})

        return box



def example_3_continuous():
    # Generating the reference
    t_domain = Interval(0, 6) # Variable à modifier pour changer le temps d'intégration
    time_step = 0.005
    a = TubeVector(t_domain, time_step, 2)
    x0 = IntervalVector([[0.5, 0.5], [0, 0]])  # Initial condition for the reference
    a.set(x0, 0.)
    f = Function("x", "y", "(-x^3-x*y^2+x-y; -y^3-x^2*y+x+y)")
    ctc_lohner = CtcLohner(f)
    ctc_lohner.contract(a)

    print("Reference generated", a)

    epsilon = 0.01 # SIVIA accuracy

    psi = Function("x10", "x20","t", "r1", "r2",
                 "((x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))-x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)))-r1;\
                 (x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))+x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(  ((1-(1/(3*exp(-2*t)+1)))/0.75)+( ((1/(3*exp(-2*t)+1))*(4+(1/0.75)))-(4/3) ) * (x10^2+x20^2)  ))-r2)")


    f_dom = Function("x10", "x20","t", "r1", "r2",
                     "((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)")


    # Création des contracteur sur nos fonctions
    ctc_psi = CtcFwdBwd(psi)
    ctc_psi_int = CtcFwdBwd(psi) | CtcFwdBwd(f_dom, Interval.NEG_REALS)

    f_circle = Function("i1","i2","( (i1-1.5)^2 + (i2-1.5)^2)")
    ctc_initial = CtcFwdBwd(f_circle, Interval(0, 0.04))
    ctc_not_initial = CtcNotIn(f_circle, Interval(0, 0.04))


    # contracteur extérieur
    cn_out = ContractorNetwork()
    box_out = IntervalVectorVar(3)
    x_init_out = cn_out.create_interm_var(IntervalVector(2))  # x0
    cn_out.add(ctc_psi, [box_out, x_init_out])
    cn_out.add(ctc_initial, [x_init_out]) # x_init in circle
    ctc_cn_out = CtcCn(cn_out, box_out)

    # contracteur intérieur
    cn_in = ContractorNetwork()
    box_in = IntervalVectorVar(3)
    x_init_in = cn_in.create_interm_var(IntervalVector(2))  # x0
    cn_in.add(ctc_psi_int, [box_in, x_init_in])
    cn_in.add(ctc_not_initial, [x_init_in]) # x_init not in circle
    ctc_cn_in = CtcCn(cn_in, box_in)

    sep = SepCtcPair(ctc_cn_in, ctc_cn_out)
    sep_proj = SepProj(sep, -t_domain, 0.01)

    mapping = IntervalVector([[-2, 2], [-2, 2]])

    beginDrawing()
    fig = VIBesFigMap('example_3_continous_py')
    fig.set_properties(50, 50, 800, 800)
    backgound_box = IntervalVector(mapping)
    fig.axis_limits(backgound_box)
    stack = deque([mapping])
    start = time.time()
    while len(stack) > 0:
        X = stack.pop()
        X_in = IntervalVector(X)
        X_out = IntervalVector(X)
        sep_proj.separate(X_in, X_out)
        if X_out[0].is_empty() and X_out[1].is_empty():
            fig.draw_box(X, "#009E73[#56B4E9]")

        elif X_in[0].is_empty() and X_in[1].is_empty():
            fig.draw_box(X, "grey[white]")

        else:
            if X.max_diam() < epsilon:
                fig.draw_box(X, "#D55E00[#CC79A7]")

            else:
                i = X.extr_diam_index(False)
                (X1, X2) = X.bisect(i)
                stack.append(X1)
                stack.append(X2)
    stop = time.time()
    print("elapsed time test-case 3 continuous: ", stop - start, "s")

    fig.draw_circle(1.5, 1.5, 0.2, "#00FF00A3[#00FF00A3]")
    fig.draw_circle(0, 0, 1, "k[]")
    fig.add_tube(a, 'reference', 0, 1)
    fig.set_tube_color(a,"k[k]")
    fig.show()
    fig.axis_limits(backgound_box)
    endDrawing()



def example_3_discrete():
    # Generating the reference
    t_domain = Interval(0, 6)  # Variable à modifier pour changer le temps d'intégration
    time_step = 0.005
    a = TubeVector(t_domain, time_step, 2)
    x0 = IntervalVector([[0.5, 0.5], [0, 0]])  # Initial condition for the reference
    a.set(x0, 0.)
    f = Function("x", "y", "(-x^3-x*y^2+x-y; -y^3-x^2*y+x+y)")
    ctc_lohner = CtcLohner(f)
    ctc_lohner.contract(a)

    print("Reference generated", a)

    epsilon = 0.01  # SIVIA accuracy

    psi = Function("x10", "x20", "t", "r1", "r2",
                   "((x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))-x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)))-r1;\
                   (x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))+x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(  ((1-(1/(3*exp(-2*t)+1)))/0.75)+( ((1/(3*exp(-2*t)+1))*(4+(1/0.75)))-(4/3) ) * (x10^2+x20^2)  ))-r2)")

    f_dom = Function("x10", "x20", "t", "r1", "r2",
                     "((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)")

    # Création des contracteur sur nos fonctions
    ctc_psi = CtcFwdBwd(psi)
    ctc_psi_int = CtcFwdBwd(psi) | CtcFwdBwd(f_dom, Interval.NEG_REALS)

    f_circle = Function("i1", "i2", "( (i1-1.5)^2 + (i2-1.5)^2)")
    ctc_initial = CtcFwdBwd(f_circle, Interval(0, 0.04))
    ctc_not_initial = CtcNotIn(f_circle, Interval(0, 0.04))

    # contracteur extérieur
    cn_out = ContractorNetwork()
    box_out = IntervalVectorVar(3)
    x_init_out = cn_out.create_interm_var(IntervalVector(2))  # x0
    cn_out.add(ctc_psi, [box_out, x_init_out])
    cn_out.add(ctc_initial, [x_init_out])  # x_init in circle
    ctc_cn_out = CtcCn(cn_out, box_out)

    # contracteur intérieur
    cn_in = ContractorNetwork()
    box_in = IntervalVectorVar(3)
    x_init_in = cn_in.create_interm_var(IntervalVector(2))  # x0
    cn_in.add(ctc_psi_int, [box_in, x_init_in])
    cn_in.add(ctc_not_initial, [x_init_in])  # x_init not in circle
    ctc_cn_in = CtcCn(cn_in, box_in)

    sep = SepCtcPair(ctc_cn_in, ctc_cn_out)

    projections = [Interval(-6),
                   Interval(-5),
                   Interval(-4),
                   Interval(-3),
                   Interval(-2),
                   Interval(-1),
                   Interval(-0.75),
                   Interval(-0.5),
                   Interval(-0.25),
                   Interval(-0.1),
                   Interval(-0.)]


    sep_proj = SepProj(sep, projections[0], 0.01)

    for i in range(1, len(projections)):
        sep_proj |= SepProj(sep, projections[i], 0.01)

    mapping = IntervalVector([[-2, 2], [-2, 2]])

    beginDrawing()
    fig = VIBesFigMap('example_3_discrete_py')
    fig.set_properties(50, 50, 800, 800)
    backgound_box = IntervalVector(mapping)
    fig.axis_limits(backgound_box)
    stack = deque([mapping])
    start = time.time()
    while len(stack) > 0:
        X = stack.pop()
        X_in = IntervalVector(X)
        X_out = IntervalVector(X)
        sep_proj.separate(X_in, X_out)
        if X_out[0].is_empty() and X_out[1].is_empty():
            fig.draw_box(X, "#009E73[#56B4E9]")

        elif X_in[0].is_empty() and X_in[1].is_empty():
            fig.draw_box(X, "grey[white]")

        else:
            if X.max_diam() < epsilon:
                fig.draw_box(X, "#D55E00[#CC79A7]")

            else:
                i = X.extr_diam_index(False)
                (X1, X2) = X.bisect(i)
                stack.append(X1)
                stack.append(X2)
    stop = time.time()
    print("elapsed time test-case 3 discrete: ", stop - start, "s")

    fig.draw_circle(1.5, 1.5, 0.2, "#00FF00A3[#00FF00A3]")
    fig.draw_circle(0, 0, 1, "k[]")
    fig.add_tube(a, 'reference', 0, 1)
    fig.set_tube_color(a,"k[k]")
    fig.show()
    fig.axis_limits(backgound_box)
    endDrawing()


if __name__ == "__main__":
    example_3_continuous()
    example_3_discrete()






