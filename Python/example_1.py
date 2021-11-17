from collections import deque

from pyibex import *
from codac import *
import time

def example_1_continuous():
    # Generate reference
    domain = Interval(0, 5)
    timestep = 0.01
    x0 = IntervalVector([[0., 0.], [1., 1.]])
    f = Function("x", "y", "(1;-y)")
    ctc_lohner = CtcLohner(f)
    a = TubeVector(domain, timestep, f.nb_arg())
    a.set(x0,0.)
    ctc_lohner.contract(a)
    print("Reference generated ", a)

    epsilon = timestep

    X0 = IntervalVector([[0,1],[2,3]])
    mapping = IntervalVector([[-0.1, 6.5], [-0.2, 3.5]])
    phi = Function("x1", "x2", "t", "(x1+t;x2/exp(t))")
    fullSep = SepFwdBwd(phi,X0)

    sep_proj = SepProj(fullSep, -domain, epsilon)

    beginDrawing()

    fig = VIBesFigMap("Example 1 continuous")
    fig.set_properties(50, 50, 800, 464)
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
    print("elapsed time test-case 1 continuous: ", stop - start, "s")

    fig.draw_box(X0,"#00FF00A3[#00FF00A3]")
    fig.add_tube(a, 'reference', 0, 1)
    fig.set_tube_color(a,"k[k]")
    fig.show()
    fig.axis_limits(backgound_box)
    endDrawing()



def example_1_discrete():
    domain = Interval(0, 5)
    timestep = 0.01
    x0 = IntervalVector([[0., 0.], [1., 1.]])
    f = Function("x", "y", "(1;-y)")
    ctc_lohner = CtcLohner(f)
    a = TubeVector(domain, timestep, f.nb_arg())
    a.set(x0,0.)
    ctc_lohner.contract(a)
    print("Reference generated ", a)

    epsilon = timestep

    X0 = IntervalVector([[0,1],[2,3]])
    mapping = IntervalVector([[-0.1, 6.5], [-0.2, 3.5]])
    phi = Function("x1", "x2", "t", "(x1+t;x2/exp(t))")
    fullSep = SepFwdBwd(phi,X0)


    projections = [Interval(0),
                   Interval(-1),
                   Interval(-2),
                   Interval(-3),
                   Interval(-4),
                   Interval(-5)]

    sep_proj = SepProj(fullSep, projections[0], epsilon)

    for i in range(1, len(projections)):
        sep_proj |= SepProj(fullSep, projections[i], 0.01)

    beginDrawing()

    fig = VIBesFigMap("Example 1 discrete")
    fig.set_properties(50,50,800,368)
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
    print("elapsed time test-case 1 discrete: ", stop - start, "s")

    fig.draw_box(X0,"#00FF00A3[#00FF00A3]")
    fig.add_tube(a, 'reference', 0, 1)
    fig.set_tube_color(a,"k[k]")
    fig.show()
    fig.axis_limits(backgound_box)
    endDrawing()


if __name__ == "__main__":

    example_1_continuous()
    example_1_discrete()