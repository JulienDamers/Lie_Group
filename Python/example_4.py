from collections import deque

from pyibex import *
from codac import *
import time
from lie_group_ex4_separator import lie_group_ex4_separator


def example_4_continuous():
    # Generate reference
    domain = Interval(0, 15)
    timestep = 0.00
    x0 = IntervalVector(4)
    x0.init(Interval(0, 0))
    f = Function("x", "y", "z", "x4", "(cos(z); sin(z); sin(0.4*x4); 1)")
    ctc_lohner = CtcLohner(f)
    a = TubeVector(domain, timestep, f.nb_arg())
    a.set(x0, 0.)
    ctc_lohner.contract(a)
    print("Reference generated ", a)

    epsilon = 0.01

    X0 = IntervalVector(3)
    X0[0] = Interval(-0.1, 0.1)
    X0[1] = Interval(-0.1, 0.1)
    X0[2] = Interval(-0.4, 0.4)

    fullSep = lie_group_ex4_separator(a, X0)

    proj = IntervalVector([[-6, 6], [-15, 0]])
    sep_proj = SepProj(fullSep, proj, epsilon)

    mapping = IntervalVector([[-4, 4], [-4, 4]])

    beginDrawing()
    fig = VIBesFigMap("example_4_continuous_py")
    fig.set_properties(50, 50, 800, 800)
    background_box = IntervalVector(mapping)
    fig.axis_limits(background_box)
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
    print("elapsed time test-case 4 continuous: ", stop - start, "s")

    fig.draw_box(X0, "#00FF00A3[#00FF00A3]")
    fig.add_tube(a, 'reference', 0, 1)
    fig.set_tube_color(a, "k[k]")
    fig.show()
    fig.axis_limits(background_box)

    endDrawing()


def example_4_discrete():
    # Generate reference
    domain = Interval(0, 15)
    timestep = 0.001
    x0 = IntervalVector(4)
    x0.init(Interval(0, 0))
    f = Function("x", "y", "z", "x4", "(cos(z); sin(z); sin(0.4*x4); 1)")
    ctc_lohner = CtcLohner(f)
    a = TubeVector(domain, timestep, f.nb_arg())
    a.set(x0, 0.)
    ctc_lohner.contract(a)
    print("Reference generated ", a)

    epsilon = 0.01

    X0 = IntervalVector(3)
    X0[0] = Interval(-0.1, 0.1)
    X0[1] = Interval(-0.1, 0.1)
    X0[2] = Interval(-0.4, 0.4)
    mapping = IntervalVector([[-4, 4], [-4, 4]])
    fullSep = lie_group_ex4_separator(a, X0)

    projections = [IntervalVector([[-6, 6], [0, 0]]),
                   IntervalVector([[-6, 6], [-1, -1]]),
                   IntervalVector([[-6, 6], [-2, -2]]),
                   IntervalVector([[-6, 6], [-3, -3]]),
                   IntervalVector([[-6, 6], [-4, -4]]),
                   IntervalVector([[-6, 6], [-14, -14]]),
                   IntervalVector([[-6, 6], [-15, -15]])
                   ]

    sep_proj = SepProj(fullSep, projections[0], epsilon)

    for i in range(1, len(projections)):
        sep_proj |= SepProj(fullSep, projections[i], epsilon)


    beginDrawing()

    fig = VIBesFigMap("example_4_discrete_py")
    fig.set_properties(50, 50, 800, 800)
    background_box = IntervalVector(mapping)
    fig.axis_limits(background_box)
    stack = deque([mapping])
    start = time.time()
    while len(stack) > 0:
        X = stack.pop()
        X_in = IntervalVector(X)
        X_out = IntervalVector(X)
        print(X)
        sep_proj.separate(X_in, X_out)
        if X_out[0].is_empty() and X_out[1].is_empty():
            print("box is out ",X)
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
    print("elapsed time test-case 4 discrete: ", stop - start, "s")

    fig.draw_box(X0,"#00FF00A3[#00FF00A3]")
    fig.add_tube(a, 'reference', 0, 1)
    fig.set_tube_color(a,"k[k]")
    fig.show()
    fig.axis_limits(background_box)

    endDrawing()

if __name__ == "__main__":
    example_4_continuous()
    example_4_discrete()
