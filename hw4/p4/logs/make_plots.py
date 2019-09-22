import math
import matplotlib.pyplot as plt

GRADIENT_DESCENT = 'gradient_descent'
CYCLIC_COORDINATE_DESCENT = 'cyclic_coordinate_descent'
RANDOM_COORDINATE_DESCENT = 'random_coordinate_descent'
METHODS = [CYCLIC_COORDINATE_DESCENT, RANDOM_COORDINATE_DESCENT, GRADIENT_DESCENT]
method2color = {
    GRADIENT_DESCENT: 'r',
    CYCLIC_COORDINATE_DESCENT: 'b',
    RANDOM_COORDINATE_DESCENT: 'g'
}
method2linestyle = {
    GRADIENT_DESCENT: '-',
    CYCLIC_COORDINATE_DESCENT: '-',
    RANDOM_COORDINATE_DESCENT: '-'
}

# Configurable param
n = 10

def parse_log(log_fn):
    xs, ys = [], []
    f = open(log_fn, 'r')
    for line in f:
        es = line.strip().split(',')
        try:
            relative_residual, work = float(es[0]), int(es[1])
            xs.append(work)
            ys.append(math.log(relative_residual, 10))
        except:
            break
    f.close()
    return xs, ys


for method in METHODS:
    xs, ys = parse_log(method + '_' + str(n) + '.txt')
    plt.plot(xs, ys, c=method2color[method], label=method.replace('_', ' '), linewidth=2)
plt.title('n={}'.format(n))
plt.ylabel('Relative Residual (Log Scale)')
plt.xlabel('Work')
plt.legend(loc='upper right')

plt.show()
