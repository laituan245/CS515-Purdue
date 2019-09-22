GRADIENT_DESCENT = 'gradient_descent'
CYCLIC_COORDINATE_DESCENT = 'cyclic_coordinate_descent'
RANDOM_COORDINATE_DESCENT = 'random_coordinate_descent'

# Configurable param
method = RANDOM_COORDINATE_DESCENT
n = 40

printed_1, printed_2, printed_3 = False, False, False
print('For method {} (n = {})'.format(method, n))
log_fn = method + '_' + str(n) + '.txt'
f = open(log_fn, 'r')
for line in f:
    es = line.strip().split(',')
    try:
        relative_residual, work = float(es[0]), int(es[1])
    except:
        break
    if 10**(-2) < relative_residual and  relative_residual <= 10**(-1) and not printed_1:
        print(line.strip())
        printed_1 = True
    if 10**(-3) < relative_residual and  relative_residual <= 10**(-2) and not printed_2:
        print(line.strip())
        printed_2 = True
    if 10**(-4) < relative_residual and  relative_residual <= 10**(-3) and not printed_3:
        print(line.strip())
        printed_3 = True
f.close()
print('{},{}'.format(relative_residual, work))
