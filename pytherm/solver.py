def minimize(min_f, bounds, fabs=1e-4):
    """

    :param function min_f:
    :param list bounds: bounds for each reaction
    :param float fabs: absolute value of sum abs min_f to complete the optimization
    :return:
    """
    
    ksi = [0] * len(bounds)
    # задание начальных кси наиболее близких к 0 из границ оптимазции
    for i in range(len(bounds)):
        m = bounds[i][0]
        for j in bounds[i]:
            if abs(j) < abs(m):
                m = j
        ksi[i] = m
    iter = 0
    while(1):
        lga = min_f(ksi)
        max_i = 0
        # поиск наиболее далеких от равновесия компонентов
        for i in range(len(lga)):
            if abs(lga[i]) > abs(lga[max_i]):
                max_i = i

        cur_bound = bounds[max_i].copy()
        while(1):
            # если > 0 то нужно уменьшать
            if lga[max_i] > 0:
                buf = ksi[max_i]
                ksi[max_i] = (cur_bound[0] + ksi[max_i]) / 2
                cur_bound[1] = buf
            else:
                buf = ksi[max_i]
                ksi[max_i] = (ksi[max_i] + cur_bound[1]) / 2
                cur_bound[0] = buf
            lga = min_f(ksi)
            if abs(lga[max_i]) < fabs / 10:
                break
        s = 0
        for i in lga:
            s += abs(i)
        if s < fabs:
            break
    return True, ksi, s


def minimize2(min_f, bounds, fabs=1e-5, n_iter=10000, n_iter2=5000):
    """

    :param function min_f:
    :param list bounds: bounds for each reaction
    :param float fabs: absolute value of sum abs min_f to complete the optimization
    :return:
    """

    ksi = [0] * len(bounds)
    # задание начальных кси наиболее близких к 0 из границ оптимазции
    for i in range(len(bounds)):
        m = bounds[i][0]
        for j in bounds[i]:
            if abs(j) < abs(m):
                m = j
        ksi[i] = m
    iter = 0
    while (1):
        lga = min_f(ksi)
        max_i = 0
        # поиск наиболее далеких от равновесия компонентов
        for i in range(len(lga)):
            if abs(lga[i]) > abs(lga[max_i]):
                max_i = i

        cur_bound = bounds[max_i].copy()
        iter_2 = 0
        while (1):
            # если > 0 то нужно уменьшать
            if lga[max_i] > 0:
                buf = ksi[max_i]
                ksi[max_i] = (cur_bound[0] + ksi[max_i]) / 2
                cur_bound[1] = buf
            else:
                buf = ksi[max_i]
                ksi[max_i] = (ksi[max_i] + cur_bound[1]) / 2
                cur_bound[0] = buf
            lga = min_f(ksi)
            if abs(lga[max_i]) < fabs / 10:
                break

            iter_2 += 1
            if iter_2 > n_iter2:
                s = 0
                for i in lga:
                    s += abs(i)
                return False, ksi, s
        s = 0
        for i in lga:
            s += abs(i)
        if s < fabs:
            break
        iter += 1
        if iter > n_iter:
            return False, ksi, s
    return True, ksi, s
