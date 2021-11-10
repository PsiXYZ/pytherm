def minimize(min_f, bounds, fabs=1e-5):
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
    return ksi, s