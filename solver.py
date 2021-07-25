def minimize(min_f, bon, fabs=1e-5):
    ksi = [0] * len(bon)
    # задание начальных кси наиболее близких к 0 из границ оптимазции
    for i in range(len(bon)):
        m = bon[i][0]
        for j in bon[i]:
            if abs(j) < abs(m):
                m = j
        ksi[i] = m

    while(1):
        lga = min_f(ksi)
        max_i = 0

        for i in range(len(lga)):
            if abs(lga[i]) > abs(lga[max_i]):
                max_i = i

        cur_bound = bon[max_i].copy()
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